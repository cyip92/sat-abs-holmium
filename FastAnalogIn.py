import pydaqmx
import numpy
import ctypes
import matplotlib.pyplot as plt
from operator import itemgetter


def sample_data(num_samples, sample_rate_in_hz, num_channels):
    """
    Read data from specified analog input channels at with a specified rate and duration.  Default channels to be
    read are 0 and 1, which should be changed via analog_in_location (within the function code) in order to suit
    particular experimental configurations.  The formatting of the taken data is such that all samples from the 1st
    channel will be contiguous, followed by all the samples in the second channel, and so on.  This also saves the
    data into a file "data.npy" in a local directory.

    :param num_samples: Number of data points to take per analog input channel
    :param sample_rate_in_hz: Number of data points taken per second
    :param num_channels: Number of channels to read data from (analog_in_location must be modified if not reading
        from the first two channels)
    :return: Analog input values from the specified channels
    """
    # This assumes the specified channels are the first ones
    analog_in_location = "Dev1/ai0:1"
    task_handle = pydaqmx.TaskHandle()
    read = ctypes.c_long(0)
    data = numpy.zeros((num_channels * num_samples,), dtype=numpy.float64)
    try:
        # Read the data
        print "Taking %d samples at %d Hz (%f seconds)" % (num_samples, sample_rate_in_hz,
                                                           1.0 * num_samples / sample_rate_in_hz)
        min_voltage = -10.0
        max_voltage = 10.0
        pydaqmx.DAQmxCreateTask("", pydaqmx.byref(task_handle))
        pydaqmx.DAQmxCreateAIVoltageChan(task_handle, analog_in_location, "", pydaqmx.DAQmx_Val_Cfg_Default,
                                         min_voltage, max_voltage, pydaqmx.DAQmx_Val_Volts, None)
        pydaqmx.DAQmxCfgSampClkTiming(task_handle, "", sample_rate_in_hz, pydaqmx.DAQmx_Val_Rising,
                                      pydaqmx.DAQmx_Val_ContSamps, num_samples)
        pydaqmx.DAQmxStartTask(task_handle)
        read_timeout = 10.0
        print "Reading data..."
        pydaqmx.DAQmxReadAnalogF64(task_handle, num_samples, read_timeout, pydaqmx.DAQmx_Val_GroupByChannel, data,
                                   num_channels * num_samples, read, None)
        print "Acquired %d points" % read.value

        f = open("data.npy", mode='wb+')
        numpy.save(f, data)
        print "Data written to file"

        return data

    except pydaqmx.DAQError as err:
        print "DAQmx Error: %s" % err
    finally:
        if task_handle:
            pydaqmx.DAQmxStopTask(task_handle)
            pydaqmx.DAQmxClearTask(task_handle)


def process_data():
    """
    Filter out invalid segments of data based on certain patterns and average together what is left.  The particular
    filtering processes are due to some unstable experimental artifacts specific to the holmium setup and may not be
    needed for general usage on other setups.

    :return: A list of (x,y) tuples containing the averaged data
    """
    # Read the raw data and split it into x/y arrays
    f = open("data.npy", mode='rb+')
    data = numpy.load(f)
    data_series = numpy.split(data, 2)
    x_raw = data_series[0]
    y_raw = data_series[1]

    """
    Due to a delay between the saturated absorption signal and the voltage ramp, there is a hysteresis effect that
    causes the signal on the ascending and descending halves of the ramp to occur at different points in ramp as well
    as having different vertical offsets.
    
    The solution here is to filter the data based on the sign of the slope and only accept data when the
    ramp slope is positive.  The ramp_jaggedness variable allows for some tolerance due to the ramp voltage not being
    strictly increasing at all points.
    """
    ramp_jaggedness = 10
    ascending_ramp_indices = [i for i in range(len(x_raw) - ramp_jaggedness) if x_raw[i] < x_raw[i + ramp_jaggedness]]
    x1 = [x_raw[i] for i in ascending_ramp_indices]
    y1 = [y_raw[i] for i in ascending_ramp_indices]

    """
    The laser sometimes has a tendency to unlock, which effectively invalidates that entire period of the ramp.
    This checks for very large jumps in the absorption signal and only adds ramps to an ongoing data array if no
    such jumps are detected.  (The loops are done weirdly, but I'm not entirely sure to make them much better)
    """
    x2 = []
    y2 = []
    x_single = []
    y_single = []
    i = 0
    num_ramps = 0
    while i < len(x1) - 1:
        curr_x = []
        curr_y = []
        ramp_is_valid = True
        while i < len(x1) - 1 and abs(x1[i] - x1[i+1]) < 1:
            ramp_is_valid = ramp_is_valid and (len(curr_y) == 0 or abs(curr_y[len(curr_y) - 1] - y1[i]) < 0.5)
            curr_x.append(x1[i])
            curr_y.append(y1[i])
            i += 1
        if ramp_is_valid:
            x2.extend(curr_x)
            y2.extend(curr_y)
            num_ramps += 1
            if len(x_single) == 0:
                x_single.append(curr_x)
                y_single.append(curr_y)
        i += 1
    print "Found %d valid sweeps" % num_ramps

    """
    Since the x-axis values are also taken from a 16-bit analog input, they are for the most part distinct.  The
    most straightforward way is to calculate a moving average, in this case with a window size of
        num_ramps * ramp_jaggedness
    """
    # Combine and sort all ramps
    xy_pairs = []
    for i in range(len(x2)):
        xy_pairs.append([x2[i], y2[i]])
    xy_pairs.sort(key=itemgetter(0))

    # Perform the moving average
    x3 = []
    y3 = []
    adj_points = num_ramps * ramp_jaggedness / 2
    for i in range(len(x2)):
        i_min = max(0, i - adj_points)
        i_max = min(len(x2), i + adj_points)
        num_points = i_max - i_min
        adj_xy = xy_pairs[i_min:i_max]
        x3.append(xy_pairs[i][0])
        y3.append(sum([adj_xy[j][1] for j in range(num_points)]) / num_points)

    plt.scatter(x_single, y_single, s=0.5)
    plt.scatter(x3, y3, s=0.5)
    plt.xlabel('Ramp Voltage (V)')
    plt.ylabel('Preamp Output (V)')
    plt.show()

    return [(x3[i], y3[i]) for i in range(len(x3))]

'''
# Read data from the analog in
numSamples = 300000
numChannels = 2
sampleRateInHz = 100000.0
sample_data(numSamples, sampleRateInHz, numChannels)'''
process_data()
