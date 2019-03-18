import pydaqmx
import numpy
import ctypes
import matplotlib.pyplot as plt

def sample_data(num_samples, sample_rate_in_hz, num_channels):
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
    # Read the raw data and split it into x/y arrays (the assumed ordering format, consistent with sample_data above,
    #   is such that all x values come first followed by all y values)
    f = open("data.npy", mode='rb+')
    data = numpy.load(f)
    data_series = numpy.split(data, 2)
    x_raw = data_series[0]
    y_raw = data_series[1]

    # Due to a delay between the saturated absorption signal and the voltage ramp, there is a hysteresis effect that
    # causes the signal on the ascending and descending halves of the ramp to occur at different points in ramp as well
    # as having different vertical offsets.
    # The chosen solution here is to filter the data based on the sign of the slope and only accept data when the
    # ramp slope is positive.  There is a 10-index offset due to jaggedness of the ramp voltage.
    ascending_ramp_indices = [i for i in range(len(x_raw) - 10) if x_raw[i] < x_raw[i+10]]
    x1 = [x_raw[i] for i in ascending_ramp_indices]
    y1 = [y_raw[i] for i in ascending_ramp_indices]

    # The laser sometimes has a tendency to unlock, which effectively invalidates that entire period of the ramp.
    # This checks for very large jumps in the absorption signal and only adds ramps to an ongoing data array if no
    # such jumps are detected.
    x2 = []
    y2 = []
    i = 0
    things = 0
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
            things += 1
        i += 1

    print things
    plt.scatter(x2, y2, s=3)
    plt.xlabel('Ramp Voltage (V)')
    plt.ylabel('Preamp Output (V)')
    plt.show()

'''numSamples = 300000
numChannels = 2
sampleRateInHz = 100000.0
sample_data(numSamples, sampleRateInHz, numChannels)'''
process_data()
