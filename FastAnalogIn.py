import pydaqmx
import numpy
import ctypes
import matplotlib.pyplot as plt
from operator import itemgetter
from tqdm import tqdm
import Queue
from sympy.physics.wigner import wigner_6j
import serial
import pyvisa as visa
import time


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
    read_segments = 10
    segment_data = numpy.zeros((num_channels * num_samples / read_segments,), dtype=numpy.float64)
    analog_data = numpy.zeros((0,), dtype=numpy.float64)
    try:
        print "Taking {} samples at {} Hz ({} seconds)".format(num_samples, sample_rate_in_hz,
                                                               1.0 * num_samples / sample_rate_in_hz)

        # Set up connection to Analog In
        min_voltage = -10.0
        max_voltage = 10.0
        pydaqmx.DAQmxCreateTask("", pydaqmx.byref(task_handle))
        pydaqmx.DAQmxCreateAIVoltageChan(task_handle, analog_in_location, "", pydaqmx.DAQmx_Val_Cfg_Default,
                                         min_voltage, max_voltage, pydaqmx.DAQmx_Val_Volts, None)
        pydaqmx.DAQmxCfgSampClkTiming(task_handle, "", sample_rate_in_hz, pydaqmx.DAQmx_Val_Rising,
                                      pydaqmx.DAQmx_Val_ContSamps, num_samples)
        pydaqmx.DAQmxStartTask(task_handle)
        read_timeout = 10.0

        # Set up GPIB communication with Frequency Counter
        rm = visa.ResourceManager()
        inst = rm.open_resource('GPIB::3::INSTR')
        initialize_frequency_counter_state(inst)

        """
        Read all Analog In data as a batch; for info on DAQmxReadAnalogF64() see
            https://www.ge.infn.it/~didomizi/documenti/daq/Ni_Docs/software/daqmxcfunc.chm/daqmxreadanalogf64.html
        The most important info is that it reads {2nd argument} samples per channel into the array {5th argument}.
        This actually interleaves reading from the Analog In and frequency counter by alternating readings between
        them a total of read_segments times.  It reads all the data from the analog in (which takes negligible time)
        and then reads from the frequency counter repeatedly until the segment time elapses.  This needs to be done
        because the GPIB port on the frequency counter can read quickly but doesn't seem to have a buffer
        """
        segment_time = num_samples / sample_rate_in_hz / read_segments
        print "Reading data, split into {} segments each {} seconds long"\
              .format(read_segments, segment_time)
        counter_data = numpy.zeros((0,), dtype=numpy.float64)
        for segment in range(read_segments):
            pydaqmx.DAQmxReadAnalogF64(task_handle, num_samples / read_segments, read_timeout,
                                       pydaqmx.DAQmx_Val_GroupByScanNumber, segment_data, num_channels * num_samples,
                                       read, None)
            analog_data = numpy.append(analog_data, segment_data)
            start_counter = time.time()
            while time.time() - start_counter < segment_time:
                curr_counter_data = inst.query("READ:FREQ?")
                next_freq = float(curr_counter_data) * 1e-6
                counter_data = numpy.append(counter_data, next_freq)

        # Save Analog In data to a file
        f = open("data.npy", mode='wb+')
        numpy.save(f, numpy.append(analog_data, counter_data))
        print "Data written to file"

        print "{} points taken from counter ({} Hz)".format(len(counter_data),
                                                            len(counter_data) / (num_samples / sample_rate_in_hz))

    except pydaqmx.DAQError as err:
        print "DAQmx Error: %s" % err
    finally:
        if task_handle:
            pydaqmx.DAQmxStopTask(task_handle)
            pydaqmx.DAQmxClearTask(task_handle)


def initialize_frequency_counter_state(counter):
    """
    Send all the relevant commands over GPIB to the frequency counter in order to optimize throughput, see
        https://literature.cdn.keysight.com/litweb/pdf/53131-90044.pdf?id=1000000328-1:epsg:man
    on section "To Optimize Throughput" (page 3-73) for details on what each command does
    :param counter: A reference to the GPIB connection to the frequency counter as given by pyvisa
    """
    counter.write("*RST")
    counter.write("*CLS")
    counter.write("*SRE 0")
    counter.write("*ESE 0")
    counter.write(":STAT:PRES")
    counter.write(":FORMAT ASCII")
    counter.write(":FUNC 'FREQ 1'")
    counter.write(":EVENT1:LEVEL 0")
    counter.write(":FREQ:ARM:STAR:SOUR IMM")
    counter.write(":FREQ:ARM:STOP:SOUR IMM")
    counter.write(":ROSC:SOUR INT")
    counter.write(":DIAG:CAL:INT:AUTO OFF")
    counter.write(":DISP:ENAB OFF")
    counter.write(":CALC:MATH:STATE OFF")
    counter.write(":CALC2:LIM:STATE OFF")
    counter.write(":CALC3:AVER:STATE OFF")
    counter.write(":HCOPY:CONT OFF")
    counter.write("*DDT #15FETC?")
    counter.write(":INIT:CONT ON")


def combine_raw_data(num_analog_samples):
    """
    Read data from the file below, assumed to be two interleaved channels from the analog input with the given total
    length followed by all of the frequency counter readings assumed to be taken over the same time span.

    :param num_analog_samples:  The total number of samples in the analog input
    :return:
    """
    f = open("data1.npy", mode='rb+')
    all_data = numpy.load(f)

    # Analog input data (index 2n+1 is spectroscopy signal)
    analog_data = all_data[:num_analog_samples]
    num_channels = 2
    samples_per_channel = len(analog_data) / num_channels
    x_analog = numpy.linspace(0, 1, samples_per_channel)
    y_analog = [analog_data[num_channels * i + 1] for i in range(samples_per_channel)]
    print "{} analog points ({} Hz)".format(len(x_analog), len(x_analog))

    # Frequency counter data
    counter_data = all_data[num_analog_samples:]
    x_counter = numpy.linspace(0, 1, len(counter_data))
    print "{} counter points ({} Hz)".format(len(x_counter), len(x_counter))

    """
    Attempt to figure out when the frequency goes negative and adjust accordingly.  This is effectively "undoing"
    an abs() operation on a noisy triangle ramp.  It does this by tracking and toggling some variables based on
    derivative changes and closeness to zero.
    
    What happens on the experiment is that the relative frequency of the two lasers sweeps past zero and is actually
    scanning from a negative value to a positive value, but the frequency counter only ever reads a positive value.
    """
    is_output_inverted = False
    index_margin = 10
    zero_margin = 50
    data_slope_increasing = counter_data[index_margin] - counter_data[0] > 0
    unfolded_counter_data = numpy.zeros((len(counter_data),), dtype=numpy.float64)
    for i in range(index_margin):
        unfolded_counter_data[i] = counter_data[i]
    for i in range(index_margin, len(counter_data)):
        curr_data_slope_increasing = counter_data[i] - counter_data[i - index_margin] > 0
        if curr_data_slope_increasing != data_slope_increasing:
            if not data_slope_increasing:
                is_output_inverted = not is_output_inverted
            data_slope_increasing = curr_data_slope_increasing
            if abs(counter_data[i]) < zero_margin:
                for j in range(index_margin / 2 + 1):
                    unfolded_counter_data[i - j] *= -1
        unfolded_counter_data[i] = counter_data[i] * (-1 if is_output_inverted else 1)

    # Combine the counter and analog data using linear interpolation from the counter data (analog_pts >> counter_pts)
    x1c = x_counter[0]
    y1c = unfolded_counter_data[0]
    x2c = x_counter[1]
    y2c = unfolded_counter_data[1]
    counter_index = 1
    for i in range(len(x_analog)):
        if x_analog[i] > x2c:
            counter_index += 1
            x1c = x2c
            y1c = y2c
            x2c = x_counter[counter_index]
            y2c = unfolded_counter_data[counter_index]
        x_analog[i] = (x_analog[i] - x1c) / (x2c - x1c) * (y2c - y1c) + y1c

    return x_analog, y_analog


def process_data(x_raw, y_raw):
    """
    Filter out invalid segments of data based on certain patterns and average together what is left.  The particular
    filtering processes are due to some unstable experimental artifacts specific to the holmium setup and may not be
    needed for general usage on other setups.

    :return: A list of (x,y) tuples containing the averaged data
    """

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
        i += 1
    print "Found %d valid sweeps" % num_ramps

    # Combine and sort all data points from the ramps in ascending x-value order
    xy_pairs = []
    for i in range(len(x2)):
        xy_pairs.append([x2[i], y2[i]])
    xy_pairs.sort(key=itemgetter(0))

    """
    The electronics on the actual saturated absorption setup cause the baseline doppler-broadened curve to be slightly
    different between different data sets, so a baseline curve is calculated by performing a moving average with a
    very large moving window (in this case 10% of the entire span) in order to effectively remove the spectroscopy
    signal from the doppler-broadened curve.
    """
    xb = []
    yb = []
    bg_avg_size = len(x2) / 10
    avg = Queue.Queue(maxsize=bg_avg_size)
    curr_total = 0
    for i in range(bg_avg_size/2):
        avg.put(xy_pairs[i][1])
        curr_total += xy_pairs[i][1]
    for i in tqdm(range(len(x2)), ascii=True):
        adjusted_index = i + bg_avg_size/2
        if avg.full() or adjusted_index > bg_avg_size:
            curr_total -= avg.get()
        if adjusted_index < len(x2):
            avg.put(xy_pairs[adjusted_index][1])
            curr_total += xy_pairs[adjusted_index][1]
        xb.append(xy_pairs[i][0])
        yb.append(curr_total / avg.qsize())

    """
    Now the actual spectroscopy data itself is constructed with a moving average, in this case with a window size of
        num_ramps * ramp_jaggedness
    with the baseline doppler-broadened curve subtracted out.
    """
    x3 = []
    y3 = []
    moving_avg_size = num_ramps * ramp_jaggedness / 2
    avg = Queue.Queue(maxsize=moving_avg_size)
    curr_total = 0
    for i in range(moving_avg_size / 2):
        avg.put(xy_pairs[i][1])
        curr_total += xy_pairs[i][1]
    for i in tqdm(range(len(x2)), ascii=True):
        adjusted_index = i + moving_avg_size / 2
        if avg.full() or adjusted_index > moving_avg_size:
            curr_total -= avg.get()
        if adjusted_index < len(x2):
            avg.put(xy_pairs[adjusted_index][1])
            curr_total += xy_pairs[adjusted_index][1]
        x3.append(xy_pairs[i][0])
        y3.append(curr_total / avg.qsize() - yb[i])

    """
    The signal tends to behave weirdly at the edges due to the way the moving averages are done; trim the edges by
    an amount equal to the BG-subtraction window size 
    """
    x4 = [x3[i] for i in range(bg_avg_size, len(x3) - bg_avg_size)]
    y4 = [y3[i] for i in range(bg_avg_size, len(y3) - bg_avg_size)]

    # Plotting code to show processed data for visual checking
    plt.scatter(x4, y4, s=0.5)
    plt.xlabel('Ramp Voltage (V)')
    plt.ylabel('Preamp Output (V)')
    plt.show()

    return [(x3[i], y3[i]) for i in range(len(x3))]


def saturated_absorption_peaks(Ae, Be, x_scale, x_offset, y_scale):
    # Atomic structure parameters (exact spin values, assumed ground state hyperfine)
    I = 7 / 2.
    Jg = 15 / 2.
    Je = 17 / 2.
    Ag = 800.58
    Bg = -1668

    transitions = []
    for Fg in range(4, 12):
        Fe = Fg + 1
        Kg = Fg*(Fg+1) - I*(I+1) - Jg*(Jg+1)
        Ke = Fe*(Fe+1) - I*(I+1) - Je*(Je+1)
        Ug = Ag*Kg / 2 + Bg * (3/2.*Kg*(Kg+1) - 2*I*(I+1)*Jg*(Jg+1)) / (4*I*(2*I-1)*Jg*(2*Jg-1))
        Ue = Ae*Ke / 2 + Be * (3/2.*Ke*(Ke+1) - 2*I*(I+1)*Je*(Je+1)) / (4*I*(2*I-1)*Je*(2*Je-1))
        peak_height = pow(wigner_6j(Jg, I, Fg, Fe, 1, Je), 2) * ((2*Fg + 1) * (2*Fe + 1))
        transitions.append((x_scale * (Ue - Ug - x_offset), peak_height * y_scale))
    return transitions


def saturated_absorption_signal(transitions, linewidth):
    low_x = min([a[0] for a in transitions]) - 5*linewidth
    high_x = max([a[0] for a in transitions]) + 5*linewidth
    x = numpy.linspace(low_x, high_x, 2000)
    y = [sum([transitions[j][1] / (1 + pow((x[i] - transitions[j][0]) / linewidth, 2))
              for j in range(len(transitions))]) for i in range(len(x))]

    plt.scatter(x, y, s=0.5)
    plt.xlabel('Detuning (MHz)')
    plt.ylabel('Absorption Signal (arb.)')
    plt.show()


def read_rs232():
    ser = serial.Serial(port='COM4', baudrate=9600, bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE,
                        stopbits=serial.STOPBITS_ONE, timeout=10)
    while ser.inWaiting():
        print float(ser.readline().split(' ')[0])


# Read data from the analog in
numSamples = 500000
numChannels = 2
sampleRateInHz = 100000.0
# sample_data(numSamples, sampleRateInHz, numChannels)
freq, signal = combine_raw_data(numSamples * numChannels)
process_data(freq, signal)
# transitions = saturated_absorption_peaks(715.8, 920, -1, 1000, 1)
# saturated_absorption_signal(transitions, 10)
# read_rs232()
