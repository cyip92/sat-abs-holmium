import pydaqmx
import numpy
import ctypes
import matplotlib.pyplot as plt
from operator import itemgetter
from tqdm import tqdm
import Queue
from sympy.physics.wigner import wigner_6j
import pyvisa as visa
import time


def sample_data(num_samples, sample_rate_in_hz, num_channels, file_name):
    """
    Read data from specified analog input channels at with a specified rate and duration.  Default channels to be
    read are 0 and 1, which should be changed via analog_in_location (within the function code) in order to suit
    particular experimental configurations.  The formatting of the data are as an array with the following segments:
        header_segment:     (header_length, analog_length, sample_rate_hz, channel_count)
        analog_segment:     An array of length num_samples*num_channels arranged such that the points from all channels
            are interleaved, eg. [ch0/t=0, ch1/t=0, ch0/t=1, ch1/t=1, ... ]
        counter_segment:    A single long array of counter values which, for the purposes of analysis, should be
            assumed to be taken over the same time span as the analog data

    :param num_samples:         Number of data points to take per analog input channel
    :param sample_rate_in_hz:   Number of data points taken per second
    :param num_channels:        Number of channels to read data from (analog_in_location must be modified if not reading
        from the first two channels)
    :param file_name:           Number of the file to save data into (assumed to be local directory)
    :return:                    None (writes data into "data.npy")
    """
    # This assumes the specified channels are the first ones
    analog_in_location = "Dev1/ai0:1"
    task_handle = pydaqmx.TaskHandle()
    read = ctypes.c_long(0)
    read_segments = 50
    header_segment = numpy.array([4, num_samples, sample_rate_in_hz, num_channels])
    segment_data = numpy.zeros((num_channels * num_samples / read_segments,), dtype=numpy.float64)
    analog_segment = numpy.zeros((0,), dtype=numpy.float64)
    try:
        print "Taking {} samples at {} Hz ({} seconds)".format(num_samples,
                                                               sample_rate_in_hz,
                                                               1.0 * num_samples / sample_rate_in_hz)

        # Set up connection to Analog In
        min_voltage = -10.0
        max_voltage = 10.0
        read_timeout = 10.0
        pydaqmx.DAQmxCreateTask("", pydaqmx.byref(task_handle))
        pydaqmx.DAQmxCreateAIVoltageChan(task_handle, analog_in_location, "", pydaqmx.DAQmx_Val_Cfg_Default,
                                         min_voltage, max_voltage, pydaqmx.DAQmx_Val_Volts, None)
        pydaqmx.DAQmxCfgSampClkTiming(task_handle, "", sample_rate_in_hz, pydaqmx.DAQmx_Val_Rising,
                                      pydaqmx.DAQmx_Val_ContSamps, num_samples)
        pydaqmx.DAQmxStartTask(task_handle)

        """
        Set up GPIB communication with Frequency Counter; for other experimental setups that try to use this, you will
        need to use pyvisa to query for what all the available connections are.  In order to do this in the console:
            import pyvisa as visa
            visa.ResourceManager().list_resources()
        
        See https://pyvisa.readthedocs.io/en/master/ for more tips if needed.
        """
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
        because the GPIB port on the frequency counter can read quickly but doesn't seem to have a buffer.
        """
        segment_time = 1.0 * num_samples / sample_rate_in_hz / read_segments
        print "Reading data, split into {} segments each {} seconds long".format(read_segments, segment_time)
        counter_segment = numpy.zeros((0,), dtype=numpy.float64)
        start_time = time.time()
        for segment in range(read_segments):
            start_time += segment_time
            # Read from the Analog Input and append it to the analog data segment
            pydaqmx.DAQmxReadAnalogF64(task_handle,
                                       num_samples / read_segments,
                                       read_timeout,
                                       pydaqmx.DAQmx_Val_GroupByScanNumber,
                                       segment_data,
                                       num_channels * num_samples,
                                       read,
                                       None)
            analog_segment = numpy.append(analog_segment, segment_data)

            # Continually read from the counter until the segment timer elapses
            while time.time() - start_time < segment_time:
                has_data = False
                total_time = 0.02
                meas_start = time.time()
                while time.time() - meas_start < total_time:
                    if not has_data:
                        curr_counter_data = inst.query("READ:FREQ?")
                        next_freq = float(curr_counter_data) * 1e-6     # Convert from Hz to MHz
                        counter_segment = numpy.append(counter_segment, next_freq)
                        has_data = True

        # Combine all data and save it to a file
        f = open(file_name, mode='wb+')
        numpy.save(f, numpy.append(numpy.append(header_segment, analog_segment), counter_segment))
        print "Data written to file"

        print "{} points taken from counter ({} Hz)".format(len(counter_segment),
                                                            len(counter_segment) / (num_samples / sample_rate_in_hz))

    # I'm not sure how helpful the exception messages are, as (surprisingly) I haven't encountered them yet.
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


def combine_raw_data(file_name):
    """
    Parse the data from the given file, reading from the array assuming the same data format convention noted in the
    docstring for sample_data().  It reads the analog data and associates it with the given timebase in the data file
    header, then assumes the frequency counter shares the same time base.  Then, assuming that the analog data is a
    function y(t) and the counter data is x(t), combines them together in order to produce data of the form y(x)

    :param file_name:   Name of the file to read data from, assumed to be in the local directory.
    :return:            x,y such that
    """
    # Read all the data out and parse the header data
    f = open(file_name, mode='rb+')
    all_data = numpy.load(f)
    header_length = int(all_data[0])
    samples_per_channel = int(all_data[1])
    num_channels = int(all_data[3])
    num_analog_samples = int(samples_per_channel * num_channels)
    time_span = 1.0 * all_data[1] / all_data[2]

    # Analog input data (index 1 is spectroscopy signal)
    channel_index = 1
    analog_data = all_data[header_length:(num_analog_samples + header_length)]
    x_analog = numpy.linspace(0, time_span, samples_per_channel)
    y_analog = [analog_data[num_channels * i + channel_index] for i in range(samples_per_channel)]
    print "Reading {} analog points".format(len(x_analog))

    # Voltage ramp data (index 0)
    channel_index = 0
    ramp_data = all_data[header_length:(num_analog_samples + header_length)]
    x_ramp = numpy.linspace(0, time_span, samples_per_channel)
    y_ramp = [ramp_data[num_channels * i + channel_index] for i in range(samples_per_channel)]
    print "Reading {} analog points".format(len(x_analog))

    """
    Parse the frequency counter data with some invalid data rejection.  First, sample the first sample_window gaps
    between successive data points and set a gap threshold such that any adjacent data points separated by more than
    1.5 times the median of these sample gaps are marked as invalid points.
    
    There is also an additional error_threshold value which will also reject data whose absolute values is higher than
    what is given.  The original purpose of this was to reject error readings which by default read something very
    around 9e32.  It can instead be used to reject clearly bogus readings such as those well outside of the counter's
    250 MHz limit.
    """
    raw_counter_data = all_data[(num_analog_samples + header_length):]
    sample_window = 40
    sample_gaps = [abs(raw_counter_data[i + 1] - raw_counter_data[i]) for i in range(sample_window)]
    sample_gaps.sort()
    gap_threshold = 1.5 * sample_gaps[sample_window / 2]
    error_threshold = 350
    x_counter_raw = numpy.linspace(0, time_span, len(raw_counter_data))
    valid_index = []
    for i in range(len(raw_counter_data) - 2):
        if (abs(raw_counter_data[i + 1] - raw_counter_data[i]) < gap_threshold
                and abs(raw_counter_data[i + 2] - raw_counter_data[i]) < 2 * gap_threshold
                and abs(raw_counter_data[i]) < error_threshold):
            valid_index.extend([True])
        else:
            valid_index.extend([False])
    valid_index.extend([False])
    valid_index.extend([False])
    x_counter = [x_counter_raw[i] for i in range(len(raw_counter_data) - 2) if valid_index[i]]
    y_counter = [raw_counter_data[i] for i in range(len(raw_counter_data) - 2) if valid_index[i]]
    print "{} counter points, filtered to {}".format(len(raw_counter_data), len(y_counter))

    # Print average frequency and standard deviation
    avg_freq = sum(raw_counter_data) / len(raw_counter_data)
    variances = [pow(raw_counter_data[i] - avg_freq, 2) for i in range(len(raw_counter_data))]
    stdev = pow(sum(variances) / len(raw_counter_data), 0.5)
    print "Avg. frequency {}, stdev {}".format(sum(raw_counter_data) / len(raw_counter_data), stdev)

    """
    Attempt to figure out when the frequency goes negative and adjust accordingly.  This is effectively "undoing"
    an abs() operation on a noisy triangle ramp.  It does this by tracking and toggling some variables based on
    derivative changes and closeness to zero.
    
    What happens on the experiment is that the relative frequency of the two lasers sweeps past zero and is actually
    scanning from a negative value to a positive value, but the frequency counter only ever reads a positive value.
    """
    is_output_inverted = False
    index_margin = 10
    data_slope_increasing = y_counter[index_margin] - y_counter[0] > 0
    unfolded_counter_data = numpy.zeros((len(y_counter),), dtype=numpy.float64)
    for i in range(index_margin):
        unfolded_counter_data[i] = y_counter[i]
    for i in range(index_margin, len(y_counter)):
        curr_data_slope_increasing = y_counter[i] - y_counter[i - index_margin] > 0
        if curr_data_slope_increasing != data_slope_increasing:
            if not data_slope_increasing:
                is_output_inverted = not is_output_inverted
            data_slope_increasing = curr_data_slope_increasing
        unfolded_counter_data[i] = y_counter[i] * (-1 if is_output_inverted else 1)
    # Go back through and invert any points missed near the turning points
    '''y_counter = unfolded_counter_data
    zero_margin = 50
    for i in range(index_margin, len(x_counter)):
        if abs(y_counter[i]) < zero_margin:
            slope = y_counter[i - 1] - y_counter[i - 3]
            estimate = y_counter[i - 1] + slope / 2.0
            if abs(-y_counter[i] - estimate) < abs(y_counter[i] - estimate):
                y_counter[i] *= -1'''

    '''
    # Combine the counter and analog data using linear interpolation from the counter data (analog_pts >> counter_pts)
    x1c = x_ramp[0]
    y1c = y_ramp[0]
    x2c = x_ramp[1]
    y2c = y_ramp[1]
    counter_index = 1
    x_combined = []
    y_combined = []
    for i in range(len(x_analog)):
        if x_analog[i] > x2c:
            counter_index += 1
            if counter_index == len(x_counter):
                break
            x1c = x2c
            y1c = y2c
            x2c = x_counter[counter_index]
            y2c = y_counter[counter_index]
        if valid_index[counter_index - 1] and valid_index[counter_index]:
            x_combined.extend([(x_analog[i] - x1c) / (x2c - x1c) * (y2c - y1c) + y1c])
            y_combined.extend([y_analog[i]])
    '''

    # Combine ramp and spectroscopy data
    x_combined = []
    y_combined = []
    for i in tqdm(range(len(x_ramp)), ascii=True):
        x_combined.extend([y_ramp[i]])
        y_combined.extend([y_analog[i]])


    # Plot raw data
    '''
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Analog In (V)', color='r')
    ax1.plot(x_analog, y_analog, color='r')
    # ax2 = ax1.twinx()
    # ax2.set_ylabel('Frequency Counter (MHz)', color='b')
    # ax2.plot(x_counter, y_counter, color='b')
    plt.show()

    plt.scatter(x_combined, y_combined, s=0.5)
    plt.xlabel('Beat Frequency (MHz)')
    plt.ylabel('Preamp Output (V)')
    plt.show()
    '''

    return x_combined, y_combined


def read_spectroscopy_data(file_name, start_index, end_index):
    """
    Temporary function for reading out raw spectroscopy data.  Assumes all x values followed by all y values.

    :param file_name:   Name of the file to read data from, assumed to be in the local directory.
    """
    # Read all the data out
    f = open(file_name, mode='rb+')
    all_data = numpy.load(f)
    x_analog = all_data[:len(all_data) / 2]
    y_analog = all_data[len(all_data) / 2:]
    x_analog = x_analog[start_index: end_index]
    y_analog = y_analog[start_index: end_index]

    return x_analog, y_analog


def process_data(x_raw, y_raw, file_name):
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
    ramp_jaggedness = 5
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
    plt.xlabel('Beat Frequency (MHz)')
    plt.ylabel('Preamp Output (V)')
    plt.show()

    # Append both data sets together and save to a file.  Further processing should assume x and y of equal lengths due
    # to the data formatting, which is a list of (x,y) points
    f = open(file_name, mode='wb+')
    numpy.save(f, numpy.append(x4, y4))
    print 'Data saved to file "{}"'.format(file_name)

    return [(x3[i], y3[i]) for i in range(len(x3))]


# Atomic structure parameters (exact spin values, assumed ground state hyperfine)
I = 7 / 2.
Jg = 15 / 2.
Je = 17 / 2.
Ag = 800.583645
Bg = -1668.00527


def get_peak_locations_and_heights(Ae, Be, x_offset, y_scale):
    """
    Returns a list of tuples (peak_freq, peak_height) corresponding to locations and relative heights of peaks in the
    saturated absorption spectrum, given the excited state hyperfine constants as input parameters.

    :param Ae:          Excited state A hyperfine coefficient in MHz (magnetic dipole)
    :param Be:          Excited state B hyperfine coefficient in MHz (electric quadrupole)
    :param x_offset:    Shift to apply to all peaks in order to appropriately fit the data
    :param y_scale:     Multiplicative factor to apply to all peak heights
    :return:            List of tuples (peak_location, peak_height)
    """
    transitions = []
    for Fg in range(4, 12):
        Fe = Fg + 1
        Kg = Fg*(Fg+1) - I*(I+1) - Jg*(Jg+1)
        Ke = Fe*(Fe+1) - I*(I+1) - Je*(Je+1)
        Ug = Ag*Kg / 2 + Bg * (3/2.*Kg*(Kg+1) - 2*I*(I+1)*Jg*(Jg+1)) / (4*I*(2*I-1)*Jg*(2*Jg-1))
        Ue = Ae*Ke / 2 + Be * (3/2.*Ke*(Ke+1) - 2*I*(I+1)*Je*(Je+1)) / (4*I*(2*I-1)*Je*(2*Je-1))
        peak_height = pow(wigner_6j(Jg, I, Fg, Fe, 1, Je), 2) * ((2*Fg + 1) * (2*Fe + 1))
        transitions.append((Ue - Ug - x_offset, peak_height * y_scale))
    return transitions


def get_peak_locations(Ae, Be):
    """
    Returns a list of floats corresponding to locations of peaks in the saturated absorption spectrum, given the excited
    state hyperfine constants as input parameters.  Note that this doesn't return peak heights.

    :param Ae:          Excited state A hyperfine coefficient in MHz (magnetic dipole)
    :param Be:          Excited state B hyperfine coefficient in MHz (electric quadrupole)
    :return:            List of floats
    """
    transitions = []
    for Fg in range(4, 12):
        Fe = Fg + 1
        Kg = Fg * (Fg + 1) - I * (I + 1) - Jg * (Jg + 1)
        Ke = Fe * (Fe + 1) - I * (I + 1) - Je * (Je + 1)
        Ug = Ag * Kg / 2 + Bg * (3 / 2. * Kg * (Kg + 1) - 2 * I * (I + 1) * Jg * (Jg + 1)) / (
                    4 * I * (2 * I - 1) * Jg * (2 * Jg - 1))
        Ue = Ae * Ke / 2 + Be * (3 / 2. * Ke * (Ke + 1) - 2 * I * (I + 1) * Je * (Je + 1)) / (
                    4 * I * (2 * I - 1) * Je * (2 * Je - 1))
        transitions.append(Ue - Ug)
    return transitions


def saturated_absorption_signal(transitions, linewidth, v_lines, x_other=None, y_other=None):
    """
    Produces a plot showing the saturated absorption signal, given a transition list from
    get_peak_locations_and_heights().

    :param transitions:     Transition list from get_peak_locations_and_heights()
    :param linewidth:       Linewidth of the main transition in MHz, scaled by the factors in the transition list for
        each individual peak
    :return:                None; produces a plot instead
    """
    low_x = min([a[0] for a in transitions]) - 5*linewidth
    high_x = max([a[0] for a in transitions]) + 5*linewidth
    x = numpy.linspace(low_x, high_x, 1000)
    y = [sum([transitions[j][1] / (1 + pow((x[i] - transitions[j][0]) / linewidth, 2))
              for j in range(len(transitions))]) for i in range(len(x))]

    for f in v_lines:
        plt.axvline(x=f)
    plt.scatter(x_other, y_other, s=0.5, color='k')
    plt.scatter(x, y, s=2, color='r')
    plt.xlabel('Detuning (MHz)')
    plt.ylabel('Absorption Signal (arb.)')
    plt.show()


def get_nearest_differences(measured, calculated):
    """
    Returns the a list of values corresponding to the absolute difference between each value in measured and the nearest
    value in calculated to that value.

    :param measured:    List of measured values to get differences for
    :param calculated:  List of reference values to compare against
    :return:
    """
    meas = measured
    calc = calculated
    meas.sort()
    calc.sort()
    return [abs(meas[i] - calc[i]) for i in range(len(meas))]


def adjust_offset_to_minimize_error(freq_measured, freq_calculated):
    """
    Adjusts the calculated frequencies with an offset in order to find a shifted frequency set that minimizes the
    mean-squared error between each measured peak and its corresponding nearest calculated peak.

    :param freq_measured:       List of measured peaks
    :param freq_calculated:     List of calculated peaks (assumed to be sorted)
    """
    # In order to get a decent first approxmiation, shift the calculated list so both have the same center
    center_measured = (min(freq_measured) + max(freq_measured)) / 2
    center_calculated = (min(freq_calculated) + max(freq_calculated)) / 2
    freq_calculated = [x - (center_calculated - center_measured) for x in freq_calculated]

    applied_offset = 0
    offset_step = 10
    num_steps = 200
    while offset_step > 0.001:
        best_offset = 0
        best_error_metric = 1e100
        for i in range(- num_steps / 2, num_steps / 2):
            curr_offset = applied_offset + i * offset_step
            freq_shifted = [x + curr_offset for x in freq_calculated]
            errors = get_nearest_differences(freq_measured, freq_shifted)
            mean_sq_error = sum(x*x for x in errors) / len(freq_measured)
            # print "offset {}, error {}".format(curr_offset, mean_sq_error)
            if mean_sq_error < best_error_metric:
                best_error_metric = mean_sq_error
                best_offset = curr_offset
        offset_step /= 3.0
        applied_offset = best_offset
    return [x + applied_offset for x in freq_calculated]


def get_fitting_error(measured_freq, Ae, Be):
    """
    Returns mean-squared errors for an optimally-shifted set of calculated frequencies from the given hyperfine
    coefficients.

    :param measured_freq:   List of measured frequencies, in MHz
    :param Ae:              Excited state A hyperfine coefficient in MHz (magnetic dipole)
    :param Be:              Excited state B hyperfine coefficient in MHz (electric quadrupole)
    :return:                Total mean-squared fitting error
    """
    transition_locations = get_peak_locations(Ae, Be)
    transition_locations.sort()
    adjusted = adjust_offset_to_minimize_error(measured_freq, transition_locations)
    diff = get_nearest_differences(measured_freq, adjusted)
    return pow(sum(x*x for x in diff) / len(measured_freq), 0.5)


def generate_contour_plot(measured_freq, A_min, A_max, B_min, B_max, grid_points):
    """
    Generates a coutour plot of the optimized mean-squared error between the measured spectrum and a spectrum generated
    from hyperfine constants within the ranges defined by the parameters.

    :param measured_freq:   List of peaks
    :param A_min:           Lower bound of hyperfine constant A (MHz)
    :param A_max:           Upper bound of hyperfine constant A (MHz)
    :param B_min:           Lower bound of hyperfine constant B (MHz)
    :param B_max:           Upper bound of hyperfine constant B (MHz)
    :param grid_points:     Number of points to divide each axis into - note that total time to generate the contour
        plot scales as this argument squared
    :return:                None
    """
    hyperfine_A = numpy.linspace(A_min, A_max, grid_points)
    hyperfine_B = numpy.linspace(B_min, B_max, grid_points)
    x, y = numpy.meshgrid(hyperfine_A, hyperfine_B)
    mean_sq_error = numpy.zeros((len(hyperfine_A), len(hyperfine_B)))
    for i in tqdm(range(grid_points), ascii=True):
        for j in range(grid_points):
            mean_sq_error[i][j] = get_fitting_error(measured_freq, x[i][j], y[i][j])
    fig, ax = plt.subplots()
    CS = ax.contour(x, y, mean_sq_error)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title("sqrt(Mean Squared Error)")
    plt.xlabel("A")
    plt.ylabel("B")
    plt.show()


# Read data from the analog in
data_time = 10
sampleRateInHz = 500000
numChannels = 2
numSamples = data_time * sampleRateInHz
# sample_data(numSamples, sampleRateInHz, numChannels, "fast_sweeps2.npy")
# freq, signal = combine_raw_data("fast_sweeps2.npy")
# process_data(freq, signal, "combined_sweeps.npy")

x_plot, plot_data = read_spectroscopy_data("combined_sweeps.npy", 8.7e5, 17.4e5)
x_plot = numpy.linspace(500, -500, len(plot_data))
plot_data = [10.5*val+0.25 for val in plot_data]

data_freq = [-221.3, -216.0, -128.5, 4.1, 92.9, 131.5, 217.7, 217.7]
data_uncertainty = [3.5, 4.9, 5.9, 2.8, 3.4, 3.6, 4.4, 4.4]
measured_freq = [-2*f for f in data_freq]
measured_uncertainty = [2*f for f in data_uncertainty]
# generate_contour_plot(measured_freq, 728.2, 729.3, 870, 930, 100)
A_range = (715.8-0.7, 716.0+0.7)
B_range = (997-64, 1017+64)
# generate_contour_plot(measured_freq, A_range[0], A_range[1], B_range[0], B_range[1], 20)

A_fitted = (A_range[0] + A_range[1]) / 2
B_fitted = (B_range[0] + B_range[1]) / 2
print "Plotting with A={} and B={}".format(A_fitted, B_fitted)
calc_transitions = get_peak_locations(A_fitted, B_fitted)
calc_transitions.sort()
adjusted = adjust_offset_to_minimize_error(measured_freq, calc_transitions)
measured_freq.sort()
'''
print measured_freq
print adjusted
chi2 = sum(pow((measured_freq[i] - adjusted[i]) / measured_uncertainty[i], 2) for i in range(len(measured_freq)))
print chi2
'''
shift = calc_transitions[0] - adjusted[0]
fitted_transitions = get_peak_locations_and_heights(A_fitted, B_fitted, shift, 1)
print fitted_transitions
saturated_absorption_signal(fitted_transitions, 17, measured_freq, x_other=x_plot, y_other=plot_data)

