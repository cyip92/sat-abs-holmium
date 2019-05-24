import pydaqmx
import numpy
import ctypes
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from operator import itemgetter
from tqdm import tqdm
import Queue
from sympy.physics.wigner import wigner_6j
import pyvisa as visa
import time

# Atomic structure parameters (exact spin values, assumed ground state hyperfine)
I = 7 / 2.
Jg = 15 / 2.
Je = 17 / 2.


def sample_data(num_samples, sample_rate_in_hz, num_channels, save_to_file_name):
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
    :param save_to_file_name:   Name of the file to save data into (assumed to be local directory)
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
        f = open(save_to_file_name, mode='wb+')
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
    on section "To Optimize Throughput" (page 3-73) for details on what each command does.  Note that many of these
    commands don't persist after power-cycling the frequency counter.

    :param counter: A reference to the GPIB connection to the frequency counter as given by pyvisa
    :return:        None
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


def combine_raw_data(read_from_file_name):
    """
    Parse the data from the given file, reading from the array assuming the same data format convention noted in the
    docstring for sample_data().  It reads the analog data and associates it with the given timebase in the data file
    header, then assumes the frequency counter shares the same time base.  Then, assuming that the analog data is a
    function y(t) and the counter data is x(t), combines them together in order to produce data of the form y(x)

    :param read_from_file_name:     Name of the file to read data from, assumed to be in the local directory.
    :return:                        Combined and filtered data in the format x,y
    """
    # Read all the data out, parsing the header data appropriately
    data_file = open(read_from_file_name, mode='rb+')
    all_data = numpy.load(data_file)
    header_length = int(all_data[0])
    samples_per_channel = int(all_data[1])
    num_channels = int(all_data[3])
    num_analog_samples = int(samples_per_channel * num_channels)
    time_span = 1.0 * all_data[1] / all_data[2]

    # Voltage ramp data (index 0)
    channel_index = 0
    ramp_data = all_data[header_length:(num_analog_samples + header_length)]
    x_ramp = numpy.linspace(0, time_span, samples_per_channel)
    y_ramp = [ramp_data[num_channels * i + channel_index] for i in range(samples_per_channel)]
    print "Reading {} analog points from voltage ramp".format(len(x_ramp))

    # Analog input data (index 1 is spectroscopy signal)
    channel_index = 1
    analog_data = all_data[header_length:(num_analog_samples + header_length)]
    x_analog = numpy.linspace(0, time_span, samples_per_channel)
    y_analog = [analog_data[num_channels * i + channel_index] for i in range(samples_per_channel)]
    print "Reading {} analog points from spectroscopy channel".format(len(x_analog))

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
    y_counter = [raw_counter_data[i] for i in range(len(raw_counter_data) - 2) if valid_index[i]]
    print "{} counter points, filtered to {}".format(len(raw_counter_data), len(y_counter))

    # Print average, standard deviation, and range of frequency readings
    # This is meant to be a sanity check to make sure filtering didn't go horribly wrong
    avg_freq = sum(raw_counter_data) / len(raw_counter_data)
    variances = [pow(raw_counter_data[i] - avg_freq, 2) for i in range(len(raw_counter_data))]
    stdev = pow(sum(variances) / len(raw_counter_data), 0.5)
    print "Avg. frequency {}, stdev {}, range ({}, {})".format(sum(raw_counter_data) / len(raw_counter_data), stdev,
                                                               min(raw_counter_data), max(raw_counter_data))

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

    # Combine ramp and spectroscopy data
    x_combined = []
    y_combined = []
    print "Combining ramp and spectroscopy data..."
    for i in tqdm(range(len(x_ramp)), ascii=True):
        x_combined.extend([y_ramp[i]])
        y_combined.extend([y_analog[i]])

    # Plot raw data
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Analog In (V)', color='r')
    ax1.plot(x_analog, y_analog, color='r')
    # ax2 = ax1.twinx()
    # ax2.set_ylabel('Frequency Counter (MHz)', color='b')
    # ax2.plot(x_counter, y_counter, color='b')
    plt.show()

    plt.scatter(x_combined, y_combined, s=0.5)
    plt.xlabel('Voltage Ramp (V)')
    plt.ylabel('Preamp Output (V)')
    plt.title('Combined analog input data')
    plt.show()

    return x_combined, y_combined


def read_spectroscopy_data(file_name, start_index, end_index):
    """
    Temporary function for reading out raw spectroscopy data from a file, assumes all x values followed by all y values.
    Arguments allow for selecting a subsequence of all the data for very quick trimming.  Since the data sets generally
    have a few hundred thousand very closely-spaced points which are mostly redundant, there is also built-in filtering
    to only take every 10th point.

    In practice this doesn't visually change when the data is immediately plotted afterwards, but it vastly speeds up
    the plotting time.

    :param file_name:   Name of the file to read data from, assumed to be in the local directory.
    :param start_index: Starting index for quick data trimming
    :param end_index:   Ending index for quick data trimming
    """
    # Read all the data out
    data_file = open(file_name, mode='rb+')
    all_data = numpy.load(data_file)
    x_analog = all_data[:len(all_data) / 2]
    y_analog = all_data[len(all_data) / 2:]
    x_analog = x_analog[start_index: end_index]
    y_analog = y_analog[start_index: end_index]

    # Only pick every 10th point to speed up plotting
    x_analog = [x_analog[i] for i in range(len(x_analog)) if i % 10 == 0]
    y_analog = [y_analog[i] for i in range(len(y_analog)) if i % 10 == 0]

    return x_analog, y_analog


def process_data(x_raw, y_raw, save_to_file_name):
    """
    Filter out invalid segments of data based on certain patterns and average together what is left.  The particular
    filtering processes are due to some unstable experimental artifacts specific to the holmium setup and may not be
    needed for general usage on other setups.

    :param x_raw:               List of x values for the data to filter
    :param y_raw:               List of y values for the data to filter, assumed to have the same dimensions of x_raw
    :param save_to_file_name:   Starting index for quick data trimming
    :return:                    A list of (x,y) tuples containing the averaged data
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
    The electronics on the actual saturated absorption setup cause a baseline curve which is a combination of slightly
    mismatched beam powers and doppler broadening, so a baseline curve is calculated by performing a moving average with
    a very large moving window (in this case 5% of the entire span) in order to effectively remove the baseline curve
    so that all which is left is the spectroscopy signal.
    
    In practice this baseline curve varies between different data sets but stays relatively constant on any particular
    run of taking data, so it needs to be calculated on a per-file basis.
    """
    xb = []
    yb = []
    bg_avg_size = len(x2) / 20
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
    
    This averaging is meant to accomplish the task of averaging together multiple sweeps, but it needs to be done with a
    moving average approach because the x values won't necessarily be identical between sweeps because they come from an
    analog input channel and not a consistent timebase.
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
    plt.title('Data after filtering and averaging')
    plt.show()

    # Append both data sets together and save to a file.  Further processing should assume x and y of equal lengths due
    # to the data formatting, which is a list of (x,y) points
    data_file = open(save_to_file_name, mode='wb+')
    numpy.save(data_file, numpy.append(x4, y4))
    print 'Data saved to file "{}"'.format(save_to_file_name)

    return [(x3[i], y3[i]) for i in range(len(x3))]


def get_peak_locations_and_heights(Ag, Bg, Ae, Be, x_offset, y_scale):
    """
    Returns a list of tuples (peak_freq, peak_height) corresponding to locations and relative heights of peaks in the
    saturated absorption spectrum, given all of the hyperfine constants as input parameters.

    :param Ag:          Ground state A hyperfine coefficient in MHz (magnetic dipole)
    :param Bg:          Ground state B hyperfine coefficient in MHz (electric quadrupole)
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


def get_peak_locations(Ag, Bg, Ae, Be):
    """
    Returns a list of floats corresponding to locations of peaks in the saturated absorption spectrum, given all of the
    hyperfine constants as input parameters.  Note that this doesn't return peak heights, and that the transitions are
    sorted in order of F and thus may not be strictly increasing or decreasing.

    :param Ag:          Ground state A hyperfine coefficient in MHz (magnetic dipole)
    :param Bg:          Ground state B hyperfine coefficient in MHz (electric quadrupole)
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
    :param v_lines:         List of vertical lines to overlap on top of the plot at specified x values
    :param x_other:         X values for a scatterplot to overlay over the calculated curve
    :param y_other:         Y values for a scatterplot to overlay over the calculated curve
    :return:                None; produces a plot instead
    """
    low_x = min([a[0] for a in transitions]) - 5*linewidth
    high_x = max([a[0] for a in transitions]) + 5*linewidth + 50
    x = numpy.linspace(low_x, high_x, 15000)
    y = [sum([transitions[j][1] / (1 + pow((x[i] - transitions[j][0]) / linewidth, 2))
              for j in range(len(transitions))]) for i in range(len(x))]

    fig, ax = plt.subplots()
    ax.grid()
    for val in v_lines:
        plt.axvline(x=val)

    # Automatically draw lines at the locations of peaks in "transitions" and label them appropriately
    for Fg in range(4, 12):
        trans = transitions[Fg - 4]
        plt.axvline(x=trans[0])
        plt.text(trans[0], 3.75 - 0.15 * Fg, "{} to {}'".format(Fg, Fg + 1), horizontalalignment='center',
                 verticalalignment='center', fontsize=16)

    black = mpatches.Patch(color='black', label='Experimental data')
    red = mpatches.Patch(color='red', label='Calculated fit')
    plt.legend(handles=[black, red], loc="upper left")
    fig.patch.set_facecolor('white')
    plt.scatter(x_other, y_other, s=0.5, color='k')
    plt.scatter(x, y, s=1, color='r')
    plt.xlabel('Beat Frequency (MHz)')
    plt.ylabel('Absorption Signal (arb. units)')
    plt.xlim(200, 1200)
    plt.ylim((-0.3, 3.3))
    plt.show()


def get_nearest_differences(measured, calculated):
    """
    Returns the a list of values corresponding to the absolute difference between each value in measured and the nearest
    value in calculated to that value.  DO NOT CHANGE THE ERROR METRIC IN THIS FUNCTION, CHANGE IT IN
    adjust_offset_to_minimize_error() INSTEAD IF NEEDED.

    :param measured:    List of measured values to get differences for
    :param calculated:  List of reference values to compare against
    :return:
    """
    meas = measured
    calc = calculated
    meas.sort()
    calc.sort()
    return [abs(meas[i] - calc[i]) for i in range(len(meas))]


def calculate_error_metric(value_list):
    """
    Given a list of values corresponding to differences between peak heights, calculate an appropriate error metric for
    the gradient descent to minimize.  This function primarily exists because the error metric is used in more than one
    place, so changing it here affects all relevant code.  In principle this should always just be the squared error.

    :param value_list:  List of values to calculate error metric for
    :return:            A float representing the error metric value
    """
    return sum(x*x for x in value_list)


def adjust_offset_to_minimize_error(freq_measured, freq_calculated):
    """
    Adjusts the calculated frequencies with an offset in order to find a shifted frequency set that minimizes the
    error metric between each measured peak and its corresponding nearest calculated peak.

    This function finds the best-fit value within a certain range of a specified center, then picks a new center based
    on the best within this range.  Then it takes this center and produces a smaller range around it before repeating
    this process, until the gap between adjacent steps in the range drops below the 1e-6 threshold.

    :param freq_measured:       List of measured peaks
    :param freq_calculated:     List of calculated peaks (assumed to be sorted)
    """
    # In order to get a decent first approxmiation, shift the calculated list so both have the same center
    center_measured = (min(freq_measured) + max(freq_measured)) / 2
    center_calculated = (min(freq_calculated) + max(freq_calculated)) / 2
    freq_calculated = [x - (center_calculated - center_measured) for x in freq_calculated]

    # Iteratively improve the offset via repeatedly picking the optimal fit within a range and shrinking the range
    applied_offset = 0
    offset_step = 10
    num_steps = 200
    while offset_step > 1e-6:
        best_offset = 0
        best_error_metric = 1e100
        for i in range(- num_steps / 2, num_steps / 2):
            curr_offset = applied_offset + i * offset_step
            freq_shifted = [x + curr_offset for x in freq_calculated]
            errors = get_nearest_differences(freq_measured, freq_shifted)
            mean_sq_error = calculate_error_metric(errors) / len(freq_measured)
            if mean_sq_error < best_error_metric:
                best_error_metric = mean_sq_error
                best_offset = curr_offset
        offset_step /= 3.0
        applied_offset = best_offset
    return [x + applied_offset for x in freq_calculated]


def get_fitting_error(measured_freq, Ag, Bg, Ae, Be):
    """
    Returns the value for the error metric for an optimally-shifted set of calculated frequencies from the given
    hyperfine coefficients.

    :param measured_freq:   List of measured frequencies, in MHz
    :param Ag:              Ground state A hyperfine coefficient in MHz (magnetic dipole)
    :param Bg:              Ground state B hyperfine coefficient in MHz (electric quadrupole)
    :param Ae:              Excited state A hyperfine coefficient in MHz (magnetic dipole)
    :param Be:              Excited state B hyperfine coefficient in MHz (electric quadrupole)
    :return:                Total value of the error metric
    """
    transition_locations = get_peak_locations(Ag, Bg, Ae, Be)
    transition_locations.sort()
    shifted = adjust_offset_to_minimize_error(measured_freq, transition_locations)
    diff = get_nearest_differences(measured_freq, shifted)
    return calculate_error_metric(diff)


def generate_contour_plot(measured_freq, A_min, A_max, B_min, B_max, grid_points):
    """
    Generates a contour plot of the minimized error between the measured spectrum and a spectrum generated from
    hyperfine constants.  The range of the contour plot is defined by the parameters.

    The ground state values are assumed to be in A_ground, B_ground from the outer scope.  These generally should not be
    changed.

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
    print "Generating Contour plot..."
    for i in tqdm(range(grid_points), ascii=True):
        for j in range(grid_points):
            mean_sq_error[i][j] = get_fitting_error(measured_freq, A_ground, B_ground, x[i][j], y[i][j])
    fig, ax = plt.subplots()
    CS = ax.contour(x, y, mean_sq_error)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title("Error Metric")
    plt.xlabel("A")
    plt.ylabel("B")
    plt.show()


def find_hyperfine_coefficients(measured_freq, init_Ag, init_Bg, init_Ae, init_Be):
    """
    Uses gradient descent to iteratively find values for the hyperfine constants which minimize the given metric for
    quantifying goodness-of-fit.  It uses the function is used in calculate_error_metric() and stops once the absolute
    error between steps drops below a set threshold.  It attempts to optimize all four hyperfine coefficients, but any
    of them can be dummied out as needed by multiplying its corresponding grad_XX by zero.

    In order to discourage large changes in the ground hyperfine constant, they are suppressed by 1/100.  There is also
    an asymmetry factor to account for any numerical problems potentially stemming from the significantly different
    sensitivity between A and B.

    This might not be fully working quite yet, but as of now it seems to work "well enough" although step_grad seems to
    need adjustment if the error metric is changed.

    :param measured_freq:   Measured frequencies to fit
    :param init_Ae:          Initial value for hyperfine A
    :param init_Be:          Initial value for hyperfine B
    :return:                Values A, B which minimize the fitting error
    """
    step_grad = 1
    curr_Ag = init_Ag
    curr_Bg = init_Bg
    curr_Ae = init_Ae
    curr_Be = init_Be
    prev_error = 1e30
    curr_error = get_fitting_error(measured_freq, curr_Ag, curr_Bg, curr_Ae, curr_Be)
    asym_factor = 100
    while abs(curr_error - prev_error) > 1e-6:
        grad_Ag = get_fitting_error(measured_freq, curr_Ag + step_grad, curr_Bg, curr_Ae, curr_Be) - curr_error
        grad_Bg = get_fitting_error(measured_freq, curr_Ag, curr_Bg + step_grad, curr_Ae, curr_Be) - curr_error
        grad_Ae = get_fitting_error(measured_freq, curr_Ag, curr_Bg, curr_Ae + step_grad, curr_Be) - curr_error
        grad_Be = get_fitting_error(measured_freq, curr_Ag, curr_Bg, curr_Ae, curr_Be + step_grad) - curr_error
        grad_Ag *= 0.01
        grad_Bg *= 0.01
        grad_Bg *= asym_factor
        grad_Be *= asym_factor
        grad_tot = pow(grad_Ag * grad_Ag + grad_Bg * grad_Bg + grad_Ae * grad_Ae + grad_Be * grad_Be, 0.5) + 1e-20
        curr_Ag -= grad_Ag / grad_tot * step_grad
        curr_Bg -= grad_Bg / grad_tot * step_grad
        curr_Ae -= grad_Ae / grad_tot * step_grad
        curr_Be -= grad_Be / grad_tot * step_grad
        step_grad = min(grad_tot / 45, step_grad)
        prev_error = curr_error
        curr_error = get_fitting_error(measured_freq, curr_Ag, curr_Bg, curr_Ae, curr_Be)
        print "Error {}, step {} ({}, {}, {}, {})".format(curr_error, step_grad,
                                                          curr_Ag, curr_Bg, curr_Ae, curr_Be)

    return curr_Ag, curr_Bg, curr_Ae, curr_Be


# Read data from the analog in
data_time = 20
sampleRateInHz = 500000
numChannels = 2
numSamples = data_time * sampleRateInHz
curr_file = "new_data.npy"
# Data acquisition and filtering/combining temporarily commented out
"""
sample_data(numSamples, sampleRateInHz, numChannels, curr_file)
freq, signal = combine_raw_data(curr_file)
process_data(freq, signal, "new_combined.npy")
"""

# Quick-and-dirty trimming and rescaling of raw data to visually fit the calculated curve
x_plot, plot_data = read_spectroscopy_data("new_combined.npy", 24.5e5, 38.9e5)
plot_data = [14 * val + 0.2 for val in plot_data]
x_plot = numpy.linspace(545+680, -485+680, len(plot_data))

# First point is duplicated since it's assumed to be two degenerate peaks.  This is roughly consistent with its height
data_freq = [142.8, 142.8, 211.6, 244.1, 325.8, 461.1, 552.1, 563.1]
# The next line is meant to force the fitting to match the taken data better
# data_freq = [142.8, 142.8, 220. , 254. , 332. , 458. , 552.1, 563.1]
data_uncertainty = [5.4, 5.4, 4.9, 3.2, 4.8, 4.1, 4.9, 4.4]
# Beat frequency is in the IR, actual spectroscopy is in the blue after SHG
measured_freq = [2*f for f in data_freq]
measured_uncertainty = [2*f for f in data_uncertainty]
A_ground = 800.583645
B_ground = -1668.00527


# A_ground, B_ground, A_fitted, B_fitted = find_hyperfine_coefficients(measured_freq, 800.583645, -1668.00527, 718, 950)
A_ground, B_ground, A_fitted, B_fitted = 800.583645, -1668.00527, 715.992, 925.975
A_div = 1.5
B_div = 120
# generate_contour_plot(measured_freq, A_fitted - A_div, A_fitted + A_div, B_fitted - B_div, B_fitted + B_div, 20)

print "Plotting with Ag={} and Bg={}".format(A_ground, B_ground)
print "              Ae={} and Be={}".format(A_fitted, B_fitted)
calc_transitions = get_peak_locations(A_ground, B_ground, A_fitted, B_fitted)
calc_transitions.sort()
adjusted = adjust_offset_to_minimize_error(measured_freq, calc_transitions)
measured_freq.sort()

# Chi-squared value calculation
print "Measured frequencies: {}".format(measured_freq)
print "Calculated shifted frequencies: {}".format(adjusted)
chi2 = sum(pow((measured_freq[i] - adjusted[i]) / measured_uncertainty[i], 2) for i in range(len(measured_freq)))
print "chi2 value = {}".format(chi2)

# Generate plot
shift = calc_transitions[0] - adjusted[0]
fitted_transitions = get_peak_locations_and_heights(A_ground, B_ground, A_fitted, B_fitted, shift, 1)
saturated_absorption_signal(fitted_transitions, 13, [], x_other=x_plot, y_other=plot_data)
