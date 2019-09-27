import pydaqmx
import numpy
import ctypes
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tqdm import tqdm
from sympy.physics.wigner import wigner_6j
import math
import time
from scipy.optimize import curve_fit

# Atomic structure parameters (exact spin values, assumed ground state hyperfine)
I = 7 / 2.
Jg = 15 / 2.
Je = 17 / 2.

def lorentzian(x, x0, base, width, h):
    """
    This is a generic Lorentzian line shape with argumenr x, center x0, vertical shift base, linewidth width, height h
    """
    return h / (1 + pow((x - x0) / width, 2)) + base


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
    analog_in_location = "Dev1/ai0:2"
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

        # Combine all data and save it to a file
        f = open(save_to_file_name, mode='wb+')
        numpy.save(f, numpy.append(header_segment, analog_segment))
        print "Data written to file"

    # I'm not sure how helpful the exception messages are, as (surprisingly) I haven't encountered them yet.
    except pydaqmx.DAQError as err:
        print "DAQmx Error: %s" % err
    finally:
        if task_handle:
            pydaqmx.DAQmxStopTask(task_handle)
            pydaqmx.DAQmxClearTask(task_handle)


def read_data_from_file(read_from_file_name, show_plots):
    """
    Parse the data from the given file, reading from the array assuming the same data format convention noted in the
    docstring for sample_data().  It reads the analog data and associates it with the given timebase in the data file
    header, then assumes the frequency counter shares the same time base.  Then, assuming that the analog data is a
    function y(t) and the counter data is x(t), combines them together in order to produce data of the form y(x)

    :param read_from_file_name:     Name of the file to read data from, assumed to be in the local directory
    :param show_plots:              Whether or not to show plots after reading data
    :return:                        Combined and filtered data in the format x,y
    """
    # Read all the data parameters out, parsing the header data appropriately
    data_file = open(read_from_file_name, mode='rb+')
    all_data = numpy.load(data_file)
    header_length = int(all_data[0])
    samples_per_channel = int(all_data[1])
    num_channels = int(all_data[3])
    num_analog_samples = int(samples_per_channel * num_channels)
    time_span = 1.0 * all_data[1] / all_data[2]
    x_data = numpy.linspace(0, time_span, samples_per_channel)
    print "Reading {} analog points from {} channels".format(len(x_data), num_channels)

    # Read the data itself
    analog_data = all_data[header_length:(num_analog_samples + header_length)]
    y_data = []
    for channel_index in range(num_channels):
        y_data.append([analog_data[num_channels * i + channel_index] for i in range(samples_per_channel)])

    # Plot raw data
    if show_plots:
        for channel_index in range(num_channels):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel("Analog In, channel {} (V)".format(channel_index), color='r')
            ax1.plot(x_data, y_data[channel_index], color='r')
            plt.show()

    return x_data, y_data


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


def process_data(x_raw, ramp, y_raw, save_to_file_name):
    """
    Filter out invalid segments of data based on certain patterns and average together what is left.  The particular
    filtering processes are due to some unstable experimental artifacts specific to the holmium setup and may not be
    needed for general usage on other setups.

    :param x_raw:               List of x values for the data to filter
    :param ramp:                List of values from the voltage ramp to use for filtering
    :param y_raw:               List of y values for the data to filter, assumed to have the same dimensions of x_raw
    :param save_to_file_name:   Starting index for quick data trimming
    :return:                    A list of x values and then a list of lists of y values
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
    ascending_ramp_indices = [i for i in range(len(x_raw) - ramp_jaggedness)
                              if ramp[i] < ramp[i + ramp_jaggedness]]
    x1 = [x_raw[i] for i in ascending_ramp_indices]
    r1 = [ramp[i] for i in ascending_ramp_indices]
    print "{} points filtered to {} on only-ascending ramp edges".format(len(x_raw), len(ascending_ramp_indices))
    y1 = []
    for channel_index in range(len(y_raw)):
        y1.append([y_raw[channel_index][i] for i in ascending_ramp_indices])
        '''
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel("Analog In, channel {} (V)".format(channel_index), color='r')
        ax1.plot(x1, y1[channel_index], color='r')
        plt.show()
        '''

    """
    The laser sometimes has a tendency to unlock, which effectively invalidates that entire period of the ramp.
    This checks for very large jumps in the absorption signal and only adds ramps to an ongoing data array if no
    such jumps are detected.  (The loops are done weirdly, but I'm not entirely sure to make them much better)
    
    This also does some horizontal shifting in order to make sure successive scans line up, using the cavity signal as
    a reference
    """
    i = 0
    num_ramps = 0
    absorption_signal_index = 0
    cavity_index = 1
    x2 = []
    y2 = [[], []]
    r2 = []
    while i < len(r1) - 1:
        curr_ramp_indices = []
        absorption_signal = []
        ramp_is_valid = True
        while i < len(r1) - 1 and abs(r1[i] - r1[i+1]) < 1:
            curr_point = y1[absorption_signal_index][i]
            ramp_is_valid = ramp_is_valid and (len(absorption_signal) == 0
                                               or abs(absorption_signal[len(absorption_signal) - 1] - curr_point) < 0.5)
            curr_ramp_indices.append(i)
            absorption_signal.append(curr_point)
            i += 1
        if ramp_is_valid:
            r2.extend([r1[i] for i in curr_ramp_indices])

            # Find the peak of the ULE transmission and shift the ramp accordingly
            low = min([y1[cavity_index][i] for i in curr_ramp_indices])
            gap = max([y1[cavity_index][i] for i in curr_ramp_indices]) - low
            above_threshold = [i for i in curr_ramp_indices if y1[cavity_index][i] > low + gap * 0.3]
            center_index = 1. * sum(above_threshold) / len(above_threshold)
            f = center_index - math.floor(center_index)
            x_avg = x1[int(math.floor(center_index))] * (1 - f) + x1[int(math.ceil(center_index))] * f
            x2.extend([x1[i] - x_avg for i in curr_ramp_indices])
            for channel_index in range(len(y_raw)):
                y2[channel_index].extend([y1[channel_index][i] for i in curr_ramp_indices])
            num_ramps += 1
            i += 1
        i += 1
    print "Found {} valid sweeps, filtered {} points to {}".format(num_ramps, len(x1), len(x2))

    # Plotting code to show processed data for visual checking
    plt.scatter(x2, y2[0], s=0.5, color='m')
    plt.scatter(x2, y2[1], s=0.5, color='b')
    plt.xlabel('Ramp Voltage (V)')
    plt.ylabel('ULE PD (V)')
    plt.title('Data after filtering')
    plt.show()

    # Append both data sets together and save to a file.  Further processing should assume x and y of equal lengths due
    # to the data formatting, which is a list of (x,y) points
    if save_to_file_name is not None:
        data_file = open(save_to_file_name, mode='wb+')
        header_info = [1 + len(y2), len(x2)]
        numpy.save(data_file, numpy.append(header_info, numpy.append(x2, y2)))
        print 'Data saved to file "{}"'.format(save_to_file_name)

    return [x2, y2]


def average_spectroscopy_signal(x_val, spectroscopy_signal):
    """
    Takes the very large number of data points from the spectroscopy signal and averages them together into something
    tractable, with error bars corresponding to the measurement spread.

    :param x_val:                   List of x values
    :param spectroscopy_signal:     List of y values corresponding to the spectroscopy signal, same length as x_val
    """
    plt.scatter(x_val, spectroscopy_signal, s=0.5, color='m')
    plt.xlabel('Ramp Voltage (V)')
    plt.ylabel('ULE PD (V)')
    plt.title('Data after filtering')
    plt.show()


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
        # plt.axvline(x=trans[0])
        plt.text(trans[0], 3.75 - 0.15 * Fg, "{} to {}'".format(Fg, Fg + 1), horizontalalignment='center',
                 verticalalignment='center', fontsize=16)

    black = mpatches.Patch(color='black', label='Older sweep data')
    red = mpatches.Patch(color='red', label='Curve from fitted A,B')
    blue = mpatches.Patch(color='blue', label='Measured frequencies')
    plt.legend(handles=[black, red, blue], loc="upper left")
    fig.patch.set_facecolor('white')
    plt.scatter(x_other, y_other, s=0.5, color='k')
    plt.scatter(x, y, s=1, color='r')
    plt.xlabel('Beat Frequency (MHz)')
    plt.ylabel('Absorption Signal (arb. units)')
    plt.xlim(-500, 600)
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


def load_and_get_error_bars(file_name, analog_channel_index, bin_size, trim_fractions):
    """
    Load data which has already been processed in such a way to line up ULE cavity peaks on successive sweeps with each
    other.  This data is assumed to have the formatting from process_data() which is
        [ num_channels, samples_per_channel, ch0[0], ch0[1], ... , ch0[n], ch1[0], ch1[1], ... ]
    The data is arranged in such a way that the x values repeatedly sweep over the range, so this function sorts them in
    ascending order.  Finally, it filters the data down using bins of the specified size and calculates error bars.

    :param file_name:               Name of the file to load the data from
    :param analog_channel_index:    This function only pulls from one analog channel, this parameter specifies which one
    :param bin_size:                Number of data points per bin for output data points
    :param trim_fractions:          A pair of numbers [left, right] which specify what proportion of data should be
        excluded on the left and right sides; for example [0.1, 0.2] will exclude the leftmost 10% and rightmost 20%
    :return:                        An array of arrays, containing [x_avg, y_avg, x_err, y_err] of all equal sizes
    """

    data_file = open(file_name, mode='rb+')
    all_data = numpy.load(data_file)
    num_channels = int(all_data[0])
    channel_length = int(all_data[1])
    all_data = all_data[2:]
    print "Loading data from {} channels with length {}".format(num_channels, channel_length)
    x_filtered = all_data[0: channel_length]
    analog_filtered = all_data[(analog_channel_index + 1) * channel_length: (analog_channel_index + 2) * channel_length]

    # Sort the data points by x value
    zipped = zip(x_filtered, analog_filtered)
    zipped.sort()
    x_filtered = [x for x, y in zipped]
    analog_filtered = [y for x, y in zipped]

    print "Binning data"
    x_bin = []
    y_bin = []
    x_error = []
    y_error = []
    for index in tqdm(range(0, channel_length, bin_size), ascii=True):
        max_index = min(channel_length, index + bin_size)
        x_window = [x_filtered[i] for i in range(index, max_index)]
        y_window = [analog_filtered[i] for i in range(index, max_index)]
        x_avg = sum(x for x in x_window) / len(x_window)
        y_avg = sum(x for x in y_window) / len(y_window)
        y_stdev = pow(sum(pow(x - y_avg, 2) for x in y_window) / len(y_window), 0.5)
        x_bin.append(x_avg)
        y_bin.append(y_avg)
        x_error.append((max(x_window) - min(x_window)) / 2)
        y_error.append(y_stdev)

    start = int(trim_fractions[0] * len(x_bin))
    end = int((1 - trim_fractions[1]) * len(x_bin))
    x_bin = [x_bin[i] for i in range(start, end)]
    y_bin = [y_bin[i] for i in range(start, end)]
    x_error = [x_error[i] for i in range(start, end)]
    y_error = [y_error[i] for i in range(start, end)]
    return x_bin, y_bin, x_error, y_error


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
    step_grad = 0.01
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
        grad_Ag *= 0
        grad_Bg *= 0
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


def get_freq_scale(x_val1, ule_val1, x_val2, ule_val2):
    """
    Takes in two sets of data containing a horizontal scale (x_val), a series of transmission data (ule_val), and a
    frequency reference taken to be where the peak is "supposed" to be.  It takes these two values together

    :param x_val1:      x values for the first peak
    :param ule_val1:    transmission values for the first peak
    :param x_val2:      x values for the second peak
    :param ule_val2:    transmission values for the second peak
    :return:            A pair of fitted [x_val1, x_val2] which correspond to where the ULE peaks fall on the x scale
    """
    max1 = max(ule_val1)
    max_index1 = ule_val1.index(max1)
    p_opt1, _ = curve_fit(lorentzian, x_val1, ule_val1,
                             bounds=([x_val1[max_index1] - 0.002, 0.00, 0.000, 0.8 * max1],
                             [x_val1[max_index1] + 0.002, 0.05, 0.0005, 1.2 * max1]))
    max2 = max(ule_val2)
    max_index2 = ule_val2.index(max2)
    p_opt2, _ = curve_fit(lorentzian, x_val2, ule_val2,
                             bounds=([x_val2[max_index2] - 0.002, 0.00, 0.000, 0.8 * max2],
                             [x_val2[max_index2] + 0.002, 0.05, 0.0005, 1.2 * max2]))
    return p_opt1[0], p_opt2[0]


# Read data from the analog in
data_time = 20
sampleRateInHz = 100000
numChannels = 3
numSamples = data_time * sampleRateInHz
curr_file = "455MHz.npy"
#sample_data(numSamples, sampleRateInHz, numChannels, curr_file)

#time_series, analog = read_data_from_file(curr_file, False)
#ramp_data, spectroscopy_data, cavity_data = analog
#x_filtered, analog_filtered = process_data(time_series, ramp_data, [spectroscopy_data, cavity_data], "pr455.npy")

freq1 = 455
freq2 = 460
file1 = "pr{}.npy".format(freq1)
file2 = "pr{}.npy".format(freq2)
trim_fractions = [0.3, 0]
x_bin1, y_bin1, x_error1, y_error1 = load_and_get_error_bars(file1, 0, 5000, trim_fractions)
u1, ule1, _, _ = load_and_get_error_bars(file1, 1, 200, trim_fractions)
x_bin2, y_bin2, x_error2, y_error2 = load_and_get_error_bars(file2, 0, 5000, trim_fractions)
u2, ule2, _, _ = load_and_get_error_bars(file2, 1, 200, trim_fractions)

# Vertical shifts to match their averages
print "Shifting data to match"
padding = len(x_bin1) / 5
avg1 = sum(y_bin1[i] for i in range(padding, len(y_bin1) - padding)) / (len(y_bin1) - 2 * padding)
v_shift = avg1 - sum(y_bin2[i] for i in range(padding, len(y_bin2) - padding)) / (len(y_bin2) - 2 * padding)
y_bin2 = [y_bin2[i] + v_shift for i in range(len(y_bin2))]

# Shift all data to have zero baseline
y_bin1 = [y_bin1[i] - min(y_bin1) for i in range(len(y_bin1))]
y_bin2 = [y_bin2[i] - min(y_bin2) for i in range(len(y_bin2))]
ule1 = [ule1[i] - min(ule1) for i in range(len(ule1))]
ule2 = [ule2[i] - min(ule2) for i in range(len(ule2))]

# Find the optimal horizontal shift to match the two data sets
max1 = y_bin1.index(max(y_bin1))
max2 = y_bin2.index(max(y_bin2))
shift = 2 * (max2 - max1) * (x_bin1[1] - x_bin1[0])
x_bin2 = [x_bin2[i] - shift for i in range(len(x_bin2))]
u2 = [u2[i] - shift for i in range(len(u2))]

# Rescale the horizontal scale to indicate frequency, based on ULE peak locations
x_peak1, x_peak2 = get_freq_scale(u1, ule1, u2, ule2)
scale = (freq2 - freq1) / (x_peak2 - x_peak1)
x_bin1 = [(x_bin1[i] - x_peak1) * scale + freq1 for i in range(len(x_bin1))]
x_bin2 = [(x_bin2[i] - x_peak2) * scale + freq2 for i in range(len(x_bin2))]

p_opt, p_cov = curve_fit(lorentzian, x_bin1 + x_bin2, y_bin1 + y_bin2, sigma=y_error1 + y_error2,
                         bounds=([freq1 - 10, 0.00,  3, 0.1],
                                 [freq2 + 10, 0.05, 15, 0.2]))
print p_opt
print numpy.sqrt(numpy.diag(p_cov))

plt.errorbar(x_bin1, y_bin1, xerr=x_error1, yerr=y_error1, color='b')
plt.errorbar(x_bin2, y_bin2, xerr=x_error2, yerr=y_error2, color='r')
plt.plot(x_bin1, lorentzian(x_bin1, *p_opt), color='k')
plt.xlabel('Frequency (MHz)')
plt.ylabel("Absorption signal (arb)")
plt.show()
