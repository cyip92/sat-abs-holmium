# Fast PCI Analog Input

The purpose of this script is to quickly read data from multiple analog input channels on a PCI card and then appropriately filter the raw data in order to reject segments of data which were invalidated by certain experimental instabilities (eg. lasers unlocking).  The code that reads in the data should be adaptable to any experiment in a relatively straightforward way, while the data filtering is more specialized to the holmium experiment and not necessarily suited for general-purpose usage on other experiments.

# Some installation notes

Some IDEs may state that the pydaqmx functions are not callable, but the script should run regardless as those functions are dynamically generated at runtime by the pydaqmx package.

The pydaqmx package by default seems to install incorrectly; the fix is to navigate to the folder with all the pydaqmx files, delete task.pyc, and rename task.py to Task.py as all the internal references seem to be to Task and not task.

# Data reading process

Currently the software actually reads from both the PCI card and the RS232 connection of a HP53131A frequency counter.  Due to needing to avoid issues with the frequency counter communication timing out between reads, the data-taking process alternates reading from the PCI card and from the frequency counter in an interleaved way.  The process of reading data appears to be fast enough to not noticeably delay or pause data acquisition from the PCI card (which is significantly faster than the frequency counter).  For more general-purpose usage, the code for reading from either device can be commented out as needed since they can operate independently of each other.  The frequency counter readings come in as strings formatted like "56.231 MHz" which requires a very minor amount of parsing before being saved to the data file as well.
