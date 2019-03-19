# Fast PCI Analog Input

The purpose of this script is to quickly read data from multiple analog input channels on a PCI card and then appropriately filter the raw data in order to reject segments of data which were invalidated by certain experimental instabilities (eg. lasers unlocking).  The code that reads in the data should be adaptable to any experiment in a relatively straightforward way, while the data filtering is more specialized to the holmium experiment and not necessarily suited for general-purpose usage on other experiments.

Some IDEs may state that the pydaqmx functions are not callable, but the script should run regardless as those functions are dynamically generated at runtime by the pydaqmx package.