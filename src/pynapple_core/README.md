# PYthon Neural Analysis Package.

**Note from forker: I've taken this lovely package and skinned it alive, because I only need it's internal organs. And the GUI + processing dependencies (notably PyQt5 and numba) are causing problems with my M1 chip since ARM-native wheels have not yet been built for these packages.**

To use core, clone this repo into your project, then do
```
from pynapple_core import core as nap
```
and use core classes and functions to your hearts content without having to install Qt or numba. Admitedly this fork is a complete hack, but it's a relatively small and useful piece of code... 


pynapple is a light-weight python library for neurophysiological data analysis.


Dependencies 

-   pandas
-   numpy
-   scipy

