# luna #

The luna library is designed to support large-scale objective studies
of sleep.  It supports EDF(+) data and XML annotations, and offers
efficient handling, filtering and manipulation of data, a range of
spectral analysis routines and spindle detection methods, and
visualization of raw data.

# dependencies #

Libraries: libfftw3, libhpdf, libz, libpbg, libsamplerate

```
sudo apt install libfftw3-dev libsamplerate0-dev libhpdf-dev libpng-dev zlib1g-dev
```



# installation #

If dependecies are not in the usual system locations, edit the `Makefile.inc` to point to the locations of headers and libraries: e.g.

```
# if dependencies are installed locally, denote here
CXXFLAGS += -I/Users/joe/src/luna-depends/depends/include
LDFLAGS += -L/Users/joe/src/luna-depends/depends/lib
```

