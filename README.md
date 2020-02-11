# pulsewaves

Python3 module for reading uncomressed lidar full waveform PulseWaves files (PLS).

## Install

    git clone https://github.com/Rheinwalt/pulsewaves.git
    cp pulsewaves/pulsewaves.py ~/your/python/path

## Usage

    from pulsewaves import PulseWaves
    
    f = PulseWaves('sample.pls')
    f.export()

This will export all Pulse Records to a HDF5 file with the same
file name (sample.hdf).

Alternatively, one can load individual Pulse Records and Waves.
This is done as in the original Python2 version of this module:

    from pulsewaves import PulseWaves
    
    f = PulseWaves('sample.pls')

    # get the first pulse record
    r = f.get_pulse(0)

    # get the corresponding wave segments
    w = f.get_waves(0)
    print(w.segments.keys())
