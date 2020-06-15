#!/usr/bin/env python
from __future__ import print_function
import pyopenms

def setCompressionOptions(opt, verbose=False, USE_ZLIB=True, USE_NUMPRESS=True):
    """
    Adds suitable compression options for an object of type
    pyopenms.PeakFileOptions
        - compresses mass / time arrays with numpress linear
        - compresses intensity with slof (log integer)
        - compresses ion mobility with slof (log integer)

    Disable Numpress (common) or zlib (rare) for mzML consumers that don't support it.
    """

    if verbose:
        print(f'MZML Compression: ZLIB {USE_ZLIB}; Numpress {USE_NUMPRESS}')

    if USE_NUMPRESS:
        cfg = pyopenms.NumpressConfig()
        cfg.estimate_fixed_point = True
        cfg.numpressErrorTolerance = -1.0 # skip check, faster
        cfg.setCompression(b"linear");
        cfg.linear_fp_mass_acc = -1; # set the desired RT accuracy in seconds
        opt.setNumpressConfigurationMassTime(cfg)
        cfg = pyopenms.NumpressConfig()
        cfg.estimate_fixed_point = True
        cfg.numpressErrorTolerance = -1.0 # skip check, faster
        cfg.setCompression(b"slof");
        opt.setNumpressConfigurationIntensity(cfg)
    if USE_ZLIB:
        opt.setCompression(True) # zlib compression

    # Now also try to compress float data arrays (this is not enabled in all
    # versions of pyOpenMS).
    if USE_NUMPRESS:
        try:
            cfg = pyopenms.NumpressConfig()
            cfg.estimate_fixed_point = True
            cfg.numpressErrorTolerance = -1.0 # skip check, faster
            cfg.setCompression(b"slof");
            opt.setNumpressConfigurationFloatDataArray(cfg)
        except Exception:
            pass

