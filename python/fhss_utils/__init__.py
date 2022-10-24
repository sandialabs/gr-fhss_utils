#
# Copyright 2008,2009 Free Software Foundation, Inc.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

# The presence of this file turns this directory into a Python package

'''
This is the GNU Radio FHSS_UTILS module. Place your Python package
description here (python/__init__.py).
'''
import os

# import pybind11 generated symbols into the fhss_utils namespace
try:
    # this might fail if the module is python-only
    from .fhss_utils_python import *
except ModuleNotFoundError:
    pass

# import any pure python here
from .fft_peak import fft_peak
from .s_and_h_detector import s_and_h_detector
from .coarse_dehopper import coarse_dehopper
from .fine_dehopper import fine_dehopper
from .fsk_burst_extractor_hier import fsk_burst_extractor_hier
from .sigmf_meta_writer import sigmf_meta_writer
#

try:
    import PIL
    from .burst_tag_debug import burst_tag_debug
except ImportError as e:
    print("Python PIL library not found, not including burst_tag_debug in import", e)
