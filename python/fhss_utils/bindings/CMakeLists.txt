# Copyright 2021 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Check if there is C++ code at all
########################################################################
if(NOT fhss_utils_sources)
    MESSAGE(STATUS "No C++ sources... skipping python bindings")
    return()
endif(NOT fhss_utils_sources)

########################################################################
# Check for pygccxml
########################################################################
GR_PYTHON_CHECK_MODULE_RAW(
    "pygccxml"
    "import pygccxml"
    PYGCCXML_FOUND
    )

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND fhss_utils_python_files
    cf_estimate_python.cc
    fft_burst_tagger_python.cc
    tagged_burst_to_pdu_python.cc
    python_bindings.cc
    constants_python.cc)

GR_PYBIND_MAKE_OOT(fhss_utils
   ../../..
   gr::fhss_utils
   "${fhss_utils_python_files}")

# copy bindings extension for use in QA test module
add_custom_command(TARGET fhss_utils_python POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:fhss_utils_python>
    ${PROJECT_BINARY_DIR}/test_modules/gnuradio/fhss_utils/
)

install(TARGETS fhss_utils_python DESTINATION ${GR_PYTHON_DIR}/gnuradio/fhss_utils COMPONENT pythonapi)
