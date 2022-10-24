#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
# retains certain rights in this software.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

try:
    from gnuradio import fhss_utils
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from gnuradio import fhss_utils

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import pmt


class qa_constants (gr_unittest.TestCase):

    def setUp(self):
        self.tb = None

    def tearDown(self):
        self.tb = None

    def test_interned_string_constants(self):
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__in(), pmt.intern("in")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__out(), pmt.intern("out")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__center_frequency(), pmt.intern("center_frequency")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.intern("relative_frequency")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__sample_rate(), pmt.intern("sample_rate")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__bandwidth(), pmt.intern("bandwidth")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__pwr_db(), pmt.intern("pwr_db")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__snr_db(), pmt.intern("snr_db")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__debug(), pmt.intern("debug")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__rx_freq(), pmt.intern("rx_freq")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__burst_id(), pmt.intern("burst_id")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__magnitude(), pmt.intern("magnitude")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__noise_density(), pmt.intern("noise_density")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__new_burst(), pmt.intern("new_burst")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__gone_burst(), pmt.intern("gone_burst")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__rx_time(), pmt.intern("rx_time")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__start_time(), pmt.intern("start_time")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__duration(), pmt.intern("duration")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__cpdus(), pmt.intern("cpdus")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__start_offset(), pmt.intern("start_offset")))
        assert(pmt.eq(fhss_utils.PMTCONSTSTR__end_offset(), pmt.intern("end_offset")))


if __name__ == '__main__':
    gr_unittest.run(qa_constants)
