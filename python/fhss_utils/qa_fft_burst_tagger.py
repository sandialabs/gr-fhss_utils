#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018-2021 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
# retains certain rights in this software.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest
from gnuradio import blocks
from gnuradio import pdu_utils
import pmt
import time
try:
    from gnuradio import fhss_utils
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    from gnuradio import fhss_utils


class qa_fft_burst_tagger (gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_instantiate(self):

        # data
        src_data = (1 + 1j, 2 + 2j, 3 + 3j)

        # blocks
        src = blocks.vector_source_c(src_data)
        dst = blocks.vector_sink_c()
        dut = fhss_utils.fft_burst_tagger(910.6e6, 1024, 1000000, 0, 0, 500000)

       # float center_freq, int fft_size, int sample_rate, int burst_pre_len, int burst_post_len, int burst_width,
       # int max_bursts, int max_burst_len, float threshold, int history_size, int lookahead, bool debug)

        self.tb.connect(src, dut)
        self.tb.connect(dut, dst)

        self.tb.start()


if __name__ == '__main__':
    gr_unittest.run(qa_fft_burst_tagger, "qa_fft_burst_tagger.xml")
