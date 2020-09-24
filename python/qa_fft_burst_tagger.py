#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
# retains certain rights in this software.
# 
# SPDX-License-Identifier: GPL-3.0-or-later
# 

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import fhss_utils_swig as fhss_utils
import pdu_utils
import pmt
import time

class qa_fft_burst_tagger (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()
        

    def tearDown (self):
        self.tb = None
        
    def test_001_instantiate (self):
        
        # data
        src_data = (1+1j,2+2j,3+3j)

        # blocks
        src = blocks.vector_source_c(src_data)
        dst = blocks.vector_sink_c()
        dut = fhss_utils.fft_burst_tagger(910.6e6, 1024, 1e6, 0.00015, 0.00015, 500e3 )
        
        self.tb.connect(src, dut)
        self.tb.connect(dut, dst)
        
        self.tb.start()
    
if __name__ == '__main__':
    gr_unittest.run(qa_fft_burst_tagger, "qa_fft_burst_tagger.xml")
