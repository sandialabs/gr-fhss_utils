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
import time

class qa_tagged_burst_to_pdu (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_instantiate (self):
        # data
        src_data = (1+1j,2+2j,3+3j)

        # blocks
        src = blocks.vector_source_c(src_data)
        self.debug = blocks.message_debug()
        dut = fhss_utils.tagged_burst_to_pdu(1, [1], 0.001, 0.5, 915e6, 30.72e6, 30.72e6, 30.72e6, 3 )
        
 
        self.tb.connect(src, dut)
        self.tb.msg_connect((dut, 'cpdus'), (self.debug, 'store'))
        
        self.tb.start()
        time.sleep(.1)
        self.tb.stop()
        self.tb.wait()


if __name__ == '__main__':
    gr_unittest.run(qa_tagged_burst_to_pdu, "qa_tagged_burst_to_pdu.xml")
