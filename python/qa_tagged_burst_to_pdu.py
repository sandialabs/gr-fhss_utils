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
import time
import pmt
import math
try:
    import fhss_utils
except ImportError:
    import os
    import sys
    dirname, filename = os.path.split(os.path.abspath(__file__))
    sys.path.append(os.path.join(dirname, "bindings"))
    import fhss_utils

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

    def test_002_simple (self):
        # This test processes a single burst and confirms that the resulting metadata is as expected

        # data
        min_data_size = 32 * 1024 # arbitrary constant in tagged_burst_to_pdu_impl.h
        src_data = (1,) * min_data_size * 8
        
        new_burst_offset = 64
        new_burst_dict = pmt.make_dict()
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1234))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(1e6/30.72e6)) # in [-1.0,1.0], not Hz
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__center_frequency(), pmt.from_float(915e6)) # of the whole signal, not the burst
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__magnitude(), pmt.from_float(40))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__sample_rate(), pmt.from_float(30.72e6))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__noise_density(), pmt.from_float(-100)) # in dBFS/Hz
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__bandwidth(), pmt.from_float(0.1e6))
        nb_tag = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        duration = 1024
        gone_burst_offset = new_burst_offset + duration
        gone_burst_dict = pmt.make_dict()
        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1234))
        gb_tag = gr.tag_utils.python_to_tag([gone_burst_offset, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])

        # blocks
        src = blocks.vector_source_c(src_data, False, 1, [nb_tag, gb_tag])
        #src.set_min_output_buffer(min_data_size*2) # not necessary, block calls set_output_multiple
        debug = blocks.message_debug()

        dec = 16
        taps = [1]
        min_time = 10e-6 # 1024/30.72e6 is about 33 usec
        max_time = 1.0
        nthreads = 3
        samp_rate = 30.72e6
        rel_span = 1
        rel_samp_rate = 1
        rel_cf = 0
        dut = fhss_utils.tagged_burst_to_pdu(dec, taps, min_time, max_time, rel_cf, rel_span, rel_samp_rate, samp_rate, nthreads )
        
        self.tb.connect(src, dut)
        self.tb.msg_connect((dut, 'cpdus'), (debug, 'store'))
        
        #self.tb.run() # blocking, vector_source will end flowgraph
        self.tb.start()
        time.sleep(0.1)
        self.tb.stop()
        time.sleep(0.1)
        self.tb.wait()
        time.sleep(0.1)

        print("test simple:")
        #print(f"how many msg? {debug.num_messages()}")
        self.assertEqual(debug.num_messages(), 1)
        #print(f"received: {pmt.car(debug.get_message(0))}")
        #print(f"received: {pmt.cdr(debug.get_message(0))}")
        rcv_meta = pmt.car(debug.get_message(0))
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("duration"), pmt.PMT_NIL)), duration/samp_rate, 6)
        self.assertEqual(pmt.to_uint64(pmt.dict_ref(rcv_meta, pmt.intern("start_offset"), pmt.PMT_NIL)), new_burst_offset)
        self.assertEqual(pmt.to_uint64(pmt.dict_ref(rcv_meta, pmt.intern("end_offset"), pmt.PMT_NIL)), gone_burst_offset)
        self.assertEqual(pmt.to_uint64(pmt.dict_ref(rcv_meta, pmt.intern("burst_id"), pmt.PMT_NIL)), 1234)
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("start_time"), pmt.PMT_NIL)), new_burst_offset/samp_rate, 6)
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("sample_rate"), pmt.PMT_NIL)), samp_rate/dec, 6)
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("bandwidth"), pmt.PMT_NIL)), 0.1e6, 1)
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("noise_density"), pmt.PMT_NIL)), -100, 1)
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("magnitude"), pmt.PMT_NIL)), 40, 1)
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("center_frequency"), pmt.PMT_NIL)), 916e6, 0) # center of burst
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("relative_frequency"), pmt.PMT_NIL)), 1e6, 0) # in Hz

    def test_003_simul_bursts (self):
        # This test forces (as much as we can) three threads to process three simultaneous bursts. We check that
        # the resulting vectors are appropriately rotated to baseband.

        # data
        min_data_size = 32 * 1024 # arbitrary constant in tagged_burst_to_pdu_impl.h
        src_data = (1,) * min_data_size * 8
        
        new_burst_offset = 64
        new_burst_dict = pmt.make_dict()
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__center_frequency(), pmt.from_float(915e6)) # of the whole signal, not the burst
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__magnitude(), pmt.from_float(40))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__sample_rate(), pmt.from_float(30.72e6))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__noise_density(), pmt.from_float(-100)) # in dBFS/Hz
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__bandwidth(), pmt.from_float(0.1e6))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1001))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(0e6/30.72e6)) # in [-1.0,1.0], not Hz
        nb_tag1 = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1002))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(1e6/30.72e6)) # in [-1.0,1.0], not Hz
        nb_tag2 = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1003))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(-1e6/30.72e6)) # in [-1.0,1.0], not Hz
        nb_tag3 = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        duration = 1024
        gone_burst_offset = new_burst_offset + duration
        gone_burst_dict = pmt.make_dict()
        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1001))
        gb_tag1 = gr.tag_utils.python_to_tag([gone_burst_offset, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])
        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1002))
        gb_tag2 = gr.tag_utils.python_to_tag([gone_burst_offset+1, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])
        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(1003))
        gb_tag3 = gr.tag_utils.python_to_tag([gone_burst_offset+2, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])

        # blocks
        src = blocks.vector_source_c(src_data, False, 1, [nb_tag1, nb_tag2, nb_tag3, gb_tag1, gb_tag2, gb_tag3])
        #src.set_min_output_buffer(min_data_size*2) # not necessary, block calls set_output_multiple
        debug = blocks.message_debug()

        dec = 256
        taps = [1]
        min_time = 10e-6 # 1024/30.72e6 is about 33 usec
        max_time = 1.0
        nthreads = 3
        samp_rate = 30.72e6
        rel_span = 1
        rel_samp_rate = 1
        rel_cf = 0
        dut = fhss_utils.tagged_burst_to_pdu(dec, taps, min_time, max_time, rel_cf, rel_span, rel_samp_rate, samp_rate, nthreads )
        
        self.tb.connect(src, dut)
        self.tb.msg_connect((dut, 'cpdus'), (debug, 'store'))
        
        #self.tb.run() # blocking, vector_source will end flowgraph
        self.tb.start()
        time.sleep(0.5)
        self.tb.stop()
        time.sleep(0.1)
        self.tb.wait()
        time.sleep(0.1)

        print("test simultaneous:")
        print(f"how many msg? {debug.num_messages()}")
        self.assertEqual(debug.num_messages(), 3)
        #print(f"received: {pmt.car(debug.get_message(0))}")
        #print(f"received: {pmt.car(debug.get_message(1))}")
        #print(f"received: {pmt.car(debug.get_message(2))}")
        #print(f"received: {pmt.cdr(debug.get_message(0))}")
        #print(f"received: {pmt.cdr(debug.get_message(1))}")
        #print(f"received: {pmt.cdr(debug.get_message(2))}")

        r1 = pmt.c32vector_elements(pmt.cdr(debug.get_message(0)))
        r2 = pmt.c32vector_elements(pmt.cdr(debug.get_message(1)))
        r3 = pmt.c32vector_elements(pmt.cdr(debug.get_message(2)))

        #### Expected Results
        # what do we expect these vectors to be if they are filtered and rotated correctly?
        # - taps are [1], so no filtering should be noticeable in the data
        # - rotation is by 0, -1MHz, 1MHz for PDUs 0,1,2 respectively
        # - decimation of 256 means only every 256 samples is kept
        # First and last elements are trivial since the decimation (256) evenly divides the number of samples (1024).

        # example: get_messages(2) is shifted +1MHz to reach baseband
        #  >>> cmath.exp(0 + 1j * 1e6/30.72e6 * 2*math.pi* 256) # second element
        #  (-0.5000000000000023+0.8660254037844373j)
        #  >>> cmath.exp(0 + 1j * 1e6/30.72e6 * 2*math.pi* 512) # third element
        #  (-0.49999999999999534-0.8660254037844414j)
        #  >>> cmath.exp(0 + 1j * 1e6/30.72e6 * 2*math.pi* 768) # fourth element
        #  (1+9.82193361864236e-16j)

        # similar math can be done for get_messages(1), resulting in similar results with different signs
        
        cp6 = math.cos(math.pi/6) # expected samples are 30 degrees off of the real-axis
        self.assertComplexTuplesAlmostEqual(r1, ((1+0j),(1+0j),(1+0j),(1+0j)), 3)
        self.assertComplexTuplesAlmostEqual(r2, ((1+0j),(-.5-cp6*1j),(-.5+cp6*1j),(1+0j)), 3)
        self.assertComplexTuplesAlmostEqual(r3, ((1+0j),(-.5+cp6*1j),(-.5-cp6*1j),(1+0j)), 3)

    def test_004_long_burst (self):
        # This test processes a single burst that exceeds the maximum burst size that will be truncated

        # data
        min_data_size = 32 * 1024 # arbitrary constant in tagged_burst_to_pdu_impl.h
        src_data = (1,) * min_data_size * 8
        
        new_burst_offset = 64
        new_burst_dict = pmt.make_dict()
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(505))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(1e6/30.72e6)) # in [-1.0,1.0], not Hz
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__center_frequency(), pmt.from_float(915e6)) # of the whole signal, not the burst
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__magnitude(), pmt.from_float(40))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__sample_rate(), pmt.from_float(30.72e6))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__noise_density(), pmt.from_float(-100)) # in dBFS/Hz
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__bandwidth(), pmt.from_float(0.1e6))
        nb_tag = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        duration = 1024
        gone_burst_offset = new_burst_offset + duration
        gone_burst_dict = pmt.make_dict()
        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(505))
        gb_tag = gr.tag_utils.python_to_tag([gone_burst_offset, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])

        # blocks
        src = blocks.vector_source_c(src_data, False, 1, [nb_tag, gb_tag])
        #src.set_min_output_buffer(min_data_size*2) # not necessary, block calls set_output_multiple
        debug = blocks.message_debug()

        dec = 16
        taps = [1]
        min_time = 1e-6
        max_time = 10e-6
        nthreads = 3
        samp_rate = 30.72e6
        rel_span = 1
        rel_samp_rate = 1
        rel_cf = 0
        dut = fhss_utils.tagged_burst_to_pdu(dec, taps, min_time, max_time, rel_cf, rel_span, rel_samp_rate, samp_rate, nthreads )
        
        self.tb.connect(src, dut)
        self.tb.msg_connect((dut, 'cpdus'), (debug, 'store'))
        
        #self.tb.run() # blocking, vector_source will end flowgraph
        self.tb.start()
        time.sleep(0.1)
        self.tb.stop()
        time.sleep(0.1)
        self.tb.wait()
        time.sleep(0.1)

        print("test long burst:")
        #print(f"how many msg? {debug.num_messages()}")
        self.assertEqual(debug.num_messages(), 1)
        #print(f"received: {pmt.car(debug.get_message(0))}")
        #print(f"received: {pmt.cdr(debug.get_message(0))}")
        #print(f"received len: {pmt.length(pmt.cdr(debug.get_message(0)))}")
        rcv_meta = pmt.car(debug.get_message(0))
        
        # we expect a duration equal to `max_time` and a vector of length that corrseponds to `max_time` samples
        self.assertAlmostEqual(pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("duration"), pmt.PMT_NIL)), max_time, 6)
        self.assertEqual(pmt.length(pmt.cdr(debug.get_message(0))), (max_time*samp_rate) // dec)
        self.assertEqual(pmt.to_uint64(pmt.dict_ref(rcv_meta, pmt.intern("burst_id"), pmt.PMT_NIL)), 505)

    def test_005_short_burst (self):
        # This test processes two bursts of which the first should be dropped for being too short

        # data
        min_data_size = 32 * 1024 # arbitrary constant in tagged_burst_to_pdu_impl.h
        src_data = (1,) * min_data_size * 8
        
        new_burst_offset = 64
        new_burst_dict = pmt.make_dict()
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__center_frequency(), pmt.from_float(915e6)) # of the whole signal, not the burst
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__magnitude(), pmt.from_float(40))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__sample_rate(), pmt.from_float(30.72e6))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__noise_density(), pmt.from_float(-100)) # in dBFS/Hz
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__bandwidth(), pmt.from_float(0.1e6))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(111))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(1e6/30.72e6)) # in [-1.0,1.0], not Hz
        nb_tag_short = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(222))
        new_burst_dict = pmt.dict_add(new_burst_dict, fhss_utils.PMTCONSTSTR__relative_frequency(), pmt.from_float(2e6/30.72e6)) # in [-1.0,1.0], not Hz
        nb_tag_normal = gr.tag_utils.python_to_tag([new_burst_offset, fhss_utils.PMTCONSTSTR__new_burst(), new_burst_dict, pmt.intern("qa_test")])

        short_duration = 128
        normal_duration = 1024
        gone_burst_dict = pmt.make_dict()
        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(111))
        gone_burst_offset = new_burst_offset + short_duration
        gb_tag_short = gr.tag_utils.python_to_tag([gone_burst_offset, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])

        gone_burst_dict = pmt.dict_add(gone_burst_dict, fhss_utils.PMTCONSTSTR__burst_id(), pmt.from_uint64(222))
        gone_burst_offset = new_burst_offset + normal_duration
        gb_tag_normal = gr.tag_utils.python_to_tag([gone_burst_offset, fhss_utils.PMTCONSTSTR__gone_burst(), gone_burst_dict, pmt.intern("qa_test")])

        # blocks
        src = blocks.vector_source_c(src_data, False, 1, [nb_tag_short, nb_tag_normal, gb_tag_short, gb_tag_normal])
        #src.set_min_output_buffer(min_data_size*2) # not necessary, block calls set_output_multiple
        debug = blocks.message_debug()

        dec = 16
        taps = [1]
        min_time = 5e-6 # 153-ish samples
        max_time = 1e-3
        nthreads = 3
        samp_rate = 30.72e6
        rel_span = 1
        rel_samp_rate = 1
        rel_cf = 0
        dut = fhss_utils.tagged_burst_to_pdu(dec, taps, min_time, max_time, rel_cf, rel_span, rel_samp_rate, samp_rate, nthreads )
        
        self.tb.connect(src, dut)
        self.tb.msg_connect((dut, 'cpdus'), (debug, 'store'))
        
        #self.tb.run() # blocking, vector_source will end flowgraph
        self.tb.start()
        time.sleep(0.1)
        self.tb.stop()
        time.sleep(0.1)
        self.tb.wait()
        time.sleep(0.1)

        print("test short burst:")
        print(f"how many msg? {debug.num_messages()}")
        self.assertEqual(debug.num_messages(), 1) # first message dropped
        #print(f"received: {pmt.car(debug.get_message(0))}")
        #print(f"received: {pmt.cdr(debug.get_message(0))}")
        rcv_meta = pmt.car(debug.get_message(0))
        
        # we expect to not receive burst 111, but to receive burst 222
        self.assertEqual(pmt.to_uint64(pmt.dict_ref(rcv_meta, pmt.intern("burst_id"), pmt.PMT_NIL)), 222)

if __name__ == '__main__':
    gr_unittest.run(qa_tagged_burst_to_pdu, "qa_tagged_burst_to_pdu.xml")
