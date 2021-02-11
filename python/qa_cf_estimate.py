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
from gnuradio.filter import firdes
import fhss_utils_swig as fhss_utils
import numpy as np
import pdu_utils
import pmt
import time
from math import pi

class qa_cf_estimate (gr_unittest.TestCase):

    class simple_modulator(gr.hier_block2):
      def __init__(self, sps=8, bt=0.5, mod_idx=0.68):
        gr.hier_block2.__init__(self,
            "simple_modulator",
            gr.io_signature(0, 0, 0),  # Input signature
            gr.io_signature(0, 0, 0)) # Output signature

        # message ports
        self.message_port_register_hier_in("in")
        self.message_port_register_hier_out("out")

        # blocks
        self.pack = pdu_utils.pack_unpack(pdu_utils.MODE_UNPACK_BYTE, pdu_utils.BIT_ORDER_LSB_FIRST)
        self.preamble = pdu_utils.pdu_preamble([], [], sps, 0)
        modulation_index = mod_idx
        sensitivity = (pi * modulation_index) / sps
        gain = 1.0
        taps = firdes.gaussian(gain, sps, bt, 5)
        self.gmsk = pdu_utils.pdu_gmsk_fc(sensitivity, taps)

        # connections
        self.msg_connect(self, "in", self.pack, "pdu_in")
        self.msg_connect(self.pack, "pdu_out", self.preamble, "pdu_in")
        self.msg_connect(self.preamble, "pdu_out", self.gmsk, "pdu_in")
        self.msg_connect(self.gmsk, "pdu_out", self, "out")

    class pdu_rotate(gr.sync_block):
      def __init__(self, rotate):
        gr.sync_block.__init__(self, "pdu_rotate", in_sig=None, out_sig=None)
        self.rotate = rotate
        self.message_port_register_in(pmt.intern("in"))
        self.set_msg_handler(pmt.intern("in"), self.msg_handler)
        self.message_port_register_out(pmt.intern("out"))

      def msg_handler(self, pdu):
        meta = pmt.car(pdu)
        data = pmt.c32vector_elements(pmt.cdr(pdu))
        data = [d*r for d,r in zip(data, np.exp(1j*np.linspace(0,2*np.pi*self.rotate*len(data),len(data))))]
        self.message_port_pub(pmt.intern("out"), pmt.cons(meta,pmt.init_c32vector(len(data), data)))

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_coerce (self):
        samp_rate = 1e6
        freq_offset = -400e3
        center_freq = 911e6

        # blocks
        self.emitter = pdu_utils.message_emitter()
        self.cf = fhss_utils.cf_estimate(fhss_utils.COERCE, [x*1e6 for x in range(900,930)])
        self.debug = blocks.message_debug()

        # connections
        self.tb.msg_connect((self.emitter, 'msg'), (self.cf, 'in'))
        self.tb.msg_connect((self.cf, 'out'), (self.debug, 'store'))
        
        # data
        in_data = (1+0j,) * 2048
        i_vec = pmt.init_c32vector(len(in_data), in_data)
        out_data = np.exp(1j*np.arange(0, 2*np.pi*(freq_offset/samp_rate * len(in_data)), 2*np.pi*(freq_offset/samp_rate), dtype=np.complex64))
        e_vec = pmt.init_c32vector(len(out_data), out_data.tolist()) # pmt doesn't play nice with numpy sometimes, convert to list

        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, pmt.intern("sample_rate"), pmt.from_float(samp_rate))
        meta = pmt.dict_add(meta, pmt.intern("center_frequency"), pmt.from_float(center_freq+freq_offset))
        in_pdu = pmt.cons(meta, i_vec)
        e_pdu = pmt.cons(meta, e_vec)

        # flowgraph
        self.tb.start()
        time.sleep(.001)
        self.emitter.emit(in_pdu)
        time.sleep(.01)
        self.tb.stop()
        self.tb.wait()

        # parse output
        #print "got ", list(pmt.to_python(pmt.cdr(self.debug.get_message(0))))
        #print "got ", self.debug.get_message(0)
        rcv = self.debug.get_message(0)
        rcv_meta = pmt.car(rcv)
        rcv_data = pmt.cdr(rcv)
        rcv_cf = pmt.dict_ref(rcv_meta, pmt.intern("center_frequency"), pmt.PMT_NIL)

        # asserts
        self.assertComplexTuplesAlmostEqual(tuple(pmt.c32vector_elements(rcv_data)), tuple(out_data), 2)
        self.assertTrue(pmt.equal(rcv_cf, pmt.from_float(911e6)))

    def test_rms_metadata (self):
        samp_rate = 1e6
        rot_offset = 1.0/16.0
        self.emitter = pdu_utils.message_emitter()
        self.cf = fhss_utils.cf_estimate(fhss_utils.RMS, [])
        self.sm = self.simple_modulator()
        self.rt = self.pdu_rotate(rot_offset)
        self.debug = blocks.message_debug()

        self.tb.msg_connect((self.emitter, 'msg'), (self.sm, 'in'))
        self.tb.msg_connect((self.sm, 'out'), (self.rt, 'in'))
        self.tb.msg_connect((self.rt, 'out'), (self.cf, 'in'))
        self.tb.msg_connect((self.cf, 'out'), (self.debug, 'store'))
        
        # original data
        in_data = [0xAA] * 10 + [0x69] * 10 + [0x55] * 10 # 30 bytes = 240 bits = 1920 samples
        i_vec = pmt.init_u8vector(len(in_data), in_data)

        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, pmt.intern("sample_rate"), pmt.from_float(1e6))
        meta = pmt.dict_add(meta, pmt.intern("center_frequency"), pmt.from_float(100.0e6))
        in_pdu = pmt.cons(meta, i_vec)

        self.tb.start()
        time.sleep(.1)
        self.emitter.emit(in_pdu)
        time.sleep(.1)
        self.tb.stop()
        self.tb.wait()

        # parse output
        #print("got ", list(pmt.to_python(pmt.cdr(self.debug.get_message(0)))))
        #print("got ", pmt.car(self.debug.get_message(0)))
        rcv = self.debug.get_message(0)
        rcv_meta = pmt.car(rcv)
        rcv_data = pmt.cdr(rcv)

        # asserts
        rcv_cf = pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("center_frequency"), pmt.PMT_NIL))
        expected_cf = 100e6 + samp_rate * rot_offset # we rotated the burst 62500 Hz off
        self.assertTrue(abs(rcv_cf - expected_cf) < 100) # less than 100 Hz off

    def test_half_power_metadata (self):
        samp_rate = 1e6
        rot_offset = 1.0/16.0
        self.emitter = pdu_utils.message_emitter()
        self.cf = fhss_utils.cf_estimate(fhss_utils.HALF_POWER, [])
        self.sm = self.simple_modulator()
        self.rt = self.pdu_rotate(rot_offset)
        self.debug = blocks.message_debug()

        self.tb.msg_connect((self.emitter, 'msg'), (self.sm, 'in'))
        self.tb.msg_connect((self.sm, 'out'), (self.rt, 'in'))
        self.tb.msg_connect((self.rt, 'out'), (self.cf, 'in'))
        self.tb.msg_connect((self.cf, 'out'), (self.debug, 'store'))
        
        # original data
        in_data = [0xAA] * 10 + [0x69] * 10 + [0x55] * 10 # 30 bytes = 240 bits = 1920 samples
        i_vec = pmt.init_u8vector(len(in_data), in_data)

        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, pmt.intern("sample_rate"), pmt.from_float(1e6))
        meta = pmt.dict_add(meta, pmt.intern("center_frequency"), pmt.from_float(100.0e6))
        in_pdu = pmt.cons(meta, i_vec)

        self.tb.start()
        time.sleep(.1)
        self.emitter.emit(in_pdu)
        time.sleep(.1)
        self.tb.stop()
        self.tb.wait()

        # parse output
        #print("got ", list(pmt.to_python(pmt.cdr(self.debug.get_message(0)))))
        #print("got ", pmt.car(self.debug.get_message(0)))
        rcv = self.debug.get_message(0)
        rcv_meta = pmt.car(rcv)
        rcv_data = pmt.cdr(rcv)

        # asserts
        rcv_cf = pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("center_frequency"), pmt.PMT_NIL))
        expected_cf = 100e6 + samp_rate * rot_offset # we rotated the burst 62500 Hz off
        self.assertTrue(abs(rcv_cf - expected_cf) < 1.0) # less than 1 Hz off (no noise)


if __name__ == '__main__':
    gr_unittest.run(qa_cf_estimate)
