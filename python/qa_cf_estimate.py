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
import numpy as np
import pdu_utils
import pmt
import time

class qa_cf_estimate (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_coerce (self):
        self.emitter = pdu_utils.message_emitter()
        self.cf = fhss_utils.cf_estimate(fhss_utils.COERCE, [x*1e6 for x in range(900,930)])
        self.debug = blocks.message_debug()
        self.tb.msg_connect((self.emitter, 'msg'), (self.cf, 'in'))
        self.tb.msg_connect((self.cf, 'out'), (self.debug, 'store'))
        
        # original data
        in_data = np.exp(1j*np.array(np.linspace(0,10*np.pi*2, 20)))
        i_vec = pmt.init_c32vector(len(in_data), in_data)

        out_data = [(1+0j), (-0.9863691-0.16454819j), (0.9458478+0.32461047j), (-0.8795409-0.47582325j), (0.7892561+0.61406416j), (-0.67745465-0.7355646j), (0.54718447+0.8370121j), (-0.40199703-0.915641j), (0.24585037+0.9693079j), (-0.083001345-0.9965495j), (-0.08211045+0.99662334j), (0.24498378-0.96952736j), (-0.4011784+0.9160001j), (0.5464361-0.83750105j), (-0.6767969+0.73617005j), (0.78870696-0.6147696j), (-0.87911534+0.47660938j), (0.94555736-0.3254559j), (-0.9862217+0.16542989j), (0.9999997-0.00089395075j)]
        e_vec = pmt.init_c32vector(len(out_data), out_data)

        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, pmt.intern("sample_rate"), pmt.from_float(1e3))
        meta = pmt.dict_add(meta, pmt.intern("center_frequency"), pmt.from_float(910.6e6))
        in_pdu = pmt.cons(meta, i_vec)
        e_pdu = pmt.cons(meta, e_vec)

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

    def test_rms (self):
        self.emitter = pdu_utils.message_emitter()
        self.cf = fhss_utils.cf_estimate(fhss_utils.RMS, [])
        self.debug = blocks.message_debug()
        self.tb.msg_connect((self.emitter, 'msg'), (self.cf, 'in'))
        self.tb.msg_connect((self.cf, 'out'), (self.debug, 'store'))
        
        # original data
        in_data = np.exp(1j*np.array(np.linspace(0,1*np.pi*.02, 20)))
        i_vec = pmt.init_c32vector(len(in_data), in_data)

        out_data = [(1+0j), (0.9999966+0.0026077442j), (0.9999864+0.0052154697j), (0.9999694+0.007823161j), (0.99994564+0.010430798j), (0.99991506+0.013038365j), (0.99987763+0.015645843j), (0.99983346+0.018253215j), (0.99978244+0.020860463j), (0.9997247+0.023467569j), (0.99966+0.026074518j), (0.99958867+0.028681284j), (0.9995105+0.03128786j), (0.9994256+0.033894222j), (0.99933374+0.03650035j), (0.99923515+0.03910623j), (0.9991298+0.04171185j), (0.99901766+0.044317182j), (0.9988987+0.046922214j), (0.9987729+0.04952693j)]
        e_vec = pmt.init_c32vector(len(out_data), out_data)

        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, pmt.intern("sample_rate"), pmt.from_float(1e6))
        meta = pmt.dict_add(meta, pmt.intern("center_frequency"), pmt.from_float(910.6e6))
        in_pdu = pmt.cons(meta, i_vec)
        e_pdu = pmt.cons(meta, e_vec)

        self.tb.start()
        time.sleep(.001)
        self.emitter.emit(in_pdu)
        time.sleep(.01)
        self.tb.stop()
        self.tb.wait()

        # parse output
        #print("got ", list(pmt.to_python(pmt.cdr(self.debug.get_message(0)))))
        #print("got ", self.debug.get_message(0))
        rcv = self.debug.get_message(0)
        rcv_meta = pmt.car(rcv)
        rcv_data = pmt.cdr(rcv)
        rcv_cf = pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("center_frequency"), pmt.PMT_NIL))

        # asserts
        self.assertComplexTuplesAlmostEqual(tuple(pmt.c32vector_elements(rcv_data)), tuple(out_data), 2)
        self.assertTrue(abs(rcv_cf - 910.6001e6) < 100)

    def test_half_power (self):
        self.emitter = pdu_utils.message_emitter()
        self.cf = fhss_utils.cf_estimate(fhss_utils.HALF_POWER, [])
        self.debug = blocks.message_debug()
        self.tb.msg_connect((self.emitter, 'msg'), (self.cf, 'in'))
        self.tb.msg_connect((self.cf, 'out'), (self.debug, 'store'))
        
        # original data
        in_data = np.exp(1j*np.array(np.linspace(0,1*np.pi*.02, 20)))
        i_vec = pmt.init_c32vector(len(in_data), in_data)

        out_data = [(1+0j), (0.9999945+0.0033069337j), (0.9999781+0.006613831j), (0.99995077+0.009920656j), (0.9999125+0.013227372j), (0.9998633+0.016533945j), (0.9998032+0.019840335j), (0.9997321+0.02314651j), (0.99965006+0.026452431j), (0.99955714+0.029758062j), (0.99945325+0.03306337j), (0.99933845+0.036368314j), (0.99921274+0.039672863j), (0.99907607+0.042976975j), (0.9989285+0.04628062j), (0.99877+0.049583755j), (0.99860054+0.05288635j), (0.9984202+0.056188367j), (0.9982289+0.059489768j), (0.9980267+0.06279052j)]
        e_vec = pmt.init_c32vector(len(out_data), out_data)

        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, pmt.intern("sample_rate"), pmt.from_float(1e6))
        meta = pmt.dict_add(meta, pmt.intern("center_frequency"), pmt.from_float(910.6e6))
        in_pdu = pmt.cons(meta, i_vec)
        e_pdu = pmt.cons(meta, e_vec)

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
        rcv_cf = pmt.to_double(pmt.dict_ref(rcv_meta, pmt.intern("center_frequency"), pmt.PMT_NIL))

        # asserts
        self.assertComplexTuplesAlmostEqual(tuple(pmt.c32vector_elements(rcv_data)), tuple(out_data), 2)
        self.assertTrue(abs(rcv_cf - 910.6e6) < 100)


if __name__ == '__main__':
    gr_unittest.run(qa_cf_estimate, "qa_cf_estimate.xml")
