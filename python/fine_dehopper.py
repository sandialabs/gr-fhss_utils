#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 <+YOU OR YOUR COMPANY+>.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

from gnuradio import gr
from gnuradio import analog
from gnuradio import blocks
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio.filter import pfb
from math import pi
import math
from s_and_h_detector import s_and_h_detector  # hier_block

class fine_dehopper(gr.hier_block2):
    """
    docstring for block fine_dehopper
    """
    def __init__(self, bias, freq_sample_delay_samps, freq_samps_to_avg, mag_samps_to_avg, resamp_rate, thresh):
        gr.hier_block2.__init__(self,
            "Fine Dehopper",
            gr.io_signature(1, 1, gr.sizeof_gr_complex*1),
            gr.io_signature(1, 1, gr.sizeof_gr_complex*1))

        ##################################################
        # Parameters
        ##################################################
        self.bias = bias
        self.freq_sample_delay_samps = freq_sample_delay_samps
        self.freq_samps_to_avg = freq_samps_to_avg
        self.mag_samps_to_avg = mag_samps_to_avg
        self.resamp_rate = resamp_rate
        self.thresh = thresh

        ##################################################
        # Blocks
        ##################################################
        self.s_and_h_detector = s_and_h_detector(
            freq_sample_delay_samps=freq_sample_delay_samps,
            freq_samps_to_avg=freq_samps_to_avg,
            mag_samps_to_avg=mag_samps_to_avg,
            thresh=thresh,
        )
        self.resamp = pfb.arb_resampler_ccf(resamp_rate * 2.0, taps=None, flt_size=32)
        self.resamp.declare_sample_delay(0)
        self.fir = filter.fir_filter_ccc(2, (firdes.low_pass_2(1,1,.25,.05,60)))
        self.fir.declare_sample_delay(0)
        self.vco = blocks.vco_c(1, 1, 1)
        self.mult_conj = blocks.multiply_conjugate_cc(1)
        self.delay = blocks.delay(gr.sizeof_gr_complex*1, int(freq_samps_to_avg) + freq_sample_delay_samps)
        self.c2mag = blocks.complex_to_mag(1)
        self.add_const = blocks.add_const_vff((-1.0 * bias * (resamp_rate), ))
        self.demod = analog.quadrature_demod_cf(1)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.demod, 0), (self.s_and_h_detector, 1))
        self.connect((self.add_const, 0), (self.vco, 0))
        self.connect((self.c2mag, 0), (self.s_and_h_detector, 0))
        self.connect((self.delay, 0), (self.mult_conj, 0))
        self.connect((self.mult_conj, 0), (self.fir, 0))
        self.connect((self.vco, 0), (self.mult_conj, 1))
        self.connect((self.fir, 0), (self.resamp, 0))
        self.connect((self, 0), (self.demod, 0))
        self.connect((self, 0), (self.c2mag, 0))
        self.connect((self, 0), (self.delay, 0))
        self.connect((self.resamp, 0), (self, 0))
        self.connect((self.s_and_h_detector, 0), (self.add_const, 0))

    def get_bias(self):
        return self.bias

    def set_bias(self, bias):
        self.bias = bias
        self.add_const.set_k((-1.0 * self.bias * (self.resamp_rate), ))

    def get_freq_sample_delay_samps(self):
        return self.freq_sample_delay_samps

    def set_freq_sample_delay_samps(self, freq_sample_delay_samps):
        self.freq_sample_delay_samps = freq_sample_delay_samps
        self.s_and_h_detector.set_freq_sample_delay_samps(self.freq_sample_delay_samps)
        self.delay.set_dly(int(self.freq_samps_to_avg) + self.freq_sample_delay_samps)

    def get_freq_samps_to_avg(self):
        return self.freq_samps_to_avg

    def set_freq_samps_to_avg(self, freq_samps_to_avg):
        self.freq_samps_to_avg = freq_samps_to_avg
        self.s_and_h_detector.set_freq_samps_to_avg(self.freq_samps_to_avg)
        self.delay.set_dly(int(self.freq_samps_to_avg) + self.freq_sample_delay_samps)

    def get_mag_samps_to_avg(self):
        return self.mag_samps_to_avg

    def set_mag_samps_to_avg(self, mag_samps_to_avg):
        self.mag_samps_to_avg = mag_samps_to_avg
        self.s_and_h_detector.set_mag_samps_to_avg(self.mag_samps_to_avg)

    def get_resamp_rate(self):
        return self.resamp_rate

    def set_resamp_rate(self, resamp_rate):
        self.resamp_rate = resamp_rate
        self.resamp.set_rate(self.resamp_rate * 2.0)
        self.add_const.set_k((-1.0 * self.bias * (self.resamp_rate), ))

    def get_thresh(self):
        return self.thresh

    def set_thresh(self, thresh):
        self.thresh = thresh
        self.s_and_h_detector.set_thresh(self.thresh)
