#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018, 2019, 2020 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
# retains certain rights in this software.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr
from gnuradio import blocks
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio.filter import pfb
from math import pi
from .fft_peak import fft_peak  # hier_block
from .s_and_h_detector import s_and_h_detector  # hier_block


class coarse_dehopper(gr.hier_block2):
    """
    Coarse FFT based dehopper
    """
    def __init__(self, fft_len, freq_sample_delay_samps, freq_samps_to_avg, mag_samps_to_avg, thresh):
        gr.hier_block2.__init__(self,
            "Coarse Dehopper",
            gr.io_signature(1, 1, gr.sizeof_gr_complex*1),
            gr.io_signature(1, 1, gr.sizeof_gr_complex*1))
        '''
        Constructor
        
        @param fft_len - 
        @param freq_sample_delay_samps - 
        @param freq_samps_to_avg -
        @param mag_samps_to_avg - 
        @param thresh - 
        '''

        ##################################################
        # Parameters
        ##################################################
        self.fft_len = fft_len
        self.freq_sample_delay_samps = freq_sample_delay_samps
        self.freq_samps_to_avg = freq_samps_to_avg
        self.mag_samps_to_avg = mag_samps_to_avg
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
        self.resamp = pfb.arb_resampler_ccf(1.0 / (fft_len / 4.0), taps=None, flt_size=32)
        self.resamp.declare_sample_delay(0)
        self.fir = filter.fir_filter_ccc(2, (firdes.low_pass_2(1,1,.30,.05,60)))
        self.fir.declare_sample_delay(0)
        self.fft_peak = fft_peak(fft_len=fft_len)
        self.vco = blocks.vco_c(1, 2.0 * pi / fft_len, 1)
        self.mult_conj = blocks.multiply_conjugate_cc(1)
        self.delay = blocks.delay(gr.sizeof_gr_complex*1, int(freq_samps_to_avg) + freq_sample_delay_samps)
        self.c2mag = blocks.complex_to_mag(1)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.c2mag, 0), (self.s_and_h_detector, 0))
        self.connect((self.delay, 0), (self.mult_conj, 0))
        self.connect((self.mult_conj, 0), (self.fir, 0))
        self.connect((self.vco, 0), (self.mult_conj, 1))
        self.connect((self.fft_peak, 0), (self.s_and_h_detector, 1))
        self.connect((self.fir, 0), (self.resamp, 0))
        self.connect((self, 0), (self.c2mag, 0))
        self.connect((self, 0), (self.delay, 0))
        self.connect((self, 0), (self.fft_peak, 0))
        self.connect((self.resamp, 0), (self, 0))
        self.connect((self.s_and_h_detector, 0), (self.vco, 0))

    def get_fft_len(self):
        return self.fft_len

    def set_fft_len(self, fft_len):
        self.fft_len = fft_len
        self.resamp.set_rate(1.0 / (self.fft_len / 4.0))
        self.fft_peak.set_fft_len(self.fft_len)

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

    def get_thresh(self):
        return self.thresh

    def set_thresh(self, thresh):
        self.thresh = thresh
        self.s_and_h_detector.set_thresh(self.thresh)
