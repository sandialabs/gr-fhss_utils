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
from gnuradio.filter import firdes
from gnuradio import timing_utils


class s_and_h_detector(gr.hier_block2):
    """
    Sample and hold detector block.
    """

    def __init__(self, freq_sample_delay_samps, freq_samps_to_avg, mag_samps_to_avg, thresh):
        gr.hier_block2.__init__(self,
                                "Sample and Hold Detector",
                                gr.io_signaturev(2, 2, [gr.sizeof_float * 1, gr.sizeof_float * 1]),
                                gr.io_signaturev(4, 4, [gr.sizeof_float * 1, gr.sizeof_float * 1, gr.sizeof_float * 1, gr.sizeof_float * 1]))
        '''
        Constructor
        
        @param freq_sample_delay_samps - 
        @param freq_samps_to_avg - 
        @param mag_samps_to_avg - 
        @param thresh - 
        
        '''

        ##################################################
        # Parameters
        ##################################################
        self.freq_sample_delay_samps = freq_sample_delay_samps
        self.freq_samps_to_avg = freq_samps_to_avg
        self.mag_samps_to_avg = mag_samps_to_avg
        self.thresh = thresh

        ##################################################
        # Blocks
        ##################################################
        self.edge_detector = timing_utils.edge_detector_bb(timing_utils.RISING_EDGE)
        self.threshold = blocks.threshold_ff(thresh / 4.0, thresh, 0)
        self.samp_hold = blocks.sample_and_hold_ff()
        self.mag_avg = blocks.moving_average_ff(int(mag_samps_to_avg), 1.0 / (mag_samps_to_avg), 4000)
        self.freq_avg = blocks.moving_average_ff(int(freq_samps_to_avg), 1.0 / (freq_samps_to_avg), 4000)
        self.f2c = blocks.float_to_char(1, 1)
        self.delay = blocks.delay(gr.sizeof_float * 1, int(freq_samps_to_avg - mag_samps_to_avg + freq_sample_delay_samps))

        ##################################################
        # Connections
        ##################################################
        self.connect((self.delay, 0), (self.mag_avg, 0))
        self.connect((self.f2c, 0), (self.edge_detector, 0))
        self.connect((self.freq_avg, 0), (self.samp_hold, 0))
        self.connect((self.freq_avg, 0), (self, 1))
        self.connect((self.mag_avg, 0), (self.threshold, 0))
        self.connect((self.mag_avg, 0), (self, 3))
        self.connect((self.samp_hold, 0), (self, 0))
        self.connect((self.threshold, 0), (self.f2c, 0))
        self.connect((self.threshold, 0), (self, 2))
        self.connect((self, 0), (self.delay, 0))
        self.connect((self, 1), (self.freq_avg, 0))
        self.connect((self.edge_detector, 0), (self.samp_hold, 1))

    def get_freq_sample_delay_samps(self):
        return self.freq_sample_delay_samps

    def set_freq_sample_delay_samps(self, freq_sample_delay_smps):
        self.freq_sample_delay_samps = freq_sample_delay_samps
        self.delay.set_dly(int(self.freq_samps_to_avg - self.mag_samps_to_avg + self.freq_sample_delay_samps))

    def get_freq_samps_to_avg(self):
        return self.freq_samps_to_avg

    def set_freq_samps_to_avg(self, freq_samps_to_avg):
        self.freq_samps_to_avg = freq_samps_to_avg
        self.freq_avg.set_length_and_scale(int(self.freq_samps_to_avg), 1.0 / (self.freq_samps_to_avg))
        self.delay.set_dly(int(self.freq_samps_to_avg - self.mag_samps_to_avg + self.freq_sample_delay_samps))

    def get_mag_samps_to_avg(self):
        return self.mag_samps_to_avg

    def set_mag_samps_to_avg(self, mag_samps_to_avg):
        self.mag_samps_to_avg = mag_samps_to_avg
        self.mag_avg.set_length_and_scale(int(self.mag_samps_to_avg), 1.0 / (self.mag_samps_to_avg))
        self.delay.set_dly(int(self.freq_samps_to_avg - self.mag_samps_to_avg + self.freq_sample_delay_samps))

    def get_thresh(self):
        return self.thresh

    def set_thresh(self, thresh):
        self.thresh = thresh
        self.threshold.set_hi(self.thresh)
        self.threshold.set_lo(self.thresh / 4.0)
