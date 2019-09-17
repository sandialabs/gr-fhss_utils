#!/usr/bin/env python
# -*- coding: utf-8 -*- #
# Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
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
from gnuradio import blocks
from gnuradio.filter import firdes
import fhss_utils
import pdu_utils


class fsk_burst_extractor_hier(gr.hier_block2):

    def __init__(self, burst_width=int(500e3), center_freq=915e6, decimation=32, fft_size=256, hist_time=0.004, lookahead_time=0.0005, max_burst_time=0.5, min_burst_time=0.001, output_attenuation=40, output_cutoff=0.26, output_trans_width=0.4, post_burst_time=0.00008, pre_burst_time=0.00008, samp_rate=int(16e6), cfo_samps_to_average=2048, cfo_bin_resolution=4000, threshold=6):
        gr.hier_block2.__init__(
            self, "FSK Burst Extractor Hier",
            gr.io_signature(1, 1, gr.sizeof_gr_complex*1),
            gr.io_signature(0, 0, 0),
        )
        self.message_port_register_hier_out("pdu_out")

        ##################################################
        # Parameters
        ##################################################
        self.burst_width = burst_width
        self.center_freq = center_freq
        self.decimation = decimation
        self.fft_size = fft_size
        self.hist_time = hist_time
        self.lookahead_time = lookahead_time
        self.max_burst_time = max_burst_time
        self.min_burst_time = min_burst_time
        self.output_attenuation = output_attenuation
        self.output_cutoff = output_cutoff
        self.output_trans_width = output_trans_width
        self.post_burst_time = post_burst_time
        self.pre_burst_time = pre_burst_time
        self.samp_rate = samp_rate
        self.threshold = threshold

        ##################################################
        # Blocks
        ##################################################
        # Low pass filter cutoff to half band.
        self.pdu_utils_pdu_fir_filter_1 = pdu_utils.pdu_fir_filter(1, (firdes.low_pass_2(1, 1, .25, .1, output_attenuation)))
        # This is a coarse filter,  Allow for transition band to alias onto itself.
        taps = firdes.low_pass_2(1, 1, output_cutoff/decimation, output_trans_width/decimation, output_attenuation)
        self.fhss_utils_tagged_burst_to_pdu_0 = fhss_utils.tagged_burst_to_pdu(decimation, taps, min_burst_time, max_burst_time, 0.0, 1.0, 1.0, samp_rate, 3)
        self.fhss_utils_fft_burst_tagger_0 = fhss_utils.fft_burst_tagger(center_freq, fft_size, samp_rate, int(round((float(samp_rate)/fft_size)*pre_burst_time)), int(round((float(samp_rate)/fft_size)*post_burst_time)), burst_width, 0, 0, threshold, int(round((float(samp_rate)/fft_size)*hist_time)), int(round((float(samp_rate)/fft_size)*lookahead_time)), False)
        (self.fhss_utils_fft_burst_tagger_0).set_min_output_buffer(2048000*2)
        self.fhss_utils_fine_burst_measure = fhss_utils.fine_burst_measure(cfo_bin_resolution, cfo_samps_to_average, .5)

        self.fine_time_measure = pdu_utils.pdu_fine_time_measure(pre_burst_time, post_burst_time, 10, 15)

        ##################################################
        # Connections
        ##################################################
        #self.msg_connect((self.pdu_utils_pdu_fir_filter_1, 'pdu_out'), (self, 'pdu_out'))
        self.msg_connect((self.fine_time_measure, 'pdu_out'), (self, 'pdu_out'))
        self.msg_connect((self.pdu_utils_pdu_fir_filter_1, 'pdu_out'), (self.fine_time_measure, 'pdu_in'))
        self.msg_connect((self.fhss_utils_fine_burst_measure, 'pdu_out'), (self.pdu_utils_pdu_fir_filter_1, 'pdu_in'))
        self.msg_connect((self.fhss_utils_tagged_burst_to_pdu_0, 'cpdus'), (self.fhss_utils_fine_burst_measure, 'pdu_in'))
        self.connect((self.fhss_utils_fft_burst_tagger_0, 0), (self.fhss_utils_tagged_burst_to_pdu_0, 0))
        self.connect((self, 0), (self.fhss_utils_fft_burst_tagger_0, 0))

    def reset(self):
        self.fhss_utils_fft_burst_tagger_0.reset()

    def get_burst_width(self):
        return self.burst_width

    def get_center_freq(self):
        return self.center_freq

    def get_decimation(self):
        return self.decimation

    def get_fft_size(self):
        return self.fft_size

    def get_hist_time(self):
        return self.hist_time

    def get_lookahead_time(self):
        return self.lookahead_time

    def get_max_burst_time(self):
        return self.max_burst_time

    def get_min_burst_time(self):
        return self.min_burst_time

    def get_output_attenuation(self):
        return self.output_attenuation

    def get_output_cutoff(self):
        return self.output_cutoff

    def get_output_trans_width(self):
        return self.output_trans_width

    def get_post_burst_time(self):
        return self.post_burst_time

    def get_pre_burst_time(self):
        return self.pre_burst_time

    def get_samp_rate(self):
        return self.samp_rate

    def get_threshold(self):
        return self.threshold
