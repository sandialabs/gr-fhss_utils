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
from gnuradio import blocks
from gnuradio import fft
from gnuradio.fft import window
from gnuradio.filter import firdes

class fft_peak(gr.hier_block2):
    """
    This block returns the peak magnitude FFT index
    """
    def __init__(self, fft_len):
        gr.hier_block2.__init__(self,
            "FFT Peak",
            gr.io_signature(1, 1, gr.sizeof_gr_complex*1),
            gr.io_signaturev(2, 2, [gr.sizeof_float*1, gr.sizeof_float*1]))
        '''
        Constructor
        
        @param fft_len -
        '''

        ##################################################
        # Parameters
        ##################################################
        self.fft_len = fft_len

        ##################################################
        # Blocks
        ##################################################
        self.fft = fft.fft_vcc(fft_len, True, (window.blackmanharris(fft_len)), False, 1)
        self.v2s = blocks.vector_to_stream(gr.sizeof_float*1, fft_len)
        self.s2v = blocks.stream_to_vector(gr.sizeof_gr_complex*1, fft_len)
        self.s2f = blocks.short_to_float(1, 1)
        self.repeat = blocks.repeat(gr.sizeof_short*1, fft_len)
        self.null_0 = blocks.null_sink(gr.sizeof_float*1)
        self.null_1 = blocks.null_sink(gr.sizeof_short*1)
        self.c2mag = blocks.complex_to_mag(fft_len)
        self.argmax = blocks.argmax_fs(fft_len)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.argmax, 1), (self.null_1, 0))
        self.connect((self.argmax, 0), (self.repeat, 0))
        self.connect((self.c2mag, 0), (self.argmax, 0))
        self.connect((self.c2mag, 0), (self.v2s, 0))
        self.connect((self.repeat, 0), (self.s2f, 0))
        self.connect((self.s2f, 0), (self, 0))
        self.connect((self.s2v, 0), (self.fft, 0))
        self.connect((self.v2s, 0), (self.null_0, 0))
        self.connect((self.v2s, 0), (self, 1))
        self.connect((self.fft, 0), (self.c2mag, 0))
        self.connect((self, 0), (self.s2v, 0))

    def get_fft_len(self):
        return self.fft_len

    def set_fft_len(self, fft_len):
        self.fft_len = fft_len
        self.repeat.set_interpolation(self.fft_len)
