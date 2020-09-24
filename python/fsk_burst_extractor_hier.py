#!/usr/bin/env python
# -*- coding: utf-8 -*- #
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
import fhss_utils
from .fhss_utils_swig import RMS
from .fhss_utils_swig import HALF_POWER
from .fhss_utils_swig import COERCE
import pdu_utils

class fsk_burst_extractor_hier(gr.hier_block2):
  """
  The burst extractor block isolates bursts in time and frequency and generates
  individual bursts that have been basedbanded and decimated.  The underlying
  block structure is:


         |---------|     |-----------|     |---------|     |---------|
         |   FFT   |     |   Tagged  |     |  Fine   |     |  Fine   |
  In---->|  Burst  |---->|  Burst To |---->|  Burst  |---->|  Time   |----> Out
         |  Tagger |     |    PDU    |     | Measure |     | Measure |
         |---------|     |-----------|     |---------|     |---------|


  The signal processing performed on each burst after a signal has been detected
  and a center frequency has been determined is

         |---------|     |------------|     |----------|     |---------|
         | CORDIC  |     |    LPF     |     |          |     |  LPF    |
  In---->|  Tune   |---->| fc =.45*fs |---->| Decimate |---->| fc = b  |----> Out
         |         |     |            |     |          |     |         |
         |---------|     |------------|     |----------|     |---------|


  where
    fs = input sampling rate
    fc = single-sided filter cutoff bandwidth
    LPF = Low-pass filter
    N = decimation factor
    b = specified output cutoff (normalized to output sampling rate of fs / N)


  Frequency Estimation:
    The initial burst tagging stage does a coarse frequency estimation based on
    bin energy.  The fine_burst_measure block is responsible for refining the
    signal center frequency estimate.  If a set of channel center frequencies is
    not specified, the frequency estimation uses a frequency histogram-based method
    to determine the signal center frequency.

    Parameters:
      * CFO Samps To Average: Block processing size (no overlap)
      * CFO Bin Resolution: Bin frequency resolution (Hz)
      * Channel Center Frequencies: Specified channel center frequencies.  If
          specified, each detected burst will be coerced to the nearest channel
          center frequency
  """
  def __init__(self, burst_width=int(500e3), center_freq=915e6, decimation=32,
               fft_size=256, hist_time=0.004, lookahead_time=0.0005,
               max_burst_time=0.5, min_burst_time=0.001,
               output_attenuation=40, output_cutoff=0.25,
               output_trans_width=0.1, post_burst_time=0.00008,
               pre_burst_time=0.00008, samp_rate=int(16e6),
               threshold=10,
               cf_method = RMS, channel_freqs = [],
               n_threads = 3):
      gr.hier_block2.__init__(
          self, "FSK Burst Extractor Hier",
          gr.io_signature(1, 1, gr.sizeof_gr_complex*1),
          gr.io_signature(0, 0, 0),
      )
      '''
      Constructor
      
      @param burst_width - max burst bandwidth, Hz
      @param center_freq - 
      @param decimation - 
      @param fft_size - 
      @param hist_time - 
      @param lookahead_time - 
      @param max_burst_time - 
      @param min_burst_time - 
      @param output_attenuation - 
      @param output_cutoff -
      @param output_trans_width - 
      @param post_burst_time - 
      @param pre_burst_time - 
      @param samp_rate - 
      @param threshold - 
      @param cf_method - Center Frequency estimation method
      @param channel_freqs - CF Coerce freq list
      @param n_threads - 
      '''
      
      
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
      self.channel_freqs = channel_freqs
      self.cf_method = cf_method
      self.n_threads = n_threads

      ##################################################
      # Blocks
      ##################################################
      # Low pass filter cutoff to half band.
      sig_taps = firdes.low_pass_2(1, 1, output_cutoff, output_trans_width, output_attenuation)
      self.pdu_utils_pdu_fir_filter_1 = pdu_utils.pdu_fir_filter(1, sig_taps)

      # This is a coarse filter,  Allow for transition band to alias onto itself.
      taps = firdes.low_pass_2(1, 1, .45 / decimation, .1 / decimation, output_attenuation)
      self.fhss_utils_tagged_burst_to_pdu_0 = fhss_utils.tagged_burst_to_pdu(decimation, taps, min_burst_time, max_burst_time, 0.0, 1.0, 1.0, samp_rate, n_threads)
      self.fhss_utils_fft_burst_tagger_0 = fhss_utils.fft_burst_tagger(center_freq, fft_size, samp_rate, int(round((float(samp_rate)/fft_size)*pre_burst_time)), int(round((float(samp_rate)/fft_size)*post_burst_time)), burst_width, 0, 0, threshold, int(round((float(samp_rate)/fft_size)*hist_time)), int(round((float(samp_rate)/fft_size)*lookahead_time)), False)
      #(self.fhss_utils_fft_burst_tagger_0).set_min_output_buffer(102400)
      self.cf_estimate = fhss_utils.cf_estimate(self.cf_method, self.channel_freqs)

      self.fine_time_measure = pdu_utils.pdu_fine_time_measure(pre_burst_time, post_burst_time, 10, 15)

      ##################################################
      # Connections
      ##################################################
      #self.msg_connect((self.pdu_utils_pdu_fir_filter_1, 'pdu_out'), (self, 'pdu_out'))
      self.msg_connect((self.fine_time_measure, 'pdu_out'), (self, 'pdu_out'))
      self.msg_connect((self.pdu_utils_pdu_fir_filter_1, 'pdu_out'), (self.fine_time_measure, 'pdu_in'))
      self.msg_connect((self.cf_estimate, 'out'), (self.pdu_utils_pdu_fir_filter_1, 'pdu_in'))
      self.msg_connect((self.fhss_utils_tagged_burst_to_pdu_0, 'cpdus'), (self.cf_estimate, 'in'))
      self.connect((self.fhss_utils_fft_burst_tagger_0, 0), (self.fhss_utils_tagged_burst_to_pdu_0, 0))
      self.connect((self, 0), (self.fhss_utils_fft_burst_tagger_0, 0))

  def reset(self):
      self.fhss_utils_fft_burst_tagger_0.reset()

  def set_channel_freqs(self, freqs):
    self.cf_estimate.set_channel_freqs(freqs)

  def set_cf_method(self, cf_method):
    self.cf_estimate.set_method(cf_method)

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
