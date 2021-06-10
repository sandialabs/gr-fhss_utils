#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2020 gr-fhss_utils author.
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


import numpy
import json
import pmt
from gnuradio import gr
from os.path import splitext
from math import isnan

class sigmf_meta_writer(gr.basic_block):
    """
    quick and dirty tool to convert detections to sigmf annotations
    """
    def __init__(self, filename, freq, rate, label, dtype):
        gr.basic_block.__init__(self,
            name="sigmf_meta_writer",
            in_sig=None,
            out_sig=None)

        self.d_filename = filename
        if not filename.endswith('.sigmf-meta'):
            pre, ext = splitext(filename)
            self.d_filename = pre + '.sigmf-meta'
            gr.log.warn("SigMF metadata filename does not end with `sigmf-meta` - using " + self.d_filename)


        self.freq = freq
        self.rate = rate
        self.soo = 0
        self.bw_min = rate/1000.0

        self.label = label

        self.initialize_sigmf_dict([{'core:sample_start': 0, 'core:frequency': freq}],
                                   {'core:datatype': dtype, 'core:sample_rate': rate, 'antenna:gain': 0})

        self.message_port_register_in(pmt.intern("in"))
        self.set_msg_handler(pmt.intern("in"), self.handler)

    def stop(self):
        try:
            f = open(self.d_filename, 'w+')
            f.write(json.dumps(self.d_dict,indent=4))
            f.close()
        except IOError as e:
            print("ERROR: could write to {}".format(self.d_filename), "because", e)
            quit()

        return True

    def initialize_sigmf_dict(self, sigmf_captures, sigmf_global, sigmf_annotations = []):
        self.d_dict = {}
        self.d_dict['captures'] = sigmf_captures
        self.d_dict['global'] = sigmf_global
        self.d_dict['annotations'] = sigmf_annotations

    def handler(self, pdu):
      if not pmt.is_pdu(pdu):
        print('input is not a PDU!, dropping')

      # there are two basic modes here: tags_to_pdu or fft burst detector
      # in either case we need to extract the following fields for the annotation:
      #     - sob:    start sample of the burst
      #     - eob:    end sample of the burst
      #     - freq:   center frequency of the burst in hz
      #     - bw:     bandwidth of the burst in hz
      #     - b_id:   unique Identifier for the burst or `None`
      #     - snr:    signal to noise ratio of annotation or 'None'
      #
      # how these are obtained differs between the two modes.

      meta = pmt.car(pdu)

      time_pmt = pmt.dict_ref(meta, pmt.intern('burst_time'), pmt.PMT_NIL)
      if pmt.is_tuple(time_pmt):
        # tags_to_pdu mode
        try:
          burst_time = pmt.to_double(pmt.tuple_ref(time_pmt, 1)) + pmt.to_uint64(pmt.tuple_ref(time_pmt, 0))
          pdu_rate = pmt.to_double(pmt.dict_ref(meta, pmt.intern('sample_rate'), pmt.from_double(self.rate)))
          freq = pmt.to_double(pmt.dict_ref(meta, pmt.intern('center_frequency'), pmt.from_double(self.freq)))
          bw = pmt.to_double(pmt.dict_ref(meta, pmt.intern('bandwidth'), pmt.from_double(self.bw_min)))

          if bw < self.bw_min:
              bw = self.bw_min

          # these can be `None` so use to_python()
          snr = pmt.to_python(pmt.dict_ref(meta, pmt.intern('snr_db'), pmt.PMT_NIL))
          b_id = pmt.to_python(pmt.dict_ref(meta, pmt.intern('pdu_num'), pmt.PMT_NIL))

          anno_len = int(pmt.length(pmt.cdr(pdu)) * (self.rate / pdu_rate))
          sob = int(self.rate * burst_time)
          eob = sob + anno_len

        except Exception as e:
          print('could not parse required data from message', pmt.car(pdu), ':',e)
          return

      else:
        # fft burst detector mode
        try:
          sob = pmt.to_uint64(pmt.dict_ref(meta, pmt.intern('start_offset'), pmt.PMT_NIL))
          eob = pmt.to_uint64(pmt.dict_ref(meta, pmt.intern('end_offset'), pmt.PMT_NIL))
          freq = pmt.to_double(pmt.dict_ref(meta, pmt.intern('center_frequency'), pmt.PMT_NIL))
          bw = pmt.to_double(pmt.dict_ref(meta, pmt.intern('bandwidth'), pmt.from_double(self.bw_min)))

          if bw < self.bw_min:
              bw = self.bw_min

          # these can be `None` so use to_python()
          snr = pmt.to_python(pmt.dict_ref(meta, pmt.intern('snr_db'), pmt.PMT_NIL))
          b_id = pmt.to_python(pmt.dict_ref(meta, pmt.intern('burst_id'), pmt.PMT_NIL))

        except Exception as e:
          print('could not parse required data from message', pmt.car(pdu), ':',e)
          return


      label = self.label
      if self.label == 'use_burst_id':
        if b_id is None:
          label = ''
        else:
          label = 'burst' + str(b_id)

      elif self.label == 'use_snr_db':
        # this probably isnt in here so it will end up blank...
        label = str(snr) + 'dB'

      # append the annotation
      try:
        if isnan(snr):
          print("Got illegal SNR value in",meta)
          self.d_dict['annotations'].append({'core:sample_start': sob-self.soo,
                      'core:sample_count': eob-sob, 'core:freq_upper_edge': int(freq+bw/2),
                      'core:freq_lower_edge': int(freq-bw/2), 'core:description': label})
        else:
          self.d_dict['annotations'].append({'core:sample_start': sob-self.soo,
                      'core:sample_count': eob-sob, 'core:freq_upper_edge': int(freq+bw/2),
                      'core:freq_lower_edge': int(freq-bw/2), 'core:description': label,
                      'capture_details:SNRdB': snr})
      except Exception as e:
        print('could not form annotation from message', pmt.car(pdu), ':', e)
