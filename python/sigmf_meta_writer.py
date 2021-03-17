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
        if filename.endswith('.sigmf-data'):
            pre, ext = splitext(filename)
            self.d_filename = pre + '.sigmf-meta'
            gr.log.warn("SigMF metadata filename ends with `sigmf-data`! Saving you from yourself and using " + self.d_filename)
        self.soo = 256
        self.bw_min = rate/1000.0

        self.label = label

        self.d_dict = {}
        self.d_dict['captures'] = [{'core:sample_start': 0, 'core:frequency': freq}]
        self.d_dict['global'] = {'core:datatype': dtype, 'core:sample_rate': rate, 'antenna:gain': 0}
        self.d_dict['annotations'] = []

        self.message_port_register_in(pmt.intern("in"))
        self.set_msg_handler(pmt.intern("in"), self.handler)

    def stop(self):
        try:
            f = open(self.d_filename, 'w+')
        except IOError:
            print("ERROR: could not open {}".format(self.d_filename))
            quit()

        f.write(json.dumps(self.d_dict,indent=4))
        f.close()

        return True

    def handler(self, msg):
      try:
        meta = pmt.car(msg)
      except:
        print('msg is not a pair')

      try:
        sob = pmt.to_uint64(pmt.dict_ref(meta, pmt.intern('start_offset'), pmt.PMT_NIL))
        eob = pmt.to_uint64(pmt.dict_ref(meta, pmt.intern('end_offset'), pmt.PMT_NIL))
        freq = pmt.to_double(pmt.dict_ref(meta, pmt.intern('center_frequency'), pmt.PMT_NIL))

        burst = self.label
        if self.label == 'use_burst_id':
            try:
              burst = 'burst' + str(pmt.to_uint64(pmt.dict_ref(meta, pmt.intern('burst_id'), pmt.PMT_NIL)))
            except:
              burst = ''
        elif self.label == 'use_snr_db':
            try:
              burst = str(round(pmt.to_double(pmt.dict_ref(meta, pmt.intern('snr_db'), pmt.PMT_NIL)),1)) + 'dB'
            except:
              burst = ''

        try:
          bw = pmt.to_double(pmt.dict_ref(meta, pmt.intern('symbol_rate'), pmt.PMT_NIL))
        except:
          pass
        try:
          bw = pmt.to_double(pmt.dict_ref(meta, pmt.intern('bandwidth'), pmt.PMT_NIL))
          if bw < self.bw_min:
            bw = self.bw_min
        except:
          pass
        self.d_dict['annotations'].append({'core:sample_start': sob-self.soo,
                    'core:sample_count': eob-sob, 'core:freq_upper_edge': int(freq+bw/2),
                    'core:freq_lower_edge': int(freq-bw/2), 'core:description': burst})
      except Exception as e:
        print('could not form annotation from message', pmt.car(msg), ':',e)
