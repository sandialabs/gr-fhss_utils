#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2020 gr-fhss_utils38 author.
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
from gnuradio import gr
import math
import numpy
import io
import pmt
from PIL import Image, ImageDraw
import time
import threading
from gnuradio import fhss_utils


class burst_tag_debug(gr.sync_block):
    """
      Build an image from fft pdus and burst metadata pdus
      - Three burst metadata ports are available for burst rectangle coloring
      - Priority is Red < Green < Blue, so Blue will overwrite a Red rect of the same burst ID
    """

    def __init__(self, filename="/tmp/btd", nrows=5000, fft_size=1024, nrows_avg=10, sample_rate=1e6, quality=75, always_on=True):
        gr.sync_block.__init__(self, name="image_raster_msg_sink", in_sig=None, out_sig=None)

        # filename should not end in '.png', will be appended. Set to "" to disable
        self.filename = filename
        if len(filename) == 0:
            self.filename = None

        # rows determine how many pdu's are stiched together to form an image
        self.nrows = nrows
        self.nrows_avg = nrows_avg

        # fft_size for pdu checking as well as overlaying bursts
        self.fft_size = fft_size

        # sample rate used for burst box drawing as well as axis (to be implemented)
        self.sample_rate = sample_rate

        # quality for jpeg conversion
        self.quality = quality

        # flag to indicate the need to publish another PDU with an image
        self.need_to_publish = True
        self.fft_ready = False
        self.always_on = always_on

        # non-jpeg image storage in a list
        self.row_idx = 0
        self.image_data = []
        self.absolute_image_count = 0
        self.starting_fft_sample = None

        self.publish_timer = None

        # burst detections stored in a list
        self.bursts = {}
        self.burst_image = []

        # common use PDUs
        self.pmt_pdu_out = pmt.intern("pdu_out")
        self.pmt_type = pmt.intern("type")
        self.pmt_detect = pmt.intern("detect")
        self.pmt_detect_image = pmt.intern("detect_image")
        self.pmt_spectral_image = pmt.intern("spectral_image")
        self.pmt_x_length = pmt.intern("x_length")
        self.pmt_y_length = pmt.intern("y_length")

        self.pmt_start_offset = fhss_utils.PMTCONSTSTR__start_offset()
        self.pmt_end_offset = fhss_utils.PMTCONSTSTR__end_offset()
        self.pmt_burst_id = fhss_utils.PMTCONSTSTR__burst_id()
        self.pmt_relative_frequency = fhss_utils.PMTCONSTSTR__relative_frequency()
        self.pmt_bandwidth = fhss_utils.PMTCONSTSTR__bandwidth()

        # register message handlers
        self.message_port_register_out(self.pmt_pdu_out)
        self.message_port_register_in(pmt.intern("fft_in"))
        self.set_msg_handler(pmt.intern("fft_in"), self.fft_msg_handler)
        # priority: red < green < blue
        self.message_port_register_in(pmt.intern("burst_in_r"))
        self.set_msg_handler(pmt.intern("burst_in_r"), self.burst_handler_r)
        self.message_port_register_in(pmt.intern("burst_in_g"))
        self.set_msg_handler(pmt.intern("burst_in_g"), self.burst_handler_g)
        self.message_port_register_in(pmt.intern("burst_in_b"))
        self.set_msg_handler(pmt.intern("burst_in_b"), self.burst_handler_b)

    #############################################
    # Control / Public Interface
    #############################################
    def stop(self):
        self.need_to_publish = False
        if self.publish_timer is not None:
            self.publish_timer.cancel()
        self.publish_timer = None
        return True

    def reset(self):
        self.image_data = []
        self.row_idx = 0
        self.starting_fft_sample = None
        if not self.always_on:
            self.need_to_publish = False
        self.bursts = {}
        self.fft_ready = False
        if self.publish_timer is not None:
            self.publish_timer.cancel()
        self.publish_timer == None

    def request_image(self):
        print("jpeg_block going to generate another spectral image")
        self.reset()
        self.need_to_publish = True

    #############################################
    # Message Handlers
    #############################################
    def fft_msg_handler(self, pdu):
        if not self.need_to_publish or self.fft_ready:
            # drop messages if we're not currently building an image
            return

        # if we are building an image, read the pdu f32 fft data
        try:
            meta = pmt.car(pdu)
            data = pmt.f32vector_elements(pmt.cdr(pdu))
        except Exception as e:
            # should we reset the image data here?
            print(f"exception in burst_tag_debug, {e}")
            return

        # check that this vector is the same size as the rest
        if not len(data) == self.fft_size:
            print(f"different length vector received in burst_tag_debug block {len(data)}, resetting")
            self.image_data = []
            self.row_idx = 0
            return

        # save off the starting sample number of this image
        if self.starting_fft_sample is None:
            self.starting_fft_sample = pmt.to_uint64(pmt.dict_ref(meta, self.pmt_start_offset, pmt.PMT_NIL))

        # looks good, add on the data to the image
        self.image_data.append(data)
        self.row_idx += 1

        # publish an event followed by an image once it is complete
        if self.row_idx >= self.nrows:
            self.fft_ready = True
            # wait N more seconds for additional burst metadata to make it in
            N = 1.0
            self.publish_timer = threading.Timer(N, self.publish_result)
            self.publish_timer.start()

    def burst_handler_r(self, pdu):
        return self.burst_handler_rgb(pdu, rgb="r")

    def burst_handler_g(self, pdu):
        return self.burst_handler_rgb(pdu, rgb="g")

    def burst_handler_b(self, pdu):
        return self.burst_handler_rgb(pdu, rgb="b")

    def burst_handler_rgb(self, pdu, rgb="r"):
        if not self.need_to_publish:
            return

        # generate burst dict if metadata is present
        burst_id, burst_dict = self.dict_from_pdu(pdu)
        if burst_id is None:
            return
        burst_dict["color"] = rgb

        # priority: red < green < blue
        if burst_id not in self.bursts and rgb == "r":
            self.bursts[burst_id] = burst_dict # red lowest priority, write if not exists
        elif burst_id in self.bursts and rgb == "b":
            self.bursts[burst_id] = burst_dict # blue highest priority, overwrite always
        elif burst_id in self.bursts and rgb == "g" and self.bursts[burst_id]["color"] == "r":
            self.bursts[burst_id] = burst_dict # green overwrite red only

    #############################################
    # Message Publishers
    #############################################
    def publish_result(self):
        meta = pmt.make_dict()
        meta = pmt.dict_add(meta, self.pmt_spectral_image, pmt.intern(f"spectral_image_{self.absolute_image_count}"))
        self.absolute_image_count += 1

        # publish event without the image so we know it's coming
        meta = pmt.dict_add(meta, self.pmt_type, self.pmt_detect)
        self.message_port_pub(self.pmt_pdu_out, pmt.cons(meta, pmt.init_u8vector(0, [])))

        # publish the event with the image
        meta = pmt.dict_add(meta, self.pmt_type, self.pmt_detect_image)
        nrows, ncols, b = self.build_jpeg(self.image_data, self.bursts)
        meta = pmt.dict_add(meta, self.pmt_x_length, pmt.from_uint64(nrows))
        meta = pmt.dict_add(meta, self.pmt_y_length, pmt.from_uint64(ncols))
        print(f"burst_tag_debug going to send a jpeg: #{self.absolute_image_count}")
        b = list(b) # bytes type no longer accepted into pmt vectors
        self.message_port_pub(self.pmt_pdu_out, pmt.cons(meta, pmt.init_u8vector(len(b), b)))

        # reset flags
        self.reset()

    #############################################
    # Helper functions
    #############################################
    def dict_from_pdu(self, pdu):
        # build a simple dictionary out of burst metadata
        meta = pmt.car(pdu)
        burst_dict = {}
        burst_id = -1
        try:
            burst_id = pmt.to_uint64(pmt.dict_ref(meta, self.pmt_burst_id, pmt.PMT_NIL))
            burst_dict["start"] = pmt.to_uint64(pmt.dict_ref(meta, self.pmt_start_offset, pmt.PMT_NIL))
            burst_dict["end"] = pmt.to_uint64(pmt.dict_ref(meta, self.pmt_end_offset, pmt.PMT_NIL))
            burst_dict["rel_cf"] = pmt.to_float(pmt.dict_ref(meta, self.pmt_relative_frequency, pmt.PMT_NIL))
            burst_dict["bw"] = pmt.to_float(pmt.dict_ref(meta, self.pmt_bandwidth, pmt.PMT_NIL))
        except Exception as e:
            print(f"malformed burst (red) in the jpeg_convertor, {e}")
            return None, {}
        return burst_id, burst_dict

    def build_jpeg(self, data_2d, burst_dicts):
        min_scale = -110
        max_scale = 0
        scaled = numpy.clip(data_2d, min_scale, max_scale)
        scaled = scaled - scaled.min()
        scaled = scaled / (max_scale - min_scale) * 255

        # get the burst rectangle coordinates as list of tuples
        burst_rects = self.build_boxes(burst_dicts)
        #print(f"b, {burst_rects[0][1]}")

        # generate image as a jpeg
        b = io.BytesIO()
        im = Image.fromarray(numpy.round(scaled).astype('uint8'), mode="L")

        # convert to color
        im = im.convert("RGB")

        # get a drawing context
        d = ImageDraw.Draw(im)
        for br in burst_rects:
            c_hex = br[3]
            d.rectangle(br[1], outline=c_hex, width=2)
            d.text((min(br[1][2] + 5, self.fft_size), br[1][3]), f"{br[0]}", fill=0x99ff33)

        # send image as jpeg
        start_time = time.time()
        im.save(b, format='JPEG', quality=self.quality)
        print(f"total time to convert to jpeg: {time.time()-start_time:02f}")
        # save local copy as png if filename not None
        if self.filename is not None:
            print(f"jpeg_block going to save an image: {self.filename}-{self.absolute_image_count}.png", flush=True)
            start_time = time.time()
            im.save(f"{self.filename}-{self.absolute_image_count}.png", format="PNG")
            print(f"total time to save to png: {time.time()-start_time:02f}")
        return (im.size[0], im.size[1], b.getvalue())

    def get_color_hex(self, rgb):
        if rgb == "r":
            return 0x3333ff
        elif rgb == "g":
            return 0x99ff33
        elif rgb == "b":
            return 0xffa31a
        else:
            return 0xf2f2f2

    def build_boxes(self, bursts_dict):
        # convert gnuradio abs indexes to fft image indexes
        #  left, bottom, right and top coordinate of a bounding box
        samples_per_image = self.fft_size * self.nrows_avg * self.nrows
        samples_per_row = self.fft_size * self.nrows_avg
        time_per_row = samples_per_row / self.sample_rate
        n_rows_per_image = self.nrows
        print("samples per image is", samples_per_image)

        # parameters that used to be settable in the jupyter notebook
        use_time_axis = False
        use_Hz_axis = False
        iq_data_center_freq = 0

        # loop through each burst and build the bounding rectangles
        burst_rects = []
        for burst_id, burst in bursts_dict.items():
            color_hex = self.get_color_hex(burst["color"])

            # which image are we on?
            image_idx = int(burst["start"] // samples_per_image)
            next_image_idx = int(burst["end"] // samples_per_image)

            # if this is the first burst of the image, start a new burst list
            while len(burst_rects) - 1 < image_idx:
                burst_rects.append([])

            # don't care which image (mod) but we do care what row (floor div)
            if use_time_axis:
                top = time_per_row * ((burst["start"] % (samples_per_image)) // samples_per_row)
                bottom = time_per_row * ((burst["end"] % (samples_per_image)) // samples_per_row)
            else:
                top = int(((burst["start"] % (samples_per_image)) // samples_per_row))
                bottom = int(((burst["end"] % (samples_per_image)) // samples_per_row))

            # some bursts span multiple images, truncate this one so it still plots
            if next_image_idx != image_idx:
                bottom = n_rows_per_image

            # change Hz width to bin width is flag is not set
            if use_Hz_axis:
                width = burst["bw"]
            else:
                width = burst["bw"] / (self.sample_rate / self.fft_size)

            # scale width HARDCODE
            width *= 4

            # center is based on the relative_frequency
            if use_Hz_axis:
                center = iq_data_center_freq + burst["rel_cf"]
            else:
                center = self.fft_size // 2 + (burst["rel_cf"] / self.sample_rate) * self.fft_size
            right = math.ceil(center + width / 2)
            left = math.floor(center - width / 2)

            # add this rectangle to the current image
            burst_rects[image_idx].append((burst_id, (left, bottom, right, top), image_idx, color_hex))

            # add it to the next_image_idx picture too if it overlaps
            if next_image_idx != image_idx:
                while len(burst_rects) - 1 < next_image_idx:
                    burst_rects.append([])
                top = 0
                if use_time_axis:
                    bottom = time_per_row * ((burst["end"] % (samples_per_image)) // samples_per_row)
                else:
                    bottom = int(((burst["end"] % (samples_per_image)) // samples_per_row))
                burst_rects[next_image_idx].append((burst_id, (left, bottom, right, top), next_image_idx, color_hex))

        try:
            #print("burst_rects is", burst_rects[0])
            print("burst_rects is len", len(burst_rects[0]))
        except IndexError:
            print("Warning: no burst squares were made", burst_rects, bursts_dict)

        # return for the first image only - this is a one-shot imaging tool
        if len(burst_rects):
            return burst_rects[0]
        else:
            return []

    #############################################
    # Misc
    #############################################
    def work(self, input_items, output_items):
        return len(output_items[0])
