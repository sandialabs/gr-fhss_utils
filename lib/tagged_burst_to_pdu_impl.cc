/* -*- c++ -*- */
/*
 * Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "tagged_burst_to_pdu_impl.h"

#include <unistd.h>
#include <inttypes.h>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/adaptor/map.hpp>
#include <omp.h>
#include <cstdlib>

namespace gr {
  namespace fhss_utils {

    tagged_burst_to_pdu::sptr
    tagged_burst_to_pdu::make(size_t decimation, const std::vector<float>& taps, float min_burst_time, float max_burst_time,
        float relative_center_frequency, float relative_span, float relative_sample_rate, float sample_rate, int num_threads)
    {
      return gnuradio::get_initial_sptr
        (new tagged_burst_to_pdu_impl(decimation, taps, min_burst_time, max_burst_time, relative_center_frequency,
            relative_span, relative_sample_rate, sample_rate, num_threads));
    }

    /*
     * The private constructor
     */
    tagged_burst_to_pdu_impl::tagged_burst_to_pdu_impl(size_t decimation, const std::vector<float>& taps, float min_burst_time, float max_burst_time, 
        float relative_center_frequency, float relative_span, float relative_sample_rate, float sample_rate, int num_threads)
      : gr::sync_block("tagged_burst_to_pdu",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0)),
          d_debug(false),
          d_relative_center_frequency(relative_center_frequency),     // May be adjusted if we are after a filter bank
          d_relative_sample_rate(relative_sample_rate),
          d_min_burst_size(min_burst_time*sample_rate/decimation),
          d_max_burst_size(max_burst_time*sample_rate/decimation),
          d_outstanding(0),
          d_max_outstanding(0),
          d_n_dropped_bursts(0),
          d_blocked(false),
          d_decimation(decimation),
          d_taps(taps),
          d_input_fir(0, taps),
          d_write_queue(6),
          d_work_queue(6),
          d_sample_rate(sample_rate),
          d_num_threads(num_threads)
    {
      d_lower_border = relative_center_frequency - relative_span / 2;
      d_upper_border = relative_center_frequency + relative_span / 2;
      message_port_register_out(PDU_OUT_PORT);

      d_current_rx_time_tag.offset = 0;
      d_current_rx_time_tag.key = pmt::mp("rx_time");
      d_current_rx_time_tag.value = pmt::make_tuple(pmt::from_uint64(0), pmt::from_double(0));
      d_rx_time_tags.push(d_current_rx_time_tag);
      d_block_increment = ((d_block_size - taps.size() + 1) / d_decimation) * d_decimation;

      process_thread = new boost::thread(boost::bind(&tagged_burst_to_pdu_impl::process_data, this));

      // Create work buffers
      for (size_t i = 0; i < 5; i++) {
        d_write_queue.bounded_push(new buffer(d_block_size));
      }

      d_current_buffer = NULL;

      // Because we are downsampling and filtering in blocks, the math is kind of funky for ensuring
      // that we account for filter tails.  In downsampling, we skip over D input points after every
      // calculation.  The filter function always starts at the first input point, so we need to
      // throw away some filter tail points to ensure that.
      printf("filter len = %zu\n", d_taps.size());

      //message_port_register_in(pmt::mp("burst_handled"));
      //set_msg_handler(pmt::mp("burst_handled"), boost::bind(&tagged_burst_to_pdu_impl::burst_handled, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    tagged_burst_to_pdu_impl::~tagged_burst_to_pdu_impl()
    {
      // Delete all of the bursts that are still around
      while(d_alloced_arrays.size() > 0) {
        two_gr_complex burst = d_alloced_arrays.front();
        d_alloced_arrays.pop();
        free(burst.one);
        free(burst.two - ROT_TMP_OFFSET);
      }
      for (auto burst : d_bursts) {
        free(burst.second.data);
        free(burst.second.rot_tmp - ROT_TMP_OFFSET);
      }

    }

    bool tagged_burst_to_pdu_impl::stop() {
      buffer* end_buffer;
      d_write_queue.pop(end_buffer);
      end_buffer->end_flag = true;
      d_work_queue.bounded_push(end_buffer);
      process_thread->join();
      delete process_thread;

      buffer* b;
      while ( d_write_queue.pop(b) ) {
        delete b;
      }
      while ( d_work_queue.pop(b) ) {
        delete b;
      }
    }


    int
    tagged_burst_to_pdu_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      size_t data_index = 0;

      // We need to fill up buffers, with overlap.
      size_t i = 0;
      for (size_t j = d_block_size; j <= noutput_items; i += d_block_increment, j = i + d_block_size) {
        // request new buffer if available
        if (!d_write_queue.pop(d_current_buffer)) {
          //printf("************We have to throw away data.  Need more threads. Bursts = %zu, data = %d\n", d_bursts.size(), noutput_items);
          //usleep(int(1e6*noutput_items/d_sample_rate/4));
          usleep(500);
          return i;
        }

        // We have to have a valid buffer at this point
        // get all the tags in this buffer, but not in the overlap region at the end
        size_t tagend = j - d_block_size + d_block_increment;
        std::vector<tag_t> new_burst_tags;
        get_tags_in_window(new_burst_tags, 0, i, tagend, NEW_BURST_TAG);
        std::vector<tag_t> rx_time_tags;
        get_tags_in_window(rx_time_tags, 0, i, tagend, RX_TIME_TAG);
        std::vector<tag_t> gone_burst_tags;
        get_tags_in_window(gone_burst_tags, 0, i, tagend, GONE_BURST_TAG);

        // Copy data and tags to buffer
        d_current_buffer->add_data(in+i, d_block_size);
        d_current_buffer->add_tags(new_burst_tags, rx_time_tags, gone_burst_tags);
        d_current_buffer->start = nitems_read(0) + i;
        d_work_queue.bounded_push(d_current_buffer);
        d_current_buffer = NULL;

      }
      if (i == 0) {
        usleep(500);
      }

      return i;
    }

    double tagged_burst_to_pdu_impl::convert_rx_time(const tag_t &rx_time_tag) {
      uint64_t seconds = pmt::to_uint64(pmt::tuple_ref(rx_time_tag.value, 0));
      double frac_seconds = pmt::to_double(pmt::tuple_ref(rx_time_tag.value, 1));
      return (double)seconds + frac_seconds;
    }

    void
    tagged_burst_to_pdu_impl::create_new_bursts(const buffer& work_buffer)
    {
      for (const tag_t &tag : work_buffer.rx_time_tags) {
        //printf("tagged_burst_to_pdu: received rx_time %f at offset %lu\n", convert_rx_time(tag), tag.offset);
        d_rx_time_tags.push(tag);
      }

      for (size_t i = 0; i < work_buffer.new_burst_tags.size(); i++) {
        tag_t tag = work_buffer.new_burst_tags[i];
        float relative_frequency = pmt::to_float(pmt::dict_ref(tag.value, META_RELATIVE_FREQ, pmt::PMT_NIL));

        if(d_lower_border < relative_frequency && relative_frequency <= d_upper_border) {
          uint64_t id = pmt::to_uint64(pmt::dict_ref(tag.value, META_ID, pmt::PMT_NIL));
          float magnitude = pmt::to_float(pmt::dict_ref(tag.value, META_MAGNITUDE, pmt::PMT_NIL));
          float center_frequency = pmt::to_float(pmt::dict_ref(tag.value, META_CENTER_FREQ, pmt::PMT_NIL));
          float sample_rate = pmt::to_float(pmt::dict_ref(tag.value, META_SAMP_RATE, pmt::PMT_NIL));

          // if we've already seen a burst with this id, throw this one away before we malloc any more memory
          if (d_bursts.count(id)){
            printf("ERROR: tagged_burst_to_pdu saw a repeated burst id (%zu), discarding...\n", id);
            break; // break to prevent memory leak
          }

          // Adjust the values based on our position behind a potential filter bank
          center_frequency += (d_relative_center_frequency + relative_frequency) * sample_rate;
          sample_rate = sample_rate * d_relative_sample_rate;
          relative_frequency = (relative_frequency - d_relative_center_frequency) / d_relative_sample_rate;
          //relative_frequency = relative_frequency * sample_rate + center_frequency;


          // pop from the d_rx_time_tags queue until d_current_rx_time_tag has the last rx_time tag received before the current new_burst tag
          while (!d_rx_time_tags.empty() && tag.offset >= d_rx_time_tags.front().offset) {
            d_current_rx_time_tag = d_rx_time_tags.front();
            d_rx_time_tags.pop();
          }
          
          uint64_t start_offset = tag.offset;
          // Adjust tag.offset by decimation, so that we get the same answer independant of block size
          start_offset = (tag.offset / d_decimation) * d_decimation;
          if (id == 1) {
            printf("start = %zu %zu\n", start_offset, tag.offset);
          }
          // It is possible (but unlikely) that the truncation will put us behind the start of the buffer.
          // In that case increment by d_decimation
          if (start_offset < work_buffer.start) start_offset += d_decimation;


          // calculate the start_time
          double start_time = convert_rx_time(d_current_rx_time_tag) + (start_offset - d_current_rx_time_tag.offset) / sample_rate;
          sample_rate /= d_decimation;
          //printf("tagged_burst_to_pdu: burst start time is %f (using rx_time %f)\n", start_time, convert_rx_time(d_current_rx_time_tag));

          tag.value = pmt::dict_delete(tag.value, META_CENTER_FREQ);
          tag.value = pmt::dict_delete(tag.value, META_SAMP_RATE);
          tag.value = pmt::dict_delete(tag.value, META_RELATIVE_FREQ);
          tag.value = pmt::dict_add(tag.value, META_CENTER_FREQ, pmt::from_float(center_frequency));
          tag.value = pmt::dict_add(tag.value, META_SAMP_RATE, pmt::from_float(sample_rate));
          tag.value = pmt::dict_add(tag.value, META_RELATIVE_FREQ, pmt::from_float(relative_frequency * sample_rate * d_decimation));
          tag.value = pmt::dict_add(tag.value, META_START_TIME, pmt::from_double(start_time));


          burst_data burst = {id, start_offset, start_offset - work_buffer.start, magnitude, relative_frequency,
            center_frequency, sample_rate, 0, 0, tag.value};
          // Get a burst if available otherwise alloc a new one
          if (!d_alloced_arrays.empty()) {
            two_gr_complex arrs = d_alloced_arrays.front();
            d_alloced_arrays.pop();
            burst.data = arrs.one;
            burst.rot_tmp = arrs.two;
          } else {
            burst.data = (gr_complex *) malloc(sizeof(gr_complex) * (d_max_burst_size + d_block_size));
            burst.rot_tmp = (gr_complex *) malloc(sizeof(gr_complex) * (d_block_size+8));
            memset(burst.rot_tmp + d_block_size, 0, 4*sizeof(gr_complex));
            memset(burst.rot_tmp, 0, 4*sizeof(gr_complex));
            burst.rot_tmp += ROT_TMP_OFFSET;
          }
          float phase_inc = 2 * M_PI * -relative_frequency;
          burst.rotate.set_phase_incr(exp(gr_complex(0, phase_inc)));
          burst.rotate.set_phase(gr_complex(1, 0));
          burst.rot_skip = 0;

          if(burst.data != NULL) {
            d_bursts[id] = burst;
            d_max_id = id;
            //append_to_burst(d_new_bursts[id], &in[relative_offset], to_copy);
            if(d_debug) {
              printf("New burst: offset=%" PRIu64 ", id=%" PRIu64 ", relative_frequency=%f, magnitude=%f\n", tag.offset, id, relative_frequency, magnitude);
            }
          } else {
            printf("Error, malloc failed\n");
          }
        }
      }
    }

    void
    tagged_burst_to_pdu_impl::publish_and_remove_old_bursts(const buffer& work_buffer)
    {
      for(tag_t tag : work_buffer.gone_burst_tags) {
        uint64_t id = pmt::to_uint64(pmt::dict_ref(tag.value, META_ID, pmt::PMT_NIL));

        if(d_bursts.count(id)) {
          burst_data &burst = d_bursts[id];
          if(d_debug) {
            printf("gone burst: %" PRIu64 " %zu\n", id, burst.len);
          }
          // Subtract off any samples from the end
          size_t offset = std::ceil((d_block_size - tag.offset + work_buffer.start) / (float)d_decimation);
          burst.len -= burst.len > offset ? offset : burst.len;
          if(burst.len >= d_min_burst_size) {
            if (id == 1) {
                printf("id %lu: len = %zu %zu\n", id, burst.len, tag.offset);
            }
            publish_burst(burst);
          }
          d_alloced_arrays.push(two_gr_complex(d_bursts[id].data, d_bursts[id].rot_tmp));
          d_bursts.erase(id);
        }
      }

      // Go through any remaining bursts and remove them if their buffer is full
      std::vector<uint64_t> removeNodes;
      for (std::pair<const uint64_t, burst_data>& burstpair : d_bursts) {
        burst_data& burst = burstpair.second;
        if (burst.len >= d_max_burst_size) {
          uint64_t id = burst.id;
          if(d_debug) {
            printf("gone burst: %" PRIu64 " %zu\n", id, burst.len);
          }
          if(burst.len >= d_min_burst_size) {
            publish_burst(burst);
          }
          printf("Long Burst: %zu, len = %zu\n", burst.id, burst.len);
          d_alloced_arrays.push(two_gr_complex(d_bursts[id].data, d_bursts[id].rot_tmp));
          removeNodes.push_back(id);
        }
      }

      for (uint64_t id : removeNodes)
        d_bursts.erase(id);

    }

    void
    tagged_burst_to_pdu_impl::publish_burst(burst_data &burst)
    {
      pmt::pmt_t d_pdu_meta = burst.dict;
      pmt::pmt_t d_pdu_vector = pmt::init_c32vector(burst.len, burst.data);

      d_pdu_meta = pmt::dict_add(d_pdu_meta, META_DURATION, pmt::from_float(burst.len / burst.sample_rate));

      pmt::pmt_t msg = pmt::cons(d_pdu_meta,
          d_pdu_vector);

      message_port_pub(PDU_OUT_PORT, msg);
    }

    void tagged_burst_to_pdu_impl::process_data() {
      buffer* work_buffer = NULL;
      d_max_id = 0;
      uint64_t max_cbursts = 0;
      // Set the number of threads we have to work with
      printf("d_num_threads = %d\n", d_num_threads);
      omp_set_num_threads(d_num_threads);

      while(true) {
        if (d_work_queue.pop(work_buffer)) {
          if(work_buffer->end_flag) break;
          // Handle all of the new tags
          create_new_bursts(*work_buffer);

          // tune, filter, decimate for each burst
          // We should be able to just do openmp here.
          //printf("Processing %zu bursts, %p, %zu\n", d_bursts.size(), work_buffer, work_buffer->data.size());
          size_t block_size =  (work_buffer->data.size() - d_taps.size() + 1) / d_decimation;
          if (d_bursts.size() > max_cbursts) {
            printf("New max bursts from %zu to %zu\n", max_cbursts, d_bursts.size());
            max_cbursts = d_bursts.size();
          }
          if (d_bursts.size() > 0) {
            std::vector<int> keys;
            keys.reserve(d_bursts.size());
            boost::copy(d_bursts | boost::adaptors::map_keys, std::back_inserter(keys));
            size_t input_size = work_buffer->data.size();
            gr_complex* input_data = &work_buffer->data[0];
            size_t useSize = std::min(d_bursts.size(), (size_t)100);
            #pragma omp parallel for
            for (size_t i = 0; i < useSize; i++) {
              burst_data& burst = d_bursts[keys[i]];
              size_t data_skip = burst.data_skip;
              size_t rot_skip = burst.rot_skip;
              size_t block_offset = burst.data_skip / d_decimation;
              size_t c_size = 0;
              if (block_size > block_offset) c_size = block_size - block_offset;
              burst.rotate.rotateN(burst.rot_tmp + rot_skip, input_data + rot_skip, input_size - rot_skip);
              d_input_fir.filterNdec(burst.data + burst.len, burst.rot_tmp + data_skip, c_size, d_decimation);
              burst.rot_skip = d_block_size - d_block_increment;
              rot_skip = burst.rot_skip;
              memcpy(burst.rot_tmp, burst.rot_tmp + input_size - rot_skip, rot_skip * sizeof(burst.rot_tmp[0]));
              burst.len += c_size;
              burst.data_skip = 0;
            }
          }
          // publish messages for any finished bursts
          publish_and_remove_old_bursts(*work_buffer);

          // return the buffer, so that it can be writen to again
          work_buffer->reset();
          d_write_queue.bounded_push(work_buffer);
          //printf("Releasing Buffer\n");
        } else {
          usleep(100);
        }
      }
      printf("Saw %zu bursts\n", d_max_id);
    }

  }

}
