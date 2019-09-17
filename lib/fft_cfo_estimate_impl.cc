/* -*- c++ -*- */
/*
 * Copyright 2018 gr-fhss_utils author.
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
#include <volk/volk.h>
#include "fft_cfo_estimate_impl.h"

namespace gr {
  namespace fhss_utils {

    fft_cfo_estimate::sptr
    fft_cfo_estimate::make(int fftsize, int offset, int num_ffts, float threshold, std::vector<float> valid_freqs)
    {
      return gnuradio::get_initial_sptr
        (new fft_cfo_estimate_impl(fftsize, offset, num_ffts, threshold, valid_freqs));
    }

    /*
     * The private constructor
     */
    fft_cfo_estimate_impl::fft_cfo_estimate_impl(int fftsize, int offset, int num_ffts, float threshold, std::vector<float> valid_freqs)
      : gr::block("fft_cfo_estimate",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
      d_fftsize(fftsize),
      d_offset(offset),
      d_nffts(num_ffts),
      d_fft(d_fftsize, true, 1)
    {
      set_freqs(valid_freqs);
      d_threshold = pow(10, threshold / 10.0f);

      d_mag     = (float*) volk_malloc(d_fftsize * sizeof(float), volk_get_alignment());
      d_mag_sum = (float*) volk_malloc(d_fftsize * sizeof(float), volk_get_alignment());

      message_port_register_in(pmt::mp("pdu_in"));
      message_port_register_out(pmt::mp("pdu_out"));
      set_msg_handler(pmt::mp("pdu_in"), boost::bind(&fft_cfo_estimate_impl::pdu_handler, this, _1));
    }

    /*
     * Our virtual destructor.
     */
    fft_cfo_estimate_impl::~fft_cfo_estimate_impl()
    {
      volk_free(d_mag);
      volk_free(d_mag_sum);
    }


    void fft_cfo_estimate_impl::set_freqs(const std::vector<float> freqs) {
      gr::thread::scoped_lock l(d_setlock);

      d_valid_freqs = freqs;
    }


    double fft_cfo_estimate_impl::get_closest_freq(double freq) {
      gr::thread::scoped_lock l(d_setlock);

      if (d_valid_freqs.size() > 0) {
        double best_freq = d_valid_freqs[0];
        double min_diff = std::abs(freq - d_valid_freqs[0]);

        for (int ii=1; ii<d_valid_freqs.size(); ii++) {
          double diff = std::abs(d_valid_freqs[ii] - freq);
          if (diff < min_diff) {
            min_diff = diff;
            best_freq = d_valid_freqs[ii];
          }
        }
        return best_freq;
      } else {
        return freq;
      }
    }

    float fft_cfo_estimate_impl::calc_center_freq(float sample_rate) {
      float power = 0;
      for (size_t i = 0; i < d_fftsize; i++) {
        power += d_mag_sum[i];
      }

      // 99% occupied bandwidth
      float bw99_thresh = power * .005;
      size_t minIndex = 0;
      size_t maxIndex = d_fftsize - 1;
      float cSum = 0;
      while(cSum < bw99_thresh && minIndex < d_fftsize) {
        cSum += d_mag_sum[minIndex];
        minIndex++;
      }
      minIndex--;
      cSum = 0;
      while(cSum < bw99_thresh && maxIndex > 0) {
        cSum += d_mag_sum[maxIndex];
        maxIndex--;
      }
      maxIndex++;

      // center freq
      // This is very sensitive to the SNR, but gives us a good first pass estimate.
      float fc1 = (float)(maxIndex + minIndex) / (2*d_fftsize) * sample_rate - sample_rate/2;

      float occPower = power * .99;
      occPower /= (maxIndex - minIndex + 1);
      // We have magnitude squared, so three db is 1/4
      float thresh = .25*occPower;
      size_t minIndex3db = 1;
      size_t maxIndex3db = d_fftsize - 2;
      while (d_mag_sum[minIndex3db] < thresh) minIndex3db++;
      minIndex3db--;
      while (d_mag_sum[maxIndex3db] < thresh) maxIndex3db--;
      maxIndex3db++;
      //printf("3db index: %zu, %zu\n", minIndex3db, maxIndex3db);

      float fc3db = (float)(maxIndex3db + minIndex3db) / (2*d_fftsize) * sample_rate - sample_rate/2;
      return fc3db;
    }

    float fft_cfo_estimate_impl::calc_bandwidth(float& power, float sample_rate) {
      // create a list of bins that pass the threshold
      float* max_ptr = std::max_element(d_mag_sum, d_mag_sum + d_fftsize);
      float max_magnitude = *max_ptr * 0.5;
      if (max_ptr != d_mag_sum) max_magnitude += .25*max_ptr[-1];
      else max_magnitude += .25*max_ptr[0];
      if (max_ptr != d_mag_sum + d_fftsize) max_magnitude += .25*max_ptr[1];
      else max_magnitude += .25*max_ptr[0];

      float threshold = max_magnitude / d_threshold;

      std::vector<int> peak_bins;
      peak_bins.reserve(d_fftsize/2);

      for (int bin = 0; bin < d_fftsize; bin++) {
        if (d_mag_sum[bin] > threshold) {
          peak_bins.push_back(bin);
        }
      }

      float median_bin;
      if (peak_bins.size() % 2 == 0) {
        median_bin = 0.5f * (peak_bins[peak_bins.size()/2 - 1] + peak_bins[peak_bins.size()/2]);
      } else {
        median_bin = peak_bins[peak_bins.size()/2];
      }

       // Search for the bandwidth.  It is roughly the first sorted bin - the last sorted bin, but we may have spurs, so look for continious bins
      size_t start_bin = peak_bins.size()/2;
      size_t stop_bin = peak_bins.size()/2;
      size_t bin_thresh = d_fftsize / 100;
      for (; stop_bin < peak_bins.size(); stop_bin++) {
        if (peak_bins[stop_bin] - peak_bins[stop_bin-1] > bin_thresh) {
          break;
        }
      }
      stop_bin = peak_bins[stop_bin-1] + 1;
      for (; start_bin > 0; start_bin--) {
        if (peak_bins[start_bin] - peak_bins[start_bin-1] > bin_thresh) {
          break;
        }
      }
      if (start_bin == 0) start_bin = peak_bins[0] - 1;
      else start_bin = peak_bins[start_bin-1] - 1;
      float bandwidth = (float)(stop_bin - start_bin + 1) / d_fftsize * sample_rate;
      power = 0;
      for (size_t i = start_bin; i <= stop_bin; i++) power += d_mag_sum[i];
      power /= (stop_bin + start_bin + 1);
      //printf("bins = %zu, %zu\n", start_bin, stop_bin);
      return bandwidth;


    }

    void fft_cfo_estimate_impl::pdu_handler(pmt::pmt_t pdu) {
      if (!pmt::is_pair(pdu)) {
        GR_LOG_WARN(d_logger, "PDU not a pair, dropping");
        return;
      }

      pmt::pmt_t metadata = pmt::car(pdu);
      pmt::pmt_t data = pmt::cdr(pdu);

      if (!pmt::is_dict(metadata)) {
        GR_LOG_WARN(d_logger, "PDU metadata not a dictionary, dropping");
        return;
      }

      float center_frequency = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("center_frequency"), pmt::PMT_NIL));
      float relative_frequency = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("relative_frequency"), pmt::PMT_NIL));
      float sample_rate = pmt::to_float(pmt::dict_ref(metadata, pmt::mp("sample_rate"), pmt::PMT_NIL));
      uint64_t bid = pmt::to_uint64(pmt::dict_ref(metadata, pmt::mp("burst_id"), pmt::PMT_NIL));

      if (!pmt::is_c32vector(data)) {
        GR_LOG_WARN(d_logger, "PDU data not complex, dropping");
        return;
      }

      size_t burst_size = pmt::length(data);
      const gr_complex* burst = (const gr_complex*) pmt::c32vector_elements(data, burst_size);

      int nffts = std::min(d_nffts, int((burst_size - d_offset)/d_fftsize));

      // set the sum to 0
      memset(d_mag_sum, 0, sizeof(float)*d_fftsize);
      size_t fftd2 = d_fftsize / 2;

      for (int fft_idx = 0; fft_idx < nffts; fft_idx++) {
        // take the next FFT
        memcpy(d_fft.get_inbuf(),                      // dest
               burst + d_offset + fft_idx * d_fftsize, // src
               sizeof(gr_complex) * d_fftsize);        // size
        d_fft.execute();

        // calculate magnitudes

        volk_32fc_magnitude_squared_32f(d_mag + fftd2, d_fft.get_outbuf(), fftd2);
        volk_32fc_magnitude_squared_32f(d_mag, d_fft.get_outbuf() + fftd2, fftd2);

        // add to the sum
        volk_32f_x2_add_32f(d_mag_sum, d_mag_sum, d_mag, d_fftsize);
      }
      float fc3db = calc_center_freq(sample_rate);
      float power;
      float bandwidth = calc_bandwidth(power, sample_rate);

      if (d_valid_freqs.size() > 0) {
        fc3db = get_closest_freq(fc3db + center_frequency) - center_frequency;
      }
      d_output.resize(burst_size);

      d_r.set_phase_incr(exp(gr_complex(0, 2 * M_PI * (-1 * fc3db/ sample_rate))));
      d_r.set_phase(gr_complex(1, 0));
      d_r.rotateN(&d_output[0], burst, burst_size);

      // build and publish the new PDU
      pmt::pmt_t pdu_vector = pmt::init_c32vector(burst_size, d_output);

      metadata = pmt::dict_add(metadata, pmt::mp("center_frequency"), pmt::mp(fc3db + center_frequency));
      metadata = pmt::dict_add(metadata, pmt::mp("relative_frequency"), pmt::mp(fc3db + relative_frequency));
      metadata = pmt::dict_add(metadata, pmt::mp("bandwidth"), pmt::mp(bandwidth));
      metadata = pmt::dict_add(metadata, pmt::mp("burst_power"), pmt::mp(power));

      pmt::pmt_t out_msg = pmt::cons(metadata, pdu_vector);
      message_port_pub(pmt::mp("pdu_out"), out_msg);

      return;


    }


  } /* namespace fhss_utils */
} /* namespace gr */
