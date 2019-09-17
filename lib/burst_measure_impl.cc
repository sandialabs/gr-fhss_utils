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
#include <gnuradio/fft/window.h>
#include <volk/volk.h>
#include "burst_measure_impl.h"

namespace gr {
  namespace fhss_utils {

    burst_measure::sptr
    burst_measure::make()
    {
      return gnuradio::get_initial_sptr
        (new burst_measure_impl());
    }

    /*
     * The private constructor
     */
    burst_measure_impl::burst_measure_impl()
      : gr::block("burst_measure",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
        d_mags(nullptr)
    {
      message_port_register_in(pmt::mp("cpdus"));
      message_port_register_out(pmt::mp("cpdus"));
      set_msg_handler(pmt::mp("cpdus"), boost::bind(&burst_measure_impl::handler, this, _1));

      fft_setup(15);
    }

    void burst_measure_impl::fft_setup(int power) {
      for (int i=d_ffts.size(); i<=power; i++) {

        // init fft
        int fftsize = pow(2, i);
        d_ffts.push_back(new gr::fft::fft_complex(fftsize, true, 1));

        // init window
        d_windows.push_back((float*)volk_malloc(sizeof(float)*fftsize, volk_get_alignment()));
        std::vector<float> window = gr::fft::window::hamming(fftsize);
        memcpy(d_windows.back(), window.data(), sizeof(float)*fftsize);

        // init d_mags
        if (d_mags != nullptr) {
          volk_free(d_mags);
        }
        d_mags = (float*) volk_malloc(sizeof(float)*fftsize, volk_get_alignment());

      }
    }

    void burst_measure_impl::fft_cleanup() {
      for (gr::fft::fft_complex* fft : d_ffts) {
        delete fft;
      }
      d_ffts.clear();
      for (float* win : d_windows) {
        volk_free(win);
      }
      d_windows.clear();
      if (d_mags != nullptr) {
        volk_free(d_mags);
      }
      d_mags = nullptr;
    }

    /*
     * Our virtual destructor.
     */
    burst_measure_impl::~burst_measure_impl()
    {
      fft_cleanup();
    }

    void burst_measure_impl::handler(pmt::pmt_t msg) {
      pmt::pmt_t meta  = pmt::car(msg);
      float sample_rate = pmt::to_float(pmt::dict_ref(meta, pmt::mp("sample_rate"), pmt::PMT_NIL));
      //float noise_power = pmt::to_float(pmt::dict_ref(meta, pmt::mp("noise_power"), pmt::PMT_NIL));

      pmt::pmt_t samples = pmt::cdr(msg);
      size_t burst_size = pmt::length(samples);
      const gr_complex* burst = (const gr_complex*)pmt::c32vector_elements(samples, burst_size);

      // make sure burst_size is divisible by 4, so that fftsize is even
      burst_size &= ~0x3lu;
      
      // only take FFT of middle half
      uint64_t burst_start = burst_size/4;
      size_t adjusted_size = burst_size/2;
      int fft_power = ceil(log2(adjusted_size));
      fft_setup(fft_power);
      int fftsize = pow(2, fft_power);
      gr_complex* inbuf = d_ffts[fft_power]->get_inbuf();

      // apply window
      volk_32fc_32f_multiply_32fc(inbuf, burst+burst_start, d_windows[fft_power], adjusted_size);
      memset(inbuf+adjusted_size, 0, sizeof(gr_complex) * (fftsize - adjusted_size));

      // take FFT
      d_ffts[fft_power]->execute();

      // calculate magnitude squared
      gr_complex* outbuf = d_ffts[fft_power]->get_outbuf();
      volk_32fc_magnitude_squared_32f(d_mags, outbuf, fftsize);

      // calculate total energy
      float energy = 0;
      volk_32f_accumulator_s32f(&energy, d_mags, fftsize);

      // calculate bandwidth
      float energy90 = energy * 0.90;
      int last_bin = 0;
      float sum = d_mags[0];
      for (int i=1; i<fftsize/2-1; i++) {
        sum += d_mags[i] + d_mags[fftsize - i];
        if (sum >= energy90) {
          last_bin = i;
          break;
        }
      }

      float bandwidth = (last_bin*2 + 1) * sample_rate / fftsize;
      float power = energy / adjusted_size;

      meta = pmt::dict_add(meta, pmt::mp("burst_power"), pmt::from_float(10*log10(power)));
      meta = pmt::dict_add(meta, pmt::mp("bandwidth"), pmt::from_float(bandwidth));

      message_port_pub(pmt::mp("cpdus"), pmt::cons(meta, samples));
    }

  } /* namespace fhss_utils */
} /* namespace gr */

// vim: ts=2:sw=2:sts=2

