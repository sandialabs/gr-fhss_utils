/* -*- c++ -*- */

#define FHSS_UTILS_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "fhss_utils_swig_doc.i"

%{
#include "fhss_utils/fft_burst_tagger.h"
#include "fhss_utils/tagged_burst_to_pdu.h"
#include "fhss_utils/burst_downmix.h"
#include "fhss_utils/pdu_quadrature_demod_cf.h"
#include "fhss_utils/burst_measure.h"
#include "fhss_utils/fine_burst_measure.h"
#include "fhss_utils/fsk_cfo_estimate.h"
#include "fhss_utils/fft_cfo_estimate.h"
#include "fhss_utils/coerce_frequency_estimate.h"
%}


%include "fhss_utils/fft_burst_tagger.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, fft_burst_tagger);
%include "fhss_utils/tagged_burst_to_pdu.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, tagged_burst_to_pdu);
%include "fhss_utils/burst_downmix.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, burst_downmix);
%include "fhss_utils/pdu_quadrature_demod_cf.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, pdu_quadrature_demod_cf);
%include "fhss_utils/burst_measure.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, burst_measure);
%include "fhss_utils/fine_burst_measure.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, fine_burst_measure);
%include "fhss_utils/fsk_cfo_estimate.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, fsk_cfo_estimate);
%include "fhss_utils/fft_cfo_estimate.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, fft_cfo_estimate);
%include "fhss_utils/coerce_frequency_estimate.h"
GR_SWIG_BLOCK_MAGIC2(fhss_utils, coerce_frequency_estimate);
