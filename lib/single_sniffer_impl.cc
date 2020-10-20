/* -*- c++ -*- */
/* 
 * Copyright 2020 Free Software Foundation, Inc.
 * 
 * This file is part of gr-bluetooth
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
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
#include "single_sniffer_impl.h"

#include <gnuradio/digital/clock_recovery_mm_ff.h>
#include <gnuradio/analog/quadrature_demod_cf.h>
#include <gnuradio/digital/binary_slicer_fb.h>
#include "gr_bluetooth/no_filter_sniffer.h"

namespace gr {
namespace bluetooth {

    single_sniffer::sptr single_sniffer::make(double sample_rate, double center_freq)
    {
        return gnuradio::get_initial_sptr (new single_sniffer_impl(sample_rate, center_freq));
    }

    /*
     * The private constructor
     */
    single_sniffer_impl::single_sniffer_impl(double sample_rate, double center_freq)
        : gr::hier_block2 ("bluetooth single sniffer block",
                gr::io_signature::make (1, 1, sizeof (gr_complex)),
                gr::io_signature::make (0, 0, 0))
    {
        float gain = 3.125 / M_PI_2;
        gr::analog::quadrature_demod_cf::sptr fm_demod = 
            gr::analog::quadrature_demod_cf::make(gain);

        float omega = 3.125;
        float gain_mu = 0.175;
        float gain_omega = .25 * gain_mu * gain_mu;
        float mu = 0.32;
        float omega_relative_limit = 0.005;
        gr::digital::clock_recovery_mm_ff::sptr mm_cr =
            gr::digital::clock_recovery_mm_ff::make(omega, gain_omega, mu, gain_mu, omega_relative_limit);

        gr::digital::binary_slicer_fb::sptr bin_slice =
            gr::digital::binary_slicer_fb::make();

        no_filter_sniffer::sptr sniffer =
            no_filter_sniffer::make(sample_rate, center_freq);

        connect(self(), 0, fm_demod, 0);
        connect(fm_demod, 0, mm_cr, 0);
        connect(mm_cr, 0, bin_slice, 0);
        connect(bin_slice, 0, sniffer, 0);
    }

    /*
     * Our virtual destructor.
     */
    single_sniffer_impl::~single_sniffer_impl()
    {
    }

} /* namespace bluetooth */
} /* namespace gr */

