/* -*- c++ -*- */
/* 
 * Copyright 2013 Christopher D. Kilgour
 * Copyright 2008, 2009 Dominic Spill, Michael Ossmann
 * Copyright 2007 Dominic Spill
 * Copyright 2005, 2006 Free Software Foundation, Inc.
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

#ifndef INCLUDED_BLUETOOTH_GR_BLUETOOTH_NO_FILTER_SNIFFER_IMPL_H
#define INCLUDED_BLUETOOTH_GR_BLUETOOTH_NO_FILTER_SNIFFER_IMPL_H

#include "gr_bluetooth/no_filter_sniffer.h"
#include "gr_bluetooth/packet.h"
#include "gr_bluetooth/piconet.h"
#include <math.h>
#include <map>

namespace gr {
namespace bluetooth {

    class no_filter_sniffer_impl : virtual public no_filter_sniffer
    {
        private:
            /* symbols per second */
            static const int SYMBOL_RATE = 1000000;
            static const int SYMBOLS_PER_BASIC_RATE_SHORTENED_ACCESS_CODE = 68; 

            /* length of time slot in symbols */
            static const int SYMBOLS_PER_BASIC_RATE_SLOT    = 625;
            static const int SYMBOLS_FOR_BASIC_RATE_HISTORY = 3125;

            /* channel 0 in Hz */
            static const uint32_t BASE_FREQUENCY = 2402000000UL;

            /* channel width in Hz */
            static const int CHANNEL_WIDTH = 1000000;

            /* total number of samples elapsed */
            uint64_t d_cumulative_count;

            /* frequency and number of the channel being decoded */
            double d_channel_freq;
            int d_channel;

            /* General Inquiry and Limited Inquiry Access Codes */
            static const uint32_t GIAC = 0x9E8B33;
            static const uint32_t LIAC = 0x9E8B00;

            /* the piconets we are monitoring */
            std::map<int, basic_rate_piconet::sptr> d_basic_rate_piconets;

            /* handle AC */
            void ac(char *symbols, int max_len, double freq, int offset);

            /* handle ID packet (no header) */
            void id(uint32_t lap);

            /* decode packets with headers */
            void decode(classic_packet::sptr pkt, basic_rate_piconet::sptr pn,
                    bool first_run);

            /* work on UAP/CLK1-6 discovery */
            void discover(classic_packet::sptr pkt, basic_rate_piconet::sptr pn);

            /* decode stored packets */
            void recall(basic_rate_piconet::sptr pn);

            /* pull information out of FHS packet */
            void fhs(classic_packet::sptr pkt);

        public:
            no_filter_sniffer_impl(double sample_rate, double center_freq);
            ~no_filter_sniffer_impl();

            // Where all the action really happens
            int work(int                        noutput_items,
                    gr_vector_const_void_star& input_items,
                    gr_vector_void_star&       output_items);
    };

} // namespace bluetooth
} // namespace gr

#endif /* INCLUDED_BLUETOOTH_GR_BLUETOOTH_NO_FILTER_SNIFFER_IMPL_H */

