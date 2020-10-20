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


#ifndef INCLUDED_GR_BLUETOOTH_NO_FILTER_SNIFFER_H
#define INCLUDED_GR_BLUETOOTH_NO_FILTER_SNIFFER_H

#include <gr_bluetooth/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
namespace bluetooth {

    /*!
     * \brief <+description of block+>
     * \ingroup bluetooth
     *
     */
    class GR_BLUETOOTH_API no_filter_sniffer : virtual public gr::sync_block
    {
        public:
            typedef boost::shared_ptr<no_filter_sniffer> sptr;

            /*!
             * \brief Return a shared_ptr to a new instance of gr::bluetooth::no_filter_sniffer.
             *
             * To avoid accidental use of raw pointers, gr::bluetooth::no_filter_sniffer's
             * constructor is in a private implementation
             * class. gr::bluetooth::no_filter_sniffer::make is the public interface for
             * creating new instances.
             */
            static sptr make(double sample_rate, double center_freq);
    };

} // namespace bluetooth
} // namespace gr

#endif /* INCLUDED_GR_BLUETOOTH_NO_FILTER_SNIFFER_H */

