# Copyright 2011,2012,2016,2018,2019 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

find_package(BTBB REQUIRED)

list(APPEND bluetooth_sources
    tun.cc
    multi_block.cc
    multi_hopper_impl.cc
    multi_LAP_impl.cc
    multi_sniffer_impl.cc
    multi_UAP_impl.cc
    no_filter_sniffer_impl.cc
    single_sniffer_impl.cc
    packet_impl.cc
    piconet_impl.cc
)

set(bluetooth_sources "${bluetooth_sources}" PARENT_SCOPE)
if(NOT bluetooth_sources)
    MESSAGE(STATUS "No C++ sources... skipping lib/")
    return()
endif(NOT bluetooth_sources)

find_package(PythonLibs)
include(GrPython)

add_library(gnuradio-bluetooth SHARED ${bluetooth_sources})
target_link_libraries(gnuradio-bluetooth gnuradio::gnuradio-runtime
  gnuradio::gnuradio-blocks
  gnuradio::gnuradio-filter
  gnuradio::gnuradio-fft
  gnuradio::gnuradio-analog
  gnuradio::gnuradio-digital
  Python::Python
  )
target_include_directories(gnuradio-bluetooth
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    PUBLIC $<INSTALL_INTERFACE:include>
  )
set_target_properties(gnuradio-bluetooth PROPERTIES DEFINE_SYMBOL "gnuradio_bluetooth_EXPORTS")

if(APPLE)
    set_target_properties(gnuradio-bluetooth PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-bluetooth
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

