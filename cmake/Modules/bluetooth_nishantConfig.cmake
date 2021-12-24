INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_bluetooth_nishant bluetooth_nishant)

FIND_PATH(
    BLUETOOTH_NISHANT_INCLUDE_DIRS
    NAMES bluetooth_nishant/api.h
    HINTS $ENV{BLUETOOTH_NISHANT_DIR}/include
        ${PC_BLUETOOTH_NISHANT_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    BLUETOOTH_NISHANT_LIBRARIES
    NAMES gnuradio-bluetooth-nishant
    HINTS $ENV{BLUETOOTH_NISHANT_DIR}/lib
        ${PC_BLUETOOTH_NISHANT_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/bluetooth_nishantTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BLUETOOTH_NISHANT DEFAULT_MSG BLUETOOTH_NISHANT_LIBRARIES BLUETOOTH_NISHANT_INCLUDE_DIRS)
MARK_AS_ADVANCED(BLUETOOTH_NISHANT_LIBRARIES BLUETOOTH_NISHANT_INCLUDE_DIRS)
