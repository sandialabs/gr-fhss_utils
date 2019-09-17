INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_FHSS_UTILS fhss_utils)

FIND_PATH(
    FHSS_UTILS_INCLUDE_DIRS
    NAMES fhss_utils/api.h
    HINTS $ENV{FHSS_UTILS_DIR}/include
        ${PC_FHSS_UTILS_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    FHSS_UTILS_LIBRARIES
    NAMES gnuradio-fhss_utils
    HINTS $ENV{FHSS_UTILS_DIR}/lib
        ${PC_FHSS_UTILS_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FHSS_UTILS DEFAULT_MSG FHSS_UTILS_LIBRARIES FHSS_UTILS_INCLUDE_DIRS)
MARK_AS_ADVANCED(FHSS_UTILS_LIBRARIES FHSS_UTILS_INCLUDE_DIRS)

