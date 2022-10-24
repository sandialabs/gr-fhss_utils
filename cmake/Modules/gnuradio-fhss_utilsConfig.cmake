find_package(PkgConfig)

PKG_CHECK_MODULES(PC_GR_FHSS_UTILS gnuradio-fhss_utils)

FIND_PATH(
    GR_FHSS_UTILS_INCLUDE_DIRS
    NAMES gnuradio/fhss_utils/api.h
    HINTS $ENV{FHSS_UTILS_DIR}/include
        ${PC_FHSS_UTILS_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GR_FHSS_UTILS_LIBRARIES
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

include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-fhss_utilsTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GR_FHSS_UTILS DEFAULT_MSG GR_FHSS_UTILS_LIBRARIES GR_FHSS_UTILS_INCLUDE_DIRS)
MARK_AS_ADVANCED(GR_FHSS_UTILS_LIBRARIES GR_FHSS_UTILS_INCLUDE_DIRS)
