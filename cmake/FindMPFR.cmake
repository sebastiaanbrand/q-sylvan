# Try to find MPFR
# Once done this will define:
# - MPFR_FOUND - True if the system has MPFR
# - MPFR_INCLUDE_DIRS - include directories for compiling
# - MPFR_LIBRARIES - libraries for linking
# - MPFR_DEFINITIONS - cflags suggested by pkg-config

find_package(PkgConfig)
pkg_check_modules(PC_MPFR QUIET mpfr)

set(MPFR_DEFINITIONS ${PC_MPFR_CFLAGS_OTHER})

# Find path of Sylvan file
find_path(MPFR_INCLUDE_DIR mpfr.h
          HINTS ${PC_MPFR_INCLUDEDIR} ${PC_MPFR_INCLUDE_DIRS})

# Find MPC library
find_library(MPFR_LIBRARIES NAMES mpfr libmpfr 
             HINTS ${PC_MPFR_LIBDIR} ${PC_MPFR_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MPFR_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_LIBRARIES MPFR_INCLUDE_DIR)

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)
