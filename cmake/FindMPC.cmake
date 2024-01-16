# Try to find MPC
# Once done this will define:
# - MPC_FOUND - True if the system has GMP
# - MPC_INCLUDE_DIRS - include directories for compiling
# - MPC_LIBRARIES - libraries for linking
# - MPC_DEFINITIONS - cflags suggested by pkg-config

find_package(PkgConfig)
pkg_check_modules(PC_MPC QUIET mpc)

set(MPC_DEFINITIONS ${PC_MPC_CFLAGS_OTHER})

# Find path of Sylvan file
find_path(MPC_INCLUDE_DIR mpc.h
          HINTS ${PC_MPC_INCLUDEDIR} ${PC_MPC_INCLUDE_DIRS})

# Find MPC library
find_library(MPC_LIBRARIES NAMES mpc libmpc 
             HINTS ${PC_MPC_LIBDIR} ${PC_MPC_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MPC_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MPC DEFAULT_MSG MPC_LIBRARIES MPC_INCLUDE_DIR)

mark_as_advanced(MPC_INCLUDE_DIR MPC_LIBRARIES)
