cmake_minimum_required (VERSION 3.5.1)

project (IBEX_HCONTROL VERSION 0.42 LANGUAGES CXX)
set (IBEX_HCONTROL_DESCRIPTION "H_infinity Control synthesis using global optimization")
set (IBEX_HCONTROL_URL "https://github.com/ibex-team/ibex-Hcontrol/")
set (IBEX_MIN_VERSION_REQ 2.8.9.1)

# If IBEX_DIR was given on the command-line and not in the environment, set it
# as an environment variable. Enable the usage of 'cmake -DIBEX_DIR=<path> ...'
if (NOT DEFINED ENV{IBEX_DIR} AND DEFINED IBEX_DIR)
  set (ENV{IBEX_DIR} "${IBEX_DIR}")
endif ()
find_package (IBEX ${IBEX_MIN_VERSION_REQ} REQUIRED)
message (STATUS "Found Ibex version ${IBEX_VERSION}")
ibex_init_common () # Ibex should have installed this function
set (CMAKE_INSTALL_CONFIGCMAKE "share/ibex_Hcontrol/cmake")

# Interval and linear programming names and versions are available
message (STATUS "Ibex was compiled with interval library "
                "${IBEX_INTERVAL_LIB_NAME} ${IBEX_INTERVAL_LIB_VERSION}")

if (IBEX_LP_LIB_NAME) # Was Ibex compiled with a linear programming library ?
  message (STATUS "Ibex was compiled with Linear Programming library "
                  "${IBEX_LP_LIB_NAME} ${IBEX_LP_LIB_VERSION}")
else ()
  message (STATUS "Ibex was compiled without any Linear Programming library")
endif ()

find_package (IBEX_AFFINE REQUIRED)
message (STATUS "Found Ibex-affine")

################################################################################
# Generate pkg-config file and install it (can be removed if no library is
# installed)
################################################################################
configure_file ("ibex-Hcontrol.pc.in" "ibex-Hcontrol.pc" @ONLY)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ibex-Hcontrol.pc
          DESTINATION ${CMAKE_INSTALL_PKGCONFIG})

################################################################################
# Compile sources
################################################################################
add_subdirectory (src)

################################################################################
# Generate CMake config file (can be removed if no library is installed)
################################################################################
set (IBEX_HCONTROL_TARGETS_CONFIG_FILE "ibex_hcontrol-targets.cmake")
install (EXPORT HcontrolExport FILE ${IBEX_HCONTROL_TARGETS_CONFIG_FILE}
                                NAMESPACE Ibex::
                                DESTINATION ${CMAKE_INSTALL_CONFIGCMAKE})
configure_file (ibex-Hcontrol-config.cmake.in ibex_Hcontrol-config.cmake @ONLY)
include (CMakePackageConfigHelpers)
set (VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/ibex_Hcontrol-config-version.cmake")
write_basic_package_version_file (${VERSION_FILE} COMPATIBILITY AnyNewerVersion)
install (FILES ${VERSION_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/ibex_Hcontrol-config.cmake
         DESTINATION ${CMAKE_INSTALL_CONFIGCMAKE})

################################################################################
# Tests
################################################################################
# The following command assumes that the tests are in a tests/ subdirectory. So
# the tests/ subdirectory must exist and must contain a CMakeLists.txt file.
# If this is not the case, comment the next line.
ibex_add_make_target_for_ctest (check)
