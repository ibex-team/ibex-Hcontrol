################################################################################
# First example: compiling a library that depends on libibex
################################################################################
# source files of libibex-Hcontrol
list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/strategy/ibex_LightLocalSolver.cpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/strategy/ibex_LightLocalSolver.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/strategy/ibex_LightOptimMinMax.cpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/strategy/ibex_LightOptimMinMax.h
                 ${CMAKE_CURRENT_SOURCE_DIR}/strategy/ibex_OptimMinMax.cpp
                 ${CMAKE_CURRENT_SOURCE_DIR}/strategy/ibex_OptimMinMax.h
                 # Add more source files here.
                 # For most plugins, the rest can be left as is.
                 )

# Create the target for libibex-myplugin
add_library (ibex-Hcontrol ${SRC})
target_include_directories (ibex-Hcontrol PUBLIC
                            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries (ibex-Hcontrol PUBLIC Ibex::ibex)

# installation of libibex-Hcontrol files
ibex_list_filter_header (HDR ${SRC}) # Ibex should have installed this function
install (TARGETS ibex-Hcontrol EXPORT HcontrolExport
                               DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (FILES ${HDR} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ibex-Hcontrol)

################################################################################
# More examples in the bin subdirectory
################################################################################
add_subdirectory (bin)
