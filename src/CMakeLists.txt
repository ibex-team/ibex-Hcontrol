################################################################################
# First example: compiling a library that depends on libibex
################################################################################
# source files of libibex-Hcontrol
list (APPEND SRC	strategy/ibex_CellMinMaxHeap.h
       				strategy/ibex_EvalMax.cpp
       				strategy/ibex_EvalMax.h
       				strategy/ibex_BxpMinMax.cpp
       				strategy/ibex_BxpMinMax.h
        # Add more source files here.
        # For most plugins, the rest can be left as is.
        )

# Create the target for libibex-Hcontrol
add_library (ibex-Hcontrol SHARED ${SRC})
target_include_directories (ibex-Hcontrol PUBLIC
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/strategy>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
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
