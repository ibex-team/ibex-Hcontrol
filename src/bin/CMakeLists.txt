################################################################################
# First example: compiling a binary that depends on libibex-Hcontrol
################################################################################
# source files of bin1
list (APPEND SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/bin1.cpp
                  # Add more source files here.
                  # For most plugins, the rest can be left as is.
                  )

# Create the target for bin1
add_executable (bin1 ${SRC1})
target_link_libraries (bin1 PUBLIC ibex-Hcontrol)

# install the binary
install (TARGETS bin1 DESTINATION ${CMAKE_INSTALL_BINDIR})

################################################################################
# Second example: compiling a binary that depends only on libibex
################################################################################
# source files of bin2
list (APPEND SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/bin2.cpp
                  # Add more source files here.
                  # For most plugins, the rest can be left as is.
                  )

# Create the target for bin2
add_executable (bin2 ${SRC2})
target_link_libraries (bin2 PUBLIC Ibex::ibex)

# install the binary
install (TARGETS bin2 DESTINATION ${CMAKE_INSTALL_BINDIR})

