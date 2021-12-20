# source files of bin2
list (APPEND SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/bin1.cpp
                  # Add more source files here.
                  # For most plugins, the rest can be left as is.
                  )

# Create the target for bin2
add_executable (bin1 ${SRC2})
target_link_libraries (bin1 PUBLIC Ibex::ibex ibex-Hcontrol)

# install the binary
install (TARGETS bin1 DESTINATION ${CMAKE_INSTALL_BINDIR})

