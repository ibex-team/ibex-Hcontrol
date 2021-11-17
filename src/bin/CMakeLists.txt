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

