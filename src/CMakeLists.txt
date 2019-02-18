cmake_minimum_required(VERSION 3.12)
project (pw85 LANGUAGES C)
set (version 1.0)
add_library (pw85 SHARED pw85.c pw85.h)
set_target_properties (pw85 PROPERTIES VERSION ${version})
install (TARGETS pw85 EXPORT pw85-targets
         RUNTIME DESTINATION bin
	 LIBRARY DESTINATION lib
	 ARCHIVE DESTINATION lib)
install (EXPORT pw85-targets
         DESTINATION lib/pw85-${version})
install (FILES pw85.h DESTINATION include)
configure_file (${pw85_SOURCE_DIR}/pw85-config.cmake.in
                ${pw85_BINARY_DIR}/pw85-config.cmake @ONLY)
configure_file (${pw85_SOURCE_DIR}/pw85-config-version.cmake.in
                ${pw85_BINARY_DIR}/pw85-config-version.cmake @ONLY)
install (FILES ${pw85_BINARY_DIR}/pw85-config.cmake
               ${pw85_BINARY_DIR}/pw85-config-version.cmake
	       DESTINATION lib/pw85-${version})
