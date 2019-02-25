cmake_minimum_required(VERSION 3.0)

#------------------------------------------------------------------------------
set (metadata_dir "../metadata")

file (STRINGS ${metadata_dir}/name.txt lines ENCODING UTF-8)
list (GET lines 0 project_name)

file (STRINGS ${metadata_dir}/version.txt lines ENCODING UTF-8)
list (GET lines 0 project_version)

file (STRINGS ${metadata_dir}/author.txt lines ENCODING UTF-8)
list (GET lines 0 project_author)

file (STRINGS ${metadata_dir}/url.txt lines ENCODING UTF-8)
list (GET lines 0 project_url)

file (STRINGS ${metadata_dir}/copyright.txt lines ENCODING UTF-8)
list (GET lines 0 project_copyright)

file (STRINGS ${metadata_dir}/short_description.txt lines ENCODING UTF-8)
list (GET lines 0 project_short_description)

file (READ ${metadata_dir}/description.txt project_description ENCODING UTF-8)

set (doc_title "Documentation of the ${project_name} library")
set (doc_basename "pw85")

#------------------------------------------------------------------------------
project (pw85 LANGUAGES C)

add_library (pw85 SHARED pw85.c pw85.h)
set_target_properties (pw85 PROPERTIES VERSION ${project_version})
configure_file (${pw85_SOURCE_DIR}/pw85-config.cmake.in
  ${pw85_BINARY_DIR}/pw85-config.cmake @ONLY)
configure_file (${pw85_SOURCE_DIR}/pw85-config-version.cmake.in
  ${pw85_BINARY_DIR}/pw85-config-version.cmake @ONLY)
configure_file (${pw85_SOURCE_DIR}/pw85.h.in
  ${pw85_BINARY_DIR}/pw85.h
  @ONLY)
include_directories (${pw85_BINARY_DIR})
# TODO: place the output file conf.py in another directory?
configure_file (${pw85_SOURCE_DIR}/../sphinx/conf.py.in
  ${pw85_SOURCE_DIR}/../sphinx/conf.py
  @ONLY)
configure_file (${pw85_SOURCE_DIR}/../sphinx/index.rst.in
  ${pw85_SOURCE_DIR}/../sphinx/index.rst
  @ONLY)
install (FILES ${pw85_BINARY_DIR}/pw85-config.cmake
  ${pw85_BINARY_DIR}/pw85-config-version.cmake
  DESTINATION lib/pw85-${project_version})
install (TARGETS pw85 EXPORT pw85-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
install (EXPORT pw85-targets
  DESTINATION lib/pw85-${project_version})
install (FILES ${pw85_BINARY_DIR}/pw85.h DESTINATION include)
