cmake_minimum_required(VERSION 3.0)

#------------------------------------------------------------------------------
set (project_name "PW85")
set (project_version 1.0)
set (project_author "Sébastien Brisard")
set (project_url "https://github.com/sbrisard/pw85.git")
set (project_copyright "2018")
set (project_short_description "Overlap test of two ellipsoids.")
set (project_description "${project_short_description}

PW85 provides an implementation of the “contact function”
defined by Perram and Wertheim (J. Comp. Phys. 58(3), 409–416) for two
ellipsoids. Given two ellipsoids, this function returns the *square*
of the common factor by which both ellipsoids must be scaled (their
centers being fixed) in order to be tangentially in contact.

PW85 is released under a BSD 3-Clause License.")

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
