find_library(MATH_LIBRARY m)

add_compile_definitions(PROJECT_VERSION_QUOTED=\"${METADATA_VERSION}\")
add_compile_definitions(PROJECT_AUTHOR_QUOTED=\"${METADATA_AUTHOR}\")

add_library(
  ${PROJECT_NAME} SHARED
  pw85.c
  ../include/${PROJECT_NAME}/pw85.h)

target_include_directories(
  ${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)

if (MATH_LIBRARY)
  target_link_libraries(${PROJECT_NAME} ${MATH_LIBRARY})
endif()

include(FindGSL)
#find_path(GSL_INCLUDE_DIR gsl_math.h)
message(${GSL_INCLUDE_DIRS})
message(${GSL_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})

install(
  TARGETS ${PROJECT_NAME}
  DESTINATION lib
  EXPORT ${PROJECT_NAME}_targets)

install(
  EXPORT ${PROJECT_NAME}_targets
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# write_basic_package_version_file(
#   "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
#   VERSION "${Tutorial_VERSION_MAJOR}.${Tutorial_VERSION_MINOR}"
#   COMPATIBILITY AnyNewerVersion
# )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})
