cmake_minimum_required(VERSION 3.13)

project("tutorial" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)

find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(${PROJECT_NAME} INTERFACE ${MATH_LIBRARY})
endif()

find_package(Boost REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Boost::headers)

target_include_directories(${PROJECT_NAME} PRIVATE "../../include")
