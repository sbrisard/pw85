cmake_minimum_required(VERSION 3.13)

project("tutorial")

find_package(GSL REQUIRED)
find_package(pw85 CONFIG REQUIRED)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.c)
target_link_libraries(${PROJECT_NAME} pw85)
