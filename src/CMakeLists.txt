cmake_minimum_required(VERSION 3.12)
project (pw85 LANGUAGES C)
add_library (pw85 SHARED pw85.c)
install (TARGETS pw85
         RUNTIME DESTINATION bin
	 LIBRARY DESTINATION lib
	 ARCHIVE DESTINATION lib/pw85)
install (FILES pw85.h DESTINATION include)
