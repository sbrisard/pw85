add_executable(test_pw85
        test_pw85.cpp)

target_include_directories(test_pw85 PUBLIC ../include)

if (MATH_LIBRARY)
    target_link_libraries(test_pw85 ${MATH_LIBRARY})
endif ()

find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(test_pw85
        pw85
        Eigen3::Eigen
        hdf5::hdf5
        hdf5::hdf5_hl)

set(PW85_REF_DATA_FILENAME "pw85_ref_data.h5")
set(PW85_REF_DATA_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME})
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PW85_REF_DATA_FILENAME})
    message("-- Downloading reference data")
    execute_process(COMMAND python download_ref_data.py --pw85_ref_data_path ${CMAKE_CURRENT_BINARY_DIR}/${PW85_REF_DATA_FILENAME} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.)
endif ()

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PW85_REF_DATA_FILENAME}
        DESTINATION ${PW85_REF_DATA_DIR})

add_compile_definitions(PW85_REF_DATA_PATH="${PW85_REF_DATA_FILENAME}")

add_test(NAME test_pw85 COMMAND test_pw85)
