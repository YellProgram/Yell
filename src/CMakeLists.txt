cmake_minimum_required (VERSION 2.6)
project (YELL)

if(MSVC)
    foreach(flag_var CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE CMAKE_CXX_FLAGS_RELEASE CMAKE_C_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
        string(REGEX REPLACE "/MDd" "/MTd" ${flag_var} "${${flag_var}}")
    endforeach(flag_var)
    SET (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING "MSVC C Debug MT flags " FORCE)
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "MSVC CXX Debug MT flags " FORCE)
    SET (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "MSVC C Release MT flags " FORCE)
    SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "MSVC CXX Release MT flags " FORCE)
    SET (CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL}" CACHE STRING "MSVC C Debug MT flags " FORCE)
    SET (CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL}" CACHE STRING "MSVC C Release MT flags " FORCE)
    SET (CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "MSVC CXX Debug MT flags " FORCE)
    SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "MSVC CXX Release MT flags " FORCE)
endif()

set(SCITBX_PATH ../lib)
set(CCTBX_PATH ../lib/cctbx)
set(LEVMAR_PATH ../lib/levmar-2.4)
set(LAPACK_PATH ../lib/clapack-3.2.1-CMAKE) #notmac
set(CXXTEST_PATH ../lib/cxxtest)
set(SZIP_PATH ../lib/szip-2.1)
set(ZLIB_PATH ../lib/zlib-1.2.5)

add_subdirectory(${LAPACK_PATH} lapack EXCLUDE_FROM_ALL)
add_subdirectory(${LEVMAR_PATH} levmar EXCLUDE_FROM_ALL)
add_subdirectory(${ZLIB_PATH} zlib EXCLUDE_FROM_ALL)
add_subdirectory(${SZIP_PATH} szip EXCLUDE_FROM_ALL)

set(BOOST_INCLUDEDIR ../lib/boost_1_45_0) #Intentionally has the same name as variable from findBoost.cmake. Later replace this definition with "FIND_PACKAGE( Boost …"

#find_package(HDF5) # As I am looser, hdf5 can not find the proper place either
set(HDF5_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/../libs/include" "${PROJECT_BINARY_DIR}/../libs/include/cpp")
set(HDF5_LIB_DIRS ${PROJECT_BINARY_DIR}/../libs/lib)
set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl)
set(LAPACK_LIBRARIES f2c lapack blas)


#for mac this will be
#set(LAPACK_LIBRARIES, LAPACK_DIRS #framework accelerate

include_directories(${SCITBX_PATH} ${CCTBX_PATH} ${LEVMAR_PATH} ${HDF5_INCLUDE_DIRS} ${BOOST_INCLUDEDIR} ${CXXTEST_PATH} ${HDF5_PATH}/c++/src ${HDF5_PATH}/src) #THe order makes difference. SCITBX_PATH shoud come before CCTBX_PATH !!
link_directories(${HDF5_LIB_DIRS} ${LAPACK_DIRS})


# Manualy get what we need from cctbx
add_library(cctbx ${CCTBX_PATH}/cctbx/sgtbx/change_of_basis_op.cpp
${CCTBX_PATH}/cctbx/uctbx/crystal_orientation.cpp
${CCTBX_PATH}/cctbx/sgtbx/rot_mx.cpp
${CCTBX_PATH}/cctbx/sgtbx/rt_mx.cpp
${CCTBX_PATH}/cctbx/uctbx/spoil_optimization.cpp
${CCTBX_PATH}/omptbx/stubs.cpp
${CCTBX_PATH}/cctbx/sgtbx/tr_vec.cpp
${CCTBX_PATH}/cctbx/uctbx/uctbx.cpp
${CCTBX_PATH}/omptbx/stubs.cpp
${CCTBX_PATH}/cctbx/eltbx/xray_scattering/wk1995.cpp
${CCTBX_PATH}/cctbx/eltbx/basic.cpp
${CCTBX_PATH}/cctbx/sgtbx/utils.cpp)



add_executable(yell basic_classes.cpp
                    diffuser_core.cpp
                    basic_io.cpp
                    FormulaParser.cpp
                    InputFileParser.cpp
                    InputFileParserI.cpp
                    InputFileParserII.cpp
                    InputFileParserIII.cpp
                    main.cpp
                    model.cpp
                    OutputHandler.cpp)

#add_executable(test-yell basic_classes.cpp
#                    diffuser_core.cpp
#                    basic_io.cpp
#                    FormulaParser.cpp
#                    InputFileParser.cpp
#                    InputFileParserI.cpp
#                    InputFileParserII.cpp
#                    InputFileParserIII.cpp
#                    test_runner.cpp
#                    model.cpp
#                    OutputHandler.cpp)

ADD_DEPENDENCIES(yell levmar cctbx zlib szip ${LAPACK_LIBRARIES})


#MS studio includes tr1 libraries by default. To disable it we use
add_definitions(-D_HAS_CPP0X=0)

target_link_libraries (yell levmar cctbx zlib szip ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries (test-yell levmar cctbx zlib szip ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})

#target_link_libraries (yell levmar cctbx ${HDF5_z_LIBRARY} ${HDF5_sz_LIBRARY} ${HDF5_hdf5_cpp_LIBRARY} ${HDF5_hdf5_cpp_hl_LIBRARY} ${HDF5_z_LIBRARY} ${HDF5_m_LIBRARY} ${HDF5_hdf5_LIBRARY} ${HDF5_hdf5_hl_LIBRARY})


#-lm