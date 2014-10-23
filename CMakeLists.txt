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
set(CCTBX_PATH lib/cctbx-v-10/cctbx_sources CACHE PATH "Path to cctbx c++ files")
set(LEVMAR_PATH lib/levmar CACHE PATH "Path to levmar library")
set(SCITBX_PATH lib/cctbx_build/include CACHE PATH "Our library depends on some of the cctbx files which are not readily available in the cctbx source bundle. Rather, they are compiled in one of the stages of cctbx build process. This variable should link to the folder ... which will be created after the command '/your/choice/bin/python ../cctbx_sources/libtbx/configure.py mmtbx' for more information on the installation process see http://cctbx.sourceforge.net/current_cvs/installation.html#manually-building-from-sources-under-unix")
#set(LAPACK_PATH lib/clapack)
#set(CXXTEST_PATH lib/cxxtest)
#set(SZIP_PATH lib/szip)
#set(ZLIB_PATH lib/zlib)

#add_subdirectory(${LAPACK_PATH} lapack EXCLUDE_FROM_ALL)
add_subdirectory(${LEVMAR_PATH} levmar EXCLUDE_FROM_ALL)
#add_subdirectory(${ZLIB_PATH} zlib EXCLUDE_FROM_ALL)
#add_subdirectory(${SZIP_PATH} szip EXCLUDE_FROM_ALL)

#set(BOOST_INCLUDEDIR ../lib/boost)

#find_package(HDF5) # As I am looser, hdf5 can not find the proper place either
#set(HDF5_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/../libs/include" "${PROJECT_BINARY_DIR}/../libs/include/cpp")
#set(HDF5_LIB_DIRS ${PROJECT_BINARY_DIR}/../libs/lib)
#set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl)
#set(LAPACK_LIBRARIES f2c lapack blas)

find_package(HDF5)
find_package(Boost)
find_package(LAPACK)
find_package(BLAS)
find_package(ZLIB)
find_package(HDF5)
set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl)
#FindBoost FindCxxTest FindBLAS FindHDF5 FindLAPACK FindZLIB FindHDF5


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

add_library(yell-lib src/basic_classes.cpp
                     src/diffuser_core.cpp
                     src/basic_io.cpp
                     src/FormulaParser.cpp
                     src/InputFileParser.cpp
                     src/InputFileParserI.cpp
                     src/InputFileParserII.cpp
                     src/InputFileParserIII.cpp
                     src/model.cpp
                     src/OutputHandler.cpp)

add_executable(yell src/main.cpp)

ADD_DEPENDENCIES(yell yell-lib levmar cctbx ${ZLIB_LIBRARY} ${LAPACK_LIBRARIES})

#MS studio includes tr1 libraries by default. To disable it we use
add_definitions(-DPATH_TO_YELL_SRC="${CMAKE_SOURCE_DIR}/src")

target_link_libraries (yell levmar cctbx ${ZLIB_LIBRARY} ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})

# cxxtest
find_package(CxxTest)
if(CXXTEST_FOUND)
    enable_testing()
    CXXTEST_ADD_TEST(test-yell test_runner.cpp ${CMAKE_SOURCE_DIR}/src/test_diffuser.h)
    target_link_libraries(test-yell yell-lib levmar cctbx ${ZLIB_LIBRARY} ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
endif()