cmake_minimum_required (VERSION 2.8.1)
project (YELL)

set(CCTBX_PATH /Users/arkadiy/ag/yell/libs/cctbx/modules/cctbx_project CACHE PATH "Path to cctbx c++ files")
set(LEVMAR_PATH lib/levmar CACHE PATH "Path to levmar library")
set(SCITBX_PATH /Users/arkadiy/ag/yell/libs/cctbx/build/include CACHE PATH "Our library depends on some of the cctbx files which are not readily available in the cctbx source bundle. Rather, they are compiled in one of the stages of cctbx build process. This variable should link to the folder ... which will be created after the command '/your/choice/bin/python ../cctbx_sources/libtbx/configure.py mmtbx' for more information on the installation process see http://cctbx.sourceforge.net/current_cvs/installation.html#manually-building-from-sources-under-unix")
set(CXXTEST_PATH lib/cxxtest CACHE PATH "Path to cxxtest")

#set(SZIP_PATH lib/szip) /Users/arkadiy/ag/yell/libs/cctbx
#set(ZLIB_PATH lib/zlib)



set(BOOST_INCLUDEDIR lib/boost)

#Unfortunately find_package(HDF5) by default finds dynamic libraries
set(HDF5_INCLUDE_DIRS lib/hdf5/include)

set(HDF5_LIB_DIRS lib/hdf5/lib)
set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl szip z)
IF(APPLE)
set(LAPACK_LIBRARIES "-framework accelerate")
ELSE(APPLE)
set(LAPACK_LIBRARIES f2c lapack blas)
ENDIF(APPLE)


#find_package(HDF5)
#find_package(Boost)
#find_package(LAPACK)
#find_package(BLAS)
#find_package(ZLIB)


#FindBoost FindCxxTest FindBLAS FindHDF5 FindLAPACK FindZLIB FindHDF5

include_directories(${SCITBX_PATH} ${CCTBX_PATH} ${LEVMAR_PATH}/include ${HDF5_INCLUDE_DIRS} ${BOOST_INCLUDEDIR} ${CXXTEST_PATH} ${HDF5_PATH}/c++/src ${HDF5_PATH}/src) #THe order makes difference. SCITBX_PATH shoud come before CCTBX_PATH !!
link_directories(${HDF5_LIB_DIRS} ${LAPACK_DIRS} ${LEVMAR_PATH}/lib)

# Manualy get what we need from cctbx
add_library(cctbx 
${CCTBX_PATH}/cctbx/sgtbx/change_of_basis_op.cpp
${CCTBX_PATH}/cctbx/sgtbx/rot_mx.cpp
${CCTBX_PATH}/cctbx/sgtbx/rt_mx.cpp
${CCTBX_PATH}/cctbx/sgtbx/tr_vec.cpp
${CCTBX_PATH}/cctbx/sgtbx/utils.cpp
${CCTBX_PATH}/cctbx/uctbx/spoil_optimization.cpp
${CCTBX_PATH}/cctbx/uctbx/crystal_orientation.cpp
${CCTBX_PATH}/cctbx/uctbx/uctbx.cpp
${CCTBX_PATH}/omptbx/stubs.cpp
${CCTBX_PATH}/cctbx/eltbx/xray_scattering/wk1995.cpp
${CCTBX_PATH}/cctbx/eltbx/basic.cpp
)

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

ADD_DEPENDENCIES(yell yell-lib cctbx) # ${ZLIB_LIBRARY} ${LAPACK_LIBRARIES})



target_link_libraries(yell yell-lib levmar cctbx  ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})

add_definitions(-DPATH_TO_YELL_SRC="${CMAKE_SOURCE_DIR}/src")

# cxxtest
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_LIST_DIR}/lib/cxxtest)
find_package(CxxTest)
if(CXXTEST_FOUND)
    enable_testing()
    CXXTEST_ADD_TEST(test-yell test_runner.cpp ${CMAKE_SOURCE_DIR}/src/test_diffuser.h)
    target_link_libraries(test-yell yell-lib levmar cctbx ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
endif()