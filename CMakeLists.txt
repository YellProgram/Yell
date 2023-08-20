cmake_minimum_required (VERSION 3.2)
project (YELL)

# Directories for dependencies
set(LEVMAR_PATH lib/levmar CACHE PATH "Path to levmar library")
set(CCTBX_PATH lib/cctbx_stubs CACHE PATH "Path to cctbx c++ files")
set(CXXTEST_PATH lib/cxxtest CACHE PATH "Path to cxxtest")

# HDF5 configurations
set(HDF5_INCLUDE_DIRS lib/hdf5/include)
set(HDF5_LIB_DIRS lib/hdf5/lib)
set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl) 

# Boost setup
set(BOOST_ROOT "${CMAKE_SOURCE_DIR}/lib/boost")
set(Boost_NO_SYSTEM_PATHS TRUE)
find_package(Boost REQUIRED)  # Add COMPONENTS <library-name> if you need specific libraries

# Platform-specific libraries and configurations
IF(APPLE)
  set(LAPACK_LIBRARIES "-framework accelerate")
ELSEIF(UNIX)
  set(LAPACK_LIBRARIES f2c lapack blas dl)
ELSEIF(MINGW)
  set(LAPACK_PATH lib/lapack)
  set(LAPACK_LIBRARIES lapack blas gfortran.a quadmath.a)
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  set(HDF5_LIBRARIES libhdf5_cpp libhdf5_hl_cpp libhdf5 libhdf5_hl libszip libzlib)
ENDIF()

include_directories(
    ${CCTBX_PATH}
    ${LEVMAR_PATH}/include
    ${HDF5_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${CXXTEST_PATH}
    ${HDF5_PATH}/c++/src
    ${HDF5_PATH}/src
)
link_directories(${HDF5_LIB_DIRS} ${LEVMAR_PATH}/lib ${LAPACK_PATH}/lib)

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
    ${CCTBX_PATH}/cctbx/eltbx/neutron.cpp
    ${CCTBX_PATH}/cctbx/eltbx/xray_scattering/wk1995.cpp
    ${CCTBX_PATH}/cctbx/eltbx/basic.cpp
)

add_library(yell-lib 
    src/basic_classes.cpp
    src/diffuser_core.cpp
    src/basic_io.cpp
    src/FormulaParser.cpp
    src/InputFileParser.cpp
    src/InputFileParserI.cpp
    src/InputFileParserII.cpp
    src/InputFileParserIII.cpp
    src/model.cpp
    src/OutputHandler.cpp 
    src/exceptions.h
)

add_executable(yell src/main.cpp)
target_link_libraries(yell yell-lib levmar cctbx ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})

add_definitions(-DPATH_TO_YELL_SRC="${CMAKE_SOURCE_DIR}/src")

# cxxtest
find_package(CxxTest)
if(CXXTEST_FOUND)
    enable_testing()
    CXXTEST_ADD_TEST(test-yell test_runner.cpp ${CMAKE_SOURCE_DIR}/src/test_diffuser.h)
    target_link_libraries(test-yell yell-lib levmar cctbx ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
endif()
