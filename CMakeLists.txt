cmake_minimum_required (VERSION 3.2)
project (YELL)
# Define Paths
set(CMAKE_SOURCE_DIR ${GITHUB_WORKSPACE} CACHE PATH "Source directory root")
set(LEVMAR_PATH ${CMAKE_SOURCE_DIR}/lib/levmar)
set(CCTBX_PATH ${CMAKE_SOURCE_DIR}/lib/cctbx)
set(CXXTEST_PATH ${CMAKE_SOURCE_DIR}/lib/cxxtest)
set(BOOST_INCLUDEDIR ${CMAKE_SOURCE_DIR}/lib/boost)
set(HDF5_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/hdf5/include)
set(HDF5_LIB_DIRS ${CMAKE_SOURCE_DIR}/lib/hdf5/lib)

# Define Libraries
set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl)

# LAPACK Configuration
if(APPLE)
  set(LAPACK_LIBRARIES "-framework accelerate")
elseif(UNIX)
  set(LAPACK_LIBRARIES f2c lapack blas dl)
elseif(MINGW)
  set(LAPACK_PATH ${CMAKE_SOURCE_DIR}/lib/lapack)
  set(LAPACK_LIBRARIES lapack blas gfortran.a quadmath.a)
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  set(HDF5_LIBRARIES libhdf5_cpp libhdf5_hl_cpp libhdf5 libhdf5_hl libszip libzlib)
endif()

# Include and Link Directories
include_directories(
  ${CCTBX_PATH}
  ${LEVMAR_PATH}/include
  ${HDF5_INCLUDE_DIRS}
  ${BOOST_INCLUDEDIR}
  ${CXXTEST_PATH}
)

link_directories(${HDF5_LIB_DIRS} ${LEVMAR_PATH}/lib)

# Add Libraries
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

# Link Libraries
target_link_libraries(yell yell-lib levmar cctbx ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})

add_definitions(-DPATH_TO_YELL_SRC="${CMAKE_SOURCE_DIR}/src")

# CxxTest Configuration
find_package(CxxTest)
if(CXXTEST_FOUND)
  enable_testing()
  CXXTEST_ADD_TEST(test-yell test_runner.cpp ${CMAKE_SOURCE_DIR}/src/test_diffuser.h)
  target_link_libraries(test-yell yell-lib levmar cctbx ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

add_library(cctbx cctbx_stubs/cctbx/sgtbx/change_of_basis_op.cpp)