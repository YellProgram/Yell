cmake_minimum_required (VERSION 3.2)
project (YELL)

include(cmake/hdf5.cmake)



set(CCTBX_PATH lib/cctbx_stubs CACHE PATH "Path to cctbx c++ files")
set(CXXTEST_PATH ../libs/cxxtest-4.4 CACHE PATH "Path to cxxtest")
set(LEVMAR_PATH lib/levmar CACHE PATH "Path to levmar library")

set(BOOST_INCLUDEDIR lib/boost/boost_1_60_0 CACHE PATH "Path to boost library")


set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl szip z)

IF(APPLE)
  set(LAPACK_LIBRARIES "-framework accelerate")
ELSE(APPLE)
IF(UNIX)
  set(LAPACK_LIBRARIES f2c lapack blas dl)
#  find_package(LAPACK)
ENDIF(UNIX)
ENDIF(APPLE)
IF(MINGW)
    set(LAPACK_PATH lib/lapack)
    set(LAPACK_LIBRARIES lapack blas gfortran.a quadmath.a)
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    set(HDF5_LIBRARIES libhdf5_cpp libhdf5_hl_cpp libhdf5 libhdf5_hl libszip libzlib)
ENDIF(MINGW)


include_directories( ${CCTBX_PATH}
        ${LEVMAR_PATH}/levmar2.6
        ${BOOST_INCLUDEDIR}
        ${CXXTEST_PATH}
)
link_directories(${LEVMAR_PATH}/lib ${LEVMAR_PATH}/levmar-build ${LAPACK_PATH}/lib)

#DUMP all include libs
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
    message(STATUS "dir='${dir}'")
endforeach()

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
${CCTBX_PATH}/cctbx/eltbx/neutron.cpp
${CCTBX_PATH}/cctbx/eltbx/xray_scattering/wk1995.cpp
${CCTBX_PATH}/cctbx/eltbx/electron_scattering/peng1996.cpp
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
                     src/OutputHandler.cpp src/exceptions.h)

DownloadHDF5For(yell-lib)

add_executable(yell src/main.cpp)

#ADD_DEPENDENCIES(yell yell-lib cctbx) # ${ZLIB_LIBRARY} ${LAPACK_LIBRARIES})

target_link_libraries(yell yell-lib levmar cctbx  #${HDF5_LIBRARIES}
        ${LAPACK_LIBRARIES})

add_definitions(-DPATH_TO_YELL_SRC="${CMAKE_SOURCE_DIR}/src")

# cxxtest
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_LIST_DIR}/lib/cxxtest)
find_package(CxxTest)
if(CXXTEST_FOUND)
    enable_testing()
    CXXTEST_ADD_TEST(test-yell test_runner.cpp ${CMAKE_SOURCE_DIR}/src/test_diffuser.h)
    target_link_libraries(test-yell yell-lib levmar cctbx ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
endif()
