cmake_minimum_required (VERSION 3.10)
project (YELL)

include(cmake/hdf5.cmake)

set (CMAKE_CXX_STANDARD 20)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CCTBX_PATH lib/cctbx_stubs CACHE PATH "Path to cctbx c++ files")
set(CXXTEST_PATH lib/cxxtest CACHE PATH "Path to cxxtestt")
set(LAPACK_PATH lib/lapack/lapack-build CACHE PATH "Path to lapack and blas")

set(CERES_PATH ../../libs/ceres-solver)

set(GLOG_PATH ../../libs/glog-0.7.1)
set(EIGEN_PATH ../../libs/eigen)

set(BOOST_INCLUDEDIR ${CMAKE_SOURCE_DIR}/../libs/boost_1_88_0 CACHE PATH "Path to boost library" FORCE)

set(HDF5_LIBRARIES hdf5_cpp hdf5_hl_cpp hdf5 hdf5_hl szip z)

include(ExternalProject)
include(FetchContent)
find_package(Threads REQUIRED)

# pocketfft — header-only, single file.  FetchContent downloads the repo and
# we expose its source dir so #include "pocketfft_hdronly.h" works everywhere.
# pocketfft — single-header, header-only FFT library.
# Use the copy in lib/ if present (committed or manually placed), otherwise
# download just the header from the canonical GitLab source.
set(POCKETFFT_HEADER_NAME pocketfft_hdronly.h)
set(POCKETFFT_LOCAL_DIR   ${CMAKE_SOURCE_DIR}/lib)
set(POCKETFFT_DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/pocketfft-download)
if(EXISTS "${POCKETFFT_LOCAL_DIR}/${POCKETFFT_HEADER_NAME}")
    set(POCKETFFT_INCLUDE_DIR ${POCKETFFT_LOCAL_DIR})
    message(STATUS "pocketfft: using local copy at ${POCKETFFT_INCLUDE_DIR}")
else()
    file(MAKE_DIRECTORY ${POCKETFFT_DOWNLOAD_DIR})
    if(NOT EXISTS "${POCKETFFT_DOWNLOAD_DIR}/${POCKETFFT_HEADER_NAME}")
        message(STATUS "pocketfft: downloading header...")
        file(DOWNLOAD
            "https://gitlab.mpcdf.mpg.de/mtr/pocketfft/-/raw/master/pocketfft_hdronly.h"
            "${POCKETFFT_DOWNLOAD_DIR}/${POCKETFFT_HEADER_NAME}"
            STATUS download_status
            TLS_VERIFY ON)
        list(GET download_status 0 download_code)
        if(NOT download_code EQUAL 0)
            message(FATAL_ERROR "pocketfft download failed: ${download_status}\n"
                "Place pocketfft_hdronly.h in ${POCKETFFT_LOCAL_DIR} manually.")
        endif()
        message(STATUS "pocketfft: downloaded to ${POCKETFFT_DOWNLOAD_DIR}")
    endif()
    set(POCKETFFT_INCLUDE_DIR ${POCKETFFT_DOWNLOAD_DIR})
endif()

ExternalProject_Add(
        eigen
        URL  https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        #SOURCE_DIR ${EIGEN_PATH}
        #BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen"
        #INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/libs/eigen-install"
        CMAKE_ARGS
        #-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/eigen-install
)

ExternalProject_Add(
        glog-lib
        URL  https://github.com/google/glog/archive/refs/tags/v0.7.1.tar.gz
        #SOURCE_DIR ${GLOG_PATH}
        #BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen"
        #INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/libs/eigen-install"
        CMAKE_ARGS
        #-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_BUILD_TYPE=RELEASE #${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_STANDARD=14
        -DBUILD_SHARED_LIBS=OFF
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/glog-install
)

ExternalProject_Add(
        extern_gflags
        URL https://github.com/gflags/gflags/archive/v2.2.2.zip
        URL_HASH SHA256=19713a36c9f32b33df59d1c79b4958434cb005b5b47dc5400a7a4b078111d9b5
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${LIB_BUILD_TYPE}
        -DGFLAGS_BUILD_STATIC_LIBS=ON
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=OFF
        -DGFLAGS_BUILD_SHARED_LIBS=OFF
        -DGFLAGS_BUILD_gflags_LIB=ON
        -DGFLAGS_BUILD_gflags_nothreads_LIB=ON
        -DGFLAGS_NAMESPACE=google
        -DGFLAGS_BUILD_TESTING=OFF
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 #fix that maybe?
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/gflags-install
)

ExternalProject_Add(ceres-solver
        PREFIX CeresSolver
        # googlesource.com generates non-deterministic tarballs; use the GitHub mirror instead.
        URL https://github.com/ceres-solver/ceres-solver/archive/76af132d02d055db2af1322e1314b03c493c6803.tar.gz
        URL_HASH SHA256=4a1862a6dda913868aeeda26058d14bca92f90235cfc118fc1ffc36beb74c243
        #SOURCE_DIR ${CERES_PATH}
        #BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libs/ceres-solver
        #INSTALL_DIR ${CMAKE_BINARY_DIR}/ceres-solver
        LOG_DIR ${CMAKE_BINARY_DIR}/ceres-solver/log
        LOG_CONFIGURE ON
        CMAKE_ARGS
        -DCERES_THREADING_MODEL=OPENMP
        -DCMAKE_CXX_STANDARD=17
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/ceres-solver
        -DUSE_CUDA=OFF

        -DBUILD_BENCHMARKS=OFF
        -DBUILD_DOCUMENTATION=OFF
        -DBUILD_EXAMPLES=OFF
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_TESTING=OFF

        -DCXSPARSE=OFF
        -DEIGENSPARSE=ON

        -DEIGEN_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}
        -DEigen3_DIR=${CMAKE_BINARY_DIR}/eigen-prefix/share/eigen3/cmake #here should be eigen-install
        -Dglog_DIR=${CMAKE_BINARY_DIR}/glog-install
        -DSUITESPARSE=ON
        -DLAPACK=ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CERES_LIB ceres-debug)
else()
    set(CERES_LIB ceres)
endif()



IF(APPLE)
  set(LAPACK_LIBRARIES "-framework accelerate")
ELSE(APPLE)
IF(UNIX)
  #set(LAPACK_LIBRARIES  lapack blas gfortran) #f2c dl

  set(LAPACK_LIBRARIES lapack blas gfortran.a quadmath.a)
  if(ENABLE_TSAN)
    # TSan requires dynamic linking — skip -Wl,-Bstatic so the sanitizer runtime loads correctly.
    set(CMAKE_EXE_LINKER_FLAGS "-fsanitize=thread")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -g")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fsanitize=thread -g")
  else()
    # -Wl,-Bstatic makes everything after it static.
    # We must ensure libpthread is linked dynamically, otherwise std::thread may be broken.
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -Wl,-Bstatic")
    set(CMAKE_THREAD_LIBS_INIT "-Wl,-Bdynamic -lpthread -Wl,-Bstatic")
  endif()
#  find_package(LAPACK)
ENDIF(UNIX)
ENDIF(APPLE)
IF(MINGW) #Windows on github actions
#    set(LAPACK_LIBRARIES lapack blas gfortran)
#    set(LAPACK_PATH lib/lapack)
    set(LAPACK_LIBRARIES lapack blas gfortran.a quadmath.a)
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -Wl,-Bstatic")
ENDIF(MINGW)
include_directories( ${CCTBX_PATH}
        ${BOOST_INCLUDEDIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
        ${POCKETFFT_INCLUDE_DIR}
        src/
        ${CXXTEST_PATH}
)
link_directories(${LAPACK_PATH}/lib)

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
                     src/model.cpp
                     src/OutputHandler.cpp src/exceptions.h
        src/Grid.cpp
        src/Grid.h
        src/IntensityMap.cpp
        src/IntensityMap.h
        src/LaueSymmetry.cpp
        src/LaueSymmetry.h
)

add_library(yell-minimizer
        src/CeresMinimizer.cpp
        src/CeresMinimizer.h )
set_target_properties(yell-minimizer
        PROPERTIES
        CXX_STANDARD 20
)

target_include_directories(yell-minimizer PUBLIC ${CMAKE_BINARY_DIR}/ceres-solver/include
        ${CMAKE_BINARY_DIR}/eigen-install/include/eigen3
        ${CMAKE_BINARY_DIR}/glog-install/include)
target_link_libraries(yell-minimizer PUBLIC Threads::Threads)

target_include_directories(yell-lib PUBLIC ${CMAKE_BINARY_DIR}/eigen-install/include/eigen3)
target_link_libraries(yell-lib PUBLIC Threads::Threads)
add_dependencies(yell-lib eigen)

DownloadHDF5For(yell-lib )

link_directories(${CMAKE_BINARY_DIR}/gflags-install/lib/ ${CMAKE_BINARY_DIR}/glog-install/lib/ ${CMAKE_BINARY_DIR}/ceres-solver/lib)

add_definitions(-DGLOG_USE_GLOG_EXPORT)
# Disable OpenMP parallelisation inside the cctbx 3D FFT.
# We provide outer parallelism (one thread per pair peak); the inner OMP team
# causes O(n_threads * omp_team_size) thread contention and completely
# serialises the forward/derivative calculation.  basic_classes.h and Grid.h
# include complex_to_complex_3d.h before IntensityMap.h has a chance to set
# the macro, so it must be a global compile definition.
add_definitions(-DSCITBX_FFTPACK_COMPLEX_TO_COMPLEX_3D_NO_PRAGMA_OMP)

add_executable(yell src/main.cpp)
add_dependencies(ceres-solver extern_gflags glog-lib)
add_dependencies(yell-minimizer eigen ceres-solver) #

#ADD_DEPENDENCIES(yell yell-lib cctbx) # ${ZLIB_LIBRARY} ${LAPACK_LIBRARIES})

target_link_libraries(yell PUBLIC
        yell-lib yell-minimizer cctbx ${CERES_LIB} glog gflags ${LAPACK_LIBRARIES})

add_definitions(-DPATH_TO_YELL_SRC="${CMAKE_SOURCE_DIR}/src")

# cxxtest
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CXXTEST_PATH})

find_package(CxxTest)
if(CXXTEST_FOUND)
    enable_testing()
    CXXTEST_ADD_TEST(test-yell test_runner.cpp
        ${CMAKE_SOURCE_DIR}/src/test_diffuser.h)
    target_link_libraries(test-yell yell-lib cctbx ${LAPACK_LIBRARIES})
    target_include_directories(test-yell PUBLIC ${CMAKE_BINARY_DIR}/eigen-install/include/eigen3)
    add_dependencies(test-yell eigen)
endif()

# Google Test
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    URL_HASH SHA256=8ad598c73ad796e0d8280b082cebd82a630d73e73cd3c70057938a6501bba5d7
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test-yell-gtest
    ${CMAKE_SOURCE_DIR}/src/test_expr_gtest.cpp
    ${CMAKE_SOURCE_DIR}/src/test_anharmonic_gtest.cpp)
target_link_libraries(test-yell-gtest
    yell-lib cctbx ${LAPACK_LIBRARIES}
    GTest::gtest_main)
target_include_directories(test-yell-gtest PUBLIC
    ${CMAKE_BINARY_DIR}/eigen-install/include/eigen3)
add_dependencies(test-yell-gtest eigen)
include(GoogleTest)
gtest_discover_tests(test-yell-gtest)

# Parallelism timing tests (slow — run separately, not part of normal CI suite)
add_executable(test-yell-parallel ${CMAKE_SOURCE_DIR}/src/test_parallel_gtest.cpp)
target_link_libraries(test-yell-parallel
    yell-lib yell-minimizer cctbx ${CERES_LIB} glog gflags ${LAPACK_LIBRARIES}
    GTest::gtest_main)
target_include_directories(test-yell-parallel PUBLIC
    ${CMAKE_BINARY_DIR}/eigen-install/include/eigen3)
add_dependencies(test-yell-parallel eigen)
