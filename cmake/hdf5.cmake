function(DownloadBitshuffleFor targetName)
    include(ExternalProject)
    message(STATUS "HDF5 bitshuffle plugin will be downloaded and compiled statically")
    ExternalProject_Add(
            download_bitshuffle
            GIT_REPOSITORY https://github.com/kiyo-masui/bitshuffle
            GIT_TAG 0.4.1
            GIT_PROGRESS TRUE
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
    )
    ExternalProject_Get_property(download_bitshuffle SOURCE_DIR)
    set(bitshuffle_source_dir ${SOURCE_DIR})
    if (NOT TARGET CopyShuffleCMakeLists)
        add_custom_target(
                CopyShuffleCMakeLists
                ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/cmake/bitshuffle.cmake ${bitshuffle_source_dir}/CMakeLists.txt
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/bitshuffle.cmake download_bitshuffle
        )
    endif ()
    ExternalProject_Add(
            build_bitshuffle
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${bitshuffle_source_dir}
            INSTALL_DIR ${CMAKE_BINARY_DIR}
            CMAKE_ARGS "${CMAKE_ARGS}"
            CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
            CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}"
            CMAKE_ARGS "-DCMAKE_C_STANDARD_INCLUDE_DIRECTORIES=${HDF5_INSTALL_INCLUDE_DIR}"
            DEPENDS CopyShuffleCMakeLists
            DEPENDS hdf5_local
            DEPENDS download_hdf5
    )
    add_library(bitshuffle STATIC IMPORTED)
    add_dependencies(bitshuffle hdf5_local build_bitshuffle)
    target_include_directories(bitshuffle INTERFACE ${HDF5_INSTALL_INCLUDE_DIR})
    set(bitshuffle_lib_name ${CMAKE_STATIC_LIBRARY_PREFIX}bitshuffle${CMAKE_STATIC_LIBRARY_SUFFIX})
    set_target_properties(bitshuffle PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/${bitshuffle_lib_name})
    add_dependencies(${targetName} bitshuffle)
    target_include_directories(${targetName} PUBLIC ${CMAKE_BINARY_DIR}/include)
    target_link_libraries(${targetName} bitshuffle)
endfunction()

function(DownloadHDF5For targetName)
    find_package(HDF5)
    if (HDF5_FOUND AND NOT CMAKE_TOTAL_STATIC)
        target_include_directories(${targetName} PUBLIC ${HDF5_INCLUDE_DIR})
        target_link_libraries(${targetName} hdf5)
        add_library(hdf5_local STATIC IMPORTED)
        add_library(download_hdf5 STATIC IMPORTED)
    else ()
        message(STATUS "HDF5 library will be downloaded and compiled statically")
        set(HDF5_EXTERNALLY_CONFIGURED 1)
        set(HDF5_PREFIX ${CMAKE_BINARY_DIR}/hdf5)
        set(HDF5_INSTALL_BIN_DIR ${HDF5_PREFIX}/install)
        set(HDF5_INSTALL_LIB_DIR ${HDF5_PREFIX}/install)
        set(HDF5_INSTALL_INCLUDE_DIR ${HDF5_PREFIX}/install/include)
        set(HDF5_INSTALL_DATA_DIR ${HDF5_PREFIX}/install)
        include(ExternalProject)
        ExternalProject_Add(
                download_hdf5
                URL https://git.3lp.cx/packages/hdf5-1.12.0.tar.bz2
                URL_MD5 1fa68c4b11b6ef7a9d72ffa55995f898
                PREFIX ${HDF5_PREFIX}
                INSTALL_DIR ${HDF5_INSTALL_BIN_DIR}
                CMAKE_ARGS "${CMAKE_ARGS}"
                CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                CMAKE_ARGS "-DCMAKE_INSTALL_LIBDIR=${HDF5_INSTALL_LIB_DIR}"
                CMAKE_ARGS "-DCMAKE_INSTALL_BINDIR=${HDF5_INSTALL_BIN_DIR}"
                CMAKE_ARGS "-DCMAKE_INSTALL_INCLUDEDIR=${HDF5_INSTALL_INCLUDE_DIR}"
                CMAKE_ARGS "-DHDF5_INSTALL_BIN_DIR=${HDF5_INSTALL_BIN_DIR}"
                CMAKE_ARGS "-DHDF5_INSTALL_LIB_DIR=${HDF5_INSTALL_LIB_DIR}"
                CMAKE_ARGS "-DHDF5_INSTALL_DATA_DIR=${HDF5_INSTALL_DATA_DIR}"
                CMAKE_ARGS "-DHDF5_INSTALL_INCLUDE_DIR=${HDF5_INSTALL_INCLUDE_DIR}"
                CMAKE_ARGS "-DHDF5_EXTERNALLY_CONFIGURED=${HDF5_EXTERNALLY_CONFIGURED}"
                CMAKE_ARGS "-DHDF5_BUILD_FORTRAN=OFF"
                CMAKE_ARGS "-DBUILD_SHARED_LIBS=OFF"
                CMAKE_ARGS "-DBUILD_TESTING=OFF"
                CMAKE_ARGS "-DHDF5_BUILD_TOOLS=OFF"
                CMAKE_ARGS "-DHDF5_BUILD_EXAMPLES=OFF"
                CMAKE_ARGS "-DHDF5_BUILD_HL_LIB=OFF"
                CMAKE_ARGS "-DHDF5_BUILD_CPP_LIB=ON"
                CMAKE_ARGS "-DDEFAULT_API_VERSION:STRING=v110"
        )
        add_library(hdf5_local STATIC IMPORTED)
        add_dependencies(hdf5_local download_hdf5)
        target_link_libraries(hdf5_local INTERFACE ${CMAKE_DL_LIBS})
        if (WIN32)
            #apparently on github actions with gnu compilers stupid prefix does not appear in hdf5 library
            #set(STUPID_PREFIX "lib")
        endif ()
        if (CMAKE_BUILD_TYPE MATCHES Release)
            set(lib_hdf5_name ${CMAKE_STATIC_LIBRARY_PREFIX}${STUPID_PREFIX}hdf5${CMAKE_STATIC_LIBRARY_SUFFIX})
            set(lib_hdf5_name_cpp ${CMAKE_STATIC_LIBRARY_PREFIX}${STUPID_PREFIX}hdf5_cpp${CMAKE_STATIC_LIBRARY_SUFFIX})
        else ()
            set(lib_hdf5_name ${CMAKE_STATIC_LIBRARY_PREFIX}${STUPID_PREFIX}hdf5_debug${CMAKE_STATIC_LIBRARY_SUFFIX})
            set(lib_hdf5_name_cpp ${CMAKE_STATIC_LIBRARY_PREFIX}${STUPID_PREFIX}hdf5_cpp_debug${CMAKE_STATIC_LIBRARY_SUFFIX})
        endif ()
        set_target_properties(hdf5_local PROPERTIES IMPORTED_LOCATION ${HDF5_INSTALL_LIB_DIR}/${lib_hdf5_name})
        add_library(hdf5_local_cpp STATIC IMPORTED)
        add_dependencies(hdf5_local_cpp download_hdf5 hdf5_local)
        target_link_libraries(hdf5_local_cpp INTERFACE ${CMAKE_DL_LIBS} hdf5_local)
        set_target_properties(hdf5_local_cpp PROPERTIES IMPORTED_LOCATION ${HDF5_INSTALL_LIB_DIR}/${lib_hdf5_name_cpp})
        add_dependencies(${targetName} download_hdf5)
        target_include_directories(${targetName} PUBLIC ${HDF5_INSTALL_INCLUDE_DIR})
        target_link_libraries(${targetName} hdf5_local hdf5_local_cpp)
        file(MAKE_DIRECTORY ${HDF5_INSTALL_INCLUDE_DIR})
    endif ()
    DownloadBitshuffleFor(${targetName})
endfunction()
