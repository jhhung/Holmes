cmake_minimum_required(VERSION 3.16)
project(Sherloc LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/json)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog)

find_package(Boost CONFIG REQUIRED COMPONENTS  
    exception
    serialization
    program_options
    iostreams
)

find_package(OpenMP REQUIRED)

# htslib
set(HTSLIB_ROOT ${CMAKE_CURRENT_LIST_DIR}/lib/htslib)
find_package(BZip2)
if(NOT BZIP2_FOUND)
  set(HTSLIB_DISABLE --disable-bz2)
endif()
find_package(LibLZMA)
if(NOT LIBLZMA_FOUND)
  set(HTSLIB_DISABLE ${HTSLIB_DISABLE} --disable-lzma)
endif()
ExternalProject_add(
    htslib
    SOURCE_DIR ${HTSLIB_ROOT}
    CONFIGURE_COMMAND cd ${HTSLIB_ROOT} && autoreconf -i && ./configure ${HTSLIB_DISABLE}
    BUILD_COMMAND cd ${HTSLIB_ROOT} && make -j
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${HTSLIB_ROOT}/libhts.so
  )
add_library(hts SHARED IMPORTED)
set_target_properties(hts PROPERTIES IMPORTED_LOCATION ${HTSLIB_ROOT}/libhts.so)
include_directories(${HTSLIB_ROOT})

add_library(holmes INTERFACE)
add_dependencies(holmes htslib)
target_include_directories(holmes
    INTERFACE ${Boost_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include ${HTSLIB_ROOT})
target_link_libraries(holmes
    INTERFACE ${Boost_LIBRARIES} nlohmann_json::nlohmann_json spdlog::spdlog hts OpenMP::OpenMP_CXX)
target_compile_definitions(holmes
    INTERFACE HOLMES_CONFIG_PATH="${CMAKE_CURRENT_LIST_DIR}/config")

set(Sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Sherloc/app/archive_compressor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Sherloc/app/database_builder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Sherloc/app/disease_json_operation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Sherloc/app/sherloc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Sherloc/app/vep_cache_builder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Sherloc/app/vep_cache_query.cpp
)
    
foreach(source ${Sources})
    get_filename_component(executable ${source} NAME_WE)
    add_executable(${executable} ${source})

    target_link_libraries(${executable} PRIVATE holmes)

    set_target_properties(${executable}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_HOME_DIRECTORY}/bin"
    )
endforeach()

# build test
option(HOLMES_TESTS "Build the tests" ON)
if(HOLMES_TESTS)
  include(CTest)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tests)
endif()