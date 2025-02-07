cmake_minimum_required(VERSION 3.16)
project(holmes-test)

add_executable(holmes-test
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/coverage.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/k.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/clinvar.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/dvd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/vcf.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/vep.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/hts.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/DB/uniprot.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/Attr/allele.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/app/sherloc/predict.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sherloc/app/sherloc/sherloc_parameter.cpp
)

set_target_properties(holmes-test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_HOME_DIRECTORY}/bin"
)

target_link_libraries(holmes-test
    PRIVATE holmes)
target_include_directories(holmes-test
    PRIVATE ${CMAKE_HOME_DIRECTORY}/lib)
target_compile_definitions(holmes-test
    PRIVATE DATA_PATH="${CMAKE_CURRENT_LIST_DIR}/data")
target_compile_options(holmes-test PRIVATE -Wno-ignored-attributes)

add_test(NAME Sherloc COMMAND holmes-test)
