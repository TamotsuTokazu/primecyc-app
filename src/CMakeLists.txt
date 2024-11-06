include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(OpenFHE CONFIG REQUIRED)
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
    message(STATUS "OpenFHE_EXE_LINKER_FLAGS: ${OpenFHE_EXE_LINKER_FLAGS}")
    message(STATUS "OpenFHE_CXX_FLAGS: ${OpenFHE_CXX_FLAGS}")
    message(STATUS "OPENMP_INCLUDES: ${OPENMP_INCLUDES}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

list(APPEND OPENFHE_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/lib/openfhe-development/src/core/include")
list(APPEND OPENFHE_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/lib/openfhe-development/src/binfhe/include")
list(APPEND OPENFHE_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/lib/openfhe-development/src/pke/include")
list(APPEND OPENFHE_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/lib/openfhe-development/third-party/cereal/include")
list(APPEND OPENFHE_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/lib/openfhe-development/src/core")

include_directories(${OPENFHE_INCLUDE_DIRS})
set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
set( CMAKE_CXX_FLAGS "${OpenFHE_CXX_FLAGS} -O3")

add_definitions(-DPARALLEL)

add_executable(bdf17 bdf17.cpp)
target_link_libraries(bdf17 PRIVATE OPENFHEcore OPENFHEbinfhe)

add_executable(bdf17-new bdf17-new.cpp)
target_link_libraries(bdf17-new PRIVATE OPENFHEcore OPENFHEbinfhe)