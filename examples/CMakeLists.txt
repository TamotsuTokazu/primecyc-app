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
set( CMAKE_CXX_FLAGS "${OpenFHE_CXX_FLAGS} -Ofast -march=native -mtune=native -fno-schedule-insns -funroll-loops -ffinite-math-only")

add_definitions(-DPARALLEL)

add_executable(tensor tensor.cpp)
target_link_libraries(tensor PRIVATE OPENFHEcore_static OPENFHEbinfhe_static)

add_executable(cyc-reduction cyc-reduction.cpp)
target_link_libraries(cyc-reduction PRIVATE OPENFHEcore_static OPENFHEbinfhe_static)

add_executable(trans trans.cpp)
target_link_libraries(trans PRIVATE OPENFHEcore_static OPENFHEbinfhe_static)

add_executable(benchpoly benchpoly.cpp)
target_link_libraries(benchpoly PRIVATE OPENFHEcore_static OPENFHEbinfhe_static)

add_executable(scratch scratch.cpp)
target_link_libraries(scratch PRIVATE OPENFHEcore_static OPENFHEbinfhe_static)

add_executable(benchtransform benchtransform.cpp)
target_link_libraries(benchtransform PRIVATE OPENFHEcore_static OPENFHEbinfhe_static)