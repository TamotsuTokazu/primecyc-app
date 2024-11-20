    cmake_minimum_required(VERSION 3.10)
    project(PrimeCycApp VERSION 1.0 LANGUAGES CXX)

    # Set C++ standard
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
    
    add_subdirectory(lib/openfhe-development)
    add_subdirectory(lib/openfhe-development/third-party/google-test)
    add_subdirectory(lib/openfhe-development/third-party/google-benchmark)

    # Define the sources directory
    add_subdirectory(src)
    add_subdirectory(examples)
    add_subdirectory(benchmark)

    # Optionally add tests
    enable_testing()
    add_subdirectory(tests)