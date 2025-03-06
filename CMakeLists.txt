
cmake_minimum_required(VERSION 3.10)

project(TEST1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set debug flags (uncomment for Debug mode)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto -ffast-math -funroll-loops -fno-math-errno")
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)



# Set source directory
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/source)

# Add your source files
add_executable(Simul
    ${SOURCE_DIR}/main.cpp
    ${SOURCE_DIR}/GF_tools.cpp
    ${SOURCE_DIR}/init.cpp
    ${SOURCE_DIR}/tools.cpp
    ${SOURCE_DIR}/Decoder_functions.cpp
)

target_compile_options(Simul PUBLIC -march=native -mtune=native -ffast-math)

# Include header files
include_directories(${CMAKE_SOURCE_DIR}/include)

