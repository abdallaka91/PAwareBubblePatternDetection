
cmake_minimum_required(VERSION 3.10)

project(TEST1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Debug build settings
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "-g -O0")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

# Release build settings
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "-O3 -march=native -flto -ffast-math -funroll-loops -fno-math-errno")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

 set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/source)



add_executable(Simul
    # ${SOURCE_DIR}/simulator_decoding.cpp
    ${SOURCE_DIR}/simulator_bubble_gen.cpp
    ${SOURCE_DIR}/GF_tools.cpp
    ${SOURCE_DIR}/init.cpp
    ${SOURCE_DIR}/tools.cpp
    ${SOURCE_DIR}/Decoder_functions.cpp
    ${SOURCE_DIR}/channel.cpp
)

# option(USE_DECODING "Build with simulator_decoding.cpp" ON) # use cmake -B build -DUSE_DECODING=ON

# if(USE_DECODING)
#     set(MAIN_FILE ${SOURCE_DIR}/simulator_decoding.cpp)
#     set(EXEC_NAME "Simul") # Executable name when using simulator_decoding.cpp
# else()
#     set(MAIN_FILE ${SOURCE_DIR}/simulator_bubble_gen.cpp)
#     set(EXEC_NAME "GenBub") # Executable name when using simulator_bubble_gen.cpp
# endif()

# add_executable(${EXEC_NAME}
#     ${MAIN_FILE}
#     ${SOURCE_DIR}/GF_tools.cpp
#     ${SOURCE_DIR}/init.cpp
#     ${SOURCE_DIR}/tools.cpp
#     ${SOURCE_DIR}/Decoder_functions.cpp
#     ${SOURCE_DIR}/channel.cpp
# )
# target_compile_options(${EXEC_NAME} PUBLIC -march=native -mtune=native -ffast-math)

include_directories(${CMAKE_SOURCE_DIR}/include)


