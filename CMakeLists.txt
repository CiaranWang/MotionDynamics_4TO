cmake_minimum_required(VERSION 3.10)
project(MotionDynamics LANGUAGES CXX)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Executable
add_executable(MotionDynamics_4TO
    src/main.cpp
    src/md.cpp
    src/ReadInput.cpp
)

# Optional: enable warnings (good practice)
if(MSVC)
    target_compile_options(MotionDynamics_4TO PRIVATE /W4 /permissive-)
else()
    target_compile_options(MotionDynamics_4TO PRIVATE -Wall -Wextra -pedantic)
endif()