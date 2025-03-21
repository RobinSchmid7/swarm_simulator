cmake_minimum_required(VERSION 3.12)

# set name of the project
project(boids)
set(CMAKE_CXX_STANDARD 17)
enable_testing()

option(BOID_BUILD_GUI "build GUI" ON)

# thirdparty code
add_subdirectory(ext)

# our code
add_subdirectory(src)
