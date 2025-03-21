  
cmake_minimum_required(VERSION 3.5)

project(boids)

add_library(${PROJECT_NAME}
    boids.h
    boids.cpp
)
target_link_libraries(${PROJECT_NAME}
    eigen
)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})