cmake_minimum_required(VERSION 3.5)

add_subdirectory(boids)
if(BOID_BUILD_GUI)
add_subdirectory(guiLib)
add_subdirectory(app)
endif(BOID_BUILD_GUI)
