cmake_minimum_required(VERSION 3.5)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(guiLib
    ${sources}
)
target_link_libraries(guiLib
    ${OPENGL_LIBRARIES} # open gl libs already linked by gl3w, linking again ...
    glad
    glfw
    nanovg
    imgui
    stb_image
    glm
)
target_include_directories(guiLib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${NANOVG_SRC_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../ext/glfw/include
)
target_compile_definitions(guiLib PUBLIC BOID_ASSETS_FOLDER="${CMAKE_CURRENT_LIST_DIR}/assets")
target_compile_definitions(guiLib PUBLIC SHADER_FOLDER="${CMAKE_CURRENT_LIST_DIR}/../gui/shaders")
target_compile_definitions(guiLib PUBLIC IMGUI_FONT_FOLDER=${BOID_IMGUI_FONT_FOLDER})
target_compile_definitions(guiLib PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD)
