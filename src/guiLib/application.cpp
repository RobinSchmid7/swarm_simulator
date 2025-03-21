#include "application.h"

#include <imgui.h>
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#define NANOVG_GL3_IMPLEMENTATION
#include <nanovg.h>
#include <nanovg_gl.h>

//#define STB_IMAGE_IMPLEMENTATION // defined in nanovg
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#if defined __APPLE__ && !defined RETINA_SCREEN
#define RETINA_SCREEN
#endif

float get_pixel_ratio() {
#ifdef RETINA_SCREEN
    return 2.f;
#endif
    GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    if(monitor == nullptr) throw "Primary monitor not found.";
    float xscale, yscale;
    glfwGetMonitorContentScale(monitor, &xscale, &yscale);
    return xscale;
}

Application::Application(const char *title, int w, int h, std::string iconPath, std::string font)
    : width(w), height(h) {

    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef SINGLE_BUFFER
    glfwWindowHint( GLFW_DOUBLEBUFFER, GL_FALSE ); // turn off framerate limit
#endif

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // fix compilation on OS X
#endif

    // get pixel ratio and adjust width/height
    pixelRatio = get_pixel_ratio();
    this->width = w * pixelRatio;
    this->height = h * pixelRatio;

#ifdef RETINA_SCREEN
    int windowWidth = w, windowHeight = h;
#else
    int windowWidth = width, windowHeight = height;
#endif
    // glfw window creation
    window = glfwCreateWindow(windowWidth, windowHeight, title, nullptr, nullptr);
    if (window == nullptr)
    {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(true);

    // app icon
    if(iconPath != ""){
        GLFWimage image;
        image.pixels = stbi_load(iconPath.c_str(), &image.width, &image.height, nullptr, 4);
        glfwSetWindowIcon(window, 1, &image);
        stbi_image_free(image.pixels);
    }

    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        throw std::runtime_error("Failed to initialize GLAD");
    }

    // Setup Dear ImGui binding
    const char* glsl_version = "#version 150";
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    // io.DisplayFramebufferScale = {2,2};
    ImFontConfig cfg;
   cfg.SizePixels = 40 * pixelRatio;
    ImFont *imFont = io.Fonts->AddFontFromFileTTF(font.c_str(), 14.0f * pixelRatio, &cfg);
////    imFont->DisplayOffset.y = pixelRatio;
    ImGuiStyle& style = ImGui::GetStyle();
#ifdef RETINA_SCREEN
   io.FontGlobalScale = 1.f/pixelRatio;
#else
    style.ScaleAllSizes(pixelRatio);
#endif

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    vg = nvgCreateGL3(NVG_ANTIALIAS);

    setCallbacks();
}

Application::~Application()
{

}

void Application::setCallbacks() {

    glfwSetWindowUserPointer(window, this);

    glfwSetErrorCallback([](int error, const char* description) {
        std::cout << "Error " << error << ": " << description << std::endl;
    });

    glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height){
        auto app = static_cast<Application*>(glfwGetWindowUserPointer(window));
        app->resizeWindow(width, height);
    });

    glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode, int action, int mods){
        if(ImGui::GetIO().WantCaptureKeyboard || ImGui::GetIO().WantTextInput){
            ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
            return;
        }

        auto app = static_cast<Application*>(glfwGetWindowUserPointer(window));
        app->keyDown[key] = (action != GLFW_RELEASE);

        if(action == GLFW_PRESS)
            app->keyPressed(key, mods);

        if(action == GLFW_RELEASE)
            app->keyReleased(key, mods);

    });

    glfwSetMouseButtonCallback(window, [](GLFWwindow* window, int button, int action, int mods){
        if(ImGui::GetIO().WantCaptureMouse){
            ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
            return;
        }

        auto app = static_cast<Application*>(glfwGetWindowUserPointer(window));

        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
#ifndef RETINA_SCREEN
        xpos /= app->pixelRatio;
        ypos /= app->pixelRatio;
#endif
        app->mouseState.onMouseClick(xpos, ypos, button, action, mods);

        if(action == GLFW_PRESS)
            app->mouseButtonPressed(button, mods);

        if(action == GLFW_RELEASE)
            app->mouseButtonReleased(button, mods);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow* window, double xpos, double ypos){
        if(ImGui::GetIO().WantSetMousePos){
            return;
        }

        auto app = static_cast<Application*>(glfwGetWindowUserPointer(window));
#ifndef RETINA_SCREEN
        xpos /= app->pixelRatio;
        ypos /= app->pixelRatio;
#endif
        app->mouseState.onMouseMove(xpos, ypos);

        app->mouseMove(xpos, ypos);
    });

    glfwSetScrollCallback(window, [](GLFWwindow *window, double xoffset, double yoffset){
        if(ImGui::GetIO().WantCaptureMouse){
            ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
            return;
        }

        auto app = static_cast<Application*>(glfwGetWindowUserPointer(window));
        app->scrollWheel( xoffset, yoffset);
    });

    glfwSetDropCallback(window, [](GLFWwindow* window, int count, const char** filenames){
        auto app = static_cast<Application*>(glfwGetWindowUserPointer(window));
        app->drop(count, filenames);
    });
}

void Application::run() {
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        process();

        draw();

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
#ifdef SINGLE_BUFFER
        glFlush();
#else
        glfwSwapBuffers(window);
#endif
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    glfwTerminate();
}

void Application::process() {

}

void Application::draw() {
    glClearColor(clearColor[0], clearColor[1], clearColor[2], 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    // nano vg
    {
        nvgBeginFrame(vg, width/pixelRatio, height/pixelRatio, pixelRatio);
        drawNanoVG();
        nvgEndFrame(vg);
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();


    // ImGui
    {
        using namespace ImGui;
        NewFrame();

        BeginMainMenuBar();
        if(BeginMenu("app")) {
            InputFloat("pixel ratio", &pixelRatio);
            InputInt("width", &width); SameLine();
            InputInt("height", &height);
            float fps = 1.f/deltaTime;
            static bool imguiFps = false;
            LabelText("fps", "fps: %.1f", (imguiFps) ? ImGui::GetIO().Framerate : fps);
            SameLine();
            Checkbox("imgui fps", &imguiFps);
            ImGui::EndMenu();
        }
        EndMainMenuBar();

        drawImGui();

        ImGui::EndFrame();
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
}

void Application::drawImGui()
{

}

void Application::drawNanoVG()
{

}

void Application::resizeWindow(int width, int height) {
    //todo: should the pixel ratio be accounted for here as well?!
//    pixelRatio = get_pixel_ratio();
//    width = pixelRatio * width;
//    height = pixelRatio * height;

    this->width = width;
    this->height = height;

// #ifdef RETINA_SCREEN
//     this->width /= 2;
//     this->height /= 2;
// #endif

    glViewport(0, 0, width, height);
}

bool Application::screenshot(const char *filename) const {
    std::vector< unsigned char > pixels( width * height * 3 );
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0]);
    stbi_flip_vertically_on_write(1);
    return (bool)stbi_write_png(filename, width, height, 3, &pixels[0], 0);
}
