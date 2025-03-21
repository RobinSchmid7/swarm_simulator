#pragma once


#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <stdexcept>
#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <nanovg.h>

#pragma warning( disable : 4244 )

inline float get_pixel_ratio();

// Keep track of the mouse state, positions, buttons pressed, etc...
class MouseState {
public:
    //keep track of the last mouse position
    double lastMouseX = 0, lastMouseY = 0;
    double mouseMoveX = 0, mouseMoveY = 0;

    bool rButtonPressed = false, lButtonPressed = false, mButtonPressed = false;
    bool dragging = false;

    int mods = 0;

public:
    MouseState() {}

    ~MouseState() {}

    void onMouseClick(double xPos, double yPos, int button, int action, int mods) {
        this->mods = mods;

        lastMouseX = xPos;
        lastMouseY = yPos;

        dragging = (action == GLFW_PRESS);

        if (button == GLFW_MOUSE_BUTTON_LEFT)
            lButtonPressed = (action != GLFW_RELEASE);
        if (button == GLFW_MOUSE_BUTTON_MIDDLE)
            mButtonPressed = (action != GLFW_RELEASE);
        if (button == GLFW_MOUSE_BUTTON_RIGHT)
            rButtonPressed = (action != GLFW_RELEASE);
    }

    void onMouseMove(double xPos, double yPos) {
        mouseMoveX = lastMouseX - xPos;
        mouseMoveY = -lastMouseY + yPos;
        lastMouseX = xPos; lastMouseY = yPos;
    }
};

// Sets up a GLFW window, its callbacks and ImGui
class Application {
public:
    GLFWwindow *window;
    int width, height;
    float pixelRatio;

    struct NVGcontext* vg;

//    float clearColor[3] = { 0.1f, 0.1f, 0.1f };
    float clearColor[3] = { .9f, .9f, .9f };

    MouseState mouseState;
    std::map<int, bool> keyDown;

    // timing
    float deltaTime = 0.0f;
    float lastFrame = 0.0f;

public:
    Application(const char *title, int width, int height, std::string iconPath = "", std::string font_path = IMGUI_FONT_FOLDER"/Cousine-Regular.ttf");

    virtual ~Application();

    virtual void setCallbacks();
    virtual void run();
    virtual void process();
    virtual void draw();
    virtual void drawImGui();
    virtual void drawNanoVG();
    virtual void resizeWindow(int width, int height);

    virtual void keyPressed(int key, int mods) { }
    virtual void keyReleased(int key, int mods) { }
    virtual void mouseButtonPressed(int button, int mods) { }
    virtual void mouseButtonReleased(int button, int mods) { }
    virtual void mouseMove(double xpos, double ypos) { }
    virtual void scrollWheel(double xoffset, double yoffset) { }
    virtual void drop(int count, const char** filenames) { }

    bool screenshot(const char *filename) const;
};
