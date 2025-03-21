#include "application.h"
#include <imgui.h>

#include <iostream>
#include <math.h>
#include <deque>
#include <chrono>

#include "../boids/boids.h"

#define T float
#define dim 2
// Helps to use different types and dimensions

class TestApp : public Application
{
#define COLOR_RED    nvgRGBA(220,50,50,255)
#define COLOR_BLUE  nvgRGBA(50,50,220,255)
#define COLOR_GREEN   nvgRGBA(50,220,50,255)
#define COLOR_LEADER nvgRGBA(50,220,50,255)
#define COLOR_BLACK nvgRGBA(0,0,0,255) // Black

typedef Matrix<T, Eigen::Dynamic, 1> VectorXT;
typedef Matrix<bool, Eigen::Dynamic, 1> VectorXBool;
typedef Matrix<T, dim, Eigen::Dynamic> TVStack;
typedef Vector<T, dim> TV;

public:

    TestApp(int w, int h, const char * title) : Application(title, w, h) {

        ImGui::StyleColorsClassic();

        const char* name = IMGUI_FONT_FOLDER"/Cousine-Regular.ttf";
        nvgCreateFont(vg, "sans", name);
    }

    void process() override {
        std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
        if(std::chrono::duration_cast<std::chrono::microseconds>(now-lastFrame).count() >= 10./60. * 1.e6)
        {
            if(keyDown[GLFW_KEY_R] || neverRun || newMethod != currentMethod) {
                currentMethod = newMethod;
                boids.initializePositionsAndGroups(currentMethod, n);
                initializeObstacles(obstacles);

                leaderPos_01 = boids.getPositions().col(0);
                leaderPos_screen = shift_01_to_screen(leaderPos_01, scale, width, height);

                neverRun = false;
            }
            if(keyDown[GLFW_KEY_SPACE])
                boids.pause();
            if(keyDown[GLFW_KEY_ESCAPE])
                exit(0);
            lastFrame = now;

            if(draggingCircle)
            {
                leaderPos_screen(0) = mouseState.lastMouseX - draggingCircleOffset(0);
                leaderPos_screen(1) = mouseState.lastMouseY - draggingCircleOffset(1);
            }
        }
    }

    void drawImGui() override {

        using namespace ImGui;

        // Method parameters
        const char* methodNames[] = {"FreeFall", "Circle", "Cohesion", "Alignment", "Separation", "Collision Avoidance", "Leader", "Collaborative"};
        SetNextWindowPos(ImVec2(0,0), ImGuiCond_Always);
        SetNextWindowSize(ImVec2(500,150), ImGuiCond_Always);
        Begin("Method Parameters");
        Combo("Boids behavior", (int*)&newMethod, methodNames, std::size(methodNames));
        const char* schemeNames[] = {"Explicit Euler", "Sympletic Euler", "Explicit Midpoint"};
        Combo("Integration scheme", (int*)&currentScheme, schemeNames, std::size(schemeNames));
        const char* velocityInitializationNames[] = {"Zero Velocity", "Random Velocity", "Circular Velocity"};
        Combo("Initial velocity", (int*)&currentInitialVelocity, velocityInitializationNames, std::size(velocityInitializationNames));
        SliderInt("# initial boids", &n, 1, 100, "%i");
        SliderFloat("Timestep", &timestep, 0.001f, 1.0f, "%.5f");
        End();

        // Simulation parameters
        SetNextWindowPos(ImVec2(500,0), ImGuiCond_Always);
        SetNextWindowSize(ImVec2(500,230), ImGuiCond_Always);
        Begin("Simulation Parameters");
        SliderFloat("Cohesion radius", &cohesionRadius, 0.0f, 1.0f, "%.5f");
        SliderFloat("Alignment radius", &alignmentRadius, 0.0f, 1.0f, "%.5f");
        SliderFloat("Separation radius", &separationRadius, 0.0f, 1.0f, "%.5f");
        SliderFloat("Boid mass", &boidMass, 0.01f, 10.0f, "%.5f");
        SliderFloat("Damping", &damping, -1.0f, 1.0f, "%.5f");
        SliderFloat("Friction", &friction, 0.0f, 10.0f, "%.5f");
        SliderFloat("Leader observation", &leaderObservationDistance, 0.0f, 1.0f, "%.5f");
        SliderFloat("Leader approach stop", &stopApproachLeaderRadius, 0.0f, 0.1f, "%.5f");
        End();

        // Collaboration parameters
        SetNextWindowPos(ImVec2(0, 780), ImGuiCond_Always);
        SetNextWindowSize(ImVec2(600,220), ImGuiCond_Always);
        Begin("Collaboration Parameters");
        const char* strategyNames[] = {"No strategy", "Go to closest of own group", "Go to center and stay", "Form groups of three and chase closest enemy", "Form groups of three and move to bounce off point"};
        Combo("Strategy red", (int*)&collaborationStrategyRed, strategyNames, std::size(strategyNames));
        Combo("Strategy blue", (int*)&collaborationStrategyBlue, strategyNames, std::size(strategyNames));
        SliderFloat("Breeding radius", &breedingRadius, 0.0f, 1.0f, "%.5f");
        SliderFloat("Elimination radius", &eliminationRadius, 0.0f, 1.0f, "%.5f");
        SliderInt("Max # boids per iteration", &maxNumNewBreeded, 0, 5, "%i");
        SliderInt("Breeding pause steps", &breedingPauseSteps, 10, 1000, "%i");
        SliderFloat("Hunter group size", &hunterGroupSize, 0.01f, 0.5f, "%.5f");
        SliderInt("Max # of boids", &maxNumParticles, 1, 1000, "%i");
        End();
    }

    void drawNanoVG() override {
        boids.setBehaviourParam(currentMethod, currentScheme, currentInitialVelocity,timestep, boidMass, cohesionRadius, alignmentRadius, separationRadius,
                                leaderObservationDistance, stopApproachLeaderRadius, friction, damping);
        boids.setCollaborationParam(breedingRadius, eliminationRadius, maxNumParticles, maxNumNewBreeded, breedingPauseSteps,
                                    collaborationStrategyRed, collaborationStrategyBlue, hunterGroupSize);
        boids.setCollisionParam(trackingVelocity, obstacles);
        boids.updateBehavior(currentMethod);

        TVStack boids_pos = boids.getPositions();
        VectorXBool boids_groups = boids.getGroups();

        switch (currentMethod) {
            case FREEFALL:
                {
                    timestep = 0.01;
                    drawBoids(boids_pos);
                }
                break;

            case CIRCLE:
                {
                    timestep = 0.05;
                    currentInitialVelocity = CIRCULAR_VELOCITY;
                    // Origin
                    if (drawOrigin) {
                        nvgBeginPath(vg);
                        TV origin_pos = shift_01_to_screen(TV(0,0), scale, width, height);
                        nvgCircle(vg, origin_pos[0], origin_pos[1], 4.f);
                        nvgFillColor(vg, COLOR_BLACK);
                        nvgFill(vg);
                    }

                    drawBoids(boids_pos);
                }
                break;

            case COHESION:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = ZERO_VELOCITY;
                    drawBoids(boids_pos);
                }
                break;

            case ALIGNMENT:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = RANDOM_VELOCITY;
                    drawBoids(boids_pos);
                }
                break;

            case SEPARATION:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = RANDOM_VELOCITY;
                    drawBoids(boids_pos);
                }
                break;

            case COLLSION_AVOIDANCE:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = ZERO_VELOCITY;

                    // Collision object
                    drawCollisionObject(obstacles);

                    drawBoids(boids_pos);
                }
                break;

            case LEADER:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = ZERO_VELOCITY;

                    // Collision object
                    drawCollisionObject(obstacles);

                    // Leader
                    leaderPos_01 = shift_screen_to_01(leaderPos_screen, scale, width, height);
                    boids_pos.col(0) = leaderPos_01;
                    boids.setLeaderPosition(leaderPos_01);
                    nvgBeginPath(vg);
                    TV leaderPos_01_shifted = shift_01_to_screen(TV(leaderPos_01(0)+leaderRadius_01,leaderPos_01(1)+leaderRadius_01), scale, width, height);
                    leaderRadius_screen = leaderPos_01_shifted(0) - leaderPos_screen(0);
                    nvgCircle(vg, leaderPos_screen[0], leaderPos_screen[1], leaderRadius_screen);
                    nvgFillColor(vg, COLOR_LEADER);
                    nvgFill(vg);

                    // Boids
                    drawBoids(boids_pos, true);
                }
                break;

            case COLLABORATION:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = RANDOM_VELOCITY;

                    for (int i = 0; i < boids.getParticleNumber(); i++) {
                        TV pos = boids_pos.col(i);

                        TV screen_pos = shift_01_to_screen(TV(pos[0], pos[1]), scale, width, height);

                        // Draw groups separately
                        if (boids_groups(i)) {
                            nvgBeginPath(vg);
                            nvgCircle(vg, screen_pos[0], screen_pos[1], 2.f);
                            nvgFillColor(vg, COLOR_RED);
                            nvgFill(vg);
                        }
                        else {
                            nvgBeginPath(vg);
                            nvgCircle(vg, screen_pos[0], screen_pos[1], 2.f);
                            nvgFillColor(vg, COLOR_BLUE);
                            nvgFill(vg);
                        }
                    }

                    // Draw playground boundary
                    TV topLeftCorner = shift_01_to_screen(TV(0, 0), scale, width, height);
                    TV bottomRightCorner = shift_01_to_screen(TV(1,1), scale, width, height);
                    nvgBeginPath(vg);
                    nvgRect(vg, topLeftCorner[0], topLeftCorner[1], bottomRightCorner[0]-topLeftCorner[0], bottomRightCorner[1]-topLeftCorner[1]);
                    nvgStrokeColor(vg, COLOR_BLACK);
                    nvgStroke(vg);
                }
                break;

            default:
                break;
        }
    }

    void drawBoids(TVStack &boids_pos, bool drawLeaderSeparately = false) {
        int startIndex = 0;
        if (drawLeaderSeparately) startIndex = 1;
        for (int i = startIndex; i < boids.getParticleNumber(); i++) {
            TV pos = boids_pos.col(i);
            nvgBeginPath(vg);

            TV screen_pos = shift_01_to_screen(TV(pos[0], pos[1]), scale, width, height);
            nvgCircle(vg, screen_pos[0], screen_pos[1], 2.f);

            nvgFillColor(vg, COLOR_RED);
            nvgFill(vg);
        }
    }

    void drawCollisionObject(Matrix<T, 3, Eigen::Dynamic> &obstacles) {
        for (int i = 0; i < obstacles.cols(); i++) {
            nvgBeginPath(vg);
            TV collision_pos = shift_01_to_screen(TV(obstacles(0,i), obstacles(1,i)), scale, width, height);
            TV collision_radius = shift_01_to_screen(TV(obstacles(0,i)+ obstacles(2,i), obstacles(1,i) + obstacles(2,i)),
                                                     scale, width, height);
            T objRadius_screen = collision_radius[0] - collision_pos[0];
            nvgCircle(vg, collision_pos[0], collision_pos[1], objRadius_screen);
            nvgFillColor(vg, COLOR_BLACK);
            nvgFill(vg);
        }
    }

    void initializeObstacles(Matrix<T, 3, Eigen::Dynamic> &obstacles) {
        for (int i = 0; i < 10; i++) {
            VectorXT randomPosAndRadius = VectorXT::Zero(3, 1).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);});
            obstacles(0,i) = randomPosAndRadius(0);
            obstacles(1,i) = randomPosAndRadius(1);
            obstacles(2,i) = 0.2*randomPosAndRadius(2);
        }
    }

    TV shift_01_to_screen(TV pos_01, T scale, T width, T height)
    {
        return TV(0.5 * (0.5 - scale) * width + scale * pos_01[0] * width, 0.5 * (0.5 - scale) * height + scale * pos_01[1] * height);
    };

    TV shift_screen_to_01(TV pos_screen, T scale, T width, T height)
    {
        return TV((pos_screen[0] - 0.5 * (0.5 - scale) * width) / (scale * width), (pos_screen[1] - 0.5 * (0.5 - scale) * height) / (scale * height));
    };

protected:
    void mouseButtonPressed(int button, int mods) override {
        TV x = {mouseState.lastMouseX, mouseState.lastMouseY};
        if(button == GLFW_MOUSE_BUTTON_LEFT && (x-leaderPos_screen).norm() < leaderRadius_screen) {
            draggingCircle = true;
            draggingCircleOffset = x - leaderPos_screen;
        }
    }

    void mouseButtonReleased(int button, int mods) override {
        draggingCircle = false;
    }

private:
    int loadFonts(NVGcontext* vg)
    {
        int font;
        font = nvgCreateFont(vg, "sans", "../example/Roboto-Regular.ttf");
        if (font == -1) {
            printf("Could not add font regular.\n");
            return -1;
        }
        font = nvgCreateFont(vg, "sans-bold", "../example/Roboto-Bold.ttf");
        if (font == -1) {
            printf("Could not add font bold.\n");
            return -1;
        }
        return 0;
    }

private:
    // Timestep
    T timestep = 0.05;

    // Number of boids
    int n = 20;

    // Behaviour parameters
    MethodTypes newMethod = FREEFALL;
    MethodTypes currentMethod = FREEFALL;
    SchemeTypes currentScheme = SYMPLETIC_EULER;
    InitialVelocityType currentInitialVelocity = ZERO_VELOCITY;
    T boidMass = 3.0;
    T cohesionRadius = 0.2;
    T alignmentRadius = 0.2;
    T separationRadius = 0.02;
    T friction = 2;
    T damping = 0.25;
    T leaderObservationDistance = 0.3;
    T stopApproachLeaderRadius = 0.05;

    // Collision object parameters
    TV trackingVelocity = {-0.1, -0.1};
    int numObstacles = 10;
    Matrix<T, 3, Eigen::Dynamic> obstacles = Matrix<T, 3, Eigen::Dynamic>::Zero(3, numObstacles);

    // Collaboration parameters
    T breedingRadius = 0.03;
    T eliminationRadius = 0.05;
    int maxNumNewBreeded = 1;
    int breedingPauseSteps = 100;
    int maxNumParticles = 300;
    CollaborativeTypes collaborationStrategyRed = GO_TO_CLOSEST;
    CollaborativeTypes collaborationStrategyBlue = NO_STRATEGY;
    T hunterGroupSize = 0.1;

    // Drawing parameters
    T scale = 0.3;
    bool drawOrigin = true;
    bool drawLeader = true;
    bool neverRun = true;
    bool draggingCircle = false;
    TV draggingCircleOffset = TV::Zero(dim);
    TV leaderPos_screen = TV::Zero(dim);
    TV leaderPos_01 = TV::Zero(dim);
    T leaderRadius_01 = 0.01;
    T leaderRadius_screen;
    std::chrono::high_resolution_clock::time_point lastFrame;

    Boids<T, dim> boids = Boids<T, dim>(n);
};

int main(int, char**)
{
    int width = 720; // Default: 720
    int height = 720; // Default: 720
    TestApp app(width, height, "Assignment 3 Boids");
    app.run();

    return 0;
}