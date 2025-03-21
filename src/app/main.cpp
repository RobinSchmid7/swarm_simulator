#include "application.h"
#include <imgui.h>

#include <iostream>
#include <math.h>
#include <deque>
#include <chrono>
#include <map>

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
                
                // Initialize both swarms
                boids1.initializePositionsAndGroups(currentMethod, n);
                boids2.initializePositionsAndGroups(currentMethod, n);
                
                initializeObstacles(obstacles);

                // Initialize leader positions - keep them clearly separated
                leaderPos1_01 = TV(0.25, 0.5); // Left side
                leaderPos1_screen = shift_01_to_screen(leaderPos1_01, scale, width, height);
                boids1.setLeaderPosition(leaderPos1_01);
                
                // Position second leader on the right side
                leaderPos2_01 = TV(0.75, 0.5); // Right side
                leaderPos2_screen = shift_01_to_screen(leaderPos2_01, scale, width, height);
                boids2.setLeaderPosition(leaderPos2_01);

                neverRun = false;
            }
            
            if(keyDown[GLFW_KEY_SPACE]) {
                boids1.pause();
                boids2.pause();
            }
            if(keyDown[GLFW_KEY_ESCAPE])
                exit(0);
            lastFrame = now;

            // Handle dragging of either leader
            if(draggingLeaderIndex == 0) {
                leaderPos1_screen(0) = mouseState.lastMouseX - draggingCircleOffset(0);
                leaderPos1_screen(1) = mouseState.lastMouseY - draggingCircleOffset(1);
            } else if(draggingLeaderIndex == 1) {
                leaderPos2_screen(0) = mouseState.lastMouseX - draggingCircleOffset(0);
                leaderPos2_screen(1) = mouseState.lastMouseY - draggingCircleOffset(1);
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
        // SetNextWindowPos(ImVec2(0, 780), ImGuiCond_Always);
        // SetNextWindowSize(ImVec2(600,220), ImGuiCond_Always);
        // Begin("Collaboration Parameters");
        // const char* strategyNames[] = {"No strategy", "Go to closest of own group", "Go to center and stay", "Form groups of three and chase closest enemy", "Form groups of three and move to bounce off point"};
        // Combo("Strategy red", (int*)&collaborationStrategyRed, strategyNames, std::size(strategyNames));
        // Combo("Strategy blue", (int*)&collaborationStrategyBlue, strategyNames, std::size(strategyNames));
        // SliderFloat("Breeding radius", &breedingRadius, 0.0f, 1.0f, "%.5f");
        // SliderFloat("Elimination radius", &eliminationRadius, 0.0f, 1.0f, "%.5f");
        // SliderInt("Max # boids per iteration", &maxNumNewBreeded, 0, 5, "%i");
        // SliderInt("Breeding pause steps", &breedingPauseSteps, 10, 1000, "%i");
        // SliderFloat("Hunter group size", &hunterGroupSize, 0.01f, 0.5f, "%.5f");
        // SliderInt("Max # of boids", &maxNumParticles, 1, 1000, "%i");
        // End();
    }

    void drawNanoVG() override {
        // Adjust boid behavior parameters to make movement smoother
        // Increase damping to reduce vibration
        T smoother_damping = 0.6; // Increased from 0.25 (assuming that was your value)
        
        // Slightly increase friction to slow down rapid changes
        T smoother_friction = 2.5; // Increased from 2.0
        
        // Update both swarms with smoother parameters
        boids1.setBehaviourParam(currentMethod, currentScheme, currentInitialVelocity, timestep, boidMass, 
                               cohesionRadius, alignmentRadius, separationRadius,
                               leaderObservationDistance, stopApproachLeaderRadius, 
                               smoother_friction, smoother_damping); // Use smoother values
        
        boids2.setBehaviourParam(currentMethod, currentScheme, currentInitialVelocity, timestep, boidMass, 
                               cohesionRadius, alignmentRadius, separationRadius,
                               leaderObservationDistance, stopApproachLeaderRadius, 
                               smoother_friction, smoother_damping); // Use smoother values
        
        boids1.setCollaborationParam(breedingRadius, eliminationRadius, maxNumParticles, maxNumNewBreeded, 
                                    breedingPauseSteps, collaborationStrategyRed, collaborationStrategyBlue, hunterGroupSize);
        boids2.setCollaborationParam(breedingRadius, eliminationRadius, maxNumParticles, maxNumNewBreeded, 
                                    breedingPauseSteps, collaborationStrategyRed, collaborationStrategyBlue, hunterGroupSize);
        
        boids1.setCollisionParam(trackingVelocity, obstacles);
        boids2.setCollisionParam(trackingVelocity, obstacles);
        
        boids1.updateBehavior(currentMethod);
        boids2.updateBehavior(currentMethod);

        // Get positions of both swarms
        TVStack boids1_pos = boids1.getPositions();
        TVStack boids2_pos = boids2.getPositions();
        
        // If needed, get groups
        VectorXBool boids1_groups = boids1.getGroups();
        VectorXBool boids2_groups = boids2.getGroups();

        // Handle the drawing based on the current method
        switch (currentMethod) {
            case FREEFALL:
                {
                    timestep = 0.01;
                    drawBoids(boids1_pos, COLOR_RED, true);
                    drawBoids(boids2_pos, COLOR_GREEN, true);
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

                    drawBoids(boids1_pos, COLOR_RED, true);
                    drawBoids(boids2_pos, COLOR_GREEN, true);
                }
                break;

            case COHESION:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = ZERO_VELOCITY;
                    drawBoids(boids1_pos, COLOR_RED, true);
                    drawBoids(boids2_pos, COLOR_GREEN, true);
                }
                break;

            case ALIGNMENT:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = RANDOM_VELOCITY;
                    drawBoids(boids1_pos, COLOR_RED, true);
                    drawBoids(boids2_pos, COLOR_GREEN, true);
                }
                break;

            case SEPARATION:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = RANDOM_VELOCITY;
                    drawBoids(boids1_pos, COLOR_RED, true);
                    drawBoids(boids2_pos, COLOR_GREEN, true);
                }
                break;

            case COLLSION_AVOIDANCE:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = ZERO_VELOCITY;

                    // Collision object
                    drawCollisionObject(obstacles);

                    drawBoids(boids1_pos, COLOR_RED, true);
                    drawBoids(boids2_pos, COLOR_GREEN, true);
                }
                break;

            case LEADER:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = ZERO_VELOCITY;

                    // Draw collision objects
                    drawCollisionObject(obstacles);

                    // Set leader 1 position
                    leaderPos1_01 = shift_screen_to_01(leaderPos1_screen, scale, width, height);
                    boids1_pos.col(0) = leaderPos1_01;
                    boids1.setLeaderPosition(leaderPos1_01);
                    
                    // Set leader 2 position
                    leaderPos2_01 = shift_screen_to_01(leaderPos2_screen, scale, width, height);
                    boids2_pos.col(0) = leaderPos2_01;
                    boids2.setLeaderPosition(leaderPos2_01);
                    
                    // Calculate leader radii for collision detection - use a larger value for visibility
                    TV leaderPos1_01_shifted = shift_01_to_screen(TV(leaderPos1_01(0)+leaderRadius_01, leaderPos1_01(1)+leaderRadius_01), 
                                                                scale, width, height);
                    leaderRadius1_screen = leaderPos1_01_shifted(0) - leaderPos1_screen(0);
                    
                    TV leaderPos2_01_shifted = shift_01_to_screen(TV(leaderPos2_01(0)+leaderRadius_01, leaderPos2_01(1)+leaderRadius_01), 
                                                                scale, width, height);
                    leaderRadius2_screen = leaderPos2_01_shifted(0) - leaderPos2_screen(0);
                    
                    // Draw leader 1 (green) with a larger size
                    nvgBeginPath(vg);
                    nvgCircle(vg, leaderPos1_screen[0], leaderPos1_screen[1], 5.f); // Increased from 2.f to 5.f
                    nvgFillColor(vg, COLOR_LEADER); // Green
                    nvgFill(vg);
                    
                    // Draw leader 2 (blue) with a larger size
                    nvgBeginPath(vg);
                    nvgCircle(vg, leaderPos2_screen[0], leaderPos2_screen[1], 5.f); // Increased from 2.f to 5.f
                    nvgFillColor(vg, COLOR_BLUE);
                    nvgFill(vg);

                    // Draw debug text to show leader positions
                    char buf[128];
                    snprintf(buf, sizeof(buf), "Leader 1: (%.2f, %.2f)", leaderPos1_screen[0], leaderPos1_screen[1]);
                    nvgFontSize(vg, 14.0f);
                    nvgFontFace(vg, "sans");
                    nvgFillColor(vg, COLOR_BLACK);
                    nvgTextAlign(vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                    nvgText(vg, 10, 500, buf, NULL);
                    
                    snprintf(buf, sizeof(buf), "Leader 2: (%.2f, %.2f)", leaderPos2_screen[0], leaderPos2_screen[1]);
                    nvgText(vg, 10, 520, buf, NULL);

                    // Draw the boids for each swarm with different colors
                    drawBoidsWithColor(boids1_pos, COLOR_RED, true);
                    drawBoidsWithColor(boids2_pos, COLOR_GREEN, true);
                }
                break;

            case COLLABORATION:
                {
                    timestep = 0.05;
                    currentScheme = SYMPLETIC_EULER;
                    currentInitialVelocity = RANDOM_VELOCITY;
                    
                    // Draw each boid according to its group
                    for (int i = 0; i < boids1.getParticleNumber(); i++) {
                        TV pos = boids1_pos.col(i);
                        TV screen_pos = shift_01_to_screen(TV(pos[0], pos[1]), scale, width, height);
                        
                        // Draw groups separately
                        if (i < boids1_groups.size() && boids1_groups(i)) {
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
                    
                    for (int i = 0; i < boids2.getParticleNumber(); i++) {
                        TV pos = boids2_pos.col(i);
                        TV screen_pos = shift_01_to_screen(TV(pos[0], pos[1]), scale, width, height);
                        
                        // Draw groups separately
                        if (i < boids2_groups.size() && boids2_groups(i)) {
                            nvgBeginPath(vg);
                            nvgCircle(vg, screen_pos[0], screen_pos[1], 2.f);
                            nvgFillColor(vg, COLOR_GREEN);
                            nvgFill(vg);
                        }
                        else {
                            nvgBeginPath(vg);
                            nvgCircle(vg, screen_pos[0], screen_pos[1], 2.f);
                            nvgFillColor(vg, COLOR_RED);
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

    void drawBoids(TVStack &boids_pos, NVGcolor boidColor, bool drawLeaderSeparately = false) {
        int startIndex = 0;
        if (drawLeaderSeparately) startIndex = 1;
        
        // Store current positions and smoothed velocities - use separate maps for each swarm
        static std::map<void*, TVStack> prev_boids_pos_map;
        static std::map<void*, std::vector<TV>> smoothed_velocities_map;
        
        // Get or create the entries for this swarm based on boids_pos pointer
        void* swarm_id = &boids_pos;
        if (prev_boids_pos_map.find(swarm_id) == prev_boids_pos_map.end()) {
            prev_boids_pos_map[swarm_id] = boids_pos;
            smoothed_velocities_map[swarm_id] = std::vector<TV>(boids_pos.cols(), TV::Zero());
            return; // Skip first frame for this swarm
        }
        
        // Reference to this swarm's data
        TVStack& prev_boids_pos = prev_boids_pos_map[swarm_id];
        std::vector<TV>& smoothed_velocities = smoothed_velocities_map[swarm_id];
        
        // Resize if needed
        if (smoothed_velocities.size() != boids_pos.cols()) {
            smoothed_velocities.resize(boids_pos.cols(), TV::Zero());
        }
        
        // Smoothing and line params
        const T smoothing_factor = static_cast<T>(0.2);
        const T lineScale = static_cast<T>(2.0);
        NVGcolor velocityLineColor = nvgRGBA(0, 180, 0, 255);
        
        for (int i = startIndex; i < boids_pos.cols(); i++) {
            TV pos = boids_pos.col(i);
            TV screen_pos = shift_01_to_screen(TV(pos[0], pos[1]), scale, width, height);
            
            // Draw boid circle with the specified color
            nvgBeginPath(vg);
            nvgCircle(vg, screen_pos[0], screen_pos[1], 2.f);
            nvgFillColor(vg, boidColor);
            nvgFill(vg);
            
            // Calculate velocity vectors (same as before)
            if (i < prev_boids_pos.cols()) {
                TV current_vel = (pos - prev_boids_pos.col(i)) / timestep;
                
                smoothed_velocities[i] = smoothing_factor * current_vel + 
                                        (1 - smoothing_factor) * smoothed_velocities[i];
                
                if (smoothed_velocities[i].norm() > 0) {
                    T maxVel = static_cast<T>(0.02);
                    TV scaled_vel = smoothed_velocities[i].normalized() * 
                                  lineScale * 
                                  std::min(smoothed_velocities[i].norm(), maxVel);
                    
                    TV end_pos_01 = pos + scaled_vel;
                    TV screen_end_pos = shift_01_to_screen(end_pos_01, scale, width, height);
                    
                    nvgBeginPath(vg);
                    nvgMoveTo(vg, screen_pos[0], screen_pos[1]);
                    nvgLineTo(vg, screen_end_pos[0], screen_end_pos[1]);
                    nvgStrokeColor(vg, velocityLineColor);
                    nvgStrokeWidth(vg, 0.7f);
                    nvgStroke(vg);
                }
            }
        }
        
        // Update previous positions for next frame
        prev_boids_pos = boids_pos;
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
            obstacles(2,i) = 0.05*randomPosAndRadius(2);
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

    void drawBoidsWithColor(TVStack &boids_pos, NVGcolor boidColor, bool drawLeaderSeparately = false) {
        int startIndex = 0;
        if (drawLeaderSeparately) startIndex = 1;
        
        // Store current positions and smoothed velocities
        static std::map<const void*, TVStack> prev_boids_pos_map;
        static std::map<const void*, std::vector<TV>> smoothed_velocities_map;
        
        // Static variables for leader position tracking
        static TV prev_leader_pos = TV::Zero();
        static float leader_direction = 0.0f;
        
        // Use the pointer to boids_pos as a unique key for this swarm
        const void* swarm_id = static_cast<const void*>(&boids_pos);
        
        // Initialize if first time for this swarm
        if (prev_boids_pos_map.find(swarm_id) == prev_boids_pos_map.end()) {
            prev_boids_pos_map[swarm_id] = boids_pos;
            smoothed_velocities_map[swarm_id].resize(boids_pos.cols(), TV::Zero());
            return; // Skip first frame for this swarm
        }
        
        // Get reference to this swarm's stored data
        TVStack& prev_pos = prev_boids_pos_map[swarm_id];
        std::vector<TV>& smoothed_vels = smoothed_velocities_map[swarm_id];
        
        // Resize smoothed velocity array if needed
        if (smoothed_vels.size() != boids_pos.cols()) {
            smoothed_vels.resize(boids_pos.cols(), TV::Zero());
        }
        
        // Smoothing parameters
        const T smoothing_factor = static_cast<T>(0.1);
        const T lineScale = static_cast<T>(1.5);
        const T maxVel = static_cast<T>(0.015);
        
        // Line color
        NVGcolor velocityLineColor = nvgRGBA(0, 180, 0, 255);
        
        // Position smoothing
        static std::map<const void*, std::vector<TV>> smoothed_positions_map;
        
        // Initialize smoothed positions if needed
        if (smoothed_positions_map.find(swarm_id) == smoothed_positions_map.end()) {
            smoothed_positions_map[swarm_id].resize(boids_pos.cols());
            for (int i = 0; i < boids_pos.cols(); i++) {
                smoothed_positions_map[swarm_id][i] = boids_pos.col(i);
            }
        }
        
        // Resize if needed
        if (smoothed_positions_map[swarm_id].size() != boids_pos.cols()) {
            smoothed_positions_map[swarm_id].resize(boids_pos.cols());
            for (int i = 0; i < boids_pos.cols(); i++) {
                smoothed_positions_map[swarm_id][i] = boids_pos.col(i);
            }
        }
        
        // Store direction for boids
        static std::map<const void*, std::vector<float>> boid_directions_map;
        
        // Initialize directions if needed
        if (boid_directions_map.find(swarm_id) == boid_directions_map.end()) {
            boid_directions_map[swarm_id].resize(boids_pos.cols(), 0.0f);
        }
        
        // Resize directions if needed
        if (boid_directions_map[swarm_id].size() != boids_pos.cols()) {
            boid_directions_map[swarm_id].resize(boids_pos.cols(), 0.0f);
        }
        
        // Reference to directions
        std::vector<float>& boid_directions = boid_directions_map[swarm_id];
        
        // Update smoothed positions
        std::vector<TV>& smoothed_positions = smoothed_positions_map[swarm_id];
        
        // Get leader position and direction (if available)
        TV leaderPos = boids_pos.col(0); // First boid is the leader
        float leaderDirection = leader_direction; // Use stored direction by default
        
        // Calculate leader direction from velocity
        if (prev_pos.cols() > 0) {
            TV leaderVel = (leaderPos - prev_pos.col(0)) / timestep;
            if (leaderVel.norm() > 0.001) {
                leaderDirection = atan2(leaderVel[1], leaderVel[0]);
                leader_direction = leaderDirection; // Store for future use
            }
        }
        
        // Update leader position for next frame
        prev_leader_pos = leaderPos;
        
        // Process each boid
        for (int i = startIndex; i < boids_pos.cols(); i++) {
            // Get the actual position from simulation
            TV pos = boids_pos.col(i);
            
            // Apply position smoothing for visualization
            const T pos_smoothing = static_cast<T>(0.2);
            smoothed_positions[i] = pos_smoothing * pos + (1 - pos_smoothing) * smoothed_positions[i];
            
            // Use smoothed position for drawing
            TV draw_pos = smoothed_positions[i];
            TV screen_pos = shift_01_to_screen(TV(draw_pos[0], draw_pos[1]), scale, width, height);
            
            // Current direction for this boid
            float dirAngle = boid_directions[i];
            
            // Calculate velocity
            TV current_vel = TV::Zero();
            if (i < prev_pos.cols()) {
                current_vel = (pos - prev_pos.col(i)) / timestep;
                
                // Apply velocity smoothing
                smoothed_vels[i] = smoothing_factor * current_vel + 
                                  (1 - smoothing_factor) * smoothed_vels[i];
            }
            
            // Distance to leader
            T distToLeader = (pos - leaderPos).norm();
            
            // Velocity magnitude
            T velMagnitude = smoothed_vels[i].norm();
            
            // Threshold for "not moving significantly"
            const T velThreshold = static_cast<T>(0.005);
            
            // Distance threshold for "close to leader"
            const T closeToLeaderThreshold = stopApproachLeaderRadius * 1.5;
            
            if (velMagnitude > velThreshold) {
                // Boid is moving - use its velocity for direction
                dirAngle = atan2(smoothed_vels[i][1], smoothed_vels[i][0]);
                boid_directions[i] = dirAngle;
                
                // Draw the velocity line
                TV scaled_vel = smoothed_vels[i].normalized() * 
                              lineScale * 
                              std::min(velMagnitude, maxVel);
                
                TV end_pos_01 = draw_pos + scaled_vel;
                TV screen_end_pos = shift_01_to_screen(end_pos_01, scale, width, height);
                
                nvgBeginPath(vg);
                nvgMoveTo(vg, screen_pos[0], screen_pos[1]);
                nvgLineTo(vg, screen_end_pos[0], screen_end_pos[1]);
                nvgStrokeColor(vg, velocityLineColor);
                nvgStrokeWidth(vg, 0.7f);
                nvgStroke(vg);
            } 
            else if (distToLeader < closeToLeaderThreshold) {
                // Boid is close to leader and not moving significantly
                // Gradually align with leader's direction
                const float alignmentSpeed = 0.05f; // How quickly to align with leader
                
                // Calculate angle difference (ensure we take the shorter path)
                float angleDiff = leaderDirection - boid_directions[i];
                if (angleDiff > M_PI) angleDiff -= 2 * M_PI;
                if (angleDiff < -M_PI) angleDiff += 2 * M_PI;
                
                // Gradually rotate towards leader direction
                boid_directions[i] += angleDiff * alignmentSpeed;
                
                // Keep angle in proper range
                if (boid_directions[i] > M_PI) boid_directions[i] -= 2 * M_PI;
                if (boid_directions[i] < -M_PI) boid_directions[i] += 2 * M_PI;
                
                // Use the updated direction
                dirAngle = boid_directions[i];
            }
            // Otherwise keep the current direction
            
            // Draw the NATO drone symbol
            const float size = 5.0f;  // Base size of the symbol
            
            // Save current transform
            nvgSave(vg);
            
            // Translate to boid position and rotate to match direction
            nvgTranslate(vg, screen_pos[0], screen_pos[1]);
            nvgRotate(vg, dirAngle);
            
            // Draw the NATO drone symbol
            // Main diamond body
            nvgBeginPath(vg);
            nvgMoveTo(vg, size, 0);        // Front
            nvgLineTo(vg, 0, size/2);      // Right
            nvgLineTo(vg, -size/2, 0);     // Back
            nvgLineTo(vg, 0, -size/2);     // Left
            nvgClosePath(vg);
            nvgFillColor(vg, boidColor);
            nvgFill(vg);
            
            // Outline
            nvgStrokeColor(vg, nvgRGBA(0, 0, 0, 100));
            nvgStrokeWidth(vg, 0.5f);
            nvgStroke(vg);
            
            // Wings
            nvgBeginPath(vg);
            nvgMoveTo(vg, 0, size*0.8f);     // Right wing tip
            nvgLineTo(vg, -size*0.3f, 0);    // Body
            nvgLineTo(vg, 0, -size*0.8f);    // Left wing tip
            nvgStrokeColor(vg, boidColor);
            nvgStrokeWidth(vg, 1.0f);
            nvgStroke(vg);
            
            // Central dot
            nvgBeginPath(vg);
            nvgCircle(vg, 0, 0, size/5);
            nvgFillColor(vg, nvgRGBA(0, 0, 0, 150));
            nvgFill(vg);
            
            // Restore transform
            nvgRestore(vg);
        }
        
        // Update stored positions for next frame
        prev_boids_pos_map[swarm_id] = boids_pos;
    }

protected:
    void mouseButtonPressed(int button, int mods) override {
        TV x = {mouseState.lastMouseX, mouseState.lastMouseY};
        
        // Make the interaction radius larger for easier selection
        const float interactionRadius = 15.0f; // Increased from leaderRadius*_screen
        
        // Check if we're clicking on leader 1
        if (button == GLFW_MOUSE_BUTTON_LEFT && (x-leaderPos1_screen).norm() < interactionRadius) {
            draggingLeaderIndex = 0;
            draggingCircleOffset = x - leaderPos1_screen;
            std::cout << "Dragging leader 1" << std::endl; // Debug output
        }
        // Check if we're clicking on leader 2
        else if (button == GLFW_MOUSE_BUTTON_LEFT && (x-leaderPos2_screen).norm() < interactionRadius) {
            draggingLeaderIndex = 1;
            draggingCircleOffset = x - leaderPos2_screen;
            std::cout << "Dragging leader 2" << std::endl; // Debug output
        }
    }

    void mouseButtonReleased(int button, int mods) override {
        draggingLeaderIndex = -1;
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
    MethodTypes newMethod = LEADER;
    MethodTypes currentMethod = LEADER;
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
    T scale = 0.8;
    bool drawOrigin = true;
    bool drawLeader = true;
    bool neverRun = true;
    bool draggingCircle = false;
    TV draggingCircleOffset = TV::Zero(dim);
    TV leaderPos1_screen = TV::Zero(dim);
    TV leaderPos1_01 = TV::Zero(dim);
    TV leaderPos2_screen = TV::Zero(dim);
    TV leaderPos2_01 = TV::Zero(dim);
    T leaderRadius_01 = 0.01;
    T leaderRadius1_screen = 0.0;
    T leaderRadius2_screen = 0.0;
    std::chrono::high_resolution_clock::time_point lastFrame;

    // Add a second boids instance for the second swarm
    Boids<T, dim> boids1 = Boids<T, dim>(n);
    Boids<T, dim> boids2 = Boids<T, dim>(n);
    
    // Track which leader is being dragged
    int draggingLeaderIndex = -1;
};

int main(int, char**)
{
    int width = 1920; // Default: 720
    int height = 1080; // Default: 720
    TestApp app(width, height, "Boids");
    app.run();

    return 0;
}