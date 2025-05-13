# Swarm Simulator
An interactive 2D boid-based swarm simulation that models complex flocking behaviors.

<video src="docs/demo.webm" width="100%" controls autoplay loop muted></video>

## Features
- Physics-based boid simulation with cohesion, alignment, and separation behaviors
- Customizable parameters for swarm behavior
- Real-time visualization with smooth animations
- Leader-following dynamics
- Collision avoidance

## Installation
```bash
mkdir build
cd build
cmake ..
make -j4
```

## Running
```bash
./build/src/app/app
```

## Controls
- `SPACE` - Pause/resume simulation
- `ESC` - Exit application
