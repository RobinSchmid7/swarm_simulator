#ifndef BOIDS_H
#define BOIDS_H
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <fstream>
template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;

// add more for yours
enum MethodTypes {
    FREEFALL=0, CIRCLE=1, COHESION=2, ALIGNMENT=3, SEPARATION=4, COLLSION_AVOIDANCE=5, LEADER=6, COLLABORATION=7
};

enum InitialVelocityType {
    ZERO_VELOCITY=0, RANDOM_VELOCITY=1, CIRCULAR_VELOCITY=2
};

enum SchemeTypes {
    EXPLICIT_EULER=0, SYMPLETIC_EULER=1, EXPLICIT_MIDPOINT=2
};

enum GroupTypes {
    BLUE=0, RED=1
};

enum CollaborativeTypes {
    NO_STRATEGY=0, GO_TO_CLOSEST=1, GO_TO_CENTER_AND_STAY=2, BUILD_GROUPS_OF_THREE_CHASE_CLOSEST=3, BUILD_GROUPS_OF_THREE_PREDICT_BOUNCE_OFF=4
};

template <class T, int dim>
class Boids
{
    typedef Matrix<T, Eigen::Dynamic, 1> VectorXT;
    typedef Matrix<bool, Eigen::Dynamic, 1> VectorXBool;
    typedef Matrix<T, dim, Eigen::Dynamic> TVStack;
    typedef Vector<T, dim> TV;
    typedef Matrix<T, dim, dim> TM;
    
private:
    // Timestep
    T h_ = 0.05;

    // Number of boids
    int n_;

    // Number of time integration steps
    int numIntegrationSteps = 0;
    bool printFlag = true;
    std::ofstream results;

    // Behaviour parameters
    MethodTypes currentMethod_ = FREEFALL;
    SchemeTypes currentScheme_ = SYMPLETIC_EULER;
    InitialVelocityType currentInitialVelocity_ = ZERO_VELOCITY;
    T initialHamiltonian_ = 0;
    T boidMass_ = 1;
    T cohesionRadius_ = 0.2;
    T alignmentRadius_ = 0.2;
    T separationRadius_ = 0.02;
    T friction_ = 0;
    T damping_ = 1;
    T leaderObservationDistance_ = 0.1;
    T stopApproachLeaderRadius_ = 0.05;

    // Collision object parameters
    TV trackingVelocity_;
    Matrix<T, 3, Eigen::Dynamic> obstacles_;

    // Boids parameters
    TVStack positions;
    TVStack velocities;
    TVStack fAttraction;

    // Initial distance from origin
    VectorXT dInit;

    // Forces
    TV fgravity = {0, 1};
    TV fcircle = TV::Zero(dim);

    // Collaboration parameters
    T breedingRadius_ = 0.01;
    T eliminationRadius_ = 0.02;
    int breedingPauseSteps_ = 100;
    int maxNumNewBreeded_ = 1;
    int maxNumParticles_ = 100;
    VectorXBool groups;
    VectorXT breedingCoolOff;
    CollaborativeTypes collaborationStrategyRed_ = GO_TO_CLOSEST;
    CollaborativeTypes collaborationStrategyBlue_ = NO_STRATEGY;
    T hunterGroupSize_ = 0.1;

    bool update = false;

public:
    Boids() :n_(1) {}
    Boids(int n) :n_(n) {
        initializePositionsAndGroups(currentMethod_,n_);
    }
    ~Boids() {}

    void setParticleNumber(int n) {n_ = n;}
    int getParticleNumber() { return n_; }
    void initializePositionsAndGroups(MethodTypes type, int n)
    {
        n_ = n;

        // Initialize random position
//        srand(1); // Used for plotting of Hamiltonian for different integration schemes
        positions = TVStack::Zero(dim, n_).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);});

        // Initial distance from origin
        dInit = VectorXT::Zero(n_, 1);
        for (int i = 0; i < n_; i++) {
            dInit(i) = positions.col(i).norm();
        }

        // Initialize velocity at rest
        velocities = TVStack::Zero(dim, n_);
        switch (currentInitialVelocity_) {
            case ZERO_VELOCITY:
                break;
            case CIRCULAR_VELOCITY:
                for (int i = 0; i < n_; i++) {
                    velocities(0,i) = -positions(1,i);
                    velocities(1,i) = positions(0,i);
                }
                break;
            case RANDOM_VELOCITY:
                sampleAllRandomVelocities(velocities, 0.2);
                break;
            default:
                break;
        }

        // Initialize groups for collaborative behavior
        groups = VectorXBool::Zero(n_,1);
        if (type == COLLABORATION) {
            // Split initial groups
            for (int i = 0; i < n_; i++) {
                groups(i) = bool(i%2);
            }
        }

        fAttraction = TVStack::Zero(dim, n_);
        breedingCoolOff = VectorXT::Zero(n_, 1);
        initialHamiltonian_ = getHamiltonian(velocities);
        numIntegrationSteps = 0;
        printFlag = true;

        if (currentMethod_ == CIRCLE) {
            switch(currentScheme_) {
                case EXPLICIT_EULER:
                    results.open("explicit_euler_"  + std::to_string(h_) + ".txt");
                    break;
                case SYMPLETIC_EULER:
                    results.open("sympletic_euler_"  + std::to_string(h_) + ".txt");
                    break;
                case EXPLICIT_MIDPOINT:
                    results.open("explicit_midpoint_"  + std::to_string(h_) + ".txt");
                    break;
                default:
                    break;
            }
        }
    }

    void updateBehavior(MethodTypes type)
    {
        if(!update)
            return;
        switch (type) {
            // Freefall
            case FREEFALL:
                {
                    boidMass_ = 1;
                    // Analysis of different integration schemes with freefall
                    // 0: explicit euler, 1: symplicit euler, 2: explicit midpoint
                    switch (currentScheme_)
                    {
                        case EXPLICIT_EULER:
                            for (int i = 0; i < n_; i++) {
                                positions.col(i) += h_ * velocities.col(i);
                                velocities.col(i) += h_ * fgravity/boidMass_;
                            }
                            break;

                        case SYMPLETIC_EULER:
                            for (int i = 0; i < n_; i++) {
                                positions.col(i) += h_ * velocities.col(i);
                                velocities.col(i) += h_ * fgravity/boidMass_;
                            }
                            break;

                        case EXPLICIT_MIDPOINT:
                        {
                            TVStack midPositions = TVStack::Zero(dim, n_);
                            TVStack midVelocities = TVStack::Zero(dim, n_);
                            for (int i = 0; i < n_; i++) {
                                midPositions.col(i) = positions.col(i) + 0.5*h_ * velocities.col(i);
                                midVelocities.col(i) = velocities.col(i) + 0.5*h_ * fgravity/boidMass_;
                                positions.col(i) += h_ * midVelocities.col(i);
                                velocities.col(i) += h_ * fgravity/boidMass_;
                            }
                            break;
                        }

                        default:
                            break;
                    }
                }
                break;

            case CIRCLE:
                {
                    boidMass_ = 1;
                    // Analysis of different integration schemes with circular motion
                    // 0: explicit euler, 1: symplicit euler, 2: explicit midpoint
                    switch (currentScheme_)
                    {
                        case EXPLICIT_EULER:
                            for (int i = 0; i < n_; i++) {
                                fcircle = -positions.col(i)*velocities.col(i).squaredNorm()/(dInit(i)*positions.col(i).norm());
                                positions.col(i) += h_ * velocities.col(i);
                                velocities.col(i) += h_ * fcircle/boidMass_;
                            }
                            results << numIntegrationSteps << ", " << getHamiltonian(velocities) << "\n";
                            break;

                        case SYMPLETIC_EULER:
                            for (int i = 0; i < n_; i++) {
                                positions.col(i) += h_ * velocities.col(i);
                                fcircle = -positions.col(i)*velocities.col(i).squaredNorm()/(dInit(i)*positions.col(i).norm());
                                velocities.col(i) += h_ * fcircle/boidMass_;
                            }
                            results << numIntegrationSteps << ", " << getHamiltonian(velocities) << "\n";
                            break;

                        case EXPLICIT_MIDPOINT:
                        {
                            TVStack midPositions = TVStack::Zero(dim, n_);
                            TVStack midVelocities = TVStack::Zero(dim, n_);
                            for (int i = 0; i < n_; i++) {
                                midPositions.col(i) = positions.col(i) + 0.5*h_ * velocities.col(i);
                                fcircle = -positions.col(i)*velocities.col(i).squaredNorm()/(dInit(i)*positions.col(i).norm());
                                midVelocities.col(i) = velocities.col(i) + 0.5*h_ * fcircle/boidMass_;
                                positions.col(i) += h_ * midVelocities.col(i);
                                fcircle = -midPositions.col(i)*midVelocities.col(i).squaredNorm()/(dInit(i)*midPositions.col(i).norm());
                                velocities.col(i) += h_ * fcircle/boidMass_;
                            }
                            results << numIntegrationSteps << ", " << getHamiltonian(velocities) << "\n";
                            break;
                        }

                        default:
                            break;
                    }
                    if (printFlag) {
                        if (getHamiltonian(velocities) > 1.1*initialHamiltonian_ || getHamiltonian(velocities) <  1/1.1*initialHamiltonian_) {
                            std::cout << "Integrations steps until 10% off from initial Hamiltonian: " << numIntegrationSteps << std::endl;
                            printFlag = false;
                        }
                        else {
                            std::cout << "Hamiltonian: " << getHamiltonian(velocities) << std::endl;
                        }
                    }
                    if (numIntegrationSteps == 200) {
                        results.close();
                        printFlag = false;
                    }
                    numIntegrationSteps++;
                }
                break;

            // Cohesion, use sympletic Euler scheme
            case COHESION:
                {
                    for (int i = 0; i < n_; i++) positions.col(i) += h_ * velocities.col(i);
                    setCohesionForce(fAttraction, cohesionRadius_, friction_);
                    for (int i = 0; i < n_; i++) velocities.col(i) += h_ * fAttraction.col(i)/boidMass_;
                }
                break;

            // Alignment, use sympletic Euler scheme
            case ALIGNMENT:
                {
                    for (int i = 0; i < n_; i++) positions.col(i) += h_ * velocities.col(i);
                    setAlignmentForce(fAttraction, alignmentRadius_);
                    for (int i = 0; i < n_; i++) velocities.col(i) += h_ * fAttraction.col(i)/boidMass_;
                }
                break;

            // Cohesion, use sympletic Euler scheme
            case SEPARATION:
                {
                    for (int i = 0; i < n_; i++) positions.col(i) += h_ * velocities.col(i);
                    setSeparationForce(fAttraction, separationRadius_, friction_, damping_);
                    for (int i = 0; i < n_; i++) velocities.col(i) += h_ * fAttraction.col(i)/boidMass_;
                }
                break;

            // Collision avoidance, use sympletic Euler scheme
            case COLLSION_AVOIDANCE:
                {
                    for (int i = 0; i < n_; i++) positions.col(i) += h_ * velocities.col(i);
                    setCollisionAvoidanceForce(fAttraction, cohesionRadius_, alignmentRadius_,separationRadius_,
                                               friction_, damping_,trackingVelocity_,obstacles_);
                    for (int i = 0; i < n_; i++) velocities.col(i) += h_ * fAttraction.col(i)/boidMass_;
                }
                break;

            // Following leader, use sympletic Euler scheme
            case LEADER:
                {
                    for (int i = 0; i < n_; i++) positions.col(i) += h_ * velocities.col(i);
                    setLeaderForce(fAttraction, cohesionRadius_, alignmentRadius_, separationRadius_,
                                   leaderObservationDistance_, stopApproachLeaderRadius_,friction_, damping_, obstacles_);
                    for (int i = 0; i < n_; i++) velocities.col(i) += h_ * fAttraction.col(i)/boidMass_;
                }
                break;

            // Collaborative and adversarial behaviour, use sympletic Euler scheme
            case COLLABORATION:
                {
                    // At each integration step reduce the time boids which have breeded need to wait for breeding a new boid
                    T breedingPauseTime = breedingPauseSteps_*h_;
                    for (int i = 0; i < n_; i++) {
                        if (breedingCoolOff(i) > 0) {
                            breedingCoolOff(i) -= h_;
                        }
                        else {
                            breedingCoolOff(i) = 0;
                        }
                    }

                    for (int i = 0; i < n_; i++) positions.col(i) += h_ * velocities.col(i);
                    TVStack fCollaboration = TVStack::Zero(dim,n_);
                    setCollaborativeForce(fCollaboration, collaborationStrategyRed_, collaborationStrategyBlue_);
                    for (int i = 0; i < n_; i++) velocities.col(i) += h_ * fCollaboration.col(i)/boidMass_;
                    applyCollaborativeRules(breedingRadius_, eliminationRadius_, breedingPauseTime, maxNumNewBreeded_);
                    keepInPlayGround(positions, velocities);
                    if (printFlag) {
                        printPopulationSize();
                        if (checkIfAllEliminated()) {
                            std::cout << "Moves until win: " << numIntegrationSteps << std::endl;
                            printFlag = false;
                        }
                    }
                    numIntegrationSteps++;
                }
                break;

            default:
                break;
        }
    }

    void pause() {
        update = !update;
    }

    void sampleAllRandomVelocities(TVStack &sampledVelocities, T velocityNorm) {
        VectorXT angles = 2.0*M_PI*VectorXT::Zero(sampledVelocities.cols(), 1).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);});
        for (int i = 0; i < sampledVelocities.cols(); i++) {
            sampledVelocities(0,i) = std::cos(angles(i))*velocityNorm;
            sampledVelocities(1,i) = std::sin(angles(i))*velocityNorm;
        }
    }

    void sampleOneRandomPosition(TV &sampledPosition) {
        sampledPosition = TV::Zero(dim).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);});
    }

    void sampleOneRandomVelocity(TV &sampledVelocity, T velocityNorm) {
        TV angle = 2.0*M_PI*TV::Zero(1, 1).unaryExpr([&](T dummy){return static_cast <T> (rand()) / static_cast <T> (RAND_MAX);});
        sampledVelocity(0) = std::cos(angle(0))*velocityNorm;
        sampledVelocity(1) = std::sin(angle(0))*velocityNorm;
    }

    void setBehaviourParam(MethodTypes &currentMethod, SchemeTypes &currentScheme, InitialVelocityType &currentInitialVelocity, T &timestep,
                           T &boidMass, T &cohesionRadius, T &alignmentRadius, T &separationRadius,
                           T &leaderObservationDistance, T &stopApproachLeaderRadius, T &friction, T &damping) {
        currentMethod_ = currentMethod;
        currentScheme_ = currentScheme;
        currentInitialVelocity_ = currentInitialVelocity;
        h_ = timestep;
        boidMass_ = boidMass;
        cohesionRadius_ = cohesionRadius;
        alignmentRadius_ = alignmentRadius;
        separationRadius_ = separationRadius;
        leaderObservationDistance_ = leaderObservationDistance;
        stopApproachLeaderRadius_ = stopApproachLeaderRadius;
        friction_ = friction;
        damping_ = damping;
    }

    void setCollaborationParam(T &breedingRadius, T &eliminationRadius, int &maxNumParticles, int &maxNumNewBreeded, int &breedingPauseSteps,
                               CollaborativeTypes &collaborationStrategyRed, CollaborativeTypes &collaborationStrategyBlue,
                               T &hunterGroupSize) {
        breedingRadius_ = breedingRadius;
        eliminationRadius_ = eliminationRadius;
        maxNumParticles_ = maxNumParticles;
        maxNumNewBreeded_ = maxNumNewBreeded;
        breedingPauseSteps_ = breedingPauseSteps;
        collaborationStrategyRed_ = collaborationStrategyRed;
        collaborationStrategyBlue_ = collaborationStrategyBlue;
        hunterGroupSize_ = hunterGroupSize;
    }

    T getHamiltonian(TVStack &velocities) {
        T hamiltonian = 0;
        for (int i = 0; i < getParticleNumber(); i++) {
            hamiltonian += 0.5*velocities.col(i).squaredNorm();
        }
        return hamiltonian;
    }

    bool checkCollision(T &objMx, T &objMy, T &objRadius, TV position) {
        return bool((objMx-position(0))*(objMx-position(0)) + (objMy-position(1))*(objMy-position(1)) < objRadius*objRadius);
    }

    void setLeaderPosition(TV &leaderPos) {
        positions.col(0) = leaderPos;
    }

    TVStack getPositions() {
        return positions;
    }

    VectorXBool getGroups() {
        return groups;
    }

    T getDistance(int idx1, int idx2, TVStack vec) {
        return (vec.col(idx2)-vec.col(idx1)).norm();
    }

    int getClosest(int &boidIdx, bool groupToApproach) {
        T closestDistance = std::numeric_limits<T>::max();
        int closestBoidIdx = -1;

        for (int i = 0; i < getParticleNumber(); i++) {
            if (i != boidIdx && groups(i) == groupToApproach) {
                T distance = getDistance(boidIdx, i, positions);
                if (distance < closestDistance) {
                    closestDistance = distance;
                    closestBoidIdx = i;
                }
            }
        }

        if (closestBoidIdx == -1) {
            // Handle the case where no boid is found
            return -1;
        }

        // No need for assertion here
        return closestBoidIdx;
    }

    void appendColumnToTVStack(TVStack &inputStack, TV &columnToAppend) {
        inputStack.conservativeResize(inputStack.rows(), inputStack.cols()+1);
        inputStack.col(inputStack.cols()-1) = columnToAppend;
    }

    void removeColumnFromTVStack(TVStack &inputStack, int &columnIdxToRemove) {
        unsigned int numRows = inputStack.rows();
        unsigned int numCols = inputStack.cols()-1;
        if (columnIdxToRemove < numCols) {
            inputStack.block(0,columnIdxToRemove,numRows,numCols-columnIdxToRemove) = inputStack.rightCols(numCols-columnIdxToRemove);
            inputStack.conservativeResize(numRows,numCols);
        }
    }

    void appendElementToVectorXT(VectorXT &inputVector, T &elementToAppend) {
        inputVector.conservativeResize(inputVector.rows()+1, inputVector.cols());
        inputVector(inputVector.rows()-1) = elementToAppend;
    }

    void appendElementToVectorXT(VectorXBool &inputVector, bool &elementToAppend) {
        inputVector.conservativeResize(inputVector.rows()+1, inputVector.cols());
        inputVector(inputVector.rows()-1) = elementToAppend;
    }

    void removeElementFromVectorXT(VectorXT &inputVector, int &elementIdxToRemove) {
        unsigned int numRows = inputVector.rows()-1;
        unsigned int numCols = inputVector.cols();
        if(elementIdxToRemove < numRows+1) {
            inputVector.block(elementIdxToRemove,0,numRows-elementIdxToRemove,numCols) = inputVector.bottomRows(numRows-elementIdxToRemove);
            inputVector.conservativeResize(numRows,numCols);
        }
    }

    void removeElementFromVectorXT(VectorXBool &inputVector, int &elementIdxToRemove) {
        unsigned int numRows = inputVector.rows()-1;
        unsigned int numCols = inputVector.cols();
        if(elementIdxToRemove < numRows+1) {
            inputVector.block(elementIdxToRemove,0,numRows-elementIdxToRemove,numCols) = inputVector.bottomRows(numRows-elementIdxToRemove);
            inputVector.conservativeResize(numRows,numCols);
        }
    }

    void updateTVStack(TVStack &inputStack, TVStack &outputStack) {
        outputStack.resize(inputStack.rows(),inputStack.cols());
        outputStack = inputStack;
    }

    void updateVectorXT(VectorXT &inputVector, VectorXT &outputVector) {
        outputVector.resize(inputVector.rows(),inputVector.cols());
        outputVector = inputVector;
    }

    void updateVectorXT(VectorXBool &inputVector, VectorXBool &outputVector) {
        outputVector.resize(inputVector.rows(),inputVector.cols());
        outputVector = inputVector;
    }

    void setCollisionParam(TV &trackingVelocity, Matrix<T, 3, Eigen::Dynamic> &obstacles) {
        trackingVelocity_ = trackingVelocity;
        obstacles_ = obstacles;
    }

    // Forces
    void setCohesionForce(TVStack &fAttraction, T &cohesionRadius, T &friction) {
        TVStack midPoint = TVStack::Zero(dim, n_);
        fAttraction = TVStack::Zero(dim, n_);

        // Go trough all particles
        for (int i = 0; i < n_; i++) {

            int c = 0;
            for (int j = 0; j < n_; j++) {
                // Particles within the attraction radius
                if (i != j && getDistance(i,j,positions) < cohesionRadius) {
                    midPoint.col(i) += positions.col(j);
                    c++;
                }
            }
            // Midpoint of particles within attraction radius
            // Go to midpoint of particles within attraction radius
            if (c) {
                midPoint.col(i) /= c;
                fAttraction.col(i) += midPoint.col(i) - positions.col(i);

                // Add some friction to slow down particles
                fAttraction.col(i) += -friction*velocities.col(i);
            }
        }
    }

    void setAlignmentForce(TVStack &fAttraction, T &alignmentRadius) {
        TVStack midPoint = TVStack::Zero(dim, n_);
        TVStack alignVelocity = TVStack::Zero(dim, n_);
        fAttraction = TVStack::Zero(dim, n_);

        // Go trough all particles
        for (int i = 0; i < n_; i++) {

            int c = 0;
            for (int j = 0; j < n_; j++) {
                // Particles within the attraction radius
                if (i != j && getDistance(i,j,positions) < alignmentRadius) {
                    alignVelocity.col(i) += velocities.col(j);
                    c++;
                }
            }
            // Midpoint and average alignment velocity of particles within attraction radius
            // Adjust velocity such that it aligns with the particles within the attraction radius
            // Do not use friction here to keep them moving
            if (c) {
                alignVelocity.col(i) /= c;
                fAttraction.col(i) += alignVelocity.col(i) - velocities.col(i);
            }
        }
    }

    void setSeparationForce(TVStack &fAttraction, T &separationRadius, T &friction, T &damping) {
        TVStack midPoint = TVStack::Zero(dim, n_);
        TVStack separationVelocity = TVStack::Zero(dim, n_);
        TVStack alignVelocity = TVStack::Zero(dim, n_);
        fAttraction = TVStack::Zero(dim, n_);

        // Go trough all particles
        for (int i = 0; i < n_; i++) {

            int c = 0;
            for (int j = 0; j < n_; j++) {
                // Particles within the separation radius
                if (i != j && getDistance(i,j,positions) < separationRadius) {
                    separationVelocity.col(i) += positions.col(i) - positions.col(j);
                    c++;
                }
            }

            // Repulsive force for particles within the repulsion radius
            // Force is proportional to distance^(-1/(damping+1)), so for bigger damping the repulsion force is less spiky
            if (c) {
                separationVelocity.col(i) /= c*std::pow(separationVelocity.col(i).norm(), 1/(damping+1));
                fAttraction.col(i) += separationVelocity.col(i) - velocities.col(i);
            }
        }
    }

    void setCollisionAvoidanceForce(TVStack &fAttraction, T &cohesionRadius, T &alignmentRadius, T &separationRadius, T &friction, T &damping,
                                    TV &trackingVelocity, Matrix<T, 3, Eigen::Dynamic> &obstacles) {
        TVStack midPoint = TVStack::Zero(dim, n_);
        TVStack separationVelocity = TVStack::Zero(dim, n_);
        TVStack alignVelocity = TVStack::Zero(dim, n_);
        fAttraction = TVStack::Zero(dim, n_);

        // Go trough all particles
        for (int i = 0; i < n_; i++) {

            int cc = 0;
            int ca = 0;
            int cs = 0;

            for (int j = 0; j < n_; j++) {
                // Particles within the cohesion radius
                if (i != j && getDistance(i,j,positions) < cohesionRadius) {
                    midPoint.col(i) += positions.col(j);
                    cc++;
                }
                // Particles within the alignment radius
                if (i != j && getDistance(i,j,positions) < alignmentRadius) {
                    alignVelocity.col(i) += velocities.col(j);
                    ca++;
                }
                // Particles within the separation radius
                if (i != j && getDistance(i,j,positions) < separationRadius) {
                    separationVelocity.col(i) += positions.col(i) - positions.col(j);
                    cs++;
                }
            }
            // Midpoint of particles within attraction radius
            // Go to midpoint of particles within attraction radius
            if (cc) {
                midPoint.col(i) /= cc;
                fAttraction.col(i) += midPoint.col(i) - positions.col(i);
            }
            // Adjust velocity such that it aligns with the particles within the attraction radius
            // Weight alignment more than finding the closest group, looks more realistic
            if (ca) {
                alignVelocity.col(i) /= ca;
                fAttraction.col(i) += alignVelocity.col(i) - velocities.col(i);
            }

            // Repulsive force for particles within the repulsion radius
            // Force is proportional to distance^(-1/(damping+1)), so for bigger damping the repulsion force is less spiky
            if (cs) {
                separationVelocity.col(i) /= cs*std::pow(separationVelocity.col(i).norm(), 1/(damping+1));
                fAttraction.col(i) += separationVelocity.col(i) - velocities.col(i);
            }

            // Add some friction to slow down particles and reduce overshooting
            fAttraction.col(i) -= friction*velocities.col(i);

            // Track desired velocity to move through obstacles
            fAttraction.col(i) += trackingVelocity - velocities.col(i);

            // Only repulsive force acts if object is within the obstacle
            bool cleanedForce = false;
            for (int k = 0; k < obstacles.cols(); k++) {
                if (checkCollision(obstacles(0,k), obstacles(1,k), obstacles(2,k), positions.col(i))) {
                    if (!cleanedForce) {
                        fAttraction.col(i) << 0, 0;
                        cleanedForce = true;
                    }
                    fAttraction(0,i) += positions(0,i) - obstacles(0,k);
                    fAttraction(1,i) += positions(1,i) - obstacles(1,k);
                    // Use big repulsive force for hard obstacle
                    fAttraction.col(i) *= 20;
                }
            }
        }
    }

    void setLeaderForce(TVStack &fAttraction, T &cohesionRadius, T &alignmentRadius, T &separationRadius, T &leaderObservationDistance, T &stopApproachLeaderRadius,
                        T &friction, T &damping, Matrix<T, 3, Eigen::Dynamic> &obstacles) {
        TVStack midPoint = TVStack::Zero(dim, n_);
        TVStack separationVelocity = TVStack::Zero(dim, n_);
        TVStack alignVelocity = TVStack::Zero(dim, n_);
        fAttraction = TVStack::Zero(dim, n_);

        // Go trough all particles
        for (int i = 0; i < n_; i++) {

            int cc = 0;
            int ca = 0;
            int cs = 0;

            for (int j = 0; j < n_; j++) {
                // Particles within the cohesion radius
                if (i != j && getDistance(i,j,positions) < cohesionRadius) {
                    midPoint.col(i) += positions.col(j);
                    cc++;
                }
                // Particles within the alignment radius
                if (i != j && getDistance(i,j,positions) < alignmentRadius) {
                    alignVelocity.col(i) += velocities.col(j);
                    ca++;
                }
                // Particles within the separation radius
                if (i != j && getDistance(i,j,positions) < separationRadius) {
                    separationVelocity.col(i) += positions.col(i) - positions.col(j);
                    cs++;
                }
            }
            // Midpoint of particles within attraction radius
            // Go to midpoint of particles within attraction radius
            if (cc) {
                midPoint.col(i) /= cc;
                fAttraction.col(i) += midPoint.col(i) - positions.col(i);
            }
            // Adjust velocity such that it aligns with the particles within the attraction radius
            // Weight alignment more than finding the closest group, looks more realistic
            if (ca) {
                alignVelocity.col(i) /= ca;
                fAttraction.col(i) += alignVelocity.col(i) - velocities.col(i);
            }

            // Repulsive force for particles within the repulsion radius
            // Force is proportional to distance^(-1/(damping+1)), so for bigger damping the repulsion force is less spiky
            if (cs) {
                separationVelocity.col(i) /= cs*std::pow(separationVelocity.col(i).norm(), 1/(damping+1));
                fAttraction.col(i) += separationVelocity.col(i) - velocities.col(i);
            }

            // Add some friction to slow down particles and reduce overshooting
            fAttraction.col(i) -= friction*velocities.col(i);

            // Particles within observation distance to leader keep moving to leader
            if (i != 0 && getDistance(0,i,positions) < leaderObservationDistance && getDistance(0,i,positions) > stopApproachLeaderRadius_) {
                fAttraction.col(i) += positions.col(0) - positions.col(i);
            }

            // Only repulsive force acts if object is within the obstacle
            bool cleanedForce = false;
            for (int k = 0; k < obstacles.cols(); k++) {
                if (checkCollision(obstacles(0,k), obstacles(1,k), obstacles(2,k), positions.col(i))) {
                    if (!cleanedForce) {
                        fAttraction.col(i) << 0, 0;
                        cleanedForce = true;
                    }
                    fAttraction(0,i) += positions(0,i) - obstacles(0,k);
                    fAttraction(1,i) += positions(1,i) - obstacles(1,k);
                    // Use big repulsive force for hard obstacle
                    fAttraction.col(i) *= 20;
                }
            }
        }
    }

    void setGoToClosestForce(TVStack &fCollaboration, int &boidIdx) {
        TVStack separationVelocity = TVStack::Zero(dim, n_);
        int closestIdx = getClosest(boidIdx, groups(boidIdx));

        if (closestIdx == -1) {
            // No closest boid found; decide on appropriate behavior
            // For example, you might skip updating this boid or set a default behavior
            return;
        }

        fCollaboration.col(boidIdx) = positions.col(closestIdx) - positions.col(boidIdx);

        int c = 0;
        for (int j = 0; j < n_; j++) {
            // Particles within the separation radius
            if (boidIdx != j && getDistance(boidIdx, j, positions) < separationRadius_) {
                separationVelocity.col(boidIdx) += positions.col(boidIdx) - positions.col(j);
                c++;
            }
        }

        // Repulsive force for particles within the repulsion radius
        if (c) {
            separationVelocity.col(boidIdx) /= c * std::pow(separationVelocity.col(boidIdx).norm(), 1 / (damping_ + 1));
            fCollaboration.col(boidIdx) += separationVelocity.col(boidIdx) - velocities.col(boidIdx);
        }

        fCollaboration.col(boidIdx) += -friction_ * velocities.col(boidIdx);
    }

    void setBuildGroupsOfThreeChaseClosestForce(TVStack &fCollaboration, int &boidIdx, T &hunterGroupSize) {

        // If a group of three is built they chase they stay in this groups and chase the closest enemy
        int boidsClose = 0;
        int closestIdx = 0;
        for (int j = 0; j < n_; j++) {
            if (getDistance(boidIdx,j,positions) < hunterGroupSize && groups(boidIdx) == groups(j)) {
                boidsClose++;
            }
        }

        // If built a group chase the closest enemy and predict its motion
        if (boidsClose > 2 && !checkIfAllEliminated()) {
            TV collisionPoint = TV::Zero(dim);
            closestIdx = getClosest(boidIdx, !groups(boidIdx));
            fCollaboration.col(boidIdx) = (positions.col(closestIdx) - positions.col(boidIdx));
        }

        // If not in a group find a group
        else {
            // Need to find a group
            closestIdx = getClosest(boidIdx, groups(boidIdx));
            fCollaboration.col(boidIdx) = (positions.col(closestIdx) - positions.col(boidIdx));
        }

        // Add some friction for better control
        fCollaboration.col(boidIdx) += -friction_*velocities.col(boidIdx);
    }

    void setBuildGroupsOfThreePredictBounceOffForce(TVStack &fCollaboration, int &boidIdx, T &hunterGroupSize) {

        // If a group of three is built they chase they stay in this groups and chase the closest enemy
        int boidsClose = 0;
        int closestIdx = 0;
        for (int j = 0; j < n_; j++) {
            if (getDistance(boidIdx,j,positions) < hunterGroupSize && groups(boidIdx) == groups(j)) {
                boidsClose++;
            }
        }

        // If built a group chase the closest enemy and predict its motion
        if (boidsClose > 2 && !checkIfAllEliminated()) {
            TV collisionPoint = TV::Zero(dim);
            T wallDistance = 0.05;
            closestIdx = getClosest(boidIdx, !groups(boidIdx));
            setSmartCatchPosition(boidIdx,closestIdx,collisionPoint, wallDistance);
            fCollaboration.col(boidIdx) = (collisionPoint - positions.col(boidIdx));
            fCollaboration.col(boidIdx) += 0.01*(positions.col(closestIdx) - positions.col(boidIdx));
        }

            // If not in a group find a group
        else {
            // Need to find a group
            closestIdx = getClosest(boidIdx, groups(boidIdx));
            fCollaboration.col(boidIdx) = (positions.col(closestIdx) - positions.col(boidIdx));
        }

        // Add some friction for better control
        fCollaboration.col(boidIdx) += -friction_*velocities.col(boidIdx);
    }

    void setGoToCenterAndStay(TVStack &fCollaboration, int &boidIdx) {
        TVStack separationVelocity = TVStack::Zero(dim, n_);
        TV center = {0.5,0.5};
        fCollaboration.col(boidIdx) = center - positions.col(boidIdx);

        int c = 0;
        for (int j = 0; j < n_; j++) {
            // Particles within the separation radius
            if (boidIdx != j && getDistance(boidIdx,j,positions) < separationRadius_) {
                separationVelocity.col(boidIdx) += positions.col(boidIdx) - positions.col(j);
                c++;
            }
        }

        // Repulsive force for particles within the repulsion radius
        // Force is proportional to distance^(-1/(damping+1)), so for bigger damping the repulsion force is less spiky
        if (c) {
            separationVelocity.col(boidIdx) /= c*std::pow(separationVelocity.col(boidIdx).norm(), 1/(damping_+1));
            fCollaboration.col(boidIdx) += separationVelocity.col(boidIdx) - velocities.col(boidIdx);
        }

        fCollaboration.col(boidIdx) += -friction_*velocities.col(boidIdx);
    }

    void setCollaborativeForce(TVStack &fCollaboration, CollaborativeTypes &collaborationStrategyRed, CollaborativeTypes &collaborationStrategyBlue) {
        for (int i = 0; i < getParticleNumber(); i++) {
            switch (int(groups(i))) {
                case RED:
                    {
                        switch(collaborationStrategyRed) {
                            case NO_STRATEGY:
                                break;
                            case GO_TO_CLOSEST:
                                setGoToClosestForce(fCollaboration, i);
                                break;
                            case GO_TO_CENTER_AND_STAY:
                                setGoToCenterAndStay(fCollaboration,i);
                                break;
                            case BUILD_GROUPS_OF_THREE_CHASE_CLOSEST:
                                setBuildGroupsOfThreeChaseClosestForce(fCollaboration, i, hunterGroupSize_);
                                break;
                            case BUILD_GROUPS_OF_THREE_PREDICT_BOUNCE_OFF:
                                setBuildGroupsOfThreePredictBounceOffForce(fCollaboration, i, hunterGroupSize_);
                                break;
                            default:
                                break;
                        }
                    }
                    break;

                case BLUE:
                    {
                        switch(collaborationStrategyBlue) {
                            case NO_STRATEGY:
                                break;
                            case GO_TO_CLOSEST:
                                setGoToClosestForce(fCollaboration, i);
                                break;
                            case GO_TO_CENTER_AND_STAY:
                                setGoToCenterAndStay(fCollaboration,i);
                                break;
                            case BUILD_GROUPS_OF_THREE_CHASE_CLOSEST:
                                setBuildGroupsOfThreeChaseClosestForce(fCollaboration, i, hunterGroupSize_);
                                break;
                            case BUILD_GROUPS_OF_THREE_PREDICT_BOUNCE_OFF:
                                setBuildGroupsOfThreePredictBounceOffForce(fCollaboration, i, hunterGroupSize_);
                                break;
                            default:
                                break;
                        }
                    }
                    break;
                default:
                    break;
            }
        }

    }

    // Methods for collaboration strategies
    void applyCollaborativeRules(T &breedingRadius, T &eliminationRadius, T &breedingPauseTime, int &maxNumNewBreeded) {
        int nCurrent = getParticleNumber();

        TV spawnPosition, spawnVelocity;
        TVStack newPositions = positions;
        TVStack newVelocities = velocities;
        VectorXBool newGroups = groups;
        VectorXT newBreedingCoolOff = breedingCoolOff;

        int numNewBreeded = 0;
        for (int i = 0; i < nCurrent; i++) {
            int enemies = 0;
            for (int j = i + 1; j < nCurrent; j++) {
                // If two boids are within breeding radius and of the same group create a new boid of the same group
                // Only spawn new boids if max number of particles is not reached and both boids are already
                // recovered from any previous breeding
                if (getDistance(i, j, positions) < breedingRadius && groups(i) == groups(j)
                && nCurrent < maxNumParticles_ && breedingCoolOff(i) <= 0 && breedingCoolOff(j) <= 0 && numNewBreeded < maxNumNewBreeded) {
                    // Sample random position and velocity
                    sampleOneRandomPosition(spawnPosition);
                    sampleOneRandomVelocity(spawnVelocity,0.2);

                    // Update new positions, velocities and groups
                    appendColumnToTVStack(newPositions, spawnPosition);
                    appendColumnToTVStack(newVelocities, spawnVelocity);
                    appendElementToVectorXT(newGroups, groups(i));

                    // Boids which created a new one and new created particle need pause from breeding
                    newBreedingCoolOff(i) = breedingPauseTime;
                    newBreedingCoolOff(j) = breedingPauseTime;
                    appendElementToVectorXT(newBreedingCoolOff, breedingPauseTime);

                    // Number of breeded boids per update step increases
                    numNewBreeded++;
                }

                // If two boids are within elimination radius and of different groups they are considered enemies
                if (getDistance(i,j,positions) < eliminationRadius && groups(i) != groups(j)) {
                    enemies++;
                }
            }

            // If there are more than or equal to 3 enemies to a boid it gets eliminated
            if (enemies > 2) {
                removeColumnFromTVStack(newPositions, i);
                removeColumnFromTVStack(newVelocities, i);
                removeElementFromVectorXT(newGroups, i);
            }
        }

        setParticleNumber(newPositions.cols());

        updateTVStack(newPositions, positions);
        updateTVStack(newVelocities, velocities);
        updateVectorXT(newGroups,groups);
        updateVectorXT(newBreedingCoolOff,breedingCoolOff);
    }

    void keepInPlayGround(TVStack &positions, TVStack &velocities) {
        // Swap velocities if outside the boundary
        for (int i = 0; i < positions.cols(); i++) {
            if (positions(0,i) + h_ * velocities(0,i) < 0) velocities(0,i) *= -1;
            if (positions(1,i) + h_ * velocities(1,i) < 0) velocities(1,i) *= -1;
            if (positions(0,i) + h_ * velocities(0,i) > 1) velocities(0,i) *= -1;
            if (positions(1,i) + h_ * velocities(1,i) > 1) velocities(1,i) *= -1;
        }
    }

    void setSmartCatchPosition(int &ownIdx, int &enemyIdx, TV &collisionPoint, T &wallDistance) {
        int closestWallIdx = 0;
        TV ownGlobalPosition = positions.col(ownIdx);
        TV ownGlobalVelocity = velocities.col(ownIdx);
        TV enemyGlobalPosition = positions.col(enemyIdx);
        TV enemyGlobalVelocity = velocities.col(enemyIdx);

        TV ownWallPosition = positions.col(ownIdx);
        TV ownWallVelocity = velocities.col(ownIdx);
        TV enemyWallPosition = positions.col(enemyIdx);
        TV enemyWallVelocity = velocities.col(enemyIdx);

        closestWallIdx = getCollisionWallIdx(enemyIdx);

        // Transform to wall frame for catch point calculation
        transformPositionFromGlobalToWallFrame(ownIdx, ownGlobalPosition, ownWallPosition, closestWallIdx);
        transformVelocityFromGlobalToWallFrame(ownIdx, ownGlobalVelocity, ownWallVelocity, closestWallIdx);
        transformPositionFromGlobalToWallFrame(enemyIdx, enemyGlobalPosition, enemyWallPosition, closestWallIdx);
        transformVelocityFromGlobalToWallFrame(enemyIdx, enemyGlobalVelocity, enemyWallVelocity, closestWallIdx);

        T collisionAngle = 0;
        collisionAngle = getAngleClosestWall(enemyIdx, enemyWallPosition);

        collisionPoint = getCollisionPointGlobalFrameWithOffset(enemyIdx, collisionAngle, wallDistance);
    }

    void transformPositionFromGlobalToWallFrame(int &boidIdx, TV &inputVector, TV &outputVector, int &wallIdx) {
        inputVector = inputVector - getCollisionPointGlobalFrame(boidIdx);
        switch(wallIdx) {
            case 0:
                outputVector(0) = -inputVector(0);
                outputVector(1) = inputVector(1);
                break;
            case 1:
                outputVector(0) = -inputVector(1);
                outputVector(1) = -inputVector(0);
                break;
            case 2:
                outputVector(0) = inputVector(0);
                outputVector(1) = -inputVector(1);
                break;
            case 3:
                outputVector(0) = inputVector(1);
                outputVector(1) = inputVector(0);
                break;
            default:
                break;
        }
    }

    void transformVelocityFromGlobalToWallFrame(int &boidIdx, TV &inputVector, TV &outputVector, int &wallIdx) {
        // Transformation from global to wall frame and vice versa
        switch(wallIdx) {
            case 0:
                outputVector(0) = -inputVector(0);
                outputVector(1) = inputVector(1);
                break;
            case 1:
                outputVector(0) = -inputVector(1);
                outputVector(1) = -inputVector(0);
                break;
            case 2:
                outputVector(0) = inputVector(0);
                outputVector(1) = -inputVector(1);
                break;
            case 3:
                outputVector(0) = inputVector(1);
                outputVector(1) = inputVector(0);
                break;
            default:
                break;
        }
    }

    T getAngleClosestWall(int &boidIdx, TV &positionWallFrame) {
        T wallCollisionLength = 0;
        T wallCollisionAngle  = 0;
        wallCollisionLength = getCollisionDistanceClosestWall(boidIdx);
        wallCollisionAngle = std::asin(positionWallFrame(0)/wallCollisionLength);
        return wallCollisionAngle;
    }

    T getCollisionDistanceClosestWall(int &boidIdx) {
        int closestWallIdx = 0;
        T wallCollisionTime = 0;
        T wallCollisionLength = 0;
        closestWallIdx = getCollisionWallIdx(boidIdx);
        wallCollisionTime = getCollisionWallTime(boidIdx, closestWallIdx);
        wallCollisionLength = wallCollisionTime*(velocities.col(boidIdx)).norm();
        return wallCollisionLength;
    }

    TV getCollisionPointGlobalFrame(int &boidIdx) {
        T wallCollisionTime = 0;
        int closestWallIdx = 0;
        closestWallIdx = getCollisionWallIdx(boidIdx);
        wallCollisionTime = getCollisionWallTime(boidIdx, closestWallIdx);
        TV collisionPoint = TV::Zero(dim);
        TV wallOffset = TV::Zero(dim);
        collisionPoint = positions.col(boidIdx) + wallCollisionTime*velocities.col(boidIdx);
        return collisionPoint;
    }

    TV getCollisionPointGlobalFrameWithOffset(int &boidIdx, T &collisionAngle, T &wallDistance) {
        T wallCollisionTime = 0;
        int closestWallIdx = 0;
        closestWallIdx = getCollisionWallIdx(boidIdx);
        wallCollisionTime = getCollisionWallTime(boidIdx, closestWallIdx);
        TV collisionPoint = TV::Zero(dim);
        TV wallOffsetWall = TV::Zero(dim);
        TV wallOffsetGlobal = TV::Zero(dim);
        wallOffsetWall << -std::sin(collisionAngle)*wallDistance, std::cos(collisionAngle)*wallDistance;
        transformVelocityFromGlobalToWallFrame(boidIdx, wallOffsetWall, wallOffsetGlobal, closestWallIdx);
        collisionPoint = positions.col(boidIdx) + wallCollisionTime*velocities.col(boidIdx) + wallOffsetGlobal;
        return collisionPoint;
    }

    int getCollisionWallIdx(int &boidIdx) {
        T wallCollisionTime = 1000;
        int closestWallIdx = 0;
        for (int i = 0; i < 4; i++) {
            if (getCollisionWallTime(boidIdx, i) > 0 && getCollisionWallTime(boidIdx, i) < wallCollisionTime) {
                wallCollisionTime = getCollisionWallTime(boidIdx, i);
                closestWallIdx = i;
            }
        }
        return closestWallIdx;
    }

    T getCollisionWallTime(int &boidIdx, int &wallIdx) {
        T wallCollisionTime = 0;
        // Wall Idx: 0: North, 1: East, 2: South, 3: West
        switch(wallIdx) {
            case 0:
                wallCollisionTime = (0.0 - positions(1,boidIdx)) / velocities(1,boidIdx);
                break;
            case 1:
                wallCollisionTime = (1.0 - positions(0,boidIdx)) / velocities(0,boidIdx);
                break;
            case 2:
                wallCollisionTime = (1.0 - positions(1,boidIdx)) / velocities(1,boidIdx);
                break;
            case 3:
                wallCollisionTime = (0.0 - positions(0,boidIdx)) / velocities(0,boidIdx);
                break;
            default:
                break;
        }
        return wallCollisionTime;
    }

    bool checkIfAllEliminated() {
        int numGroup1 = 0;
        int numGroup2 = 0;
        for (int i = 0; i < getParticleNumber(); i++) {
            if (groups(i))
                numGroup1++;
            else
                numGroup2++;
        }
        return (numGroup1 == 0 || numGroup2 == 0);
    }

    void printPopulationSize() {
        int numGroup1 = 0;
        int numGroup2 = 0;
        for (int i = 0; i < getParticleNumber(); i++) {
            if (groups(i))
                numGroup1++;
            else
                numGroup2++;
        }
        std::cout << "Group 1: " << numGroup1 << ", Group 2: " << numGroup2 << std::endl;
    }
};
#endif