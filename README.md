# 🌊 2DFluidSimulation
This is a university project for the 'Advanced Animation and Physics Simulation' course, part of the Master's Degree in Computer Graphics, Games and Virtual Reality at Rey Juan Carlos University (URJC). 

The simulation contains an **Eulerian (Grid-Based) Semi-Lagrangian solver** and hybrids **Lagrangian-Eulerian solvers**, specifically **FLIP (Fluid-Implicit Particle) 0.95** and **APIC (Affine Particle-in-Cell)**.

The fluid solvers were independently developed and can be found in `src/Fluid2Exercise.cpp`, whereas the engine framework was provided as part of the course materials.

# 🧮 Fluid Simulation Methodology
We use a fractional step method to solve Navier-Stokes equations.
1. Advection. This step represents how the fluid transports itself.  
&emsp; &emsp; &emsp;
$$\frac{\partial \mathbf{u}}{\partial t} = -(\mathbf{u} \cdot \nabla)\mathbf{u}$$  
&emsp; **Implementation:** Depending on the selected mode, advection is calculated using:  
&emsp;&emsp;-Semi-Lagrangian grid advection  
&emsp;&emsp;-Particle-based advection  
&emsp; _Note:_ more detail on each method on the subsections below.  

3. Viscosity. It models the internal friction of the fluid, smoothing out velocity differences between neighboring cells.    
&emsp; &emsp; &emsp;
$$\frac{\partial \mathbf{u}}{\partial t} = \nu \nabla^2 \mathbf{u}$$  
&emsp; **Implementation:** Since the viscosity in this simulation is small, we can use an aproximation with explicit finite differences.  

3. Volume forces. This are the forces that act over the whole fluid volume, such as gravity and buoyancy.  
&emsp; &emsp; &emsp;
$$\frac{\partial \mathbf{u}}{\partial t} = \mathbf{f}_{external}$$  
&emsp; **Implementation:** Euler explicit.  

5. Pressure gradient. Since we don't care about energy conservation on a graphic simulation, the pressure gradient is solely responsible for enforcing mass conservation. An incompressible fluid must meet the condition of zero divergence.  
&emsp; &emsp; &emsp; $$\nabla \cdot \mathbf{u} = 0$$  
&emsp; Pressure acts as the corrective scalar force field per unit area that opposes compression and stretching in the continuum, keeping a uniform density.
&emsp; **Implementation:** to solve we create a system with a matrix that represents the connections between cells, the pressure incognita and the velocity divergence. Once the pressure is resolved, the velociy is calculated with finite differences.

# 🔬 Solvers

## Semi-lagrangian advection (Eulerian smoke simulation)
When the variable at `src/Fluid2.h` `flipEnabled` is inactive, the engine utilizes a pure Eulerian approach using unconditionally stable backward-tracing advection.
* **Concept:** the algorithm can guess how a particle would move on a velocity grid integrating the trayectory of time backwards. The intuition for this is that, if we know the velocity of a cell at the current time, we can integrate backwards to find where that originated from in the previous step. The properties of the fluid (ink and velocity) at that previous location on the current time are then sampled and transferred to the current cell for the next step.
  
It uses bilinear interpolation to sample values from the previous time step. This guarantees unconditional stability but introduces noticeable numerical diffusion (mass and velocity dissipation), which is why it is used for the smoke simulation.

<img width="359" height="360" alt="adv" src="https://github.com/user-attachments/assets/38bcbad5-64bb-46f8-a4b8-766ced3184ff" />

## FLIP 0.95 (Hybrid liquid simulation)
Fluid-Implicit-Particle (FLIP) and Particle-in-Cell (PIC) are hybrid formulations. They use a grid to calculate viscosity, external forces and pressure projection, but use particles to translate the fluid properties, such as the ink color and velocity, on the advection step. At the start of a frame, particle properties are transferred to the grid. After the grid-based physics step solves for forces and incompressibility, the updated grid velocities are transferred back to the particles.

* **PIC**: This method simply interpolates the velocities that surround a particle on the grid to update the velocity. This is highly stable but suffers from energy dissipation.
* **FLIP**: This method fixes the dissipation of PIC adding a small delta velocity based on the difference of the velocity calculated for the particle in the advection step and the grid velocity at the end of the simulation step. This preserves energy and vorticity, but can become unstable and noisy.
* **FLIP 0.95**: This is the chosen method for the solver. To balance the two previous methods, it merges both solutions, for a more stable solver with less dissipation.  
&emsp; &emsp; &emsp; $$\mathbf{v}_{\text{particle}} = 0.95 \cdot \mathbf{v}_{\text{FLIP}} + 0.05 \cdot \mathbf{v}_{\text{PIC}}$$  

<img width="359" height="360" alt="flii" src="https://github.com/user-attachments/assets/91aac883-35a2-4b1d-b86f-36ed5b26276d" />

## APIC (Hybrid liquid simulation)
In the code Affine-Particle in Cell (APIC) is an upgrade from flip, so it can't be enabled if both xx and xx aren't active.  
APIC follows the same concept as the methods seen before, but fixing both the dissipation of PIC and the instability of FLIP, at the cost of computational power. By storing not just the velocity of a particle, but also the local velocity gradient, it can make a more informed calculation, to keep information of vortices, shearing, and complex rotational behaviors that FLIP and PIC typically lose or distort. 

* During the Particle-to-Grid phase, this gradient allows particles to make a much more informed, structurally accurate velocity contribution to the grid nodes.
* During the Grid-to-Particle phase, the spatial derivatives of the grid's interpolation weights are calculated to update the particle's gradient matrices.

<img width="359" height="360" alt="api" src="https://github.com/user-attachments/assets/e2cc4dd9-4cd5-476b-abaf-353dc8190e0a" />

# 🛠️ Build instructions 
```
# 1. Clone repository
git clone https://github.com/aurogr/2DFluidSimulation
cd 2DFluidSimulation

# 2. Configure project and create build folder
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release

# 3. Build
cmake --build build --config Release

# 4. Run executable
# Linux/macOS:
./build/2DFluidSimulation
# Windows (Command Prompt/PowerShell):
.\build\Release\2DFluidSimulation.exe
```
