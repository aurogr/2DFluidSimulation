//////////////////////////////////////////////////////////////
// Master en Informatica Grafica, Juegos y Realidad virtual //
// Animacion y simulacion avanzada II                       //
// Fluidos hibridos                                         //
// Aurora Garcia Raigon                                     //
//////////////////////////////////////////////////////////////

// Changes to Fluid2Exercise.cpp to add the necessary code for FLIP 0.95 and APIC methods
// To try FLIP 0.95 enable flipEnabled on Fluid.h and disable apicEnabled
// To try APIC enable both, flipEnabled and apicEnabled on Fluid.h

#include "Scene.h"

#include "Numeric/PCGSolver.h"

namespace asa
{
namespace
{
////////////////////////////////////////////////
// Add any reusable classes or functions HERE //
////////////////////////////////////////////////

template <typename T>
T interpolate(const Vector2 &index, const Array2<T> &array)
{
    int x0 = floor(index.x);
    int y0 = floor(index.y);

    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float tx = index.x - x0;
    float ty = index.y - y0;

    // clamp borders, left/bottom is 0 and right/top is size - 1
    x0 = clamp(x0, 0, array.getSize().x - 1);
    x1 = clamp(x1, 0, array.getSize().x - 1);
    y0 = clamp(y0, 0, array.getSize().y - 1);
    y1 = clamp(y1, 0, array.getSize().y - 1);

    return bilerp<T>(
        array.getValue(x0, y0), array.getValue(x1, y0), array.getValue(x0, y1), array.getValue(x1, y1), tx, ty);
}

float getRandomFloat(float min, float max)
{
    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    return min + r * (max - min);
}

float partial_cx(float inv_dx, float ty, float u00, float u10, float u01, float u11)
{
    // derivatives of weights
    float dw00dx = -inv_dx * (1 - ty);
    float dw10dx = inv_dx * (1 - ty);
    float dw01dx = -inv_dx * ty;
    float dw11dx = inv_dx * ty;

    // sumatory of velocities * derivatives
    return u00 * dw00dx + u10 * dw10dx + u01 * dw01dx + u11 * dw11dx;
}

float partial_cy(float inv_dy, float tx, float u00, float u10, float u01, float u11)
{
    // derivatives of weights
    float dw00dy = -inv_dy * (1 - tx);
    float dw10dy = -inv_dy * tx;
    float dw01dy = inv_dy * (1 - tx);
    float dw11dy = inv_dy * tx;
    
    // sumatory of velocities * derivatives
    return u00 * dw00dy + u10 * dw10dy + u01 * dw01dy + u11 * dw11dy;
}


}  // namespace

// init particles
void Fluid2::initParticles()
{
    // particle creation HERE

    // THEORY: 
    // Particles are not just dummy indicators of if a cell has fluid in it, they transport the ink and velocity properties
    // We need to init the grid, populating it with particles (with 0 velocity and black ink to start)
    // Create 4 particles for each grid cell, distributed along the cell

    float dx = grid.getDx().x;
    float dy = grid.getDx().y;

    for (uint i = 0; i < grid.getSize().x; i++) 
    {
        for (uint j = 0; j < grid.getSize().y; j++) 
        {
            Vector2 cellLeftCorner = grid.getCellPos(Index2(i, j)) - Vector2(dx * 0.5f, dy * 0.5f);

            for (uint h = 0; h < 2; h++) 
            {
                for (uint v = 0; v < 2; v++) 
                {
                    float paddingh = dx/10;
                    float minXPercent = (0.5f * h) + paddingh;
                    float maxXPercent = (0.5f * (h + 1)) - paddingh;

                    float paddingv = dy / 10;
                    float minYPercent = (0.5f * v) + paddingv;
                    float maxYPercent = (0.5f * (v + 1)) - paddingv;

                    float randomXPercent = getRandomFloat(minXPercent, maxXPercent);
                    float randomYPercent = getRandomFloat(minYPercent, maxYPercent);

                    Vector2 finalPos = cellLeftCorner + Vector2(randomXPercent * dx, randomYPercent * dy);

                    particles.addParticle(finalPos, Vector2(0.f, 0.f), Vector3(0.f, 0.f, 0.f), Vector2(0.f, 0.f), Vector2(0.f, 0.f));
                }
            }
        }
    }
}

// advection
void Fluid2::fluidAdvection(const float dt)
{
    // THEORY:
    // Represents how the fluid transports itself
    // The fluid moves based on its own velocity (u* = velocity gradient * u)

    uint Nx = grid.getSize().x;
    uint Ny = grid.getSize().y;

    if (flipEnabled) {
        
        // ---- Move particles with RK2 with grid velocities HERE ---
        for (uint i = 0; i < particles.getSize(); i++) 
        {
            Vector2 pk0 = particles.getPosition(i);

            // vk1 = interpolate(Pk0)
            Vector2 idx0_x = grid.getFaceIndexX(pk0);
            Vector2 idx0_y = grid.getFaceIndexY(pk0);

            Vector2 vk0 = Vector2(interpolate(idx0_x, velocityX), interpolate(idx0_y, velocityY));

            // Pk1 = Pk0 + dt/2 * vk0
            Vector2 pk1 = pk0 + dt * 0.5f * vk0;

            // vk1 = interpolate(Pk1)
            Vector2 idx1_x = grid.getFaceIndexX(pk1);
            Vector2 idx1_y = grid.getFaceIndexY(pk1);

            Vector2 vk1 = Vector2(interpolate(idx1_x, velocityX), interpolate(idx1_y, velocityY));

            // Pk2 = Pk0 + dt * vk1;
            Vector2 pk2 = pk0 + dt * vk1;
            
            // --- Ensure particle remains inside the domain HERE ---
            AABox2 domain = grid.getDomain();

            pk2.x = std::max(pk2.x, domain.minPosition.x);
            pk2.x = std::min(pk2.x, domain.maxPosition.x);
            pk2.y = std::max(pk2.y, domain.minPosition.y);
            pk2.y = std::min(pk2.y, domain.maxPosition.y);

            particles.setPosition(i, pk2);
        }

        // --- Create ink grid from particles HERE ---
        // empty ink array
        inkRGB.clear();
        inkRGB.resize(grid.getSize());
        
        // create a temp array to store how many particles are in a cell (weight)
        Array2<float> inkWeights(grid.getSize());

        for (uint i = 0; i < particles.getSize(); i++) 
        {
            Vector2 pos = particles.getPosition(i);

            Vector2 idx = grid.getCellIndex(pos);
            
            // get the 4 nearest points on grid
            int x0 = floor(idx.x);
            int y0 = floor(idx.y);
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            float tx = idx.x - x0;
            float ty = idx.y - y0;

            x0 = clamp(x0, 0, Nx - 1);
            x1 = clamp(x1, 0, Nx - 1);
            y0 = clamp(y0, 0, Ny - 1);
            y1 = clamp(y1, 0, Ny - 1);

            float w00 = (1.0f - tx) * (1.0f - ty);
            float w10 = tx * (1.0f - ty);
            float w01 = (1.0f - tx) * ty;
            float w11 = tx * ty;

            // update ink and weight of every point
            Vector3 ink = particles.getInk(i);

            inkRGB.setValue(x0, y0, inkRGB.getValue(x0, y0) + ink * w00);
            inkWeights.setValue(x0, y0, inkWeights.getValue(x0, y0) + w00);

            inkRGB.setValue(x1, y0, inkRGB.getValue(x1, y0) + ink * w10);
            inkWeights.setValue(x1, y0, inkWeights.getValue(x1, y0) + w10);

            inkRGB.setValue(x0, y1, inkRGB.getValue(x0, y1) + ink * w01);
            inkWeights.setValue(x0, y1, inkWeights.getValue(x0, y1) + w01);

            inkRGB.setValue(x1, y1, inkRGB.getValue(x1, y1) + ink * w11);
            inkWeights.setValue(x1, y1, inkWeights.getValue(x1, y1) + w11);
        }

        // normalize
        for (uint x = 0; x < grid.getSize().x; x++) 
        {
            for (uint y = 0; y < grid.getSize().y; y++) 
            {
                float w = inkWeights.getValue(x, y);

                if (w > 0.0f) 
                {
                    Vector3 accumulatedInk = inkRGB.getValue(x, y);
                    inkRGB.setValue(x, y, accumulatedInk / Vector3(w,w,w));
                }
            }
        }

        
        // --- create velocityX grid from particles HERE ---
        velocityX.clear();
        velocityX.resize(grid.getSizeFacesX());

        Array2<float> velocityXWeights(grid.getSizeFacesX());

        for (uint p_idx = 0; p_idx < particles.getSize(); p_idx++) {
            Vector2 p_pos = particles.getPosition(p_idx);

            Vector2 grid_idx = grid.getFaceIndexX(p_pos);

            int x0 = floor(grid_idx.x);
            int y0 = floor(grid_idx.y);
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            float tx = grid_idx.x - x0;
            float ty = grid_idx.y - y0;

            x0 = clamp(x0, 0, grid.getSizeFacesX().x - 1);
            x1 = clamp(x1, 0, grid.getSizeFacesX().x - 1);
            y0 = clamp(y0, 0, grid.getSizeFacesX().y - 1);
            y1 = clamp(y1, 0, grid.getSizeFacesX().y - 1);

            float w00 = (1.0f - tx) * (1.0f - ty);
            float w10 = tx * (1.0f - ty);
            float w01 = (1.0f - tx) * ty;
            float w11 = tx * ty;

            //update velocity
            float p_vel = particles.getVelocity(p_idx).x;

            if (apicEnabled) 
            {
                // THEORY:
                // PIC transfers velocity using simple distance weights sum(v_particle * w_node). This dissipates a lot of energy
                // FLIP makes up for it with an added delta velocity, but it is unstable
                
                // APIC uses information of the fluid around the particle to correctly calculate the velocity value for the grid
                // for that, it stores at the end of the grid to particle transfer phase the velocity gradient
                // sum(v_particle + c_particle.dot(x_node - xparticle) * w_node)

                Vector2 cx = particles.getCx(p_idx);

                Vector2 pos00 = grid.getFacePosX(Index2(x0, y0));
                float vel00 = p_vel + cx.dot(p_pos - pos00);
                velocityX.setValue(x0, y0, velocityX.getValue(x0, y0) + vel00 * w00);
                velocityXWeights.setValue(x0, y0, velocityXWeights.getValue(x0, y0) + w00);

                Vector2 pos10 = grid.getFacePosX(Index2(x1, y0));
                float vel10 = p_vel + cx.dot(p_pos - pos10);
                velocityX.setValue(x1, y0, velocityX.getValue(x1, y0) + vel10 * w10);
                velocityXWeights.setValue(x1, y0, velocityXWeights.getValue(x1, y0) + w10);

                Vector2 pos01 = grid.getFacePosX(Index2(x0, y1));
                float vel01 = p_vel + cx.dot(p_pos - pos01);
                velocityX.setValue(x0, y1, velocityX.getValue(x0, y1) + vel01 * w01);
                velocityXWeights.setValue(x0, y1, velocityXWeights.getValue(x0, y1) + w01);

                Vector2 pos11 = grid.getFacePosX(Index2(x1, y1));
                float vel11 = p_vel + cx.dot(p_pos - pos11);
                velocityX.setValue(x1, y1, velocityX.getValue(x1, y1) + vel11 * w11);
                velocityXWeights.setValue(x1, y1, velocityXWeights.getValue(x1, y1) + w11);
            } 
            else 
            {
                velocityX.setValue(x0, y0, velocityX.getValue(x0, y0) + p_vel * w00);
                velocityXWeights.setValue(x0, y0, velocityXWeights.getValue(x0, y0) + w00);

                velocityX.setValue(x1, y0, velocityX.getValue(x1, y0) + p_vel * w10);
                velocityXWeights.setValue(x1, y0, velocityXWeights.getValue(x1, y0) + w10);

                velocityX.setValue(x0, y1, velocityX.getValue(x0, y1) + p_vel * w01);
                velocityXWeights.setValue(x0, y1, velocityXWeights.getValue(x0, y1) + w01);

                velocityX.setValue(x1, y1, velocityX.getValue(x1, y1) + p_vel * w11);
                velocityXWeights.setValue(x1, y1, velocityXWeights.getValue(x1, y1) + w11);
            }
        }

        // normalize
        for (uint x = 0; x < grid.getSizeFacesX().x; x++) {
            for (uint y = 0; y < grid.getSizeFacesX().y; y++) {
                float w = velocityXWeights.getValue(x, y);

                if (w > 0.0f) {
                    float accumulatedVelocity = velocityX.getValue(x, y);
                    velocityX.setValue(x, y, accumulatedVelocity / w);
                }
            }
        }

        // create velocityY grid from particles HERE
        velocityY.clear();
        velocityY.resize(grid.getSizeFacesY());

        Array2<float> velocityYWeights(grid.getSizeFacesY());

        for (uint p_idx = 0; p_idx < particles.getSize(); p_idx++) {
            Vector2 p_pos = particles.getPosition(p_idx);

            Vector2 grid_idx = grid.getFaceIndexY(p_pos);

            int x0 = floor(grid_idx.x);
            int y0 = floor(grid_idx.y);
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            float tx = grid_idx.x - x0;
            float ty = grid_idx.y - y0;

            x0 = clamp(x0, 0, grid.getSizeFacesY().x - 1);
            x1 = clamp(x1, 0, grid.getSizeFacesY().x - 1);
            y0 = clamp(y0, 0, grid.getSizeFacesY().y - 1);
            y1 = clamp(y1, 0, grid.getSizeFacesY().y - 1);

            float w00 = (1.0f - tx) * (1.0f - ty);
            float w10 = tx * (1.0f - ty);
            float w01 = (1.0f - tx) * ty;
            float w11 = tx * ty;

            // update velocity
            float p_vel = particles.getVelocity(p_idx).y;

            if (apicEnabled) 
            {
                // same theory as velocity x
                Vector2 cy = particles.getCy(p_idx);

                Vector2 pos00 = grid.getFacePosY(Index2(x0, y0));
                float vel00 = p_vel + cy.dot(p_pos - pos00);
                velocityY.setValue(x0, y0, velocityY.getValue(x0, y0) + vel00 * w00);
                velocityYWeights.setValue(x0, y0, velocityYWeights.getValue(x0, y0) + w00);

                Vector2 pos10 = grid.getFacePosY(Index2(x1, y0));
                float vel10 = p_vel + cy.dot(p_pos - pos10);
                velocityY.setValue(x1, y0, velocityY.getValue(x1, y0) + vel10 * w10);
                velocityYWeights.setValue(x1, y0, velocityYWeights.getValue(x1, y0) + w10);

                Vector2 pos01 = grid.getFacePosY(Index2(x0, y1));
                float vel01 = p_vel + cy.dot(p_pos - pos01);
                velocityY.setValue(x0, y1, velocityY.getValue(x0, y1) + vel01 * w01);
                velocityYWeights.setValue(x0, y1, velocityYWeights.getValue(x0, y1) + w01);

                Vector2 pos11 = grid.getFacePosY(Index2(x1, y1));
                float vel11 = p_vel + cy.dot(p_pos - pos11);
                velocityY.setValue(x1, y1, velocityY.getValue(x1, y1) + vel11 * w11);
                velocityYWeights.setValue(x1, y1, velocityYWeights.getValue(x1, y1) + w11);
            } 
            else 
            {
                velocityY.setValue(x0, y0, velocityY.getValue(x0, y0) + p_vel * w00);
                velocityYWeights.setValue(x0, y0, velocityYWeights.getValue(x0, y0) + w00);

                velocityY.setValue(x1, y0, velocityY.getValue(x1, y0) + p_vel * w10);
                velocityYWeights.setValue(x1, y0, velocityYWeights.getValue(x1, y0) + w10);

                velocityY.setValue(x0, y1, velocityY.getValue(x0, y1) + p_vel * w01);
                velocityYWeights.setValue(x0, y1, velocityYWeights.getValue(x0, y1) + w01);

                velocityY.setValue(x1, y1, velocityY.getValue(x1, y1) + p_vel * w11);
                velocityYWeights.setValue(x1, y1, velocityYWeights.getValue(x1, y1) + w11);
            }
        }

        // normalize
        for (uint x = 0; x < grid.getSizeFacesY().x; x++) {
            for (uint y = 0; y < grid.getSizeFacesY().y; y++) {
                float w = velocityYWeights.getValue(x, y);

                if (w > 0.0f) {
                    float accumulatedVelocity = velocityY.getValue(x, y);
                    velocityY.setValue(x, y, accumulatedVelocity / w);
                }
            }
        }

        // save current state velocities HERE
        if (!apicEnabled) // only for flip, apic doesn't use the delta velocity
        {
            oldVelocityX = velocityX;
            oldVelocityY = velocityY;
        }
        
    } else {

        // THEORY: Semi-lagrangian
        // Guess how a particle would move on a velocity grid, but for each cell integrate the trayectory of time backwards 
        // (We have a particle p1 on cell c1 on t, if we can guess where this particle was on t-1 (c2) it means that the particle now on c2 will be on c1 on t+1)

        {
            // Ink SL advection HERE
            static Array2<Vector3> inkRGB_Copy;
            inkRGB_Copy = inkRGB;

            for (uint i = 0; i < grid.getSize().x; i++) 
            {
                for (uint j = 0; j < grid.getSize().y; j++) 
                {
                    // get current cell pos
                    Vector2 x = grid.getCellPos(Index2(i, j));

                    // get velocity of current cell (interpolated because velocity is on faces and we want it for the position on the center)
                    Vector2 u = Vector2(interpolate(grid.getFaceIndexX(x), velocityX),
                                        interpolate(grid.getFaceIndexY(x), velocityY));

                    // integrate position backwards
                    Vector2 xPrev = x - u * dt;

                    // interpolate to find value on center of the cell that holds the xPrev position
                    asa::Vector3 ink_xPrev_t = interpolate(grid.getCellIndex(xPrev), inkRGB_Copy);

                    // update ink with advection
                    inkRGB.setValue(i, j, ink_xPrev_t);
                }
            }
        }

        {
            // Velocity SL advection term HERE
            static Array2<float> u_X_Copy;
            static Array2<float> u_Y_Copy;
            u_X_Copy = velocityX;
            u_Y_Copy = velocityY;

            // velocity x
            for (uint i = 0; i < grid.getSizeFacesX().x; i++) 
            {
                for (uint j = 0; j < grid.getSizeFacesX().y; j++) 
                {
                    // get current face pos
                    Vector2 x = grid.getFacePosX(Index2(i, j));

                    // get velocity of current face (interpolate because velocity Y is on a different face)
                    Vector2 u = Vector2(u_X_Copy.getValue(i, j),
                                        interpolate(grid.getFaceIndexY(x), u_Y_Copy));

                    // integrate position backwards
                    Vector2 xPrev = x - u * dt;

                    // interpolate to find value on the cell face that holds the xPrev position
                    float u_xPrev_t = interpolate(grid.getFaceIndexX(xPrev), u_X_Copy);

                    // update velocity with advection
                    velocityX.setValue(i, j, u_xPrev_t);
                }
            }

            // velocity y
            for (uint i = 0; i < grid.getSizeFacesY().x; i++) 
            {
                for (uint j = 0; j < grid.getSizeFacesY().y; j++) 
                {
                    // get current face pos
                    Vector2 y = grid.getFacePosY({i, j});

                    // get velocity of current face (interpolate because velocity X is on a different face)
                    Vector2 u = Vector2(interpolate(grid.getFaceIndexX(y), u_X_Copy), 
                                        u_Y_Copy.getValue(i, j));

                    // integrate position backwards
                    Vector2 yPrev = y - u * dt;

                    // interpolate to find value on the cell face that holds the yPrev position
                    float u_yPrev_t = interpolate(grid.getFaceIndexY(yPrev), u_Y_Copy);

                    // update velocity with advection
                    velocityY.setValue(i, j, u_yPrev_t);
                }
            }
        }
    }
}

void Fluid2::fluidEmission()
{
    float minArea1X = -0.2f;
    float maxArea1X = -0.1f;
    float minArea2X = maxArea1X;
    float maxArea2X = 0.1f;
    float minArea3X = maxArea2X;
    float maxArea3X = 0.2f;

    float minAreaY = -1.9f;
    float maxAreaY = -1.75f;

    if (flipEnabled) {
        // Emitters contribution to particles HERE

        // THEORY:
        // we check to see which particles are in range of the ink and velocity injection domain and update them
        for (uint i = 0; i < particles.getSize(); i++)
        {
            Vector2 pos = particles.getPosition(i);

            if (pos.y >= minAreaY && pos.y <= maxAreaY) 
            {
                if (pos.x >= minArea1X && pos.x < maxArea1X) 
                {
                    particles.setInk(i, Vector3(1.0f, 0.0f, 1.0f));
                } 
                else if (pos.x >= minArea2X && pos.x <= maxArea2X) {

                    particles.setInk(i, Vector3(1.0f, 1.0f, 0.0f));
                } 
                else if (pos.x > maxArea2X && pos.x <= maxArea3X) 
                {
                    particles.setInk(i, Vector3(0.0f, 1.0f, 1.0f));
                }

                if (pos.x >= minArea1X && pos.x <= maxArea3X) {
                    particles.setVelocity(i, Vector2(0.0f, 8.0f));
                }
            }
        }

    } else {
        // Emitters contribution to grid HERE

        // THEORY:
        // emitters are just the velocity injected to the grid system so that it has some movement
        // just go through every cell, and if the cell is inside a certain emissive area set velocity and ink initial values

        // Ink injection (center of cells)
        for (uint i = 0; i < grid.getSize().x; ++i) 
        {
            for (uint j = 0; j < grid.getSize().y; ++j) 
            {
                auto pos = grid.getCellPos(Index2(i, j));

                if (pos.y >= minAreaY && pos.y <= maxAreaY) 
                {
                    if (pos.x >= minArea1X && pos.x < maxArea1X) 
                    {
                        inkRGB.setValue(i, j, Vector3(1.0f, 0.0f, 1.0f));
                    } 
                    else if (pos.x >= minArea2X && pos.x <= maxArea2X) 
                    {
                        inkRGB.setValue(i, j, Vector3(1.0f, 1.0f, 0.0f));
                    } 
                    else if (pos.x > maxArea2X && pos.x <= maxArea3X) 
                    {
                        inkRGB.setValue(i, j, Vector3(0.0f, 1.0f, 1.0f));
                    }
                }
            }
        }

        // Velocity injection (faces of cells)
        
        for (uint i = 0; i < grid.getSizeFacesY().x; ++i) 
        {
            for (uint j = 0; j < grid.getSizeFacesY().y; ++j) 
            {
                auto pos = grid.getFacePosY(Index2(i, j));

                if (pos.y >= minAreaY && pos.y <= maxAreaY) 
                {
                    if (pos.x >= minArea1X && pos.x <= maxArea3X) 
                    {
                        velocityY.setValue(i, j, 8.0f);
                    }
                }
            }
        }

    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    // THEORY:
    // euler explicit: u* = u + dt/rho * f_ext * rho
    
    // our only external force (at the moment) is gravity
    //  u* = u + dt * G

    // Apply gravity only to velocity Y
    for (uint i = 0; i < grid.getSizeFacesY().x; i++) 
    {
        for (uint j = 0; j < grid.getSizeFacesY().y; j++) 
        {
            float new_u = velocityY.getValue(i, j) + dt * Scene::kGravity;

            velocityY.setValue(i, j, new_u);
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    // THEORY: 
    // viscosity is small so we can aproximate with explicit finite difference
    // u* = u + dt/rho * mu * laplacian

    static Array2<float> u_X_Copy;
    static Array2<float> u_Y_Copy;
    u_X_Copy = velocityX;
    u_Y_Copy = velocityY;

    const float dx2_x = grid.getDx().x * grid.getDx().x;
    const float dx2_y = grid.getDx().y * grid.getDx().y;

    const auto mu = Scene::kViscosity;
    const auto rho = Scene::kDensity;


    // Apply viscosity to velocity X
    for (uint i = 0; i < grid.getSizeFacesX().x; i++) 
    {
        for (uint j = 0; j < grid.getSizeFacesX().y; j++) 
        {
            // clamp indices to edges
            int iPrev = std::max((int)i - 1, 0);
            int iNext = std::min((int)i + 1, (int)grid.getSizeFacesX().x - 1);
            int jPrev = std::max((int)j - 1, 0);
            int jNext = std::min((int)j + 1, (int)grid.getSizeFacesX().y - 1);

            // get velocity values
            const auto u_ij = u_X_Copy.getValue(i, j);
            const auto u_iPrevj = u_X_Copy.getValue(iPrev, j);
            const auto u_iNextj = u_X_Copy.getValue(iNext, j);
            const auto u_ijPrev = u_X_Copy.getValue(i, jPrev);
            const auto u_ijNext = u_X_Copy.getValue(i, jNext);

            // calculate second spatial derivatives with finite difference
            float dudx2 = (u_iNextj - 2 * u_ij + u_iPrevj) / dx2_x;
            float dudy2 = (u_ijNext - 2 * u_ij + u_ijPrev) / dx2_y;

            // laplacian is the sum of second derivatives
            float laplacian = dudx2 + dudy2;

            // update velocity with viscosity term
            velocityX.setValue(i, j, u_ij + (dt / rho * mu) * laplacian);
        }
    }

    // Apply viscosity to velocity Y
    for (uint i = 0; i < grid.getSizeFacesY().x; i++) 
    {
        for (uint j = 0; j < grid.getSizeFacesY().y; j++) 
        {
            // clamp indices to edges
            int iPrev = std::max((int)i - 1, 0);
            int iNext = std::min((int)i + 1, (int)grid.getSizeFacesY().x - 1);
            int jPrev = std::max((int)j - 1, 0);
            int jNext = std::min((int)j + 1, (int)grid.getSizeFacesY().y - 1);
            
            // get velocity values
            const auto u_ij = u_Y_Copy.getValue(i, j);
            const auto u_iPrevj = u_Y_Copy.getValue(iPrev, j);
            const auto u_iNextj = u_Y_Copy.getValue(iNext, j);
            const auto u_ijPrev = u_Y_Copy.getValue(i, jPrev);
            const auto u_ijNext = u_Y_Copy.getValue(i, jNext);
            
            // calculate second spatial derivatives with finite difference
            float dudx2 = (u_iNextj - 2 * u_ij + u_iPrevj) / dx2_x;
            float dudy2 = (u_ijNext - 2 * u_ij + u_ijPrev) / dx2_y;
            
            // laplacian is the sum of second derivatives
            float laplacian = dudx2 + dudy2;
            
            // update velocity with viscosity term
            velocityY.setValue(i, j, u_ij + (dt / rho * mu) * laplacian);
        }
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    // Incompressibility / Pressure term HERE
    
    // THEORY:
    // Pressure opposes compressions in the continuum, keeping a uniform density
    // Its tied to energy conservation and other properties such as temperature but we are going to ignore because its not essencial for this grapphics simulations
    // To ensure the incomprensibility condition (mass conservation) we need zero divergence

    // Formula: pressure laplacian (because it searches for changes against its neighbors) = density/deltaTime * velocity divergence (directions where velocity fluctuates)
    // this calculates a pressure that pushes back with the exact ammount of force needed to stop fluctuation
    
    // Steps: we need to rearrange formula to solve as a system Ax = b where:
    // A is a matrix that encodes the spatial grid geometry and boundary conditions, defining how a cell's pressure interacts with its neighbors
    // x is the pressure we want to calculate
    // b is the velocity divergence as a vector, scaled by density and time

    // Variables
    float dx = grid.getDx().x;
    float dx2 = dx * dx;
    float dy = grid.getDx().y;
    float dy2 = dy * dy;

    uint Nx = grid.getSize().x;
    uint Ny = grid.getSize().y;
        
    const float cellsNumber = Nx * Ny;

    // Set normal velocity components in all boundaries to 0
    for (uint i = 0; i < Nx; i++) {
        velocityY.setValue(i, 0, 0);
        velocityY.setValue(i,Ny, 0);
    }

    for (uint j = 0; j < Ny; j++) {
        velocityX.setValue(0, j, 0);
        velocityX.setValue(Nx, j, 0);
    }

    // Calc pressure (build system Ax = b)
    PCGSolver<float> solver = PCGSolver<float>();
    solver.set_solver_parameters(1e-3, 200);

    // Fill RHS (b)
    static std::vector<float> RHS;

    if (RHS.size() != cellsNumber) {
        RHS.resize(cellsNumber);
    }

    for (uint j = 0; j < Ny; j++) 
    {
        int rowOffset = j * Nx;

        for (uint i = 0; i < Nx; i++) 
        {
            float index = i + rowOffset;

            const float u_ij_x = velocityX.getValue(i, j);
            const float u_ij_y = velocityY.getValue(i, j);
            const float u_iNextj = velocityX.getValue(i + 1, j);
            const float u_ijNext = velocityY.getValue(i, j + 1);

            float dudx = (u_iNextj - u_ij_x) / dx;
            float dudy = (u_ijNext - u_ij_y) / dy;

            float divergence = dudx + dudy;

            RHS[index] = -Scene::kDensity * divergence / dt;
        }
    }

    // Needed to solve mathematically
    RHS[0] = 0.0f;

    // Fill A (we only need to fill it once)
    static SparseMatrix<float> A;

    if (A.n != cellsNumber) {
        A = SparseMatrix<float>(cellsNumber, 5); // second parameter is max number of non zero elements
        const float hValue = -1.0f / dx2;
        const float vValue = -1.0f / dy2;

        for (int j = 0; j < Ny; ++j) {
            int rowOffset = j * Nx;

            for (int i = 0; i < Nx; ++i) {
                int row = i + rowOffset;

                float diag = 0.0f;

                if (i > 0) {  // has left neighbor
                    A.set_element(row, row - 1, hValue);
                    diag -= hValue;
                }
                if (i < Nx - 1) {  // has right neighbor
                    A.set_element(row, row + 1, hValue);
                    diag -= hValue;
                }
                if (j > 0) {  // has bottom neighbor
                    A.set_element(row, row - Nx, vValue);
                    diag -= vValue;
                }
                if (j < Ny - 1) {  // has top neighbor
                    A.set_element(row, row + Nx, vValue);
                    diag -= vValue;
                }

                A.set_element(row, row, diag);
            }
        }

        // To solve we need the first row to be [1 0 0 0 0 ...]
        A.set_element(0, 0, 1.0f);
        A.set_element(0, 1, 0.0f);
        A.set_element(0, Nx, 0.0f);
    }        

    // Initialize P to 0 before solving
    static std::vector<float> P;

    if (P.size() != cellsNumber) {
        P.resize(cellsNumber);
    }

    // Solve system for pressure
    float residual;
    int iterations;
    solver.solve(A, RHS, P, residual, iterations);

    // Apply pressure to velocities
    const auto rho = Scene::kDensity;

    for (int j = 0; j < Ny; j++)
    {
        int rowOffset = j * Nx;

        for (int i = 1; i < Nx; i++) 
        {
            int idx_left = (i - 1) + rowOffset;
            int idx_right = i + rowOffset;

            float gradient = (P[idx_right] - P[idx_left]) / dx;

            velocityX.setValue(i, j, velocityX.getValue(i, j) - dt / rho * gradient);
        }
    }

    for (int j = 1; j < Ny; j++) 
    {
        int rowOffset_Bottom = (j - 1) * Nx;
        int rowOffset_Top = j * Nx;
        for (int i = 0; i < Nx; i++) 
        {
            int idx_bottom = i + rowOffset_Bottom;
            int idx_top = i + rowOffset_Top;

            float gradient = (P[idx_top] - P[idx_bottom]) / dy;

            velocityY.setValue(i, j, velocityY.getValue(i, j) - dt / rho * gradient);
        }
    }

    if (flipEnabled) {
        if (apicEnabled) 
        {
            for (uint p_idx = 0; p_idx < particles.getSize(); p_idx++) {
                Vector2 p_pos = particles.getPosition(p_idx);

                Vector2 grid_idx_x = grid.getFaceIndexX(p_pos);
                Vector2 grid_idx_y = grid.getFaceIndexY(p_pos);

                // velocity is just interpolation of current velocity (like in pic)
                Vector2 picVel = Vector2(interpolate(grid_idx_x, velocityX), interpolate(grid_idx_y, velocityY));
                particles.setVelocity(p_idx, picVel);

                // THEORY:
                // APIC stores particle velocity AND information about the state of the fluid around it, with partial derivatives
                // that way it reduces energy loss, specially rotational energy, so that it can capture vorticity, shearing and divergence
                // C is the weighted average of the velocity changes in space, how the fluid velocity changes in the grid points around the particle
                // for that we need the spatial derivates of the interpolation weights 

                // calculate velocity gradient with interpolation weights
                int x0_x = floor(grid_idx_x.x);
                int y0_x = floor(grid_idx_x.y);
                int x1_x = x0_x + 1;
                int y1_x = y0_x + 1;

                float tx_x = grid_idx_x.x - x0_x;
                float ty_x = grid_idx_x.y - y0_x;

                x0_x = clamp(x0_x, 0, grid.getSizeFacesX().x - 1);
                x1_x = clamp(x1_x, 0, grid.getSizeFacesX().x - 1);
                y0_x = clamp(y0_x, 0, grid.getSizeFacesX().y - 1);
                y1_x = clamp(y1_x, 0, grid.getSizeFacesX().y - 1);

                float u00 = velocityX.getValue(x0_x, y0_x);
                float u10 = velocityX.getValue(x1_x, y0_x);
                float u01 = velocityX.getValue(x0_x, y1_x);
                float u11 = velocityX.getValue(x1_x, y1_x);

                int x0_y = floor(grid_idx_y.x);
                int y0_y = floor(grid_idx_y.y);
                int x1_y = x0_y + 1;
                int y1_y = y0_y + 1;

                float tx_y = grid_idx_y.x - x0_y;
                float ty_y = grid_idx_y.y - y0_y;
                
                x0_y = clamp(x0_y, 0, grid.getSizeFacesY().x - 1);
                x1_y = clamp(x1_y, 0, grid.getSizeFacesY().x - 1);
                y0_y = clamp(y0_y, 0, grid.getSizeFacesY().y - 1);
                y1_y = clamp(y1_y, 0, grid.getSizeFacesY().y - 1);

                float v00 = velocityY.getValue(x0_y, y0_y);
                float v10 = velocityY.getValue(x1_y, y0_y);
                float v01 = velocityY.getValue(x0_y, y1_y);
                float v11 = velocityY.getValue(x1_y, y1_y);

                float dudx = partial_cx(1.0f / dx, ty_x, u00, u10, u01, u11);
                float dvdx = partial_cx(1.0f / dx, ty_y, v00, v10, v01, v11);

                float dudy = partial_cy(1.0f / dy, tx_x, u00, u10, u01, u11);
                float dvdy = partial_cy(1.0f / dy, tx_y, v00, v10, v01, v11);

                particles.setCx(p_idx, dudx, dvdx);
                particles.setCy(p_idx, dudy, dvdy);
            }
        
        } 
        else // FLIP 0.95
        { 
            // calculate FLIP velocity delta HERE
            Array2<float> deltaVelocityX(grid.getSizeFacesX());

            for (uint i = 0; i < grid.getSizeFacesX().x; i++)
                for (uint j = 0; j < grid.getSizeFacesX().y; j++)
                    deltaVelocityX.setValue(i, j, velocityX.getValue(i, j) - oldVelocityX.getValue(i, j));

            Array2<float> deltaVelocityY(grid.getSizeFacesY());

            for (uint i = 0; i < grid.getSizeFacesY().x; i++)
                for (uint j = 0; j < grid.getSizeFacesY().y; j++)
                    deltaVelocityY.setValue(i, j, velocityY.getValue(i, j) - oldVelocityY.getValue(i, j));

            // apply PIC/FLIP to update particles velocities HERE
            for (uint p_idx = 0; p_idx < particles.getSize(); p_idx++) {
                Vector2 p_pos = particles.getPosition(p_idx);

                Vector2 grid_idx_x = grid.getFaceIndexX(p_pos);
                Vector2 grid_idx_y = grid.getFaceIndexY(p_pos);

                // FLIP velocity is calculated adding the interpolation of the delta velocity to current velocity
                Vector2 deltaVel = Vector2(interpolate(grid_idx_x, deltaVelocityX), interpolate(grid_idx_y, deltaVelocityY));

                Vector2 flipVel = particles.getVelocity(p_idx) + deltaVel;

                // PIC velocity is just interpolation of current velocity
                Vector2 picVel = Vector2(interpolate(grid_idx_x, velocityX), interpolate(grid_idx_y, velocityY));

                // FLIP 0.95 is a mix of the two velocities
                particles.setVelocity(p_idx, flipVel * 0.95 + picVel * 0.05);
            }
        }
    }
}
}  // namespace asa
