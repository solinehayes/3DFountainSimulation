#pragma once

#include "vcl/vcl.hpp"



// SPH Particle
struct particle_element
{
    vcl::vec3 p; // Position
    vcl::vec3 v; // Speed
    vcl::vec3 f; // Force

    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    int lifetime; // add lifetime to remove particle after a certain time (to avoid slowing the simulation)

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},rho(0),pressure(0),lifetime(0) {}
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.12f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

     // Total mass of a particle (consider rho0 h^2)
    float m = rho0*h*h;

    // viscosity parameter
    float nu = 0.02f;
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
};


void simulate(float dt, vcl::buffer<particle_element>& particles, sph_parameters_structure const& sph_parameters);
