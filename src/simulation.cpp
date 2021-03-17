#include "simulation.hpp"

using namespace vcl;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
	float const r = norm(p_i-p_j);
	assert_vcl_no_msg(r<=h);
	return 45/(3.14159f*std::pow(h,6.0f))*(h-r);
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
	float const r = norm(p_i-p_j);
	assert_vcl_no_msg(r<=h);
	return -45/(3.14159f*std::pow(h,6.0f))*std::pow(h-r,2)*(p_i-p_j)/r;
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    assert_vcl_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}


void update_density(buffer<particle_element>& particles, float h, float m)
{
    size_t const N = particles.size();

    for(size_t i=0; i<N; ++i)
        particles[i].rho = 0.0f;

    for(size_t i=0; i<N; ++i)
    {
        for(size_t j=0; j<N; ++j)
        {
            vec3 const& pi=particles[i].p;
            vec3 const& pj=particles[j].p;

            float const r = norm(pi-pj);
            if(r<h)
                particles[i].rho += m * W_density(pi,pj,h);
        }
    }
}

// Convert the particle density to pressure
void update_pressure(buffer<particle_element>& particles, float rho0, float stiffness)
{
	const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
void update_force(buffer<particle_element>& particles, float h, float m, float nu)
{
	// gravity
    const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
        particles[i].f = m * vec3{0,-9.81f,0};

    for(size_t i=0; i<N; ++i)
    {
        for(size_t j=0; j<N; ++j)
        {
            if (i==j)
                continue;

            const vec3& pi = particles[i].p;
            const vec3& pj = particles[j].p;
            float r = norm(pi-pj);

            if(r<h)
            {
                const vec3& vi = particles[i].v;
                const vec3& vj = particles[j].v;

                const float pressure_i = particles[i].pressure;
                const float pressure_j = particles[j].pressure;

                const float rho_i = particles[i].rho;
                const float rho_j = particles[j].rho;

				vec3 force_pressure = {0,0,0};
				vec3 force_viscosity = {0,0,0};

				force_pressure = - m/rho_i * m*(pressure_i+pressure_j)/(2*rho_j)*W_gradient_pressure(pi,pj,h);
				force_viscosity = nu * m * m * (vj-vi)/rho_j * W_laplacian_viscosity(pi,pj,h);

				
                particles[i].f += force_pressure + force_viscosity;
            }

        }
    }

}

void simulate(float dt, buffer<particle_element>& particles, sph_parameters_structure const& sph_parameters)
{

	// Update values
    update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces

	// Numerical integration
	float const damping = 0.005f;
	size_t const N = particles.size();
	float const m = sph_parameters.m;
	for(size_t k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}


	// Collision
    float const epsilon = 1e-3f;
    for(size_t k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_interval();  v.y *= -0.5f;}
        if( p.x<-0.5 ) {p.x = -0.5+epsilon*rand_interval();  v.x *= -0.5f;}
        if( p.x>0.5 )  {p.x =  0.5-epsilon*rand_interval();  v.x *= -0.5f;}
        if( p.z<-0.5 ) {p.z = -0.5+epsilon*rand_interval();  v.z *= -0.5f;}
        if( p.z>0.5 )  {p.z =  0.5-epsilon*rand_interval();  v.z *= -0.5f;}
    }

}