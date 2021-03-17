#include "vcl/vcl.hpp"
#include <iostream>
#include <list>

#include "simulation.hpp"


using namespace vcl;


struct gui_parameters {
	bool display_color     = true;
    bool display_particles = true;
	bool display_radius    = false;
	bool display_frame		= true;
};

struct user_interaction_parameters {
	vec2 mouse_prev;
	timer_fps fps_record;
	gui_parameters gui;
	bool cursor_on_gui;
};
user_interaction_parameters user;

struct scene_environment
{
	camera_around_center camera;
	mat4 projection;
	vec3 light;
};
scene_environment scene;


void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_click_callback(GLFWwindow* window, int button, int action, int mods);
void window_size_callback(GLFWwindow* window, int width, int height);


void initialize_data();
void update_particles();
void display_scene();
void display_interface();
void update_field_color(grid_3D<vec3>& field, vcl::buffer<particle_element> const& particles);


timer_basic timer;
int flow;
int compt;
int maxLifetime;
int fountain_water_storey;
float fountain_angle;


sph_parameters_structure sph_parameters; // Physical parameter related to SPH
buffer<particle_element> particles;      // Storage of the particles
mesh_drawable sphere_particle; // Sphere used to display a particle
curve_drawable curve_visual;   // Circle used to display the radius h of influence

grid_3D<vec3> field;      // grid used to represent the volume of the fluid under the particles
mesh_drawable field_quad; // quad used to display this field color

mesh_drawable ground;
mesh_drawable fountain_borders;
mesh_drawable fountain_center;



int main(int, char* argv[])
{
	std::cout << "Run " << argv[0] << std::endl;

	GLFWwindow* window = create_window(1280,1024);
	window_size_callback(window, 1280, 1024);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window);
	glfwSetCursorPosCallback(window, mouse_move_callback);
	glfwSetWindowSizeCallback(window, window_size_callback);
	glfwSetMouseButtonCallback(window, mouse_click_callback);
	
	std::cout<<"Initialize data ..."<<std::endl;
	initialize_data();


	std::cout<<"Start animation loop ..."<<std::endl;
	user.fps_record.start();
	timer.start();
	flow= 4;
	compt = 0;
	maxLifetime = 150;
	fountain_water_storey = 0;
	fountain_angle = 0.0f;
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{
		scene.light = scene.camera.position();
		user.fps_record.update();
		timer.update();
		
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);
		imgui_create_frame();
		if(user.fps_record.event) {
			std::string const title = "VCL Display - "+str(user.fps_record.fps)+" fps";
			glfwSetWindowTitle(window, title.c_str());
		}
		user.cursor_on_gui = ImGui::IsAnyWindowFocused();

		ImGui::Begin("GUI",NULL,ImGuiWindowFlags_AlwaysAutoResize);

		float const dt = 0.005f * timer.scale;
        update_particles();
		simulate(dt, particles, sph_parameters);

		display_interface();
		display_scene();
		
		ImGui::End();
		imgui_render_frame(window);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	imgui_cleanup();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}






void initialize_sph()
{
    // Initial particle spacing (relative to h)
    float const c = 0.7f;
	float const h = sph_parameters.h;


	// Fill a cube with particles
    particles.clear();
    /*float const epsilon = 1e-3f;
    for(float x=h; x<1.0f-h; x=x+c*h)
    {
        for(float y=-1.0f+h; y<1.0f-h; y=y+c*h)
        {
			
			particle_element particle;
            particle.p = {x+h/8.0*rand_interval(),y+h/8.0*rand_interval(),-1+h/8.0*rand_interval()}; // a zero value in z position will lead to a 2D simulation
            particle.p = {h/8*rand_interval(),h/8*rand_interval(),-1+h/8*rand_interval()}; // Same initial position
            particle.v = {0.0f,5.0f,0.0f}; // Initial velocity
            particles.push_back(particle);
        }
    }*/

}

void update_particles(){
    float const h = sph_parameters.h;
    std::list<unsigned int> particles_to_recycle = {};
    // Update lifetime
    for(size_t i=0; i<particles.size(); ++i){
        particles[i].lifetime += 1;
        if (particles[i].lifetime>maxLifetime){
            particles_to_recycle.push_back(i);
            particles[i].lifetime = 0;
			particles[i].trajectory.clear();
        }
    }
    // Add or reclycle particle
	if(compt>flow){
		float height;
		if(fountain_water_storey == 0) height = 4.0f;
		else if (fountain_water_storey == 1) height = 7.0f;

		if (particles_to_recycle.size() != 0){
			int first = particles_to_recycle.front();
			particles[first].p = {0.0f,-0.5f,0.0f};
			particles[first].v = {cos(fountain_angle)/50,height,sin(fountain_angle)/50};
			particles_to_recycle.pop_front();
			
		}
		else{
			particle_element particle;
			particle.p = {0.0f,-0.5f,0.0f};
			particle.v = {cos(fountain_angle)/50,height,sin(fountain_angle)/50};
			particle.lifetime = 0;
			particles.push_back(particle);
		}
		compt = 0;
		fountain_water_storey+=1;
		fountain_water_storey = fountain_water_storey % 2;
		fountain_angle+= 2*3.14f/10;
		if(fountain_angle>= 2*3.14) fountain_angle = 0.0f;
	}
	compt ++;
}

void initialize_data()
{
	GLuint const shader_mesh = opengl_create_shader_program(opengl_shader_preset("mesh_vertex"), opengl_shader_preset("mesh_fragment"));
	GLuint const shader_uniform_color = opengl_create_shader_program(opengl_shader_preset("single_color_vertex"), opengl_shader_preset("single_color_fragment"));
	GLuint const texture_white = opengl_texture_to_gpu(image_raw{1,1,image_color_type::rgba,{255,255,255,255}});
	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;
	curve_drawable::default_shader = shader_uniform_color;
	segments_drawable::default_shader = shader_uniform_color;

	scene.camera.look_at({0,0,1.0f}, {0,0,0}, {0,1,0});
	float billboard_size = 0.03f;
    field.resize(10,10,10);
    field_quad = mesh_drawable( mesh_primitive_quadrangle({-billboard_size,-billboard_size,0},
															{billboard_size,-billboard_size,0},
															{billboard_size,billboard_size,0},
															{-billboard_size,billboard_size,0}) );
    field_quad.shading.phong = {1.0f,0,0};
    field_quad.texture = opengl_texture_to_gpu(image_load_png("assets/field.png"));

	user.gui.display_frame = false;


	initialize_sph();
	sphere_particle = mesh_drawable(mesh_primitive_sphere());
	sphere_particle.transform.scale = 0.01f;
	curve_visual.color = {1,0,0};
	curve_visual = curve_drawable(curve_primitive_circle());

    ground = mesh_drawable(mesh_primitive_quadrangle({-1.0f,-1.02f,-1.0f},{-1.0f,-1.02f,1.0f},{1.0f,-1.02f,1.0f},{1.0f,-1.02f,-1.0f}));

	fountain_borders = mesh_drawable(mesh_primitive_cubic_grid({-0.5f,-1.0f,-0.5f},{0.5f,-1.0f,-0.5f},{0.5f,-0.5f,-0.5f},{-0.5f,-0.5f,-0.5f},
																{-0.5f,-1.0f,0.5f},{0.5f,-1.0f,0.5f},{0.5f,-0.5f,0.5f},{-0.5f,-0.5f,0.5f}));
	fountain_center = mesh_drawable(mesh_primitive_cylinder(0.1f,{0,-1,0},{0,-0.5,0},10,20,true));

	
	fountain_borders.shading.alpha = 0.2f;
}


void display_scene()
{
	
	draw(ground, scene);
	draw(fountain_center,scene);
	draw(fountain_borders,scene);

	if(user.gui.display_particles){
		for (size_t k = 0; k < particles.size(); ++k) {
			vec3 const& p = particles[k].p;
			sphere_particle.transform.translate = p;
			draw(sphere_particle, scene);

			//Display trajectories
			particles[k].trajectory.add(p, timer.t);
			trajectory_drawable trajectory;
			 particles[k].trajectory.visual.color = {0.18f,0.82f,0.85f};
			particles[k].trajectory.visual.alpha = 0.4;
			draw(particles[k].trajectory, scene);
		}
	}

	if(user.gui.display_radius){
		curve_visual.transform.scale = sph_parameters.h;
		for (size_t k = 0; k < particles.size(); k+=10) {
			curve_visual.transform.translate = particles[k].p;
			draw(curve_visual, scene);
		}
	}

	if(user.gui.display_color){
        glEnable(GL_BLEND);
        glDepthMask(false);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // Billboard in a grid
        update_field_color(field, particles);

        // Billboard that follows particle
        /*for(size_t k = 0; k < particles.size(); ++k)
        {
            field_quad.transform.translate = particles[k].p;
            field_quad.transform.rotate = scene.camera.orientation();
            field_quad.shading.alpha = 0.1f;
            draw(field_quad, scene);
        }*/
	}
	
		

    

}
void display_interface()
{
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");
	ImGui::SliderInt("Maximum lifetime of particules", &maxLifetime, 100, 800);
	ImGui::SliderInt("Particule creation step", &flow, 1, 10);

	bool const restart = ImGui::Button("Restart");
	if(restart)
		initialize_sph();

	ImGui::Checkbox("Color", &user.gui.display_color);
	ImGui::Checkbox("Particles", &user.gui.display_particles);
	ImGui::Checkbox("Radius", &user.gui.display_radius);
	ImGui::Checkbox("Frame", &user.gui.display_frame);


}


void window_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	float const aspect = width / static_cast<float>(height);
	scene.projection = projection_perspective(50.0f*pi/180.0f, aspect, 2.0f, 200.0f);
}

void mouse_click_callback(GLFWwindow* window, int button, int action, int mods)
{
	ImGui::SetWindowFocus(nullptr);
}

void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
	vec2 const& p0 = user.mouse_prev;
	glfw_state state = glfw_current_state(window);

	auto& camera = scene.camera;
	if(!user.cursor_on_gui){
		if(state.mouse_click_left && !state.key_ctrl)
			scene.camera.manipulator_rotate_trackball(p0, p1);
		if(state.mouse_click_left && state.key_ctrl)
			camera.manipulator_translate_in_plane(p1-p0);
		if(state.mouse_click_right)
			camera.manipulator_scale_distance_to_center( (p1-p0).y );
	}

	user.mouse_prev = p1;
}

void opengl_uniform(GLuint shader, scene_environment const& current_scene)
{
	opengl_uniform(shader, "projection", current_scene.projection);
	opengl_uniform(shader, "view", scene.camera.matrix_view());
	opengl_uniform(shader, "light", scene.light, false);
}

void update_field_color(grid_3D<vec3>& field, vcl::buffer<particle_element> const& particles)
{
    // field.fill({1,1,1});
	// float const d = 0.1f;
	// int const Nf = int(field.dimension.x);
	// for (int kx = 0; kx < Nf; ++kx) {
	// 	for (int ky = 0; ky < Nf; ++ky) {
	// 		for (int kz = 0; kz < Nf; ++kz) {
	// 			float f = 0.0f;
	// 			vec3 const p0 = { 2.0f*(kx/(Nf-1.0f)-0.5f), 2.0f*(ky/(Nf-1.0f)-0.5f), 2.0f*(kz/(Nf-1.0f)-0.5f)};
	// 			for (size_t k = 0; k < particles.size(); ++k) {
	// 				vec3 const& pi = particles[k].p;
	// 				float const r = norm(pi-p0)/d;
	// 				f += 0.25f*std::exp(-r*r);
	// 			}
    //             //field(kx,Nf-1-ky,kz) = vec3(clamp(1-f,0,1),clamp(1-f,0,1),1);

    //             // Add billboard
    //             field_quad.transform.translate = { 2.0f*(kx/(Nf-1.0f)-0.5f), 2.0f*(ky/(Nf-1.0f)-0.5f), 2.0f*(kz/(Nf-1.0f)-0.5f)};
    //             field_quad.transform.rotate = scene.camera.orientation();
    //             field_quad.shading.alpha = clamp(f,0,1);
    //             draw(field_quad, scene);

	// 		}
	// 	}
	// }
	for (int k = 0 ; k< particles.size(); k++){
		trajectory_drawable trajectory = particles[k].trajectory;
		// Add billboard
		for(int k2 = 0 ; k2< trajectory.position_record.size() ; k2+=2){ 
			vec3 position = trajectory.position_record[k2];
			if(position.x == 0 && position.y == 0 && position.z == 0) continue;

			field_quad.transform.translate = trajectory.position_record[k2];
			field_quad.transform.rotate = scene.camera.orientation();
			field_quad.shading.alpha = clamp(0.2f,0,1);
			draw(field_quad, scene);
		}
	}


}

