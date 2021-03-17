std::string s = R"(
#version 330 core

layout(location=0) out vec4 FragColor;
uniform vec3 color = vec3(1.0, 1.0, 1.0); // Unifor color of the object
uniform float alpha = 1.0f; // alpha coefficient

void main()
{
	FragColor = vec4(color, alpha);
}
)";