#version 400

uniform vec3 d;

in vec3 pbrt_sh0f;
in vec3 pbrt_sh1f;
in vec3 pbrt_sh2f;
in vec3 pbrt_sh3f;
in vec3 pbrt_sh4f;
in vec3 pbrt_sh5f;
in vec3 pbrt_sh6f;
in vec3 pbrt_sh7f;
in vec3 pbrt_sh8f;

//in vec3 basis0;
//in vec3 basis1;
//in vec3 basis2;

out vec4 frag_colour;

void main () {

	vec3 response0 = vec3(0.282095, 0.488603 * d[1], 0.488603 * d[2]);
	vec3 response1 = vec3(0.488603 * d[0], 1.092548 * d[0] * d[1], 1.092548 * d[1] * d[2]);
	vec3 response2 = vec3(0.315392 * (3 * d[2] * d[2] - 1), 1.092548 * d[0] * d[2], 0.546274 * (d[0] * d[0] - d[1] * d[1]));
		
   float intensity_r = dot(pbrt_sh0f, response0) + dot(pbrt_sh1f, response1) + dot(pbrt_sh2f, response2);
   float intensity_g = dot(pbrt_sh3f, response0) + dot(pbrt_sh4f, response1) + dot(pbrt_sh5f, response2);
   float intensity_b = dot(pbrt_sh6f, response0) + dot(pbrt_sh7f, response1) + dot(pbrt_sh8f, response2);
   
   frag_colour = vec4(intensity_r, intensity_g, intensity_b, 1);
   //frag_colour = vec4(pbrt_sh0f, 1);
}