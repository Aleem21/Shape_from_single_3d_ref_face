#version 400

uniform vec3 d;
uniform float envmap0[9];
uniform float envmap1[9];
uniform float envmap2[9];

uniform float diffusiveness;
uniform float glossiness;

in vec3 pbrt_d_sh0f;
in vec3 pbrt_d_sh1f;
in vec3 pbrt_d_sh2f;

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

    float intensity_r =  max(dot(pbrt_d_sh0f,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                        dot(pbrt_d_sh1f,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                        dot(pbrt_d_sh2f,vec3(envmap0[6],envmap0[7],envmap0[8])), 0) * diffusiveness;

    float intensity_g =  max(dot(pbrt_d_sh0f,vec3(envmap1[0],envmap1[1],envmap1[2])) +
                        dot(pbrt_d_sh1f,vec3(envmap1[3],envmap1[4],envmap1[5])) +
                        dot(pbrt_d_sh2f,vec3(envmap1[6],envmap1[7],envmap1[8])), 0) * diffusiveness;

    float intensity_b =  max(dot(pbrt_d_sh0f,vec3(envmap2[0],envmap2[1],envmap2[2])) +
                        dot(pbrt_d_sh1f,vec3(envmap2[3],envmap2[4],envmap2[5])) +
                        dot(pbrt_d_sh2f,vec3(envmap2[6],envmap2[7],envmap2[8])), 0) * diffusiveness;

    //colorV = vec4(intensity0*3, intensity1*3, intensity2*3, 1);

	vec3 response0 = vec3(0.282095, 0.488603 * d[1], 0.488603 * d[2]);
	vec3 response1 = vec3(0.488603 * d[0], 1.092548 * d[0] * d[1], 1.092548 * d[1] * d[2]);
	vec3 response2 = vec3(0.315392 * (3 * d[2] * d[2] - 1), 1.092548 * d[0] * d[2], 0.546274 * (d[0] * d[0] - d[1] * d[1]));
		
   intensity_r += max(0.5*(dot(pbrt_sh0f, response0) + dot(pbrt_sh1f, response1) + dot(pbrt_sh2f, response2)), 0) * glossiness;
   intensity_g += max(0.5*(dot(pbrt_sh3f, response0) + dot(pbrt_sh4f, response1) + dot(pbrt_sh5f, response2)), 0) * glossiness;
   intensity_b += max(0.5*(dot(pbrt_sh6f, response0) + dot(pbrt_sh7f, response1) + dot(pbrt_sh8f, response2)), 0) * glossiness;
   
   frag_colour = vec4(3*intensity_r, 3*intensity_g, 3*intensity_b, 1);
   //frag_colour = vec4(pbrt_sh0f[0], 0, 0, 1);
   //frag_colour = vec4(pbrt_sh0f[0], pbrt_sh5f[0], pbrt_sh4f[2], 1);
}