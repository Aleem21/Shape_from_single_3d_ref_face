#version 400
uniform mat3 rot;
uniform mat3 mvp;
uniform float envmap0[9];
uniform float envmap1[9];
uniform float envmap2[9];
in vec3 p;
in vec3 n;

in float visualizer;

in vec3 pbrt_sh0;
in vec3 pbrt_sh1;
in vec3 pbrt_sh2;
out vec4 colorV;
void main () {
    gl_Position = vec4(mvp*(p + vec3(0.025, -0.1, 0)), 1.0);
    float intensity0 =  dot(pbrt_sh0,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                        dot(pbrt_sh1,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                        dot(pbrt_sh2,vec3(envmap0[6],envmap0[7],envmap0[8]));

	//float intensity0 =  dot(pbrt_sh0,vec3(1,0,0)) +
    //                    dot(pbrt_sh1,vec3(0,0,0)) +
    //                    dot(pbrt_sh2,vec3(0,0,0));

    float intensity1 =  dot(pbrt_sh0,vec3(envmap1[0],envmap1[1],envmap1[2])) +
                        dot(pbrt_sh1,vec3(envmap1[3],envmap1[4],envmap1[5])) +
                        dot(pbrt_sh2,vec3(envmap1[6],envmap1[7],envmap1[8]));

    float intensity2 =  dot(pbrt_sh0,vec3(envmap2[0],envmap2[1],envmap2[2])) +
                        dot(pbrt_sh1,vec3(envmap2[3],envmap2[4],envmap2[5])) +
                        dot(pbrt_sh2,vec3(envmap2[6],envmap2[7],envmap2[8]));
    colorV = vec4(intensity0*3, intensity1*3, intensity2*3, 1);
    //colorV = vec4(envmap0[0],envmap0[1],envmap0[2], 1);

	//colorV = vec4(pbrt_sh0[1], pbrt_sh0[1], pbrt_sh0[1], 1);
	//colorV = vec4(visualizer, pbrt_sh0[1], pbrt_sh0[1], 1);
}