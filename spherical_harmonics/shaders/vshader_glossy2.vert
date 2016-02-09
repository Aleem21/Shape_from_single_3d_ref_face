#version 400
//uniform mat3 rot;
uniform mat3 mvp;
//uniform vec3 d;

in vec3 p;
in vec3 pbrt_sh0;
in vec3 pbrt_sh1;
in vec3 pbrt_sh2;
in vec3 pbrt_sh3;
in vec3 pbrt_sh4;
in vec3 pbrt_sh5;
in vec3 pbrt_sh6;
in vec3 pbrt_sh7;
in vec3 pbrt_sh8;

out vec3 pbrt_sh0f;
out vec3 pbrt_sh1f;
out vec3 pbrt_sh2f;
out vec3 pbrt_sh3f;
out vec3 pbrt_sh4f;
out vec3 pbrt_sh5f;
out vec3 pbrt_sh6f;
out vec3 pbrt_sh7f;
out vec3 pbrt_sh8f;

//in vec3 n;
//out vec3 dir;

void main () {
    gl_Position = vec4(mvp*(p + vec3(0.025, -0.1, 0)), 1.0);
	pbrt_sh0f = pbrt_sh0;
	pbrt_sh1f = pbrt_sh1;
	pbrt_sh2f = pbrt_sh2;
	pbrt_sh3f = pbrt_sh3;
	pbrt_sh4f = pbrt_sh4;
	pbrt_sh5f = pbrt_sh5;
	pbrt_sh6f = pbrt_sh6;
	pbrt_sh7f = pbrt_sh7;
	pbrt_sh8f = pbrt_sh8;

	//gl_Position = vec4(p, 1);
	//dir = d;
}