#version 400
//uniform mat3 rot;
uniform mat3 mvp;
//uniform vec3 d;

in vec3 p;
in vec3 pbrt_d_sh0;
in vec3 pbrt_d_sh1;
in vec3 pbrt_d_sh2;


out vec3 pbrt_d_sh0f;
out vec3 pbrt_d_sh1f;
out vec3 pbrt_d_sh2f;


//in vec3 n;
//out vec3 dir;

void main () {
	vec3 pos_3d = vec3(p);
	//vec3 pos_3d = vec3(mvp*p);
	//vec3 pos_3d = mvp*p;
    //gl_Position = vec4(pos_3d, 1.0);
    //gl_Position = vec4(-0.5+0.7071*pos_3d.x+0.7071*pos_3d.z, -1+pos_3d.y,-0.7071*pos_3d.x+0.7071*pos_3d.z,1.0);
    gl_Position = vec4(pos_3d.x, pos_3d.y,pos_3d.z/3,1.0);
	//nf = n;
	pbrt_d_sh0f = pbrt_d_sh0;
	pbrt_d_sh1f = pbrt_d_sh1;
	pbrt_d_sh2f = pbrt_d_sh2;


	//gl_Position = vec4(p, 1);
	//dir = d;
}