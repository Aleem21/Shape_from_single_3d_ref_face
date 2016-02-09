#version 400
uniform mat3 rot;
uniform mat3 mvp;
uniform float envmap0[9];
uniform float envmap1[9];
uniform float envmap2[9];
in vec3 p;
in vec3 n;

in float visualizer;

in vec3 pbrt_sh00;
in vec3 pbrt_sh01;
in vec3 pbrt_sh02;
in vec3 pbrt_sh03;
in vec3 pbrt_sh04;
in vec3 pbrt_sh05;
in vec3 pbrt_sh06;
in vec3 pbrt_sh07;
in vec3 pbrt_sh08;
in vec3 pbrt_sh09;
in vec3 pbrt_sh10;
in vec3 pbrt_sh11;
in vec3 pbrt_sh12;
in vec3 pbrt_sh13;
in vec3 pbrt_sh14;
in vec3 pbrt_sh15;
in vec3 pbrt_sh16;
in vec3 pbrt_sh17;
in vec3 pbrt_sh18;
in vec3 pbrt_sh19;
in vec3 pbrt_sh20;
in vec3 pbrt_sh21;
in vec3 pbrt_sh22;
in vec3 pbrt_sh23;
in vec3 pbrt_sh24;
in vec3 pbrt_sh25;
in vec3 pbrt_sh26;

out vec3 basis0;
out vec3 basis1;
out vec3 basis2;

void main () {
    gl_Position = vec4(mvp*(p + vec3(0.025, -0.1, 0)), 1.0);

	basis0 = vec3(0, 0, 0);
	basis1 = vec3(0, 0, 0);
	basis2 = vec3(0, 0, 0);

    basis0[0] =  dot(pbrt_sh00,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                       dot(pbrt_sh01,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                       dot(pbrt_sh02,vec3(envmap0[6],envmap0[7],envmap0[8]));

	basis0[1] =  dot(pbrt_sh03,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                        dot(pbrt_sh04,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                       dot(pbrt_sh05,vec3(envmap0[6],envmap0[7],envmap0[8]));

    basis0[2] =  dot(pbrt_sh06,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                       dot(pbrt_sh07,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                      dot(pbrt_sh08,vec3(envmap0[6],envmap0[7],envmap0[8]));

	basis1[0] =  dot(pbrt_sh09,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                       dot(pbrt_sh10,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                       dot(pbrt_sh11,vec3(envmap0[6],envmap0[7],envmap0[8]));

	//basis1[1] =  dot(pbrt_sh12,vec3(envmap0[0],envmap0[1],envmap0[2])) +
     //                  dot(pbrt_sh13,vec3(envmap0[3],envmap0[4],envmap0[5])) +
       //               dot(pbrt_sh14,vec3(envmap0[6],envmap0[7],envmap0[8]));

	basis1[2] =  dot(pbrt_sh15,vec3(envmap0[0],envmap0[1],envmap0[2])) +
                      dot(pbrt_sh16,vec3(envmap0[3],envmap0[4],envmap0[5])) +
                     dot(pbrt_sh17,vec3(envmap0[6],envmap0[7],envmap0[8]));

	//basis2[0] =  dot(pbrt_sh18,vec3(envmap0[0],envmap0[1],envmap0[2])) +
      //                  dot(pbrt_sh19,vec3(envmap0[3],envmap0[4],envmap0[5])) +
        //                dot(pbrt_sh20,vec3(envmap0[6],envmap0[7],envmap0[8]));

	//basis2[1] =  dot(pbrt_sh21,vec3(envmap0[0],envmap0[1],envmap0[2])) +
      //                  dot(pbrt_sh22,vec3(envmap0[3],envmap0[4],envmap0[5])) +
        //                dot(pbrt_sh23,vec3(envmap0[6],envmap0[7],envmap0[8]));

	//basis2[2] =  dot(pbrt_sh24,vec3(envmap0[0],envmap0[1],envmap0[2])) +
      //                  dot(pbrt_sh25,vec3(envmap0[3],envmap0[4],envmap0[5])) +
        //                dot(pbrt_sh26,vec3(envmap0[6],envmap0[7],envmap0[8]));

}