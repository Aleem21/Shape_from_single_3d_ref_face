#version 400
out vec4 frag_colour;
in vec4 colorV;
void main () {
   frag_colour  = colorV;
}