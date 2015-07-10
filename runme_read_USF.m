clear variables
close all
[pts,tri,rgb] = read_USF_eko('D:\Research UCSD\Ravi\Sony SFS\datasets\USF 3D Face Data\USF Raw 3D Face Data Set\data_files\test',512,512,1);


[rendered,z] = render_rgb_USF( pts,tri,rgb,512,512);

figure;imshow(rendered);
figure;surf(z);
