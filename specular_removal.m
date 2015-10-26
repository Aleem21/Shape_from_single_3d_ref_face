function [ SUV ] = specular_removal( input,source_color )
%SPECULAR_REMOVAL Summary of this function goes here
%   Detailed explanation goes here
Rot(1,:) = source_color;
temp = rand(1,3);temp=temp/norm(temp);
temp = cross(source_color,temp);
Rot(2,:) = temp/norm(temp);
Rot(3,:) = cross(Rot(1,:),Rot(2,:));

R = input(:,:,1);
G = input(:,:,2);
B = input(:,:,3);
pixRGB = [R(:)';G(:)';B(:)'];
pixSUV = Rot*pixRGB;
SUV(:,:,1) = reshape(pixSUV(1,:),size(input,1),size(input,2));
SUV(:,:,2) = reshape(pixSUV(2,:),size(input,1),size(input,2));
SUV(:,:,3) = reshape(pixSUV(3,:),size(input,1),size(input,2));
end

