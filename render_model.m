function [ c ] = render_model( ply_path, sh_coeff, talk,Rpose, cRes, rRes )
%RENDER_MODEL Summary of this function goes here
%   Detailed explanation goes here
import plyread.plyread
a0 = pi;
a1 = 2*pi/sqrt(3);
a2= 2*pi/sqrt(8);
if nargin < 4
    Rpose = eye(4);
end
if nargin < 5
    cRes = 400;
end
if nargin < 6
    rRes = 400;
end
%% read model
[tri,pts,output] = plyread(ply_path,'tri');

%% data conditioning

%pose correction
pts = [pts ones(size(pts,1),1)]';
% Rpose = makehgtform('scale',0.005);

pts = Rpose * pts;
pts = pts(1:3,:);

tri = tri';
pts = pts;

n = [output.vertex.ny';-output.vertex.nx';output.vertex.nz'];
n = n./repmat(sum(n.^2).^0.5,3,1);
% pre-mex conditioning 
ptsn = pts;
% ptsn(1,:) = (ptsn(1,:)/cRes-0.5)*2;
% ptsn(2,:) = (ptsn(2,:)/rRes-0.5)*2;
% ptsn(1,:) = ptsn(1,:);
% ptsn(2,:) = ptsn(2,:);
ptsn(3,:) = -ptsn(3,:);

tri = tri-1;

% nx = n(1,:); ny = n(2,:); nz = n(3,:);
% Y = [ones(size(nx)); nx; ny; nz; nx.*ny; nx.*nz; ny.*nz; nx.^2-ny.^2; 3*nz.^2-1];
% p = sh_coeff'*Y;
%% render

c = mex_render_sh(ptsn, tri, n, sh_coeff, cRes, rRes)*2;
c = im2double(c(:,end:-1:1)');
if talk
    figure; imshow(c);
end
end

