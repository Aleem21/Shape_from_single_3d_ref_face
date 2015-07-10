function [ c ] = render_model_general( ply_path, sh_coeff,Rpose, rRes, cRes, xrange, yrange, talk )
%RENDER_MODEL Assumes the input sh_coeff are actually pre-multiplied with
%an's and ci's such that R = dot(sh_coeff, [1,nx,ny,nz,...]')
import plyread.plyread

if nargin < 8
    talk = 0;
end
if nargin < 7
    yrange = [-1 1];
end
if nargin < 6
    xrange = [-1 1];
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

% read normals
n = [output.vertex.nx';output.vertex.ny';output.vertex.nz'];

% rotate the normals according to the Rpose
n = Rpose(1:3,1:3) * n;
n = -n./repmat(sum(n.^2).^0.5,3,1);
n(2,:) = -n(2,:);
n(3,:) = n(3,:);
% pre-mex conditioning 
pts_scaled = pts;
% scale xrange space of model into [-1,1]
pts_scaled(1,:) = (pts_scaled(1,:)-xrange(1))/diff(xrange)*2-1;
% scale yrange space of model into [-1,1]
pts_scaled(2,:) = (pts_scaled(2,:)-yrange(1))/diff(yrange)*2-1;
% scale depth space of model into [-1,1]
minz = min(pts_scaled(3,:)); maxz = max(pts_scaled(3,:));
pts_scaled(3,:) = (pts_scaled(3,:)-minz)/(maxz-minz)*2-1;
pts_scaled(3,:) = -pts_scaled(3,:);

tri = tri-1;

% nx = n(1,:); ny = n(2,:); nz = n(3,:);=
% Y = [ones(size(nx)); nx; ny; nz; nx.*ny; nx.*nz; ny.*nz; nx.^2-ny.^2; 3*nz.^2-1];
% p = sh_coeff'*Y;
%% render

c = double(mex_render_sh(pts_scaled, tri, n, sh_coeff, cRes, rRes))*2;
c = im2double(c(:,end:-1:1)');
if talk
    figure; imshow(c);
    colormap('gray')
    title('rendered, GL')
end
end

