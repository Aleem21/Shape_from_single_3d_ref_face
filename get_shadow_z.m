function [ s_ratio,shadow_o,no_shadow_o ] = get_shadow_z( z,im,sh_coeff,compute)
%GET_SHADOW Summary of this function goes here
%   Detailed explanation goes here
% [tri,ptsx,ptsy] = generate_tri_USF(1:size(dmap,1),1:size(dmap,2));
% tri = tri(:,:);
if nargin<4
    compute = 1;
end
a0 = pi;
a1 = 2*pi/3;
a2= 2*pi/8;

c0 =    sqrt(1  /(4*pi));
c1 =    sqrt(3  /(4*pi));
c2 = 3* sqrt(5  /(12*pi));

xrange = [1 size(z,2)];
yrange = [1 size(z,1)];

tri = generate_tri_USF(1:size(z,1),1:size(z,2));
tri = tri(:,:)';
[x,y] = meshgrid(1:size(im,2),1:size(im,1));
y = yrange(2)-y;
z = z(:);
valid = sum(isnan(z(tri)))==0;
tri(:,~valid) = [];
z(isnan(z)) = 0;

tri = tri-1;

pts(1,:) = x(:);
pts(2,:) = y(:);
pts(3,:) = z(:);
% 
FV.vertices = pts';
FV.faces = tri'+1;

addpath(genpath('./toolbox_graph'))
running = 1;
while(running)
    try
        [FV2.vertices,FV2.faces]=perform_mesh_simplification(FV.vertices,FV.faces,20000);
        running = 0;
    catch Err
        pause(1)
        disp('Error writing file. Attempting again')
    end
end
% FV2 = FV;
% FV2 = reducepatch(FV,10000);
pts = FV2.vertices';
tri = FV2.faces'-1;


pts(1:2,:) = pts(1:2,:)-0.5;

%%
tri2 = tri+1;
trin = tri2(:)';
pts(3,:) = pts(3,:);
if compute
    plyhelper.ply_write(pts,tri2,'./spherical_harmonics/temp.ply','ascii')
    system('cse291_exe.exe');
end
%%
n = [];
% for i=1:size(tri2,2)
%     pt = pts(:,tri2(:,i));
%     n((i-1)*3 + (1:3)) = cross(pt(:,2)-pt(:,1),pt(:,3)-pt(:,1));
%     n((i-1)*3 + (1:3)) = n((i-1)*3 + (1:3))/norm(n((i-1)*3 + (1:3)));
% end
%
% model2.FV.facecolor = 'interp';
% model2.FV.facevertexcdata = zeros(size(model2.FV.vertices));
% try
% for i=1:size(tri2,2)
%     model2.FV.facevertexcdata(tri2(1,i)+1,:) = n(tri2(1,i)*3 + (1:3));
%     model2.FV.facevertexcdata(tri2(1,i)+1,:) = model2.FV.facevertexcdata(tri2(1,i)+1,:)/norm(model2.FV.facevertexcdata(tri2(1,i)+1,:));
% end
% catch err
%     disp err
% end
% model2.FV.facevertexcdata = model2.FV.facevertexcdata(:,1);
% normals = zeros(size(model2.FV.vertices));
% for i=1:size(model2.FV.faces,1)
%     is= model2.FV.faces(i,:);
%     norm_cur = [n((is(1)-1)*3 + (0:2) + 1) n((is(2)-1)*3 + (0:2) + 1) n((is(3)-1)*3 + (0:2) + 1)];
%     normals(is,:) = normals(is,:) + reshape(norm_cur,3,3);
% end
% [FV2]=plyhelper.plyread('./spherical_harmonics/temp.ply');
% pts = FV2.vertex;
% tri = FV2.face;
% tri = cell2mat(FV2.face.vertex_indices)';
% trin = tri(:)'-1;
% pts = [FV2.vertex.x FV2.vertex.y FV2.vertex.z]';

% scale xrange space of model into [-1,1]
% pts(3,:) = -pts(3,:);
pts(1,:) = (pts(1,:))/diff(xrange)*2-1;
% % scale yrange space of model into [-1,1]
pts(2,:) = (pts(2,:))/diff(yrange)*2-1;
% % % scale depth space of model into [-1,1]
% pts(3,:) = -pts(3,:);
% % 
% pts([1 3],:) = pts([3 1],:);
minz = min(pts(3,:)); maxz = max(pts(3,:));
pts(3,:) = (pts(3,:)-minz)/(maxz-minz);

ptsn = pts(:)';

cRes = 600;rRes = 600;
rRes = yrange(2);
cRes = xrange(2);


if rRes<200
    rRes_new = 200;
else
    rRes_new = rRes;
end
if cRes<200
    cRes_new = 200;
else
    cRes_new = cRes;
end

sh_coeff2 = sh_coeff' ./[a0*c0  a1*c1  a1*c1  a1*c1  a2*c2  a2*c2  a2*c2  a2*c2/2  a2*c2/2/sqrt(3) ];
% sh_coeff2 = sh_coeff' ./ [c0 c1 c1 c1 c2 c2 c2 c2/2  c2/2/sqrt(3)];
% sh_coeff2 = sh_coeff';
% sh_coeff2(2:4) = sh_coeff2([3 4 2]);
sh_coeff2([2 3 4 5 6 7 8 9]) = sh_coeff2([3 4 2 5 7 9 6 8]);
sh_coeff2([4 5 8]) = -sh_coeff2([4 5 8]);

sh_coeff2 = reshape(sh_coeff2,1,[]);
sh_coeff2(28) = 0;
sh_coeff2(28) = [];
sh_coeff2(10:18) = sh_coeff2(1:9);
sh_coeff2(19:27) = sh_coeff2(1:9);
s_ratio_mean = zeros([rRes_new,cRes_new]);
den = zeros([rRes_new,cRes_new]);
clear valid
for i=-8:15
[r,g,b] = cse291(ptsn,trin-1,sh_coeff2*2^i,n,floor(numel(ptsn)/3),floor(numel(trin)/3),cRes_new,rRes_new,logical(1));img = cat(3,r,g,b);
img = permute(img,[2 1 3]);
img = img(end:-1:1,:,:);
shadow = im2double(img(:,:,1));
no_shadow = im2double(img(:,:,2));
s_ratio = shadow./no_shadow;
valid_i = no_shadow>0.05 & no_shadow<0.95;
den_i = den(valid_i);
s_ratio_mean(valid_i) = s_ratio_mean(valid_i).*den_i./(den_i+1) + s_ratio(valid_i)./(den_i+1);
valid{i+9} = valid_i;
den(valid_i) = den_i+1;
if i==0
    shadow_o = imresize(shadow,[rRes, cRes],'nearest');
    no_shadow_o = imresize(no_shadow,[rRes, cRes],'nearest');
end
end

% 
% shadow_o = shadow;
% no_shadow_o = no_shadow;
shadow = no_shadow - shadow;
shadow(no_shadow<=0) = 1;
shadow(shadow<=0) = 1;
no_shadow(no_shadow<=0) = 1;
s_ratio = shadow./no_shadow;
% figure;imshow(s_ratio);
s_ratio = min(imresize(s_ratio_mean,[rRes, cRes],'nearest'),1);
clear mex;
end
