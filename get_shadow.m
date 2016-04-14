function [ s_ratio ] = get_shadow( model,im,sh_coeff)
%GET_SHADOW Summary of this function goes here
%   Detailed explanation goes here
% [tri,ptsx,ptsy] = generate_tri_USF(1:size(dmap,1),1:size(dmap,2));
% tri = tri(:,:);


a0 = pi;
a1 = 2*pi/3;
a2= 2*pi/8;

c0 =    sqrt(1  /(4*pi));
c1 =    sqrt(3  /(4*pi));
c2 = 3* sqrt(5  /(12*pi));

load('mask')

xrange = [1 size(im,2)];
yrange = [1 size(im,1)];
% pts = [ptsx(:) ptsy(:)];
% pts(:,3)  = dmap(:);
valid = sum(mask(model.FV.faces),2)==0;
model.FV.faces(~valid,:)=[];

% FV = reducepatch(model.FV,0.25);
model2 = model;
% model2.FV = FV;

[pts,tri] = model2GL(model2,mean(xrange),mean(yrange));




tri2 = tri+1;
valid = sum(mask(tri2))==0;
tri2 = tri2(:,valid)-1;
trin = tri2(:)';
pts(3,:) = -pts(3,:);
plyhelper.ply_write(pts,tri2+1,'./spherical_harmonics/temp.ply','ascii')
system('cse291_exe.exe');
%%
n = [];
for i=1:size(tri2,2)
    pt = pts(:,tri2(:,i)+1);
    n((i-1)*3 + (1:3)) = cross(pt(:,2)-pt(:,1),pt(:,3)-pt(:,1));
    n((i-1)*3 + (1:3)) = n((i-1)*3 + (1:3))/norm(n((i-1)*3 + (1:3)));
end
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
% trin = tri(:)';
% pts = [FV2.vertex.x FV2.vertex.y FV2.vertex.z]';

pts(3,:) = -pts(3,:);
% scale xrange space of model into [-1,1]
pts(1,:) = (pts(1,:)-xrange(1))/diff(xrange)*2-1;
% scale yrange space of model into [-1,1]
pts(2,:) = (pts(2,:)-yrange(1))/diff(yrange)*2-1;
% scale depth space of model into [-1,1]
minz = min(pts(3,:)); maxz = max(pts(3,:));
pts(3,:) = (pts(3,:)-minz)/(maxz-minz);
%pts(3,:) = -pts(3,:);
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

% sh_coeff2 = sh_coeff' ./ [a0*c0  a1*c1  a1*c1  a1*c1  a2*c2  a2*c2  a2*c2  a2*c2/2  a2*c2/2/sqrt(3) ];
sh_coeff2 = sh_coeff' ./ [c0 c1 c1 c1 c2 c2 c2 c2/2  c2/2/sqrt(3)];
%sh_coeff2 = sh_coeff';
sh_coeff2(2:4) = sh_coeff2([3 4 2]);
sh_coeff2(3) = -sh_coeff2(3);
sh_coeff2(5) = -sh_coeff2(5);
sh_coeff2([6 7 8 9]) = sh_coeff2([7 9 6 8]);

sh_coeff2 = reshape(sh_coeff2,1,[]);
sh_coeff2(28) = 0;
sh_coeff2(28) = [];

[r,g,b] = cse291(ptsn,trin,sh_coeff2,n,floor(numel(ptsn)/3),floor(numel(trin)/3),cRes_new,rRes_new,logical(1));img = cat(3,r,g,b);
img = permute(img,[2 1 3]);
img = img(end:-1:1,end:-1:1,:);
shadow = im2double(imresize(img(:,:,1),[rRes cRes]));
no_shadow = im2double(imresize(img(:,:,2),[rRes cRes]));
figure;imshow(shadow);

shadow(no_shadow==0) = 1;
no_shadow(no_shadow==0) = 1;
s_ratio = shadow./no_shadow;
end
