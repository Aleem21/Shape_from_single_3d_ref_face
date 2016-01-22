function [ costfun,jacobianPattern ] = get_albedo_costfun( depth, im,alb_ref, sh_coeff, eye_mask,lambda2)
%GET_COSTFUN Summary of this function goes here
%   Detailed explanation goes here
%% Pre processing
[ face,face_inds, inface_inds,in_face,r_face,c_face,r_inface,c_inface,...
    b_out_full,b_in_full ] = preprocess_estimate_depth( depth );



rho_ref = alb_ref(face);

% r,c to index number in face Z's
[fr,fc] = meshgrid(1:size(face,1),1:size(face,2));
sub2ind_face = find_elements([r_face c_face],fr',fc');

%% cost function
% Data term
xp = sub2ind_face(sub2ind(size(im),r_face,max(c_face-1,1)));
xn = sub2ind_face(sub2ind(size(im),r_face,c_face));

yp = sub2ind_face(sub2ind(size(im),max(r_face-1,1),c_face));
yn = sub2ind_face(sub2ind(size(im),r_face,c_face));

ind = find(xp==0);
for i = 1:numel(ind)
    xp(ind(i))=  xn(ind(i));
    xn(ind(i))=sub2ind_face(sub2ind(size(im),r_face(ind(i)),c_face(ind(i))+1));
end
ind = find(yp==0);
for i = 1:numel(ind)
    yp(ind(i))=  yn(ind(i));
    yn(ind(i))= sub2ind_face(sub2ind(size(im),r_face(ind(i))+1,c_face(ind(i))));
end
% regularization term
% in_inface = (in_face-b_in_full);
in_inface = in_face;

[r_innface,c_innface] = find(in_inface);
innface_inds = sub2ind(size(face),r_innface,c_innface);


sz = 3; dev = 2;
gauss = fspecial('gaussian',sz,dev);
rhs_reg_mat = lambda2*(alb_ref - conv2(alb_ref,gauss,'same'));

rhs_reg = rhs_reg_mat(innface_inds);
f_w = floor(sz/2);
[boxc, boxr] = meshgrid(-f_w:f_w,-f_w:f_w);
modified_gauss = diag([0 1 0])-gauss; % 1- gauss



gaussVec = lambda2*modified_gauss(:);
iz_reg = zeros(numel(r_innface),numel(boxc));
for i=1:numel(r_innface)
    elems3x3 = sub2ind_face(sub2ind(size(face),boxr(:)+r_innface(i),boxc(:)+c_innface(i)));
    iz_reg(i,:) = elems3x3;
end

if nargout > 1
    % Jacobian Pattern
        %data terms
    nR = sum(face(:))+sum(in_inface(:));
    nC = sum(face(:));
    nOnes = sum(in_face(:))*9 + sum(face(:))*3;
    jacobianPattern = sparse([],[],[],nR,nC,nOnes);
    n_data = sum(face(:));
    jacobianPattern(sub2ind([nR nC],1:n_data,1:n_data)) = 1;
        %regularization terms
    offset = n_data;
    constNumber = repmat(1:size(iz_reg,1),9,1)' + offset;
    jacobianPattern(sub2ind([nR nC],constNumber,...
        iz_reg)) = 1;
end
z = depth(face);
p = z(xp)-z(xn);
q = z(yp)-z(yn);
costfun=@(albedo)cost_nonlin_albedo(albedo,p,q,iz_reg,im(face),rhs_reg,...
    sh_coeff,gaussVec,eye_mask(face));
end

