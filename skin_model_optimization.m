function [ z_o,alb_o,n_o,rho_o,L_sh_o,chrom_o,delta_a_o,delta_s1_o,P_o ] =...
    skin_model_optimization( is_face,L_sh_mean,L_sh_var,z_ref,chrom_p,...
    alb_init,n_init,rho_init,L_sh_init,D,S,w_o,labels,im,alb_mean,alb_var,...
    n_mean,n_var,rho_mean,rho_var,delta_a_init,delta_s1_init,P_init,rad,type )
%SKIN_MODEL_OPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here
is_face3 = repmat(is_face,1,1,3);
labels3 = repmat(labels,1,1,3);
sz = [sum(is_face(:)) numel(L_sh_init)];


 




theta_init = [n_init; rho_init];
chrom_init = chrom_p;
z_init = z_ref(is_face);
params = [z_init(:); alb_init; theta_init(:); L_sh_init(:); chrom_init(:);...
    delta_a_init; delta_s1_init; P_init];
labels3 = labels3.*is_face3;






% wr,wc weights albedo smoothness
wr{1} = (1 + (labels(:,1:end-1)~=labels(:,2:end))*0.5).*is_face(:,1:end-1).*is_face(:,2:end);
wr{1} = wr{1}(:);
wr{2} = (1 + (labels(2:end,:)~=labels(1:end-1,:))*0.5).*is_face(1:end-1,:).*is_face(2:end,:);
wr{2} = wr{2}(:);

chrom = reshape(im,[],3);
chrom = chrom./repmat(sqrt(sum(chrom.^2,2)),1,3);
chrom = reshape(chrom,size(im));
chroml = reshape(chrom(:,1:end-1,:),[],3);
chromr = reshape(chrom(:,2:end,:),[],3);
wc{1} = (1 - (dot(chroml',chromr')<0.9)*0.6);
chromt = reshape(chrom(1:end-1,:,:),[],3);
chromb = reshape(chrom(2:end,:,:),[],3);
wc{2} = (1 - (dot(chromt',chromb')<0.9)*0.6);





% Jacobian pattern
sz_fc = sum(is_face(:));
[r_face,c_face] = find(is_face);
[boxc, boxr] = meshgrid(-2:1,-2:1);

[boxc3, boxr3] = meshgrid(-1:1,-1:1);
sub2ind_face = zeros(size(is_face));
sub2ind_face(is_face) = 1:sum(is_face(:));
iz_diff = zeros(15,sum(is_face(:)));
iz_spec = zeros(3,sum(is_face(:)));
iz_alb = zeros(9,sum(is_face(:)));
inds3x3 = zeros(9,sum(is_face(:)));
for i=1:sum(is_face(:))
    %     belems4x4 = sub2ind_face(sub2ind(size(is_face),min(max(bboxr(2:end)+r_face(i),1),size(is_face,1)),min(max(bboxc(2:end)+c_face(i),1),size(is_face,2))));
    elems4x4 = sub2ind_face(sub2ind(size(is_face),min(max(boxr(2:end)+r_face(i),1),size(is_face,1)),min(max(boxc(2:end)+c_face(i),1),size(is_face,2))));
    
    elems3x3 = sub2ind_face(sub2ind(size(is_face),boxr3(:)+r_face(i),boxc3(:)+c_face(i)));
    iz_spec(:,i) = elems4x4([6 9 10]) + elems4x4([14 11 10]).*(elems4x4([6 9 10])==0);
    elems4x4(elems4x4==0) = repmat(elems4x4(10),sum(elems4x4==0),1);
    
    iz_diff(:,i) = elems4x4;
    elems3x3nan = elems3x3;
    elems3x3nan(elems3x3nan==0) = nan;
    inds3x3(:,i) = elems3x3nan;
    elems3x3(elems3x3==0) = repmat(elems3x3(5),sum(elems3x3==0),1);
    iz_alb(:,i) = elems3x3;
end

% Diffuse data term
P_offset = numel(params);
chrom_offset = 6*sz_fc  + numel(L_sh_mean);
sh_offset = 6*sz_fc;

alb_offset = sz_fc;
r_1_1 = repmat(1:size(iz_diff,2),15+9+1+3+numel(L_sh_init),1);
c_1_1 = [iz_diff; alb_offset + iz_alb ; ones(1,size(iz_diff,2))*P_offset;repmat((1:3)',1,size(iz_diff,2))+chrom_offset; repmat((1:numel(L_sh_init))',1,size(iz_diff,2))+sh_offset; ];
r_1_2 = r_1_1(end) + repmat(1:size(iz_diff,2),15+9+1+3+numel(L_sh_init),1);
c_1_2 = [iz_diff; alb_offset + iz_alb + sz_fc ; ones(1,size(iz_diff,2))*P_offset;repmat((1:3)',1,size(iz_diff,2))+chrom_offset; repmat((1:numel(L_sh_init))',1,size(iz_diff,2))+sh_offset;];
r_1_3 = r_1_2(end) + repmat(1:size(iz_diff,2),15+9+1+3+numel(L_sh_init),1);
c_1_3 = [iz_diff; alb_offset + iz_alb + sz_fc*2 ; ones(1,size(iz_diff,2))*P_offset;repmat((1:3)',1,size(iz_diff,2))+chrom_offset; repmat((1:numel(L_sh_init))',1,size(iz_diff,2))+sh_offset;];

% specular data term
r_2 = repmat(1:size(iz_spec,2),3+1+1+numel(L_sh_mean)+1,1) + r_1_3(end);
c_2 = [iz_spec; iz_spec(3,:)+sz_fc*4; iz_spec(3,:)+sz_fc*5; repmat((1:numel(L_sh_mean))'+sz_fc*6,1,sz_fc); ones(1,sz_fc)*P_offset];

% absorption term
abs_offset = 6*sz_fc + numel(L_sh_mean) + 3;
r_3 = r_2(end) + ones(1,6);
c_3 = abs_offset + (1:6);

% chrom prior
r_4 = r_3(end) + (1:3);
c_4 = chrom_offset + (1:3);

if type
    J = sparse([r_1_1(:);r_1_2(:);r_1_3(:);r_2(:);r_3(:);r_4(:)],...
        [c_1_1(:);c_1_2(:);c_1_3(:);c_2(:);c_3(:);c_4(:)],...
        1,r_4(end),numel(params));
else
    J = sparse([r_3(:);r_4(:)]-r_2(end),...
        [c_3(:);c_4(:)],...
        1,4,numel(params));
end
J(J>0) = 1;

% algo =  {'levenberg-marquardt',.005};
algo = 'trust-region-reflective';
options = optimset('MaxIter',15,'JacobPattern',J,'Display','iter-detailed','Algorithm',algo,'Jacobian','on');

costfn_opti = @(params)costfn_skin_model(params, z_ref(is_face),chrom_p,D(is_face3),S(is_face),w_o,...
    labels3,wr,wc,alb_mean,alb_var,L_sh_mean,L_sh_var,n_mean,n_var,rho_var,rho_mean,iz_diff,sz,is_face,is_face3,rad,type,inds3x3,options );
% contsfn = inline('0');
% options = optimoptions('fminunc','GradObj','on','Display','iter-detailed');
% optimum = fmincon(costfn_opti,params,ones(1,numel(params)),1,[],[],[],[],[],options);
% optimum = fminunc(costfn_opti,params,options);

optimum = lsqnonlin(costfn_opti,params,[],[],options);

figure;plot(params-optimum);
z_o = nan(size(z_ref));
z_o(is_face) = optimum(1:sz_fc);
alb_o = zeros(size(im));
alb_o(is_face3) = optimum(sz_fc+1:4*sz_fc);
n_o = nan(size(z_ref));
n_o(is_face) = optimum(4*sz_fc+1:5*sz_fc);
rho_o = nan(size(z_ref));
rho_o(is_face) = optimum(5*sz_fc+1:6*sz_fc);
L_sh_o = optimum(6*sz_fc+[1:numel(L_sh_mean)]);
chrom_o = optimum(6*sz_fc+numel(L_sh_mean)+[1:3]);
delta_a_o = optimum(end-6:end-4);
delta_s1_o = optimum(end-3:end-1);
P_o = optimum(end);
% params = [z_init(:); alb_init; theta_init(:); L_sh_init(:); chrom_init(:); delta_a; delta_s1];

% options = optimset('MaxIter',2);
figure;
subplot(2,1,1);
ci = costfn_opti(params).^2;
plot(ci);
% plot(1:r_1_3(end),ci(1:r_1_3(end)),r_2(1):r_2(end),ci(r_2(1):r_2(end)),...
%     r_3(1):r_3(end),ci(r_3(1):r_3(end)),r_4(1):r_8(end),ci(r_4(1):r_8(end)),...
%     r_9(1):r_9(end),ci(r_9(1):r_9(end)),r_9(end):r_11(1),ci(r_9(end):r_11(1)),r_11(1):r_15(end),ci(r_11(1):r_15(end)));
subplot(2,1,2);
co = costfn_opti(optimum).^2;
plot(co)
% plot(1:r_1_3(end),co(1:r_1_3(end)),r_2(1):r_2(end),co(r_2(1):r_2(end)),...
%     r_3(1):r_3(end),co(r_3(1):r_3(end)),r_4(1):r_8(end),co(r_4(1):r_8(end)),...
%     r_9(1):r_9(end),co(r_9(1):r_9(end)),r_9(end):r_11(1),co(r_9(end):r_11(1)),r_11(1):r_15(end),co(r_11(1):r_15(end)));

end

