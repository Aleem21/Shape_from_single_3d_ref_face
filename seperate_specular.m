function [ D,S,chrom_light ] = seperate_specular( im,mask,labels )
%SEPERATE_SPECULAR Summary of this function goes here
%   Detailed explanation goes here
mask(isnan(im(:,:,1))) = 0;
mask3 = repmat(mask,1,1,3);
clrs = im(mask3);
sz = numel(clrs)/3;
clrs2 = reshape(clrs,sz,3);
[~,~,V] = svd(clrs2, 0);
pln = [V(:,end)];
% figure;plot(clrs2*pln);

im_g = rgb2gray(im);
intensity = im_g(mask);
strt = min(intensity);
span = range(intensity);
valid_range = strt+[0.3 0.5]*span;
D_30_50_i = repmat(im_g>valid_range(1) & im_g<valid_range(2),1,1,3);
D_30_50 = zeros(size(im));
D_30_50(D_30_50_i) = im(D_30_50_i);

mean_dif_chrom = [];

for n_i = 1:10
    if n_i==5
        continue
    end
    mask_n_i = repmat((labels==n_i),1,1,3);
%     mask_n_i = repmat(mask,1,1,3);
    clrs_r = reshape(im(mask_n_i),[],3);
    if numel(clrs_r)==0
        continue
    end
    [~,~,V] = svd(clrs_r, 0);
%     no=3;%smallest number of points required
%     k=50;%number of iterations
%     t=1;%threshold used to id a point that fits well
%     d=1000;%number of nearby points required
%     [p_best,n_best,ro_best,X_best,Y_best,Z_best,error_best]=ransac_tim(clrs_r,3,k,t,d);
    
    N(:,n_i) = V(:,end)/norm(V(:,end))*sum(mask_n_i(:))/sum(mask(:));
    
    
    
    dif_colors_r = reshape(im(D_30_50_i.*repmat(labels==n_i,1,1,3)>0),[],3);
    dif_mag_r = sum(dif_colors_r.^2,2).^0.5;
    dif_colors_r(repmat(dif_mag_r==0,1,3))=[];
    dif_mag_r(repmat(dif_mag_r==0,1,3))=[];
    dif_chrom_r = dif_colors_r./repmat(dif_mag_r,1,3);
    mean_dif_chrom = [mean_dif_chrom;dif_chrom_r];
    mean_dif_chrom_r(:,n_i) = mean(dif_chrom_r)';
end
mean_dif_chrom = mean(mean_dif_chrom)';
costfn = @(chr)costfn_plankian(chr,N);
chrom_init = rand(3,1);
chrom_light = fminunc(costfn,chrom_init);
chrom_light = [1;1;1];
chrom_light = chrom_light/norm(chrom_light);

alb_r = mean_dif_chrom_r./repmat(chrom_light,1,10);
alb_mag_r = sum(alb_r.^2).^0.5;
alb_mag_r(alb_mag_r==0) = 1;
alb_r = alb_r./repmat(alb_mag_r,3,1);
S = zeros(size(im,1),size(im,2));
D = zeros(size(im,1),size(im,2));

for r = 1:10
    if r==5
        continue
    end
    B = [chrom_light mean_dif_chrom/norm(mean_dif_chrom)];
    dif_colors_r = reshape(im(repmat(labels==r,1,1,3)>0),[],3);
    S_D = B\dif_colors_r';
    S(labels==r) = S_D(1,:);
    D(labels==r) = S_D(2,:);
end
S = S * rgb2gray(reshape(chrom_light,1,1,3));
D = im - repmat(S,1,1,3);
gauss = fspecial('gaussian',3,1);
D = convn(D,gauss,'same');
end
% im_v = im(:);
% im_v = reshape(im_v,[],3);
% S = im_v*chrom_light;
% S = reshape(S,size(im,1),size(im,2));
% % 
% xLim = [0 1];
% zLim = [0 1];
% [X,Z] = meshgrid(xLim,zLim);
% A = V(1,end);
% B = V(2,end);
% C = V(3,end);
% % A = n_best(1);
% % B = n_best(2);
% % C = n_best(3);
% 
% D = 0;
% Y = (A * X + C * Z + D)/ (-B);
% reOrder = [1 2  4 3];
% figure();patch(X(reOrder),Y(reOrder),Z(reOrder),'b');
% grid on;
% alpha(0.3);
% hold on;
% plot3(clrs_r(:,1),clrs_r(:,2),clrs_r(:,3),'.')