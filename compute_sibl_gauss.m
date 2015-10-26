function [ sh_mean,sigma ] = compute_sibl_gauss( impath,order,is_gray )
%COMPUTE_SIBL_GAUSS Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    is_gray = 0;
end
try
    load([impath '/mean_var_cache.mat' ])
catch
    impaths = getAllFiles(impath,1);
    len = (order+1)^2;
    channels = 1+2*(is_gray<1);
    sh_vec = zeros(len*channels,numel(impaths));
    for i=1:numel(impaths)
        disp(i)
        i_path = impaths{i};
        if strcmp(i_path(end-2:end),'mat')
            continue
        end
        i_env = imresize(imread(i_path),1);
        if is_gray
            i_env = rgb2gray(i_env);
        end
        
        sh=env2sh(i_env,order,0);
        sh_vec(:,i) = sh(:);
    end
    
    
    sh_mean = mean(sh_vec,2);
    centered = sh_vec-repmat(sh_mean,1,size(sh_vec,2));
    sigma = centered*centered'/size(sh_vec,2);
    
    sh_mean = reshape(sh_mean,len,channels);
    save([impath '/mean_var_cache.mat' ],'sh_mean','sigma')
    
end