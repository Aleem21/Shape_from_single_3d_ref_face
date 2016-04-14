function [ cost ] = costfn_plankian( chrom,N )
%COSTFN_PLANKIAN Summary of this function goes here
%   Detailed explanation goes here
% for sun, i.e. 5000K plus temperature black body
chrom = chrom/norm(chrom);
chrom_xy = rgb2xy(chrom);

predicted_y = 3.081758*chrom_xy(1)^3 - 5.8733867*chrom_xy(1)^2 + 3.75112997*chrom_xy(1) - 0.37001483;
N_wt = sqrt(sum(N.^2));
N_wt(N_wt==0)=1;
N = N./repmat(N_wt,3,1);
cost = (predicted_y-chrom_xy(2))^2 + sum((N'* chrom).^2.*N_wt');
% cost = sum((N'* chrom).^2.*N_wt');
end

