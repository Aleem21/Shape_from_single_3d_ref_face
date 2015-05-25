function [ gx,gy,g ] = find_countour_normal( boundary )
%FIND_COUNTOUR_NORMAL Summary of this function goes here
%   Detailed explanation goes here

sq2 = 1/sqrt(2);
sq2 = 1;
%derivative helping masks
maskx = [ sq2    0   sq2;
           1     0    1;
          sq2    0   sq2];
masky1= [-sq2    1   sq2;
           0     0    0;
           sq2   1   -sq2];
masky2= [-sq2   -1   sq2;
           0     0    0;
           sq2  -1   -sq2];


       
maskcnt = [1 1 1;
           1 0 1;
           1 1 1];
       
g_cnt = conv2(boundary, maskcnt,'same');
gx = (conv2(boundary, maskx,'same')./g_cnt);
gy(:,:,1) = (conv2(boundary, masky1,'same')./g_cnt);
gy(:,:,2) = (conv2(boundary, masky2,'same')./g_cnt);

[~,i] = max(abs(gy),[],3);
gy = gy(:,:,1).*(i==1) + gy(:,:,2).*(i==2);
gx(~boundary) = NaN;
gy(~boundary) = NaN;

g(:,:,1) = gx;
g(:,:,2) = gy;
g(:,:,3) = 0;
g = g./repmat(sum(g.^2,3).^0.5,[1 1 3]);
gx = g(:,:,1);
gy = g(:,:,2);

end

