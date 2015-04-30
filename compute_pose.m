function [ sRT,S ] = compute_pose( pp,Pp )
%COMPUTE_POSE Given point correspondences between 2D(p) and 3D(P), computes
%rough pose of the 3D. Assumptions are orthographic projection.
TP = eye(4);
TP(1:3,4) = -mean(Pp,2);
Tp = eye(4);
Tp(1:2,4) = -mean(pp,2);

p = pp - repmat(mean(pp,2),1,size(pp,2));
P = Pp - repmat(mean(Pp,2),1,size(Pp,2));

sA = p*P'/(P*P');
S =sqrt(sum(sA.^2,2));
A = sA./repmat(S,1,size(sA,2));
R = A;
R(3,:) = cross(R(1,:),R(2,:));
R(4,4)=1;
t = mean(pp,2)-A*mean(Pp,2);
t(3) = 0;
T = [[eye(3) t]; [ 0 0 0 1]];
sRT = Tp\diag([S; 1; 1]) * R * TP;

end

