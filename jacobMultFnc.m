function [ W ] = jacobMultFnc( J,Y,flag )
%JACOBMULTFNC Summary of this function goes here
%   Detailed explanation goes here
if flag == 0
    W = J'*(J*Y);
elseif flag > 0
    W = J*Y;
else
    W = J'*Y;
end

end

