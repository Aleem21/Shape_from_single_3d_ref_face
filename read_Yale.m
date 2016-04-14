function [ faces ] = read_Yale( yale_path )
%READ_YALE Summary of this function goes here
%   Detailed explanation goes here
load(yale_path);


for i=1:165
    face{i} = reshape(fea(i,:),64,64);
end
end

