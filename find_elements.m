function [ out ] = find_elements( list,r,c )
%FIND_ELEMENTS Summary of this function goes here
%   Detailed explanation goes here

if size(list,2)>2
    list = list';
end

[~, inds] = ismember(list,[r(:) c(:)],'rows');
out = zeros(numel(r),1);
try
for i=1:numel(r)
    t = find(inds==i);
    if ~isempty(t)
        out(i) = t;
    end
end
out = reshape(out,size(r));
catch err
    disp(1)
    
end
