function [ out ] = cell2im( inp,labels,is_face )
%CELL2IM Summary of this function goes here
%   Detailed explanation goes here

out = zeros(size(is_face));
for i=1:size(is_face,3)
    temp = zeros(size(is_face,1),size(is_face,2));
    for r=1:10
        temp(labels(:,:,i)==r) = inp{r}(i);
    end
    out(:,:,i) = temp;
end
out = out.*is_face;
out = out(is_face);

end

