function [ landmarks ] = stasm_tracker( img,talk )
%STASM_TRACKER Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    talk = 0;
end
imwrite(img,'.\data\temp.jpg');
[~, output] = system(['stasm4\minimal  .\data\temp.jpg .\stasm4']);
delete('.\data\temp.jpg')
landmarks = str2num(output);
if isempty(landmarks)
    if talk
        disp(output)
    end
    landmarks = [];
    return
end
landmarks = reshape(landmarks,2,77);

end

