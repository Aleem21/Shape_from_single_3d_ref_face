function [ landmarks ] = stasm_tracker_YaleB( fused,talk )
%STASM_TRACKER_YALEB Summary of this function goes here
%   Detailed explanation goes here
if nargin <2
    talk = 0;
end
landmarks = stasm_tracker(fused,talk);
if isempty(landmarks)
    return
end
landmarks = landmarks - repmat([67;111],1,77);
end

