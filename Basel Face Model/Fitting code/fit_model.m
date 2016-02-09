function [ FV,b,R_est,t_est,s_est ] = fit_model( im,talk )
if nargin <2
    talk = 0;
end
%% Setup basic parameters:
landmarks=[8333,7301,6011,9365,10655,8320,8311,8286,8275,5959,4675,3642,4922,3631,2088,27391,39839,40091,40351,6713,10603,11641,12673,11244,12661,14472,27778,41804,41578,41310,9806,8345,7442,6936,5392,7335,7851,8354,9248,9398,11326,9399,9129,9406,9128,8890,8367,7858,7580,7471,8374,23002,32033]';
load('01_MorphableModel.mat')
load('BFMedgestruct.mat')


%% Landmark Detection)
xp= RamananDetector(im);

%% Display fitted model using xp :
numsd=1;    %Number of standard deviations

[b,R_est,t_est,s_est] = FitEdges( im,xp,landmarks,shapePC,shapeMU,shapeEV,Ef,Ev,tl );
ndims = length(b);
FV.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,53490)';
FV.faces = tl;
if talk
    figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong','AmbientStrength',0,'DiffuseStrength',1,'SpecularStrength',0); axis equal; axis off; light('Position', [0 0 1], 'Style', 'infinite')
end

end

