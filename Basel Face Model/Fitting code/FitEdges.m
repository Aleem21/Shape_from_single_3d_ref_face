function [b,R,t,s] = FitEdges( im,xp,landmarks,shapePC,shapeMU,shapeEV,Ef,Ev,tri )
%FITEDGES Perform morphable model fitting using edges
%   This function initialises by fitting to a sparse set of landmarks, then
%   iteratively fits to edges by finding nearest neighbour between
%   projected model edges and image edges
%
% Inputs:
%   im        - Image to fit to
%   xp        - 2 by nlandmarks matrix containing locations of 2D landmarks
%   landmarks - Indices of vertices corresponding to points in xp
%   shapePC   - 3DMM principal components
%   shapeMU   - 3DMM average shape
%   shapeEV   - standard deviations of each 3DMM dimension
%   Ef        - nedges by 2 matrix storing faces adjacent to each edge
%   Ev        - nedges by 2 matrix storing vertices adjacent to each edge
%   tri       - Face structure of 3DMM mesh
%
% Outputs:
%   b         - Estimated shape parameters
%   R,t,s     - Estimated rotation, translation and scale of face

FV.faces = tri;

% PARAMETERS
% Number of model dimensions used for initial landmark fit
initialndims = 5;
% Number of edge fitting iterations (may want to replace with a convergence
% test)
niter = 10;
% Proportion of nearest-neighbour edge matches used at each iteration, may
% be more sensible to threshold on distance to nearest neighbour
percentile = [0.95 0.95 0.95 0.95 0.95 0.99 0.99 0.99 0.99 0.99];
% Number of model dimensions used at each iteration of edge-based fitting
ndims = [50 50 50 50 50 75 75 75 75 75];

% Perform initial landmark-only fit
[b,R,t,s] = FitSingleSOP( xp,shapePC,shapeMU,shapeEV,initialndims,landmarks );
FV.vertices = reshape(shapePC(:,1:initialndims)*b+shapeMU,3,53490)';

% Display input image
figure; imshow(im)

% Extract image edges
edges = edge(rgb2gray(im));
% Probably worth thinking about best edge detection algorithm, parameters
% etc

figure; imshow(edges)

% Store pixel locations of edges
[r,c]=find(edges);
r = size(edges,1)+1-r;

% Open figure for showing development of edges
figure;

for iter=1:niter
    % Compute vertices lying on occluding boundary
    [ occludingVertices ] = occludingBoundaryVertices( FV,Ef,Ev,R );
    
    % Display occluding vertices on top of mesh
    %figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong','AmbientStrength',0,'DiffuseStrength',1,'SpecularStrength',0,'BackFaceLighting','lit'); axis equal; axis off; light('Position',[0 0 -1],'Style','infinite');
    %hold
    %plot3(FV.vertices(occludingVertices,1),FV.vertices(occludingVertices,2),FV.vertices(occludingVertices,3),'.')
    
    % Project occluding boundary vertices
    x2 = R*FV.vertices(occludingVertices,:)';
    x2 = x2(1:2,:);
    x2(1,:) = x2(1,:)+t(1);
    x2(2,:) = x2(2,:)+t(2);
    x2 = x2.*s;
    
    % Find edge correspondences
    [idx,d] = knnsearch([c r],x2');
    
    % Show matches
    %figure; axis equal
    %hold
    %for i=1:length(idx)
    %    plot([c(idx(i)) x2(1,i)],[r(idx(i)) x2(2,i)])
    %end

    % Filter edge matches - probably want to do something better here
    sortedd=sort(d);
    threshold = sortedd(round(percentile(iter)*length(sortedd)));
    idx = idx(d<threshold);
    occludingVertices = occludingVertices(d<threshold);

    % Display projected occluding vertices and edge pixels
    hold off
    plot(c,r,'.'); axis equal
    hold on
    plot(x2(1,d<threshold),x2(2,d<threshold),'.g');
    plot(x2(1,d>=threshold),x2(2,d>=threshold),'.r');
    drawnow
    
    % Refit to new edge landmarks
    % Increase number of model dimensions now that data is denser
    %ndims = 50;
    [b,R,t,s] = FitSingleSOP( [c(idx)' xp(1,:); r(idx)' xp(2,:)],shapePC,shapeMU,shapeEV,ndims(iter),[occludingVertices; landmarks] );
    % Note: this is completely refitting from scratch. We could just
    % re-start the nonlinear optimisation using previous estimates as
    % initialisation but this doesn't seem to work well (perhaps prone to
    % local minima whereas coordinate ascent moves further through the
    % search space?)
    
    % Display current mesh
    %FV.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,53490)';
    %FV.vertices = (R*FV.vertices')';
    %figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong','AmbientStrength',0,'DiffuseStrength',1,'SpecularStrength',0,'BackFaceLighting','lit'); axis equal; axis off; light('Position',[0 0 -1],'Style','infinite');
    
    FV.vertices = reshape(shapePC(:,1:ndims(iter))*b+shapeMU,3,53490)';

end

figure; imshow(im)
hold
plot(x2(1,:),size(im,1)+1-x2(2,:),'.')

%FV.vertices = reshape(shapePC(:,1:ndims)*b+shapeMU,3,53490)';
FV.vertices = (R*FV.vertices')';
figure; patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong','AmbientStrength',0,'DiffuseStrength',1,'SpecularStrength',0,'BackFaceLighting','lit'); 
axis equal; axis off; 
light('Position',[0 0 1],'Style','infinite');
end

