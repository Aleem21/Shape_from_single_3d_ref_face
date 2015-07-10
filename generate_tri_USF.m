function [ tri ] = generate_tri_USF( i,j )
%GENERATE_TRI_USF makes triangulation using delauny tri for USF dataset

[i,j] = meshgrid(i,j);
i = i'; 
j = j'; 
i=i(:);
j=j(:);
tri = delaunayTriangulation(i,j);

end

