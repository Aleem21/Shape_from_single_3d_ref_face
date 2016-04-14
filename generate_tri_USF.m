function [ tri,i,j ] = generate_tri_USF( i,j )
%GENERATE_TRI_USF makes triangulation using delauny tri for USF dataset

[i,j] = meshgrid(i,j);
i = i'; 
j = j'; 
i2=i(:);
j2=j(:);
tri = delaunayTriangulation(i2,j2);

end

