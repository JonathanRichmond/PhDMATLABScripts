%%% plotColorLine
%%% Jonathan Richmond
%%% C: 10 June 2025

function [p] = plotColorLine(x, y, theta4, w)

dx = gradient(x);
dy = gradient(y);
l = sqrt(dx.^2+dy.^2);
nx = -dy./l;
ny = dx./l;
v1 = [x+w*nx, y+w*ny];
v2 = [x-w*nx, y-w*ny];
n = length(x);
vertices = zeros(2*n, 2);
for j = 1:n
    vertices(2*j-1,:) = v1(j,:);
    vertices(2*j,:) = v2(j,:);
end
faces = zeros(n-1, 4);
for j = 1:n-1
    faces(j,:) = [2*j-1, 2*j, 2*j+2, 2*j+1];
end
c = angleColor(theta4);
cVertex = zeros(2*n, 3);
for j = 1:n
    cVertex(2*j-1,:) = c(j,:);
    cVertex(2*j,:) = c(j,:);
end

p = patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', cVertex, 'FaceColor', 'interp', 'EdgeColor', 'none', 'HandleVisibility', 'off');