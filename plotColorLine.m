%%% plotColorLine
%%% Jonathan Richmond
%%% C: 10 June 2025
%%% U: 12 June 2025

function [p] = plotColorLine(x, y, z, theta4, w)

dx = gradient(x);
dy = gradient(y);
dz = gradient(z);
l = sqrt(dx.^2+dy.^2+dz.^2);
dx = dx./l;
dy = dy./l;
dz = dz./l;
vx = zeros(size(dx));
vy = zeros(size(dy));
vz = ones(size(dz));
nx = dy.*vz-dz.*vy;
ny = dz.*vx-dx.*vz;
nz = dx.*vy-dy.*vx;
N = sqrt(nx.^2 + ny.^2 +nz.^2);
nx = nx./N;
ny = ny./N;
nz = nz./N;
x1 = x+w*nx;
y1 = y+w*ny;
z1 = z+w*nz;
x2 = x-w*nx;
y2 = y-w*ny;
z2 = z-w*nz;
n = length(x);
vertices = [x1 y1 z1; x2 y2 z2];
faces = zeros(n-1, 4);
for j = 1:n-1
    faces(j,:) = [j, j+1, n+j+1, n+j];
end
c = angleColor(theta4);
cVertex = [c; c];

p = patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', cVertex, 'FaceColor', 'interp', 'EdgeColor', 'none', 'HandleVisibility', 'off');