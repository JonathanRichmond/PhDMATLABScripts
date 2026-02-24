%%% plot3WithArrows
%%% Jonathan LeFevre Richmond
%%% C: 18 February 2026

function [plotHandle] = plot3WithArrows(xdata, ydata, zdata, varargin)
%% Set up Parsing
defLineWidth = 2;
defNumArr = 5;
defFlipDir = 'Off';
defArrowScale = 1;
defLineType = '-';
defColor = 'b';

p = inputParser;

validColor = @(x) (isnumeric(x) && length(x) == 3) || ischar(x);
notNegNum = @(x) isnumeric(x) && x >= 0;

addRequired(p, 'xdata', @isnumeric);
addRequired(p, 'ydata', @isnumeric);
addRequired(p, 'zdata', @isnumeric);

addOptional(p, 'LineType', defLineType, @ischar);

addParameter(p, 'Color', defColor, validColor);
addParameter(p, 'LineWidth', defLineWidth, notNegNum);
addParameter(p, 'NumArrows', defNumArr, notNegNum);
addParameter(p, 'FlipDir', defFlipDir, @ischar);
addParameter(p, 'ArrowScale', defArrowScale, notNegNum);

%% Parse
% Check to see if the first optional input could be a LineType
if(~isempty(varargin) && length(varargin{1}) < 5)
    % It is, proceed normally
    parse(p, xdata, ydata, zdata, varargin{:});
else
    % It isn't, throw in the default to avoid errors and parse the optional
    % inputs as param-value pairs
    parse(p, xdata, ydata, zdata, defLineType, varargin{:});
end

color = p.Results.Color;
lineType = p.Results.LineType;
lineWidth = p.Results.LineWidth;
numArrows = p.Results.NumArrows;
flipDir = strcmpi(p.Results.FlipDir, 'On');
scale = p.Results.ArrowScale;

% Error checking
if(length(xdata) ~= length(ydata) || length(xdata) ~= length(zdata) || length(ydata) ~= length(zdata)) 
    error('Data vectors are different lengths!');
end

if(numArrows > length(xdata))
    error('Number of arrows exceeds number of data points! Cannot create arrows...');
end

% Make the color match a color specified in LineType
% letters = ischarprop(lineType, 'alpha'); %logical array: whether or not each char is a letter
letters = isletter(lineType);
if(sum(letters) > 0)
    color = lineType(letters);
end

%% Compute arrow directions and locations
numBasePts = 8;
arrows = struct('tip', cell(1, numArrows), 'basePts', cell(1, numArrows));
stepSize = round(length(xdata)/(numArrows+1));

% Range of x and y; use to choose size of arrowhead
xExtent = abs(max(xdata) - min(xdata));
yExtent = abs(max(ydata) - min(ydata));
zExtent = abs(max(zdata) - min(zdata));
avgExt = mean([xExtent, yExtent, zExtent]);

% Compute dimensions
l = 0.04*avgExt*scale;  % Length of arrowhead
w = l;                              % Width of arrowhead
s = -0.5*l;                         % Distance from base point to bottom (flat edge) of arrowhead

for n = 1:numArrows
    ix = (n)*stepSize+1;
    
    if(ix > length(xdata))
        break;
    end
    
    loc = [xdata(ix), ydata(ix), zdata(ix)];
    
    if(ix < length(xdata))
        dir = [xdata(ix+1), ydata(ix+1), zdata(ix+1)] - loc;
    else
        dir = loc - [xdata(ix-1), ydata(ix-1), zdata(ix-1)];
    end
    
    % Normalize length of dir and flip it if desired
    dir = 0.1*(-1)^flipDir * dir/norm(dir);

    % Build local coordinate frame
    u = dir/norm(dir);

    % Create a vector not parallel to u
    if abs(dot(u, [0 0 1])) < 0.9
        temp = [0 0 1];
    else
        temp = [0 1 0];
    end
    v = cross(u, temp);
    v = v/norm(v);
    wvec = cross(u, v);
    
    tip = loc+l*u;
    base = loc+s*u;

    theta = linspace(0, 2*pi, numBasePts+1);
    theta(end) = [];
    basePts = zeros(numBasePts, 3);
    for k = 1:numBasePts
        basePts(k,:) = base+(w/2)*cos(theta(k))*v+(w/2)*sin(theta(k))*wvec;
    end

    arrows(n).tip = tip;
    arrows(n).basePts = basePts;
end

%% Plotting

% are we holding?
wasHolding = ishold;

if(~ishold)
    hold on;
end

figure(gcf); hold on;
plotHandle = plot3(xdata, ydata, zdata, lineType, 'Color', color, 'LineWidth', lineWidth);
for n = 1:numArrows
    tip = arrows(n).tip;
    basePts = arrows(n).basePts;

    N = size(basePts, 1);
    verts = [tip; basePts];
    faces = zeros(N,3);
    for k = 1:N
        next = mod(k, 8)+1;
        faces(k,:) = [1, k+1, next+1];
    end
    patch('Faces', faces, 'Vertices', verts, 'FaceColor', color, 'EdgeColor', color, 'HandleVisibility', 'off');
    patch('Faces', 2:N+1, 'Vertices', verts, 'FaceColor', color, 'EdgeColor', color, 'HandleVisibility', 'off');
end

if(wasHolding)
    hold on;
else
    hold off;
end