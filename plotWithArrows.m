function [plotHandle] = plotWithArrows(xdata, ydata, varargin)
%PLOTWITHArrows Plot data with arrows along the curve showing the direction
%   This function is meant to extend the abilities of the PLOT command by
%   adding arrowheads to the curve.
%
%   H = plotWithArrows(xdata, ydata) - plot the specified data with the default
%   options. The function returns the handle to the plot
%   
%   H = plotWithArrows(xdata, ydata, LineType) - plot the specified data with
%   the specified line type (e.g. 'b*' or '--r'). See PLOT help for more
%   details
%
%   H = plotWithArrows(xdata, ydata, LineType, ...) - Allows you to specify
%   additional options for the plot. See the OPTIONS section below for
%   details.
%
%   OPTIONS:
%
%       'Color'         -   color string | rgb array - The color of the
%                           data and arrowheads. Default is 'b'
%
%       'LineWidth'     -   a number specifying the width of the plotted
%                           line. See PLOT help for more details. Default
%                           is 1.25
%
%       'NumArrows'     -   The number of arrowheads to place on the data.
%                           Default is 10
%
%       'FlipDir'       -   'On' | 'Off' - Whether or not to flip the
%                           direction of the arrowheads (e.g. if the data
%                           is listed in reverse time). Default is 'Off'
%
%       'ArrowScale'    -   Scale factor for the arrowheads. Default is 1
%
%
%   See also PLOT
%
%   Author: Andrew Cox
%   Version: January 30, 2014


%% Set up Parsing
defLineWidth = 2;
defNumArr = 10;
defFlipDir = 'Off';
defArrowScale = 1;
defLineType = '-';
defColor = 'b';

p = inputParser;

validColor = @(x) (isnumeric(x) && length(x) == 3) || ischar(x);
notNegNum = @(x) isnumeric(x) && x >= 0;

addRequired(p, 'xdata', @isnumeric);
addRequired(p, 'ydata', @isnumeric);

addOptional(p, 'LineType', defLineType, @ischar);

addParameter(p, 'Color', defColor, validColor);
addParameter(p, 'LineWidth', defLineWidth, notNegNum);
addParameter(p, 'NumArrows', defNumArr, notNegNum);
addParameter(p, 'FlipDir', defFlipDir, @ischar);
addParameter(p, 'ArrowScale', defArrowScale, notNegNum);

%% Parse
% Check to see if the first optional input could be a LineType
if(length(varargin) > 0 && length(varargin{1}) < 5)
    % It is, proceed normally
    parse(p, xdata, ydata, varargin{:});
else
    % It isn't, throw in the default to avoid errors and parse the optional
    % inputs as param-value pairs
    parse(p, xdata, ydata, defLineType, varargin{:});
end

color = p.Results.Color;
lineType = p.Results.LineType;
lineWidth = p.Results.LineWidth;
numArrows = p.Results.NumArrows;
flipDir = strcmpi(p.Results.FlipDir, 'On');
scale = p.Results.ArrowScale;

% Error checking
if(length(xdata) ~= length(ydata))
    error('xdata and ydata are different lengths!');
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
arrows = zeros(numArrows, 4, 2);
stepSize = round(length(xdata)/numArrows);

% Range of x and y; use to choose size of arrowhead
xExtent = abs(max(xdata) - min(xdata));
yExtent = abs(max(ydata) - min(ydata));
avgExt = mean([xExtent, yExtent]);

% Compute dimensions
l = 0.04*avgExt*scale;      % Length of arrowhead
w = l;                      % Width of arrowhead
s = -0.5*l;                 % Distance from base point to bottom (flat edge) of arrowhead
m = 0.33*l;                 % Indent distance from bottom (flat edge) of arrowhead

for n = 1:numArrows
    ix = (n-1)*stepSize+1;
    
    if(ix > length(xdata))
        break;
    end
    
    loc = [xdata(ix), ydata(ix)];
    
    if(ix < length(xdata))
        dir = [xdata(ix+1), ydata(ix+1)] - loc;
    else
        dir = loc - [xdata(ix-1), ydata(ix-1)];
    end
    
    % Normalize length of dir and flip it if desired
    dir = 0.1*(-1)^flipDir * dir/norm(dir);
    
    % Angle between x-axis and direction vector
    phi = atan2(dir(2), dir(1));
    
    % Four points of arrow head; use patch() to fill these points later
    arrows(n,:,:) = [loc(1) + (s+l)*cos(phi), loc(2) + (s+l)*sin(phi);...
        loc(1) + (s*cos(phi) + w/2*sin(phi)), loc(2) + (s*sin(phi)-w/2*cos(phi));...
        loc(1) + (s+m)*cos(phi), loc(2) + (s+m)*sin(phi);...
        loc(1) + (s*cos(phi) - w/2*sin(phi)), loc(2) + (s*sin(phi)+w/2*cos(phi))];
end

%% Plotting

% are we holding?
wasHolding = ishold;

if(~ishold)
    hold on;
end

figure(gcf); hold on;
plotHandle = plot(xdata, ydata, lineType, 'Color', color, 'LineWidth', lineWidth);
for n = 1:size(arrows,1)
    patch(arrows(n,:,1), arrows(n,:,2), 'r', 'FaceColor', color, 'EdgeColor', color, 'HandleVisibility', 'off');
end

if(wasHolding)
    hold on;
else
    hold off;
end