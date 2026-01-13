function phi = circleMeanPlot(B, r, doPlot, trueColor)
% CIRCLEMEANPLOT  Mean of a boolean slice inside a centered circle + plot overlay.
%
%   phi = circleMeanPlot(B, r)
%   phi = circleMeanPlot(B, r, doPlot)
%   phi = circleMeanPlot(B, r, doPlot, trueColor)
%
% Inputs:
%   B         : 2D logical (or numeric 0/1) matrix
%   r         : circle radius in pixels
%   doPlot    : logical, overlay plot (default = true)
%   trueColor : 1x3 RGB for TRUE pixels (default = [0 0 0] = black)
%
% Output:
%   phi       : mean(B) inside the circle (area fraction of TRUEs)

if nargin < 3 || isempty(doPlot)
    doPlot = true;
end
if nargin < 4 || isempty(trueColor)
    trueColor = [0 0 0]; % black
end

% ensure logical
if ~islogical(B)
    B = B ~= 0;
end

[ny,nx] = size(B);
cx = (nx + 1)/2;
cy = (ny + 1)/2;

% disk mask
[x,y] = meshgrid(1:nx, 1:ny);
mask = (x - cx).^2 + (y - cy).^2 <= r^2;

% mean inside disk
phi = mean(B(mask));

% plot with colored TRUE pixels on white background
if doPlot
    % Build RGB image: background white, TRUE pixels = trueColor
    RGB = ones(ny, nx, 3); % white
    for c = 1:3
        tmp = RGB(:,:,c);
        tmp(B) = trueColor(c);
        RGB(:,:,c) = tmp;
    end

    image(RGB);
    axis image off;
    hold on;

    th = linspace(0, 2*pi, 600);
    xp = cx + r*cos(th);
    yp = cy + r*sin(th);
    plot(xp, yp, 'r', 'LineWidth', 2);

    title(sprintf('\\phi(r)=%.5g  (r=%g px)', phi, r));
    hold off;
end

end
