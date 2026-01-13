function [rvals, C_r] = slice2cov_bool(B, doCircleCrop, fracR, lagStep)
% SLICE2COV_BOOL  Direct 2-D covariance C(r) from a single boolean slice.
%
%   [rvals, C_r] = slice2cov_bool(B)
%   [rvals, C_r] = slice2cov_bool(B, doCircleCrop, fracR, lagStep)
%
% Inputs:
%   B            : 2-D logical (or numeric 0/1) matrix (one slice)
%   doCircleCrop : logical, if true keep central circle (default: true)
%   fracR        : radius fraction of inscribed circle (default: 1.0)
%                  fracR=1 uses the inscribed circle of the rectangle.
%                  fracR<1 uses a smaller concentric circle.
%   lagStep      : lag step in pixels (default: 1)
%
% Outputs:
%   rvals : lags r (pixels)
%   C_r   : covariance averaged over x and y directions (2-D)
%
% Notes:
%   - The function crops by masking outside the circle and then restricts
%     computations to the bounding box of that circle.
%   - Mean is computed only over the valid circle region.
%   - Covariance is computed using products of shifted arrays, but only
%     counting pairs where BOTH pixels are inside the circle.

if nargin < 2 || isempty(doCircleCrop), doCircleCrop = true; end
if nargin < 3 || isempty(fracR),       fracR = 1.0;          end
if nargin < 4 || isempty(lagStep),     lagStep = 1;          end

% --- validate / cast ---
if ~islogical(B)
    B = B ~= 0;
end
B = single(B);

[Ny, Nx] = size(B);

% --- build circular mask (centered) ---
cx = (Nx + 1)/2;
cy = (Ny + 1)/2;
R0 = min(Nx, Ny)/2;          % inscribed circle radius
R  = fracR * R0;

[x, y] = meshgrid(1:Nx, 1:Ny);
mask = ((x - cx).^2 + (y - cy).^2) <= R^2;

if doCircleCrop
    % crop to bounding box of circle to reduce work
    cols = find(any(mask,1));
    rows = find(any(mask,2));
    x1 = cols(1); x2 = cols(end);
    y1 = rows(1); y2 = rows(end);

    B    = B(y1:y2, x1:x2);
    mask = mask(y1:y2, x1:x2);
    [Ny, Nx] = size(B);
else
    % if no crop, still use full mask if desired:
    % (here: if doCircleCrop=false, use all pixels)
    mask(:) = true;
end

% --- mean subtraction on valid region only ---
m = mean(B(mask));
V = B - m;

% --- maximum lag limited by circle bbox ---
maxLag = floor(min(Nx, Ny)/2);
lagList = 0:lagStep:maxLag;
nLags = numel(lagList);

Cx = zeros(nLags,1,'double');
Cy = zeros(nLags,1,'double');

% lag 0 variance (valid region)
var0 = mean(V(mask).^2);

for iLag = 1:nLags
    h = lagList(iLag);

    if h == 0
        Cx(iLag) = var0;
        Cy(iLag) = var0;
    else
        % ---- X direction ----
        A1 = V(:,1:end-h);
        A2 = V(:,1+h:end);
        M1 = mask(:,1:end-h);
        M2 = mask(:,1+h:end);
        Mp = M1 & M2;
        if any(Mp(:))
            Cx(iLag) = mean( (A1(Mp) .* A2(Mp)) );
        else
            Cx(iLag) = NaN;
        end

        % ---- Y direction ----
        A1 = V(1:end-h,:);
        A2 = V(1+h:end,:);
        M1 = mask(1:end-h,:);
        M2 = mask(1+h:end,:);
        Mp = M1 & M2;
        if any(Mp(:))
            Cy(iLag) = mean( (A1(Mp) .* A2(Mp)) );
        else
            Cy(iLag) = NaN;
        end
    end
end

C_r = (Cx + Cy) / 2;
rvals = lagList(:);

end
