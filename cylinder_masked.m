function [Bout, mask3D] = cylinder_masked(B, r_ero)
%CYLINDER_MASKED Apply a centered cylindrical mask (axis along z) to 3D volume.
%   Bout = cylinder_masked(B)
%   Bout = cylinder_masked(B, r_ero)
%   [Bout, mask3D] = cylinder_masked(...)
%
% Inputs:
%   B     : nx-by-ny-by-nz logical or numeric (0/1 or grayscale; nonzero treated as true)
%   r_ero : erosion radius in pixels (default 0)
%
% Outputs:
%   Bout  : masked logical volume (outside cylinder is false)
%   mask3D: nx-by-ny-by-nz logical cylindrical mask

if nargin < 2 || isempty(r_ero), r_ero = 0; end
r_ero = max(0, r_ero);

fprintf('--- Block 1: setup ---\n'); tic
B = logical(B);
[nx, ny, nz] = size(B);

cx = (ny + 1)/2;      % columns
cy = (nx + 1)/2;      % rows
R  = 0.5*min(nx, ny) - 1;
Re = max(0, R - r_ero);
toc; fprintf('--- Block 1 done: R=%.2f, r_ero=%.2f, Re=%.2f ---\n', R, r_ero, Re);

fprintf('--- Block 2: build 3D cylindrical mask ---\n'); tic
[X, Y]  = meshgrid(1:ny, 1:nx);
mask2D  = (X - cx).^2 + (Y - cy).^2 <= Re^2;
mask3D  = repmat(mask2D, 1, 1, nz);     % boolean 3D mask, same size as B
toc; fprintf('--- Block 2 done ---\n');

fprintf('--- Block 3: apply mask ---\n'); tic
Bout = B & mask3D;
toc; fprintf('--- Block 3 done ---\n');
end
