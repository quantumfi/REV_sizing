function [k, Chat] = radialSpectrum2D(r, C, Nk)
% RADIALSPECTRUM2D  Approximate 2D isotropic FT of a radial covariance C(r).
%   [k, Chat] = radialSpectrum2D(r, C)
%   [k, Chat] = radialSpectrum2D(r, C, Nk)
%
% Inputs:
%   r  : vector of radii (monotone increasing, in pixels or physical units)
%   C  : covariance values C(r)
%   Nk : number of k points (default 200)
%
% Outputs:
%   k    : radial wavenumbers (inverse units of r)
%   Chat : 2D radial spectrum:
%          Chat(k) ≈ 2*pi * ∫_0^∞ C(r) * r * J0(k r) dr
%
% Notes:
%   - Removes large-r plateau so C0(r) -> 0 as r -> rmax.
%   - Uses trapz integration (no periodic continuation).
%   - This is the 2D isotropic (Hankel) transform of order 0.

if nargin < 3 || isempty(Nk)
    Nk = 200;
end

fprintf('[radialSpectrum2D] Preparing data... ');
t0 = tic;

r = r(:);
C = C(:);

% enforce sorting by r
[rs, idx] = sort(r);
Cs = C(idx);

n = numel(rs);
if n < 5
    error('Need at least ~5 samples in r to compute spectrum.');
end

% estimate plateau (average over last 20% of points)
tailIdx = max(1, round(0.8*n)) : n;
C_inf   = mean(Cs(tailIdx));
C0      = Cs - C_inf;

rmin = rs(1);
rmax = rs(end);

fprintf('done (%.3fs). r in [%.3g, %.3g], C_inf = %.3g\n', toc(t0), rmin, rmax, C_inf);

% choose k-grid
dr   = min(diff(rs));
kmax = pi / max(dr, eps);          % conservative; you can lower if needed
k    = linspace(0, kmax, Nk)';     % k=0 is fine in 2D (J0(0)=1)

Chat = zeros(Nk,1);

fprintf('[radialSpectrum2D] Integrating for %d k values... ', Nk);
t1 = tic;

for i = 1:Nk
    ki = k(i);
    kr = ki * rs;

    % 2D isotropic kernel is J0(kr)
    kernel = besselj(0, kr);

    integrand = C0 .* rs .* kernel;

    Chat(i) = 2*pi * trapz(rs, integrand);
end

fprintf('done (%.3fs)\n', toc(t1));

end
