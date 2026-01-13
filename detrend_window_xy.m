function Kex = detrend_window_xy(x,phi,window_size)
% detrend_window_xphi
% Detrend phi(x) using a rolling-mean trend in a window of size window_size,
% plot: (i) signal + trend, (ii) detrended noise.
%
% Kurtosis details (computed on detrended noise z):
%   For a finite sample z_i, i=1..M:
%     mu  = (1/M) * sum_i z_i
%     m2  = (1/M) * sum_i (z_i - mu)^2
%     m4  = (1/M) * sum_i (z_i - mu)^4
%     kurtosis        = m4 / (m2^2)
%     excess kurtosis = kurtosis - 3

x = x(:);
phi = phi(:);
N = numel(phi);

if numel(x) ~= N
    error('x and phi must have same length');
end
if window_size < 3 || window_size > N
    error('window_size must be between 3 and length(phi)');
end

trend = movmean(phi, window_size, 'Endpoints','shrink');
phi_detrended = phi - trend;

% ---- Kurtosis calculation on detrended noise ----
z  = phi_detrended(isfinite(phi_detrended));
M  = numel(z);

mu = mean(z);
m2 = mean( (z - mu).^2 );
m4 = mean( (z - mu).^4 );

K   = m4 / (m2^2);   % Pearson kurtosis
Kex = K - 3;         % excess kurtosis

fprintf('Detrended noise kurtosis details:\n');
fprintf('M  = %d\n', M);
fprintf('mu = %.6g\n', mu);
fprintf('m2 = %.6g\n', m2);
fprintf('m4 = %.6g\n', m4);
fprintf('kurtosis        = %.6g\n', K);
fprintf('excess kurtosis = %.6g\n', Kex);

% ---- Plots ----
clf
hold on

subplot(2,1,1)
plot(x, phi, 'k', 'LineWidth', 1); hold on
plot(x, trend, 'r', 'LineWidth', 2);
grid on
xlabel('k-slice'); ylabel('\phi');
legend('Field $\phi$','Trend ($\bar\phi)$','Location','southwest','interpreter','latex');
title(sprintf('Field and trend (window = %d slices )', window_size));
set(gca,'FontName','Times new roman','fontsize',20)

subplot(2,1,2)
plot(x, phi_detrended, 'b', 'LineWidth', 1);
grid on
xlabel('k-slice'); ylabel('$\phi - \bar{\phi}$','Interpreter','latex');
title(sprintf('Detrended noise (excess kurtosis = %.4g)', Kex));
set(gca,'FontName','Times new roman','fontsize',20)
axis tight
end
