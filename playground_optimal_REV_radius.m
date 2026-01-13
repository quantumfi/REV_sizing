% calculate the correlation of Thalassinoides rocks (single 2x2 tiled figure)
clear


rockname = 'Thalassinoides180';
filename = rockname + ".tif";

% --- Read + binarize + mask ---
V = tiffreadVolume(filename);
B = imbinarize(mat2gray(V));
B = cylinder_masked(B);

% --- Pick slice ---
slice   = 10;
B_slice = B(:,:,slice);
if ~islogical(B_slice), B_slice = B_slice ~= 0; end

[ny,nx] = size(B_slice);
Rmax = floor(min(nx,ny)/2);

% --- 1) <phi>(r) vs radius ---
rList = (5:2:Rmax)';             % radii (pixels)
phi_r = zeros(size(rList));
for i = 1:numel(rList)
    phi_r(i) = circleMeanPlot(B_slice, rList(i), false); % no plots here
end

% choose one radius for overlay in panel (a)
r_overlay = min(100, Rmax);


% --- 2) covariance C(r) vs r ---
doCircleCrop = true;
fracR        = 1.0;
lagStep      = 1;
[rvals, C_r] = slice2cov_bool(B_slice, doCircleCrop, fracR, lagStep);

% --- 3) radial spectrum Chat(k) vs k_r ---
Nk = 200;
[k, Chat] = radialSpectrum2D(rvals, C_r, Nk);

% --- Single 2x2 tiled figure ---
fig = figure;

set(fig,'Renderer','painters');

t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) Slice image with circumference
nexttile
% Plot B_slice with TRUE pixels black and background white
RGB = ones(ny, nx, 3); % white background
for c = 1:3
    tmp = RGB(:,:,c);
    tmp(B_slice) = 0;  % TRUE -> black
    RGB(:,:,c) = tmp;
end
image(RGB); axis image off; hold on; 

cx = (nx + 1)/2; cy = (ny + 1)/2;
th = linspace(0,2*pi,800);
plot(cx + r_overlay*cos(th), cy + r_overlay*sin(th), 'r', 'LineWidth', 2);
plot(cx + Rmax*cos(th), cy + Rmax*sin(th), 'k', 'LineWidth', 2);
title(sprintf('Slice %d with r=%d px', slice, r_overlay));
box on

% panel label
text(0.05,0.98,'(a)','Units','normalized','FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top','Color','k');

% (b) <phi>(r)
nexttile
plot(rList, phi_r, 'o-','LineWidth',2,'MarkerSize',5);
grid on;
xlabel('radius $r$ (pixels)','Interpreter','latex');
ylabel('$\langle \phi \rangle(r)=\left\langle B\right\rangle_{|\mathbf{x}|\le r}$', ...
      'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',16);
text(0.1,0.98,'(b)','Units','normalized','FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top');
axis tight
% (c) C(r)
nexttile
plot(rvals, C_r, 'LineWidth',2);
grid on;
xlabel('lag $r$ (pixels)','Interpreter','latex');
ylabel('$C(r)$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',16);
text(0.05,0.98,'(c)','Units','normalized','FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top');
axis tight
% (d) Chat(k)
nexttile
loglog(k(2:end), abs(Chat(2:end)), 'LineWidth',2); % skip k=0 for log
grid on;
xlabel('$k_r$ (1/pixels)','Interpreter','latex');
ylabel('$\hat{C}(k_r)$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',16);
text(0.5,0.98,'(d)','Units','normalized','FontSize',16,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top');

% overall title
%title(t, sprintf('Thalassinoides: slice-wise statistics (slice %d)', slice), ...
%    'FontName','Times New Roman','FontSize',18);

% export
print(fig, 'REV_r.pdf', '-dpdf', '-painters');
