% stationarity test (axial mean signal)

rockname  = "Thalassinoides180";

A = tiffreadVolume(rockname + ".tif");

% --- Robust global binarization (avoid slice-dependent behavior) ---
Ad = mat2gray(A);
T  = graythresh(Ad(:));
B  = Ad > T;                 % if pores are dark, use: B = Ad < T;

% --- Apply cylindrical mask (outside cylinder -> false) ---
B = cylinder_masked(B);

% --- Physical scan size [mm] ---
Dx = 190;
Dy = 180;
Dz = 310;

% --- Voxel size [mm/voxel] (corrected) ---
vx_mm = Dx / size(B,2);      % X = columns
vy_mm = Dy / size(B,1);      % Y = rows
vz_mm = Dz / size(B,3);      % Z = slices

% --- Slice-wise mean (phase fraction per slice) ---
ns  = size(B,3);
x   = (1:ns)';               % slice index (or convert to mm using vz_mm)
phi = zeros(ns,1);

for k = 1:ns
    slice = B(:,:,k);
    phi(k)= mean(slice(:));
end

% --- Detrending analysis (do NOT change per your instruction) ---
tws = 3:100;
nw  = length(tws);
kexs = zeros(nw,1);

for i = 1:nw
    tw = tws(i);
    kexs(i) = detrend_window_xy(x, phi, tw);
    pause(0.1)
end

% --- Plot kurtosis vs window ---
clf
plot(tws, kexs, 'o-','MarkerSize',10,'LineWidth',2)
xlabel('size of the MA window (slices)')
ylabel("excess kurtosis ($k-3$)",'Interpreter','latex')
grid on
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf,'Renderer','painters');
print(gcf,'kurtosis.pdf','-dpdf','-painters');

% --- Optimal rolling window (fixed bug: do not overwrite kexs(end)) ---
pause(5)
clf
tw_opt = 43;
kex_opt = detrend_window_xy(x, phi, tw_opt);
set(gcf,'Renderer','painters');
print(gcf,'MA.pdf','-dpdf','-painters');
