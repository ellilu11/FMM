dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcs = readmatrix(strcat(dir,"config\part3D\uniform_dip.txt"));
phi = readmatrix(strcat(dir,"out\phi.txt"));
fld = readmatrix(strcat(dir,"out\fld.txt"));

phiAnl = readmatrix(strcat(dir,"out\phiAnl.txt"));
fldAnl = readmatrix(strcat(dir,"out\fldAnl.txt"));

nvec = 1:size(phi,1);
close all;
%%
phiSort = sortrows(phi);
phiAnlSort = sortrows(phiAnl);

figure(1);
% [phiAnl, phi]
hold on;
plot(nvec, phiAnlSort, nvec, phiSort);

% hold on;
figure(2);
relErr = abs(phiSort-phiAnlSort)./abs(phiAnlSort);
semilogy(nvec, relErr, '-o');

%%
ele = 1; 
fldSort = sortrows(fld,ele);
fldAnlSort = sortrows(fldAnl,ele);

figure(3);
plot(nvec, fldAnlSort(:,ele), nvec, fldSort(:,ele));

figure(4);
relErrFld = abs(fldSort(:,ele)-fldAnlSort(:,ele))./abs(fldAnlSort(:,ele));
semilogy(nvec, relErrFld, '-o');

%%
X = srcs(:,1);
Y = srcs(:,2);
Z = srcs(:,3);
% npts = 100;
% [Xg, Yg, Zg] = meshgrid(linspace(min(X), max(X), npts), linspace(min(Y), max(Y), npts), linspace(min(Z), max(Z), npts));
% phig = griddata(X, Y, Z, phi, Xg, Yg, Zg, 'natural');
% 

figure(5);
% hold on;
% contour3(Xg, Yg, Zg, phig, 100)
scale = 2.5;
quiver3InLogScale(X, Y, Z, fld(:,1), fld(:,2), fld(:,3), scale)
% quiver3(X, Y, Z, fldAnl(:,1), fldAnl(:,2), fldAnl(:,3), scale)
% hold off;