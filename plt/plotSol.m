dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcs = readmatrix(strcat(dir,"config\uniform_dip.txt"));
phi = readmatrix(strcat(dir,"out\phi.txt"));
fld = readmatrix(strcat(dir,"out\fld.txt"));

phiAnl = readmatrix(strcat(dir,"out\phiAnl.txt"));
fldAnl = readmatrix(strcat(dir,"out\fldAnl.txt"));

nvec = 1:size(phi,1);
pvec = [7];
% Nobs = 1000;
close all;
%%
phiSort = sortrows(phi);
phiAnlSort = sortrows(phiAnl);

figure(1);
% [phiAnl, phi]
plot(nvec, phiAnlSort, nvec, phiSort);
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErr = abs(phiSort-phiAnlSort)./abs(phiAnlSort);
semilogy(nvec, relErr, '-o');

%%
ele = 2; 
fldSort = sortrows(fld,ele);
fldAnlSort = sortrows(fldAnl,ele);

figure(1);
plot(nvec, fldAnlSort(:,ele), nvec, fldSort(:,ele));
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErrFld = abs(fldSort-fldAnlSort)./abs(fldAnlSort);
semilogy(nvec, relErrFld, '-o');

%%
X = srcs(:,1);
Y = srcs(:,2);
npts = 1000;
[Xg, Yg] = meshgrid(linspace(min(X), max(X), npts), linspace(min(Y), max(Y), npts));
phig = griddata(X, Y, phi, Xg, Yg, 'natural');

scale = 2.5;
figure(3);
hold on;
contour(Xg, Yg, phig, 100)
quiver(X, Y, fld(:,1), fld(:,2), scale)
hold off;