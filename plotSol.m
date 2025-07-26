dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcfile = strcat(dir,"config\uniform.txt");
srcs = readmatrix(srcfile);

phiFile = strcat(dir,"out\phi.txt");
phi = readmatrix(phiFile);

fldFile = strcat(dir,"out\fld.txt");
fld = readmatrix(fldFile);

phiAnlFile = strcat(dir,"out\phiAnl.txt");
phiAnl = readmatrix(phiAnlFile);

fldAnlFile = strcat(dir,"out\fldAnl.txt");
fldAnl = readmatrix(fldAnlFile);

nvec = 1:size(phi,1);
% Nobs = 1000;
%%
phi = sortrows(phi);
phiAnl = sortrows(phiAnl);

close all;
figure(1);
% [phiAnl, phi]
plot(nvec, phiAnl, nvec, phi);
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(nvec, relErr, '-o');

%%
ele = 2; 
fld = sortrows(fld,ele);
fldAnl = sortrows(fldAnl,ele);

figure(1);
plot(nvec, fldAnl(:,ele), nvec, fld(:,ele));
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErrFld = abs(fld-fldAnl)./abs(fldAnl);
semilogy(nvec, relErrFld, '-o');

%%
X = srcs(:,1);
Y = srcs(:,2);
npts = 1000;
[Xg, Yg] = meshgrid(linspace(min(X), max(X), npts), linspace(min(Y), max(Y), npts));
phig = griddata(X, Y, phi, Xg, Yg, 'natural');

figure(3);
hold on;
contour(Xg, Yg, phig, 100)
quiver(X, Y, fldAnl(:,1), fldAnl(:,2))
hold off;