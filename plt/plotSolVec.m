%%
dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";

phiFile = strcat(dir,"out\ff.txt");
phiAnlFile = strcat(dir,"out\ffAnl.txt");
% phiFile = strcat(dir,"out\nf.txt");
% phiAnlFile = strcat(dir,"out\nfAnl.txt");
% phiFile = strcat(dir,"out\local.txt");
% phiAnlFile = strcat(dir,"out\localAnl.txt");

phi = readmatrix(phiFile)';
phiAnl = readmatrix(phiAnlFile)';

% phi = sortrows(phi);
% phiAnl = sortrows(phiAnl);

nvec = 1:size(phi,1);
% nvec = linspace(0,2*pi,size(phi,1));
% nobs = 100;
% nvec = 0:(2*pi/nobs):(2*pi*(nobs-1)/nobs);
numP = size(phi,2);

%%
pmax = 10;
pvec = pmax-numP+1:pmax;

close all;
figure(1);
% hold on;
[phiAnl,phi(:,numP)]
plot(nvec, phiAnl, nvec, phi(:,1:numP));
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');
% hold off;

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(nvec, relErr(:,1:numP), '-o');

legend(strcat(" p = ", arrayfun(@num2str,pvec,...
    'UniformOutput',false)) );

figure(3);
meanRelErr = mean(relErr,1);
semilogy(pvec,meanRelErr(1:numP), '-o');