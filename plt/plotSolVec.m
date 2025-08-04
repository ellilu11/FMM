%%
dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
phiFile = strcat(dir,"out\ff.txt");
phiAnlFile = strcat(dir,"out\ffAnl.txt");

phi = readmatrix(phiFile)';
phiAnl = readmatrix(phiAnlFile)';

% phi = sortrows(phi,1);
% phiAnl = sortrows(phiAnl,1);
nvec = 1:size(phi,1);
numP = size(phi,2);

%%
pmin = 1;
pmax = pmin+numP-1;
pvec = pmin:pmax;

close all;
figure(1);
[phi(:,1), phiAnl(:,1)]
plot(nvec, phiAnl, nvec, phi(:,pvec));
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','northeast');

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(nvec, relErr(:,pvec), '-o');

legend(strcat(" p = ", arrayfun(@num2str,pvec,...
    'UniformOutput',false)) );

% figure(3);
% meanRelErr = mean(relErr,1);
% semilogy(pvec,meanRelErr(pvec), '-o');

