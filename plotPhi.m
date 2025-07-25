dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\out\";
phiFile = strcat(dir,"phi.txt");
phiAnlFile = strcat(dir,"phiAnl.txt");

phi = readmatrix(phiFile).';
phiAnl = readmatrix(phiAnlFile).';
numP = size(phi,2);

phi = sortrows(phi);
phiAnl = sortrows(phiAnl);

% Nobs = 1000;
% phi = reshape(phi,Nobs,[],2);

%%
nvec = 1:size(phi,1);
% nvec = linspace(0,2*pi,size(phi,1));
pmax = 4;
pvec = pmax-numP+1:pmax;

close all;
figure(1);
[phiAnl, phi]
plot(nvec, phiAnl, nvec, phi);
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(nvec, relErr, '-o');

legend(strcat(" p = ", arrayfun(@num2str,pvec,...
    'UniformOutput',false)) );

% figure(3);
% meanRelErr = mean(relErr,1);
% semilogy(pvec,meanRelErr, '-o');
