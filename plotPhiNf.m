%%
phiFile = strcat(dir,"nf.txt");

phiAnlFile = strcat(dir,"nfAnl.txt");

phi = readmatrix(phiFile)';
phiAnl = readmatrix(phiAnlFile);

phi = sortrows(phi,1);
phiAnl = sortrows(phiAnl,1);
nvec = 1:size(phi,1);
pmax = size(phi,2);

%%
pmin = 2;
pvec = 1:5; % pmin:pmax;

figure(1);
[phi(:,1), phiAnl(:,1)]
plot(nvec, phiAnl, nvec, phi(:,pvec));
% plot(nvec, phiAnl(:,2), nvec, phi(:,2));
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(nvec, relErr(:,pvec), '-o');

legend(strcat(" p = ", arrayfun(@num2str,pvec,...
    'UniformOutput',false)) );

figure(3);
meanRelErr = mean(relErr,1);
semilogy(pvec,meanRelErr(pvec), '-o');

