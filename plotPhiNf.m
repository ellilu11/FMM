%%
phiFile = strcat(dir,"nf.txt");
phiAnlFile = strcat(dir,"nfAnl.txt");

phi = readmatrix(phiFile);
phiAnl = readmatrix(phiAnlFile);

phi = sortrows(phi,1);
phiAnl = sortrows(phiAnl,1);
nvec = 1:length(phi);

%%
pvec = [10];

figure(1);
[phi(:,1), phiAnl(:,1)]
plot(nvec, phiAnl(:,1), nvec, phi(:,1));
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))],...
     'Location','southeast');

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(nvec, relErr, '-o');

legend(strcat(" p = ", arrayfun(@num2str,pvec,...
    'UniformOutput',false)) );
%%
% [phiAnl(:,1), phi(:,3)]
% plot(nvec, phiAnl(:,1), nvec, phi(:,3))
% legend(' Analytic ', ' FMM (far) + Direct (near) ', ...
%        'Location','southeast');
