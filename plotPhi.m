dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\out\";
%%
% phiFile = strcat(dir,"ff.txt");
% phiAnlFile = strcat(dir,"ffAnl.txt");
phiFile = strcat(dir,"local.txt");
phiAnlFile = strcat(dir,"localAnl.txt");

phi = readmatrix(phiFile).';
phiAnl = readmatrix(phiAnlFile).';
pmax = size(phi,2);

% Nobs = 1000;
% phi = reshape(phi,Nobs,[],2);

%%
theta = linspace(0,2*pi,length(phi));
pmin = 1;
pvec = pmin:pmax;

figure(1);
plot(theta, phiAnl, theta, phi(:,pvec), '-o');
        % 'LineWidth',[2,ones(1,pmax-pmin+1)]
% ylim([-4620,-4600])
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pvec,...
     'UniformOutput',false))] );

figure(2);
relErr = abs(phi-phiAnl)./abs(phiAnl);
semilogy(theta, relErr, '-o');

legend(strcat(" p = ", arrayfun(@num2str,pvec,...
    'UniformOutput',false)) );

figure(3);
semilogy(pvec,mean(relErr,1), '-o');
