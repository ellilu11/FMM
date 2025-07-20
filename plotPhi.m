dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\out\";
ffFile = strcat(dir,"ff.txt");
ffAnlFile = strcat(dir,"ffAnl.txt");

phi = readmatrix(ffFile).';
phiAnl = readmatrix(ffAnlFile).';

% Nobs = 1000;
% phi = reshape(phi,Nobs,[],2);

%%
theta = linspace(0,2*pi,length(phi));

figure(1);
plot(theta, phi(:,1), theta, phiAnl);
ylim([-4620,-4600])

figure(2);
semilogy(theta, abs(phi-phiAnl)./abs(phiAnl));

legend(strcat(" p = ", arrayfun(@num2str,1:size(phi,2),...
    'UniformOutput',false)) );

% figure(2);
% plot(theta, phi(:,2), theta, phi(:,4));
