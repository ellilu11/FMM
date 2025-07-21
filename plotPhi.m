dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\out\";
%%
ffFile = strcat(dir,"ff.txt");
ffAnlFile = strcat(dir,"ffAnl.txt");

ff = readmatrix(ffFile).';
ffAnl = readmatrix(ffAnlFile).';
pmax = size(ff,2);

% Nobs = 1000;
% ff = reshape(phi,Nobs,[],2);

%%
theta = linspace(0,2*pi,length(phi));

figure(1);
pmin = 1;
plot(theta, ffAnl, theta, ff(:,pmin:pmax) );
% ylim([-4620,-4600])
legend([' Analytic',strcat(" p = ", arrayfun(@num2str,pmin:pmax,...
     'UniformOutput',false))] );

figure(2);
semilogy(theta, abs(ff-ffAnl)./abs(ffAnl));

legend(strcat(" p = ", arrayfun(@num2str,pmin:pmax,...
    'UniformOutput',false)) );

%%
nfFile = strcat(dir,"nf.txt");
nfAnlFile = strcat(dir,"nfAnl.txt");

nf = readmatrix(nfFile);
nfAnl = readmatrix(nfAnlFile);
% nfJoin = join(nfAnl,nf,'Keys',[])
pmax = size(nf,2);

figure(1);
plot(nf(:,3));

figure(2);
scatter(nfAnl(:,1),nfAnl(:,2));