%%
phiFile = strcat(dir,"nf.txt");
phiAnlFile = strcat(dir,"nfAnl.txt");

phi = readmatrix(phiFile);
phiAnl = readmatrix(phiAnlFile);

figure(1);
plot(nf(:,3));

figure(2);
scatter(nfAnl(:,1),nfAnl(:,2));