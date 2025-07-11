dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\"
posfile = strcat(dir,"positions.txt");
nodefile = strcat(dir,"nodes.txt");

positions = readmatrix(posfile);
nodes = readmatrix(nodefile);

%%
figure();
scatter(positions(:,1),positions(:,2));
for node = 1:length(nodes)
    rectangle()
end