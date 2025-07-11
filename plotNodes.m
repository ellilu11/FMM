dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\"
posfile = strcat(dir,"positions.txt");
nodefile = strcat(dir,"nodes.txt");

positions = readmatrix(posfile);
nodes = readmatrix(nodefile);

%%
figure(1);
scatter(positions(:,1),positions(:,2));
nodepos = [nodes,nodes(:,3)];
nodepos(:,1:2) = nodepos(:,1:2) -nodepos(:,3:4)/2

for node = 1:length(nodes)
    rectangle('Position',nodepos(node,:))
end