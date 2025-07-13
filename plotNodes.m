dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
posfile = strcat(dir,"positions.txt");
nodefile = strcat(dir,"nodes.txt");

positions = readmatrix(posfile);
nodes = readmatrix(nodefile);
% nodes = sortrows(nodes,3,"descend");
%%
figure(1);
hold on;
nodepos = [nodes(:,1:3),nodes(:,3)];
nodepos(:,1:2) = nodepos(:,1:2) - nodepos(:,3:4)/2;

for node = 1:length(nodes)
    nborFlag = nodes(node,4);
    if (nodes(node,3) == 0.3125 && ~nborFlag)
        rgb = 'none';
    else
        rgb = [1 1-nborFlag/2 1];
    end
    rectangle('Position',nodepos(node,:),...
        'FaceColor',rgb,...
        'LineWidth',nodes(node,3)*2)
end
scatter(positions(:,1),positions(:,2));