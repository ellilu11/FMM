dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcs = readmatrix(strcat(dir,"config\uniform_dip.txt"));
% obss = readmatrix(strcat(dir,"config\obss.txt"));
nodes = readmatrix(strcat(dir,"out\nodes.txt"));

lengs = nodes(:,3);
nodeStats = nodes(:,4);
nodeVec = [nodes(:,1:3),lengs];
nodeVec(:,1:2) = nodeVec(:,1:2) - nodeVec(:,3:4)/2;
% nodeLeng = lengs(isNode);

fprintf('# List 1 nodes: %d\n', nnz(nodeStats == 2));
fprintf('# List 2 nodes: %d\n', nnz(nodeStats == 3));
fprintf('# List 3 leaf nodes: %d\n', nnz(nodeStats == 4));
fprintf('# List 3 stem nodes: %d\n', nnz(nodeStats == 5));
fprintf('# List 4 nodes: %d\n', nnz(nodeStats == 6));
fprintf('# Overlapping nodes: %d\n', nnz(nodeStats > 6));

%%
% nodeVec = [nodes(:,1:3),lengs];
% nodeVec(:,1:2) = nodeVec(:,1:2) - nodeVec(:,3:4)/2;
% 
% close all;
% figure(1);
% hold on;
% 
% for node = 1:length(nodes)
%     rectangle('Position',nodeVec(node,:),...
%         'LineWidth',2*lengs(node))
% end
% 
% scatter(srcs(:,1),srcs(:,2));
% % scatter(obss(:,1),obss(:,2));
% 
% hold off;

%%
close all;
figure(1);

for node = 1:length(nodes)
    nborStat = nodes(node,4);
    if (lengs(node) < nodeLeng && ~nborStat)
        rgb = 'none';
    else
        rgb = stat2rgb(nborStat);
    end
    rectangle('Position',nodeVec(node,:),...
        'FaceColor',rgb,...
        'LineWidth',2*lengs(node))
end

%%
function rgb = stat2rgb(stat)
    switch stat
        case 0
            rgb = "none";
        case 1
            rgb = "black"; % self
        case 2
            rgb = "green"; % list 1
        case 3 
            rgb = "cyan";   % list 2
        case 4 
            rgb = "magenta"; % list 3 (leaf)
        case 5 
            rgb = "red";  % list 3 (stem)
        case 6 
            rgb = "blue"; % list 4
        case 7 
            rgb = "magenta"; 
        case 8 
            rgb = "yellow"; 
    end
end