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