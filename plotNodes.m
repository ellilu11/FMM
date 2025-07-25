dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcfile = strcat(dir,"config\srcs.txt");
% obsfile = strcat(dir,"config\obss.txt");
nodefile = strcat(dir,"out\nodes.txt");

srcs = readmatrix(srcfile);
% obss = readmatrix(obsfile);
nodes = readmatrix(nodefile);

lengs = nodes(:,3);
% nodes = sortrows(nodes,3,"descend");

%%
nodeVec = [nodes(:,1:3),lengs];
nodeVec(:,1:2) = nodeVec(:,1:2) - nodeVec(:,3:4)/2;

close all;
figure(1);
hold on;

for node = 1:length(nodes)
    rectangle('Position',nodeVec(node,:),...
        'LineWidth',2*lengs(node))
end

scatter(srcs(:,1),srcs(:,2));
% scatter(obss(:,1),obss(:,2));

hold off;

%%
% nodeVec = [nodes(:,1:3),lengs];
% nodeVec(:,1:2) = nodeVec(:,1:2) - nodeVec(:,3:4)/2;
% isNode = nodes(:,4) == 3;
% % assert( nnz(isNode) == 1 );
% nodeLeng = lengs(isNode);
% 
% figure();
% hold on;
% 
% for node = 1:length(nodes)
%     nborStat = nodes(node,4);
%     if (lengs(node) < nodeLeng && ~nborStat)
%         rgb = 'none';
%     else
%         rgb = [1 1-nborStat/nodes(isNode,4) 1];
%     end
%     rectangle('Position',nodeVec(node,:),...
%         'FaceColor',rgb,...
%         'LineWidth',2*lengs(node))
% end
% 
% scatter(srcs(:,1),srcs(:,2));
% 
% hold off;