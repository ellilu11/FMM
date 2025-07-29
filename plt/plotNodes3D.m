dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcs = readmatrix(strcat(dir,"config\part3D\uniform_plus.txt"));
% obss = readmatrix(strcat(dir,"config\obss.txt"));
nodes = readmatrix(strcat(dir,"out\nodes.txt"));

lengs = nodes(:,4);
% nodes = sortrows(nodes,3,"descend");

%%
faces = [1 2; 1 3; 2 3];
nodeVec = [nodes,lengs,lengs];
nodeVec(:,1:3) = nodeVec(:,1:3) - nodeVec(:,4:6)/2;
close all;

for i=1:3
    face = faces(i,:);
    nodeVec2D = [nodeVec(:,face), nodeVec(:,4:5)];
    figure(i);
    hold on;
    for node = 1:length(nodes)
        rectangle('Position',nodeVec2D(node,:),...
            'LineWidth',2*lengs(node))
    end
    scatter(srcs(:,face(1)),srcs(:,face(2)));
    hold off;
end

figure(4)
scatter3(srcs(:,1),srcs(:,2),srcs(:,3));
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