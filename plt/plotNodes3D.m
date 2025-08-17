dir = "C:\Users\ellil\Documents\WORK\FMM\FMM\out\build\x64-debug\";
srcs = readmatrix(strcat(dir,"config\part3D\uniform_plus.txt"));
% obss = readmatrix(strcat(dir,"config\obss.txt"));
nodes = readmatrix(strcat(dir,"out\nodes.txt"));

nodeLengs = nodes(:,4);
nodeStats = nodes(:,5);
nodeVec = [nodes(:,1:4),nodeLengs,nodeLengs];
nodeVec(:,1:3) = nodeVec(:,1:3) - nodeVec(:,4:6)/2;

clc;
isNode = (nodeStats == 1);
assert( nnz(isNode) == 1 );
nodeLeng = nodeLengs(isNode);

fprintf('# Neighbor nodes: %d\n', nnz(nodeStats == 2));
fprintf('# Interaction nodes: %d\n', nnz(nodeStats >= 3));

%%
% faces = [1 2; 1 3; 2 3];
% close all;
% 
% for i=1:3
%     face = faces(i,:);
%     nodeVec2D = [nodeVec(:,face), nodeVec(:,4:5)];
%     figure(i);
%     for node = 1:length(nodes)
%         nborStat = nodeStats(node);
% 
%         rectangle('Position',nodeVec2D(node,:),...
%                   'LineWidth',2*nodeLengs(node))
%     end
%     hold on;
%     scatter(srcs(:,face(1)),srcs(:,face(2)));
%     hold off;
% end

%%
% figure(1);
% scatter3(srcs(:,1),srcs(:,2),srcs(:,3));
% hold on;
% scatter3(obss(:,1), obss(:,2), obss(:,3));
% hold off;
%%
rootLeng = 1.0;
lim = [-rootLeng/2 rootLeng/2];

figure(4)
for stat=1:8
    nodes = nodeVec(nodeStats == stat,:);
    if (stat==2)
        scatter3(nodes(:,1),nodes(:,2),nodes(:,3),stat2rgb(stat))
    else
        scatter3(nodes(:,1),nodes(:,2),nodes(:,3),stat2rgb(stat),'filled')
    end
    hold on
end
hold off

xlim(lim); ylim(lim); zlim(lim);
xlabel('x'); ylabel('y'); zlabel('z'); 

%%
function rgb = stat2rgb(stat)
    switch stat
        case 1
            rgb = "black";
        case 2
            rgb = "black";
        case 3 
            rgb = "red"; % uplist
        case 4 
            rgb = "green"; % downlist
        case 5 
            rgb = "blue"; % northlist
        case 6 
            rgb = "cyan"; % southlist
        case 7 
            rgb = "magenta"; % eastlist
        case 8 
            rgb = "yellow"; % westlist
    end
end