%% example code for edge PC calculation
% load efc
load('/data/users/tfang/ePC_master/data/efc.mat');
% load edge community assignment
load('/data/users/tfang/ePC_master/data/group_idx.mat');

efc(efc<0)=0;  % using only positive connectivity

%Calculate the edge PC
ePC_pos = participation_coef(efc,idx,0); % BCT function

%convert to matrix form
mat = zeros(90);
mat(triu(ones(90),1) > 0) = ePC_pos; 
mat_ePC = mat + mat';
[gx,gy,order] = grid_communities(Net7s); 
load('/data/users/tfang/ePC_master/data/colorbar_blue.mat');
figure,
imagesc(mat_ePC(order,order),[0.83 0.87]); axis square; colormap(colorbar_blue)
hold on; plot(gx,gy,'k','linewidth',1.5); title('Edge participant coefficient',FontSize=22)
set(gcf,'color','white')
