%% Fig.1 plot
% load group-representative edge community labels
load('/data/users/tfang/ePC_master/Data_source/group_idx.mat')
load('/data/users/tfang/ePC_master/Data_source/eFC_colormap.mat')
load('/data/users/tfang/ePC_master/Data_source/efc.mat')
%% Plot eFC orderred by community label
ci = idx;
[gx1,gy1,order] = grid_communities(ci); % BCT function %oreder the martix by different communities or networks
[~,idx_] = sort(ci); dffidx = find(diff(ci(idx_)));
figure, imagesc(efc(order,order),[-.5,.5]); colormap(eFC_colormap);
axis square;
hold on; 
plot(gx1,gy1,'k','linewidth',1); 
title('Edge Functional Connectivity','FontSize',24)

%% Plot edge community and the similarity of edge community assignments associated with the two connected regions.
N = 90;
mat_FT = zeros(N);
mat_FT(triu(ones(N),1) > 0) = idx; % this step gives you the indices of each edge
mat_FT = mat_FT + mat_FT';
% calculate community similarity
s_FT = fcn_profilesim(mat_FT);
load('/data/users/tfang/ePC_master/Data_source/Net7s_alter.mat')
load('/data/users/tfang/ePC_master/Data_source/Gray2red.mat')
[gx,gy,idx_] = grid_communities(Net7s_alter); % BCT function %oreder the martix by different communities or networks
figure, 
imagesc(mat_FT(idx_,idx_),[0 20]);colormap("colorcube");
axis square; hold on; plot(gx,gy,'k','linewidth',1); 
title('Edge communities','FontSize',24)
set(gcf,'Color','white')
box off
set(gca,'Visible','off');
figure,
imagesc(s_FT(idx_,idx_),[0 1]);axis square; colormap((slanCM(9)));
 hold on; plot(gx,gy,'k','linewidth',1); 
title('Edge community similarity','FontSize',24)
set(gcf,'Color','white')
box off
set(gca,'box','on');
%% polt within-system similarity
S=s_FT(idx_,idx_);
Net7s_sort=sort(Net7s_alter);
S_within=[];
for i=1:8
    row=find(Net7s_sort==i);
    mask=triu(true(length(row),length(row)),1);
    O=S(row,row);
    varName = ['Group', num2str(i)];  % Create variable name dynamically
    eval([varName, ' = O(mask);']);

    S_within=O(mask);
    L(i)=length(S_within);
    data(i)=mean(S_within);
    err(i)=std(S_within);
end
[data_rank,net_order]=sort(data,'descend');
err_rank=err(net_order);
x=1:8;
figure,
% bar(x,(data+data1)/.2,0.7);

bar(x,data_rank,0.7);
hold on
er=errorbar(x,data_rank,err_rank);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold off
set(gca,'YTick',0:0.2:1);
set(gca,'Box','off');
set(gca,'XTickLabelRotation',90);
set(gca,'YTickLabelRotation',90);
lab={'Vis','SMN','VAN','Sub','DAN','Lim','FP','DMN'};
lab_rank=lab(net_order);
xticklabels(lab_rank);
set(gcf,'color','white')

% Perform one-way ANOVA
data_anova = [Group1; Group2; Group3; Group4; Group5; Group6; Group7; Group8];
% Create a grouping variable
group_labels = [repmat({'Group1'}, length(Group1), 1); ...
                repmat({'Group2'}, length(Group2), 1); ...
                repmat({'Group3'}, length(Group3), 1); ...
                repmat({'Group4'}, length(Group4), 1); ...
                repmat({'Group5'}, length(Group5), 1); ...
                repmat({'Group6'}, length(Group6), 1); ...
                repmat({'Group7'}, length(Group7), 1); ...
                repmat({'Group8'}, length(Group8), 1)];
[p, tbl, stats] = anova1(data_anova, group_labels,"off");
% Display F-value and p-value
disp(['F-value: ', num2str(tbl{2,5})]);
disp(['p-value: ', num2str(p)]);
%% Community overlapping 
n=90;[u,v] = find(triu(ones(n),1));
% ci=lab(:,6);
ci=idx;
h = zeros(n,max(ci));
num=length(unique(ci));
for i = 1:n 
    idx1 = u == i | v == i;
    h(i,:) = hist(ci(idx1),1:max(ci));%conut the times of reigon i appears in each cluster
end
p = bsxfun(@rdivide,h,sum(h,2));  % Proportion of all edges assigned to each community that included a given node as one of its endpoints.
e = -nansum(p.*log2(p),2); %entropy= pi*log2(pi)
enorm = e/log2(max(ci));% normalized entropy e/log2(k)
P = p;
for c = 1:num
    P(:,c) = P(:,c)+c-1;
end
load('Net7s_aal_new.mat');
[~,index] = sort(Net7s_alter);
dffindex = find(diff(Net7s_alter(index)));
%% Plot the community overlapping
load('/data/users/tfang/ePC_master/Data_source/Community_color.mat')
figure, 
imagesc(P(index,:), [0 num]);colormap(CommunityColormap)  % replace with an user-defined colormap
hold on;
for i = 1:length(dffindex)
    plot([0.5,num + 0.5],dffindex(i)*ones(1,2),'k','LineWidth',0.7);
end
hold on;
y = 1.5:1:(num-0.5);
for i = 1:length(y)
    plot(y(i)*ones(1,2),[0.5,0.5+90],'k','LineWidth',0.7);
end
set(gca,'TickDir','none')
set(gca,'YTickLabel',[])
title('Overlapping communities','FontSize',22)
set(gcf,"Color",'white');
