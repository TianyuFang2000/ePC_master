load('/data/users/tfang/ePC_master/Data_source/PMA_effect_glm.mat')
load('/data/users/tfang/ePC_master/Data_source/Net7s_aal.mat')
load('/data/users/tfang/ePC_master/Data_source/blue2red.mat')
[~,order]=sort(Net7s);
P_thres=0.01;
T1=zeros(4005,1);
T1(find(P_PT1<P_thres))=T_val_PT1(find(P_PT1<P_thres));
T2=zeros(4005,1);
T2(find(P_PT2<P_thres))=T_val_PT2(find(P_PT2<P_thres));
T3=zeros(4005,1);
T3(find(P_FT<P_thres))=T_val_FT(find(P_FT<P_thres));
T4=zeros(4005,1);
T4(find(P_PT1<P_thres/4005))=T_val_PT1(find(P_PT1<P_thres/4005));
T5=zeros(4005,1);
T5(find(P_PT2<P_thres/4005))=T_val_PT2(find(P_PT2<P_thres/4005));
T6=zeros(4005,1);
T6(find(P_FT<P_thres/4005))=T_val_FT(find(P_FT<P_thres/4005));
figure,
for i=1:6
    subplot(2,3,i)
    M = zeros(90);
    M(triu(ones(90),1) > 0) = eval(['T' num2str(i)]);
    M = M + M';
    imagesc(M(order,order),[-6 6]); 
    axis square;
    colormap(blue2red); 
    hold on; 
    plot(gx,gy,'k','linewidth',0.45); 
    box off
    set(gca,'Visible','off')
    set(gca,'TickDir','none')
    set(gca,'xtick',0:0:0);
    set(gca,'ytick',0:0:0);
end
set(gcf,'color','white')
%% Plot proportion by network
PP(:,1) = P_PT1;
PP(:,2) = P_PT2;
PP(:,3) = P_FT;
figure,
for g = 1:3
    P = PP(:,g);
    T = zeros(4005,1);
    index = find(P<0.01);
    index2 = find(P<0.01/4005);
    T(index2) = 1;
    M = zeros(90);
    M(triu(ones(90),1) > 0) = T;
    M = M + M';
    MM = M(order,order);
    for i = 1:8
        row = find(Net7s_sort==i);
        for j = 1:8
            col = find(Net7s_sort==j);
            Q(i,j) = length(find(MM(row,col)))/length(row)/length(col);
        end
    end
    subplot(1,3,g),
    imagesc(Q,[0 0.4]); axis square; colormap(blue2red)
    set(gcf,'color','white')
    x = 1:1:8;
    y = {'Vis','SMN','DAN','VAN','Lim','FP','DMN','Sub'};
    set(gca,'ytick',x, 'yticklabel',y)
    set(gca,'xtick',x, 'xticklabel',y)
    set(gca,'XTickLabelRotation',90)
    set(gca,'Fontsize',24);
    set(gca,'TickDir','none');
end
