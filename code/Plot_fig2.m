load('/data/users/tfang/ePC_master/data/Net7s_aal.mat')
load('/data/users/tfang/ePC_master/data/eFC_PC_dhcp.mat')
load('/data/users/tfang/ePC_master/data/colorbar_blue.mat')
%% Plot ePC matrix for each group
[gx,gy,order] = grid_communities(Net7s);
ePC_mean(:,1) = mean(ePC_PT1,2);
ePC_mean(:,2) = mean(ePC_PT2,2);
ePC_mean(:,3) = mean(ePC_FT,2);
tit = {'PT1' 'PT2' 'FT'};
figure,
for s = 1:3
    mat_ePC = zeros(90);
    mat_ePC(triu(ones(90),1) > 0) = ePC_mean(:,s); % this step gives you the indices of each edge
    mat_ePC = mat_ePC + mat_ePC';
    subplot(1,3,s),
    imagesc(mat_ePC(order,order),[0.83 0.87]); axis square; colormap(colorbar_blue)
    hold on; plot(gx,gy,'k','linewidth',2); title(tit{1,s},FontSize=22)
    set(gcf,'color','white')
end
%% Plot ePC distributions of 3 groups

ePC1 = mean(ePC_PT1,2);
ePC2 = mean(ePC_PT2,2);
ePC3 = mean(ePC_FT,2);
color_bar = colorbar_blue;
colorList = color_bar([220 160 100],:);
group_name = {'' 'PT1' 'PT2' 'FT' ''};
fig = figure('Units','normalized');
ax = axes('Parent',fig);hold on;
set(ax,'LineWidth',1.1,'Box','off','TickDir','out',...
'XMinorTick','on','YMinorTick','on','XGrid','off','YGrid','off',...
'FontName','Times New Roman','FontSize',12,'GridAlpha',.09,'Tickdir','none')
ax.YLabel.String = 'Edge PC';
ax.YLabel.FontSize = 14;
ax.XTickLabel = group_name;
for i = 1:3
    T_val = eval(['ePC' num2str(i)]);
    Data=T_val;fullData = Data;
    [f,xi] = ksdensity(Data);
    f = f./max(f)/2;
    fill(ax,[0,f,0].*1.2+i+.15,[xi(1),xi,xi(end)],colorList(i,:),'FaceAlpha',0.7,...
            'EdgeColor',colorList(i,:),'LineWidth',1.2,'EdgeAlpha',1);
    scatter(ax,Data.*0+i-.1+.08.*(rand(size(Data,1),1)-.5).*2,Data,45,'filled',...
            'CData',colorList(i,:),'LineWidth',1,'MarkerFaceAlpha',0.05)
    outliBool=isoutlier(Data,'quartiles');
    Data(outliBool) = [];
    qt25 = quantile(fullData,0.25);
    qt75 = quantile(fullData,0.75);
    med = median(fullData);
    plot(ax,[i,i]-.1,[max(Data),qt75],'LineWidth',1.2,'Color','k')
    plot(ax,[i,i]-.1,[min(Data),qt25],'LineWidth',1.2,'Color','k')
    fill(ax,i-.1+.1.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'FaceColor','none','EdgeColor','k','LineWidth',1.2);
    plot(ax,[-.1,.1]+i-.1,[med,med],'LineWidth',1.2,'Color','k')
    ax.XLim = [0.5,4];
    ax.YLim = [0.8,0.9];
    set(gca,'XTick',0:1:4,'Fontsize',20);
    set(gca,'YTick',0.8:0.05:0.9);
    set(gcf,'color','white');
end
axis square 

%% ePC by node
PC1 = mean(PC_bynode_PT1,2);
PC2 = mean(PC_bynode_PT2,2);
PC3 = mean(PC_bynode_FT,2);

fig = figure('Units','normalized');
ax = axes('Parent',fig);hold on;
set(ax,'LineWidth',1.1,'Box','off','TickDir','out',...
'XMinorTick','on','YMinorTick','on','XGrid','off','YGrid','off',...
'FontName','Times New Roman','FontSize',12,'GridAlpha',.09,'Tickdir','none')
ax.YLabel.String = 'Edge PC by node';
ax.YLabel.FontSize = 14;
ax.XTickLabel = group_name;
for i = 1:3
    T_val = eval(['PC' num2str(i)]);
    Data = T_val;fullData=Data;
    [f,xi] = ksdensity(Data);
    f = f./max(f)/2;
    fill(ax,[0,f,0].*1.2+i+.15,[xi(1),xi,xi(end)],colorList(i,:),'FaceAlpha',0.7,...
            'EdgeColor',colorList(i,:),'LineWidth',1.2,'EdgeAlpha',1);
    scatter(ax,Data.*0+i-.1+.08.*(rand(size(Data,1),1)-.5).*2,Data,45,'filled',...
            'CData',colorList(i,:),'LineWidth',1,'MarkerFaceAlpha',0.5)
    outliBool = isoutlier(Data,'quartiles');
    Data(outliBool) = [];
    qt25 = quantile(fullData,0.25);
    qt75 = quantile(fullData,0.75);
    med = median(fullData);
    plot(ax,[i,i]-.1,[max(Data),qt75],'LineWidth',1.2,'Color','k')
    plot(ax,[i,i]-.1,[min(Data),qt25],'LineWidth',1.2,'Color','k')
    fill(ax,i-.1+.1.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'FaceColor','none','EdgeColor','k','LineWidth',1.2);
    plot(ax,[-.1,.1]+i-.1,[med,med],'LineWidth',1.2,'Color','k')
    ax.XLim = [0.5,4];
    ax.YLim = [0.82,0.86];
    set(gca,'XTick',0:1:4,'Fontsize',20);
    set(gca,'YTick',0.82:0.02:0.86);
    set(gcf,'color','white');
end
axis square 