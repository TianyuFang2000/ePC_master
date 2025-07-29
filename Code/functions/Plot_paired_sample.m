function Plot_paired_sample(EE1,EE2)
% input: EE1 and EE2 are paired vector of size [N*1]

% specify colors you prefer
box1_color = [232 132 130]/255;
box2_color = [142 139 254]/255;
dur = 0.15;
% EE1 = eFC_PC_FT_net;
% EE2 = eFC_PC_FT_net_betw;
tit1 = {' ' 'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'Sub' ' '};

fig = figure('Units','normalized');
ax = axes('Parent',fig);hold on;
set(ax,'LineWidth',1.1,'Box','off','TickDir','in',...
'XMinorTick','off','YMinorTick','off','XGrid','off','YGrid','off','GridLineStyle','-',...
'FontName','Times New Roman','FontSize',12,'GridAlpha',.29,'TickLength',[.006 .015])

ax.YLabel.String = 'Edge PC';
ax.YLabel.FontSize = 14;

ax.YLabel.FontSize = 14;

ax.XTickLabel = tit1;

for i = 1:8
    Data1 = EE1(i,:)';fullData1 = Data1;
    scatter(ax,Data1.*0+i-dur+.07.*(rand(size(Data1,1),1)-.5).*2,Data1,40,'filled',...
            'CData',box1_color,'LineWidth',1.2,'MarkerFaceAlpha',0.3);
    outliBool = isoutlier(Data1,'quartiles');
    Data1(outliBool) = [];

    hold on
    Data2 = EE2(i,:)';fullData2=Data2;
    scatter(ax,Data2.*0+i+dur+.07.*(rand(size(Data2,1),1)-.5).*2,Data2,40,'filled',...
            'CData',box2_color,'LineWidth',1.2,'MarkerFaceAlpha',0.3);
    outliBool = isoutlier(Data2,'quartiles');
    Data2(outliBool) = [];

    qt125 = quantile(fullData1,0.25);
    qt175 = quantile(fullData1,0.75);
    med1 = median(fullData1);

    H1 = plot(ax,[i,i]-dur,[max(Data1),qt175],'LineWidth',1.2,'Color',box1_color);
    plot(ax,[i,i]-dur,[min(Data1),qt125],'LineWidth',1.2,'Color',box1_color)
    fill(ax,i-dur+.1.*1.*[-1 1 1 -1],[qt125,qt125,qt175,qt175],[1,1,1],...
        'FaceAlpha',0.95,'FaceColor','none','EdgeColor',box1_color,'LineWidth',1.2);
    plot(ax,[-.1,.1]+i-dur,[med1,med1],'LineWidth',1.2,'Color',box1_color)

    qt25 = quantile(fullData2,0.25);
    qt75 = quantile(fullData2,0.75);
    med = median(fullData2);

    H2 = plot(ax,[i,i]+dur,[max(Data2),qt75],'LineWidth',1.2,'Color',box2_color);
    plot(ax,[i,i]+dur,[min(Data2),qt25],'LineWidth',1.2,'Color',box2_color)
    fill(ax,i+dur+.1.*1.*[-1 1 1 -1],[qt25,qt25,qt75,qt75],[1,1,1],...
        'FaceAlpha',0.95,'FaceColor','none','EdgeColor',box2_color,'LineWidth',1.2);
    plot(ax,[-.1,.1]+i+dur,[med,med],'LineWidth',1.2,'Color',box2_color)
    ax.XLim = [0,9];
    ax.YLim = [0.8,0.9];
    set(gca,'YTick',0.8:0.05:0.9);
    set(gcf,'color','white');
    legend([H1,H2],{'Within','Between'},FontSize=14,Orientation="horizontal",Box="off");
end
