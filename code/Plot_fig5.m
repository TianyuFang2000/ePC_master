%% Plot Fig. 5
% load data
load('/data/users/tfang/ePC_master/data/FT_behavior.mat')
load('/data/users/tfang/ePC_master/data/ePC_VIS_High-order.mat')
% calculate 95% CI
R = [];
p = [];
y=FT_lang;
[R(1),p(1)]=corr(y1',y);
[R(2),p(2)]=corr(y2',y);
y=FT_cog;
[R(3),p(3)]=corr(y1',y);
[R(4),p(4)]=corr(y2',y);
y=FT_mot;
[R(5),p(5)]=corr(y1',y);
[R(6),p(6)]=corr(y2',y);
fdr = mafdr(p,'BHFDR',1);
n = length(y1);        % Sample size
alpha = 0.05;  % For 95% CI
for i=1:6
    [r_lo, r_hi] = corrCI(R(i), n, alpha);
    R_lo(i) = r_lo;
    R_hi(i) = r_hi;
    fprintf('95%% CI for r = %.3f: [%.3f, %.3f]\n', R(i), r_lo, r_hi);
end
%% Plot correlation between ePC and BSID-III score

x_lab = {'Edge PC between Vis' 'Edge PC between LIM / FPN /DMN'};
figure,
for i = 1:2
    y = FT_lang;
    x = eval(['ePC' num2str(i)]);
    x = x';
    [r,~]=corr(x,y);
    p_fdr = fdr(i);
    [X, I] = sort(x);
    Y = y(I);
    x = X;
    y = Y;
    colorbar = slanCM(102);
    [p, s] = polyfit(x, y, 1);
    Y1 = polyval(p, x);
    [yfit, dy] = polyconf(p, x, s, 'predopt', 'curve');

    cc = colorbar(20,:);
    subplot(1,2,i)
    patch([x; flipud(x)], [yfit - dy; flipud(yfit + dy)], 'k', 'FaceA', 0.1, 'EdgeA', 0);
    hold on;
    s = scatter(x, y, 30, 'filled','CData',[0.15 0.15 0.15],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none','Marker','o');
    s.LineWidth = 0.5;
    L1=plot(x, Y1, 'Color', cc, 'linewidth', 1.5);
    text(max(x)*0.9,20*0.15,['\it r = ' num2str(r,'%.4f')],'FontSize',16,'Color','k','Interpreter','tex');

    if p_fdr>0.0001
        txt=['\it p = ' num2str(p_fdr,'%.4f')];
    else
        txt='\it p < 0.0001';
    end
    text(max(x)*0.9,20*0.05,txt,'FontSize',16,'Color','k','Interpreter','tex');
    set(gca,'ytick',0:10:30,'XTick',0.76:0.04:0.88,'XTickLabel',0.76:0.04:0.88);
    ylim([0,30])
    % set(gca,'ytick',0:5:20,'XTick',0.76:0.04:0.88,'XTickLabel',0.76:0.04:0.88);
    xlim([0.76 0.88])
    % ylim([0,20])
    set(gca, 'LineWidth',1.5,'TickDir','out','box','off','Fontsize',14);
    xlabel(x_lab(i));
    % ylabel('BSID-III Cognitive')
    ylabel('BSID-III Language')
    axis square
end

