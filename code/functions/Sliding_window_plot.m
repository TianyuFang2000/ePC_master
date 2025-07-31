function Sliding_window_plot(age,ePC_net,ePC_net_betw,x_tick,x_ticklabel)
% sliding window settings
win_leng = 20;
win_step = 10;

Num_subj = length(age);
[~,index] = sort(age);
FT_std1 = [];
FT_med1 = [];
E_val_FT1 = [];            
NumWin = round((Num_subj-win_leng)/win_step);
for i = 1:NumWin
    k = i-1;
    l = 10*k+1;
    u = 10*k+win_leng;
    age_ind = index(l:u);
    for n = 1:8
        E_val_FT1(:,n) = ePC_net(n,age_ind);
        E_val_FT2(:,n) = ePC_net_betw(n,age_ind);
    end
    FT_med1(:,i) = mean(E_val_FT1);
    FT_std1(:,i) = std(E_val_FT1);
    FT_med2(:,i) = mean(E_val_FT2);
    FT_std2(:,i) = std(E_val_FT2);
end

% Plot 
bar_colors = {'k' 'c' 'm' 'r' 'g' 'b' 'y'};
% bar_colors{8} = [0.5 0.6 0.1];

% x = [1 4 10 21 32 39 47];
% y = {'38' '39' '40' '41' '42' '43' '44'};
tit = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'Sub'};
figure,
for i = 1:8
    subplot(2,4,i)
    h1 = shadedErrorBar(1:1:NumWin,FT_med1(i,:),(FT_std1(i,:)./sqrt(20)),'lineprops',{['-s',bar_colors{4} ],'markerfacecolor',bar_colors{4},'markersize',4},'patchSaturation',0.2);
    hold on
    h2 = shadedErrorBar(1:1:NumWin,FT_med2(i,:),(FT_std2(i,:)./sqrt(20)),'lineprops',{['-^',bar_colors{6} ],'markerfacecolor',bar_colors{6},'markersize',4},'patchSaturation',0.2);
   
    set(gca,'YLim',[0.82 0.88]);
    xlabel('PMA/(week)','FontSize',14);
    ylabel('Edge PC','FontSize',14);
    title(tit(i),'FontSize',14);
    set(gca,'YTick',0.82:0.02:0.88);
    set(gca,'xtick',x_tick, 'xticklabel',x_ticklabel)    
end
lgd = legend([h1.mainLine,h2.mainLine],{'Within','Between'},FontSize=14,Orientation="horizontal",Box="on");
set(gcf,'color','white')