% loding ePC values at system level
load('/data/users/tfang/ePC_master/Data_source/ePC_network.mat')
%% Developmental curve at system level using sliding window analysis
% Full-term cohort
FT_age = table2array(Info(indFT,7));
x_tick = [1 4 10 21 32 39 47];
x_ticklabel = {'38' '39' '40' '41' '42' '43' '44'};
Sliding_window_plot(FT_age,eFC_PC_FT_net,eFC_PC_FT_net_betw,x_tick,x_ticklabel)

% Preterm cohort (merging two PT groups)
PT1_age = table2array(Info(indPT1,7));
PT2_age = table2array(Info(indPT2,7));
PT_age = cat(1,PT1_age,PT2_age);
eFC_PC_PT_net = [eFC_PC_PT1_net,eFC_PC_PT2_net];
eFC_PC_PT_net_betw = [eFC_PC_PT1_net_betw,eFC_PC_PT2_net_betw];
x_tick = [1.5 5 9 15 19.5 24.4];
x_ticklabel = {'30' '33' '35' '37' '40' '42'};
Sliding_window_plot(PT_age,eFC_PC_PT_net,eFC_PC_PT_net_betw,x_tick,x_ticklabel)

%% plot the within-network and between-network ePC simultaneously
Plot_paired_sample(eFC_PC_PT1_net,eFC_PC_PT1_net_betw)
Plot_paired_sample(eFC_PC_PT2_net,eFC_PC_PT2_net_betw)
Plot_paired_sample(eFC_PC_FT_net,eFC_PC_FT_net_betw)