clear;close all;clc; 
% load required data
load('tFang_dHCPts_0314.mat')
amount_of_FT = 494;
ts_FT = []; % initialize variable
for subj = 1:amount_of_FT  
    ts_toload_FT = zscore(FT_ts(:,:,subj),0,2); % z-score time series of each subject
    ts_FT = [ts_FT ts_toload_FT]; % accumulate all the time series for next analysis
end
ts_FT = ts_FT';
ets_FT = fcn_edgets(ts_FT,1);

%% Perform k-means on edge timeseries merged across full-term infants
k = 10; 
[idx,ci_FT]= kmeans(ets_FT',k,...
    'distance','sqeuclidean',...
    'maxiter',100000);

%% Subject-specific clustering for the FT group
for subj = 1:amount_of_FT 
    ts1_toload_FT = zscore(FT_ts(:,:,1),0,2); % z-score time series of each subject
    ts1_FT_entropy =  ts1_toload_FT'; 
    ets1_FT = fcn_edgets(ts1_FT_entropy,1); % this step does the z-scoring and element-wise products
    TT=find(triu(mat_FT1,1));
    cluster=mat_FT1(TT);
    for i=1:nc
        start_FT(i,:)=mean(ets1_FT(:,find(idx==i)),2)'; % averaging each subject's eTS according to the group-level community labels
    end
  
    idx_FT= kmeans(ets1_FT',k,...
        'distance','sqeuclidean',...
        'maxiter',100000,'start',start_FT); % averaged time series in each cluster were used as the initial centroids 
end