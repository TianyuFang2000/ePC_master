% load gene expression matrix 
% BrainSpan Atlas is freely available at https://www.BrainSpan.org/static/download.html/
G = readmatrix('/data/users/tfang/ePC_master/Data_source/genes_matrix_csv/expression_matrix.csv');
Gene_expr = G(:,[223 234 231 232 236 227 224 233 235 226 237 225 238 228]);  % gene expression profiles from a donor aged 37 postmenstrual weeks

%  excluding genes with missing values in any region
ind_delet = [];
for i = 1:length(Gene_expr(:,1))
    if ~isempty(find(Gene_expr(i,:)==0))
        ind_delet = [ind_delet i];
    end
end
Gene_expr(ind_delet,:) = [];  % 22494*14 [Num_gene*Num_region]

% co-expression amplitude for each gene is computed as the product of its expression levels in paired regions
Gene_coexpr=fcn_edgets(Gene_expr,1);  % a function available at https://github.com/brain-networks (Faskowitz, J., Esfahlani, F. Z., Jo, Y., Sporns, O. & Betzel, R. F. 
% Edge-centric functional network representations of human cerebral cortex reveal overlapping system-level architecture. Nat Neurosci 23, 1644-1654 (2020). )

%% plot heatmap of T_val
load('/data/users/tfang/ePC_master/Data_source/T_values_GLM_91edges.mat')
load('/data/users/tfang/ePC_master/Data_source/colorbar_blue.mat')
ROI = string(['A1C'; 'DFC'; 'IPC'; 'ITC'; 'M1C'; 'MFC'; 'OFC'; 'S1C'; 'STC'; 'V1C'; 'VFC'; 'AMY'; 'HIP'; 'MDT']);
Mat = zeros(14,14);
mask = tril(true(14),-1);
Mat(mask) = T_val;
figure,imagesc(Mat,[-7,0]),axis square;colorbar,colormap(flip(colorbar_blue))
set(gca,'XTick',1:14,'XTickLabel',ROI,'XTickLabelRotation',45)
set(gca,'YTick',1:14,'YTickLabel',ROI,'FontSize',18,'TickDir','none')
box off
%% Plot gene expression levels
figure,imagesc((zscore(Gene_coexpr'))),colormap(flip(colorbar_blue))
set(gca,"Visible",'off')
%% Generate null models
addpath('/data/users/tfang/eFC_entropy/eFC_entropy/Entropy_developing/Gene')
parfor n = 1:50000
    % Permute the edge labels
    permuted_T_val = permute_edge_labels(T_val, 14);
    r_null(n,:) = corr(permuted_T_val,Gene_coexpr','type','Pearson');
end
[r_observe,p_observe] = corr(T_val,Gene_coexpr','type','Pearson');
%% Calculate p relative to null models
for i=1:length(r_observe)
    P_val(i)= mean(abs(r_observe(i)) < r_null(:,i));
end

%% Obtain the ranked list of significant genes
% load gene ID and label
rowsmetadata = readtable('/data/users/tfang/ePC_master/Data_source/genes_matrix_csv/rows_metadata.csv');
Gene_label = table2array(rowsmetadata(:,"gene_symbol"));
Gene_ID = table2array(rowsmetadata(:,"entrez_id"));
Gene_label(ind_delet) = [];
Gene_ID(ind_delet) = [];
index = find(P_val<(0.001));
selectedgenes = Gene_label(index);
[sort_r,sort_ind] = sort(r_observe(index));
orderedgenes = selectedgenes(sort_ind);
Weighted = sort_r'; 
orderedtable = table((string(orderedgenes)),Weighted);
writetable(orderedtable,'/data/users/tfang/ePC_master/Data_source/OrderedSigGenes_50kpermutations_P001.csv')
% The gene list was then subjected to Gene Ontology (GO) enrichment analysis using Metascape (https://metascape.org)
%% Plot gene list
% load gene list
Gene_list = readtable('/data/users/tfang/ePC_master/Data_source/OrderedSigGenes_50kpermutations_P001.csv');
weight = table2array(Gene_list(:,2));
plot_gene_list(weight)
%% Plot enrichment results
Path='/data/users/tfang/ePC_master/Data_source';
Category={'GO Biological Processes';'GO Molecular Functions';'GO Cellular Components'};
data_path= [Path '/Enrichment_results_developing'];
GO_csv = [data_path '/Enrichment_GO/GO_AllLists.csv' ];
data= readtable(GO_csv);
Term_category = string(table2array(data(:,1)));
figure,
for c=1:3
    index = find(Term_category == Category{c,1});
    Term_description = string(table2array(data(index,4)));
    x = table2array(data(index,7));   % enrichment score
    y = -(table2array(data(index,17)));     % -log10(q) value
    Num_gene = table2array(data(index,12)); % number of genes hit
    nan_term = find(y<-log10(0.05));  % only terms with P_fdr > 0.05 are remained
    x(nan_term) = [];
    y(nan_term) = [];
    Num_gene(nan_term) = [];
    Term_description(nan_term) = [];

    subplot(1,3,c),
    scatter(x,y,Num_gene.*1.5,y,'filled'),colormap(colorbar_blue);
    axis square
    title(Category{c,1});
    set(gca,"fontsize",16,'LineWidth',1.5);
    xlabel('Enrichment score');
    ylabel('-log10(q)')
end
set(gcf,"color",'w');
%% Cell type-specific analysis
GO_size_min = 0;

gene_all_num = length(Gene_label);
load('/data/users/tfang/ePC_master/Data_source/gene_cell_type.mat')
% get GO table
GOTable = struct;
celltypesPSP = readtable('/data/users/tfang/ePC_master/Data_source/celltypes_PSP.csv');
GOname = unique(table2array(celltypesPSP(:,2)));
cell_num = 7;
% calculate the real GO score
GO_results = struct;
gene_select = table2array(Gene_list(:,1));

for iGO = 1:length(GOname)
    GO_results(iGO).GOname = string(GOname(iGO));
    cell_type = table2array(celltypesPSP(:,2));
    ind_ = find (ismember(cell_type,GOname(iGO)));
    if iGO==4 || iGO==5
        A=char(GOname(iGO));A(6)=[];
        GOname(iGO) = ((cellstr(A)));
    end
    
    GOTable.(string(GOname(iGO))) = table2array(celltypesPSP(ind_,1));
    GO_results(iGO).GOsize = sum(ismember(GOTable.(string(GOname(iGO))),gene_select));
    GO_results(iGO).GOscore = (gene_all_num/length(GOTable.(string(GOname(iGO)))))/(length(gene_select)/GO_results(iGO).GOsize);

    if isnan(GO_results(iGO).GOscore), GO_results(iGO).GOscore = 0; end
end
% get the permutated GO score
for iper = 1:10000
    X = Gene_label; % Example vector (replace with your data)
    % Randomly shuffle X
    X_shuffled = X(randperm(length(X)));
    % Assign elements to groups based on specified sizes
    groups = cell(1, 7);
    start_idx = 1;
    for i = 1:7
        group_sizes = length(GOTable.(string(GOname(i))));
        end_idx = start_idx + group_sizes - 1;
        groups{i} = X_shuffled(start_idx:end_idx);
        start_idx = end_idx + 1;
    end

    for iGO = 1:length(GOname)
        GOscore_rand = (gene_all_num/length(GOTable.(string(GOname(iGO)))))/(length(gene_select)/sum(ismember(groups{iGO},gene_select)));

        if isnan(GOscore_rand) 
            GOscore_rand = 0; 
        end

        GO_results(iGO).GOscore_rand(iper,1) = GOscore_rand;
    end
end

% calculated the p-value of each GO
for iGO = 1:length(GOname)
    GO_results(iGO).GOscore_per_p = mean(GO_results(iGO).GOscore_rand >= GO_results(iGO).GOscore);

    GO_results(iGO).GOscore_Z_p = 1 - normcdf(GO_results(iGO).GOscore,mean(GO_results(iGO).GOscore_rand),std(GO_results(iGO).GOscore_rand));
end

all_size = cat(2,GO_results(:).GOsize);
GO_results(all_size < GO_size_min) = [];

% FDR correction using BH method
per_p_fdr = mafdr(cat(2,GO_results(1:cell_num).GOscore_per_p),'BHFDR',true,'showPlot',false);
z_p_fdr = mafdr(cat(2,GO_results(1:cell_num).GOscore_Z_p),'BHFDR',true,'showPlot',false);
for iGO = 1:cell_num
    GO_results(iGO).GOscore_per_p_FDR = per_p_fdr(iGO);
    GO_results(iGO).GOscore_Z_p_FDR = z_p_fdr(iGO);
end
save('/data/users/tfang/ePC_master/Data_source/Cell_specific.mat',"GO_results",'-mat')
%% PLot cell-specific 
load('/data/users/tfang/ePC_master/Data_source/Cell_specific.mat');
Cell_category = [];
Cell_genes = [];
Cell_p_val = [];
for i =1:7
    Cell_genes = [Cell_genes; GO_results_pos(i).GOsize ];
    Cell_p_val = [Cell_p_val; GO_results_pos(i).GOscore_Z_p_FDR ];
end
Cell_category = {'Astrocytes';'Endothelial' ; 'Microglia' ; 'Excitatory Neurons' ;'Inhibitory Neurons' ; ...
     'OPCs'; 'Oligodendrocytes'};
[~,order] = sort(Cell_genes);
Cell_category = flip(Cell_category(order));
% Data for the bars
values = Cell_genes(order);  % Heights of the bars
values = flip(values);
[sorted_Values,order2] = sort(values);
Cell_category_sorted = Cell_category(order2);
% Define custom colors (one color per bar)
randomColors = rand(length(values), 3);
cm = colorbar_blue;
% Create the horizontal bar chart
figure;
bars = barh(sorted_Values, 'FaceColor', 'flat','EdgeColor','none'); % Use 'flat' for individual coloring

% Assign colors to each bar
for i = 1:length(values)
    bars.CData(i, :) = cm(round(i/9*256), :);
end
set(gca, 'YTick', 1:length(Cell_category_sorted), 'YTickLabel', Cell_category_sorted,'FontSize',18);
set(gca, 'LineWidth',1.5,'TickDir','out','box','off');
xlabel('Number of Genes overlapped');
ylabel('Cell types');
grid off;
%% GO enrichment results
GO_table = readtable('/data/users/tfang/ePC_master/Data_source/Enrichment_results_developing/Enrichment_GO/_FINAL_GO.csv');
Enrichment_table = readtable('/data/users/tfang/ePC_master/Data_source/Enrichment_results_developing/metascape_result.xlsx','Sheet','Enrichment');
Enrich_term  = table2array(Enrichment_table(:,1));
Enrich_GOterm = table2array(GO_table(:,9));
k = 0;
termlist = [];
for i = 1:length(Enrich_term)
    if contains(Enrich_term(i),'Summary')
        k = k + 1;
        ind = find(string(Enrich_GOterm) == string(table2array(Enrichment_table(i,4))));
        termlist = [termlist ind];
        
    end
end
Enrich_term_BP = GO_table(termlist,[9 11:13 17 22]);
Cate = GO_table(termlist,7);
Term_table = Enrich_term_BP;
ind = 1:20;
X = table2array(Term_table(ind,6));
X = flip(X); 
Y_lab = flip(table2array(Term_table(ind,1)));
Enrichment_score = table2array(Term_table(ind,3));
Enrichment_score = flip(Enrichment_score);
%% plot top 20 GO terms
color_bar = colorbar_blue;
figure,
bars = barh(abs(X), 'EdgeColor','none','BarWidth',0.75); % Use 'flat' for individual coloring
bars.FaceColor = 'flat';  % Allow individual coloring of bars
for i = 1:20
    
    if (Enrichment_score(i)-2)/2 > 1
        Color_list(i,:) = color_bar(240,:);
    elseif (Enrichment_score(i)-2)/2<0
        Color_list(i,:) = color_bar(1,:);
    else
        Color_list(i,:) = color_bar(round((Enrichment_score(i)-2)./2*length(color_bar(:,1))),:);
    end
end
% Assign colors to each bar
bars.CData = Color_list;
set(gca,'YTick',1:1:20,'YTickLabel',Y_lab,'TickDir','out','TickLength',[0.0075 0.0075],'Fontsize',16)
set(gcf,'Color',"white")
set(gca, 'LineWidth',1.5,'TickDir','out','box','off');
