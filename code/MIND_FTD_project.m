%% dependencies
% 你需要安装


%% spin test null model
%sphere_coords = table2array(readtable('your_pathway/MIND/code/subParcellation-master/500mm_parcellation_308/308_regions_coordinates.txt'));
%perm_id = rotate_parcellation_single(sphere_coords, 10000);
%coord_l = sphere_coords(1:152,:);
%perm_id_lh = rotate_parcellation_single(coord_l, 10000);
%save('spin_null.mat', 'perm_id', 'perm_id_lh');
load('spin_null.mat');
%% Create MIND-degree: Take the average value for each subjects
subjects = {'DIY001', 'DIY002', '...'};
MIND_name = 'MIND_DK308.csv';

result_table = table();

% Process each subject in a loop
for i = 1:numel(subjects)
    subject_name = subjects{i};
    file_path = ['your clinica outcome/CAPS/subjects/sub-', subject_name, '/ses-001/t1/freesurfer_cross_sectional/sub-', subject_name, '_ses-001/scripts/', MIND_name];
    
    My_map = readtable(file_path);
    VariableName = My_map.Properties.VariableNames;
    My_map = table2array(My_map);
    [row_num, col_num] = size(My_map); 
    
    % Sum each column of the table to calculate MIND_degree
    sum_per_col = sum(My_map, 1);
    MIND_degree = sum_per_col ./ (row_num - 1);
    
    temp_table = array2table(MIND_degree, 'VariableNames', VariableName);
    temp_table.Properties.RowNames = {subject_name};
    result_table = [result_table; temp_table];
end

writetable(result_table, 'your pathway/MIND/outcome/result_table_DK308.csv', 'WriteRowNames', true);


%% Obtain the group average of MIND_degree
% 
file_path = 'your pathway/MIND/outcome/Wu-FTD_MIND_degree.xlsx'; % This file contains cohort characteristic information and the results from the result_table_DK308
data_sheet1 = readtable(file_path, 'Sheet', 1); % characteristic information
data_sheet2 = readtable(file_path, 'Sheet', 2); % result_table_DK308


data_sheet1.Group = categorical(data_sheet1.Group);
unique_groups = unique(data_sheet1.group); 
group_means = table(); 

for i = 1:length(unique_groups)
    group = unique_groups(i);
    group_data = data_sheet2(data_sheet1.group == group, 2:end); 
    group_mean = array2table(mean(table2array(group_data), 1, 'omitnan'), 'VariableNames', group_data.Properties.VariableNames);
    group_mean.Group = repmat(group, height(group_mean), 1); 
    group_means = [group_means; group_mean]; 
end
FTD = table2array(group_means(2,1:308))';
NC = table2array(group_means(1,1:308))';

%Cortical surface visualization (Fig2a)
MINDmap = parcel_to_surface(mymap, 'aparcDK308_fsa5');
f = figure,
    plot_cortical(MINDmap, 'surface_name', 'fsa5', 'color_range', ...
                  [0.09 0.16], 'cmap', 'RdBu_r')


%% group information and setting covariates
MIND_degree_compared = table2array(data_sheet2(ismember(data_sheet1.Group, {'NC','bvFTD','SD'}), 2:end));
Demo_compared = data_sheet1(ismember(data_sheet1.Group, {'NC','bvFTD','SD'}), 2:end);
% covariates
age = Demo_compared.age;
sex = Demo_compared.sex;
TIV = Demo_compared.TIV;
education = Demo_compared.Education_y;
group = Demo_compared.group;
MMSE = Demo_compared.MMSE;
CDR = Demo_compared.CDRGlobal;
FBI = Demo_compared.FBI_total;
FBI_apathy = Demo_compared.FBI_apathy;
FBI_disin = Demo_compared.FBI_disinhibition;
dummy=MIND_degree_compared;


%% global average MIND (FigS2)
global_MIND = mean(MIND_degree_compared, 2);

clear mytstat mypval

  tbl = table(age,sex,group,TIV,education,global_MIND);
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'global_MIND~age+sex+TIV+education+group');
  mytstat=lm.Coefficients{4,3};
  mypval=lm.Coefficients{4,4};


%% case-control regional MIND differences
nregs = 308;
clear mytstat mypval

for region=1:nregs
  tbl = table(age,sex,group,TIV,education,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var6~age+sex+TIV+education+group');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end

mytstat=transpose(mytstat);  % mytstat is the FTD MIND network alterations map
mypval=transpose(mypval);


% BF correction
sigregs_BF=find(mypval<0.05/308); % list of the statistically significant regions
filtered_tstat = mytstat;
filtered_tstat(find(mypval>0.05/308)) = 0;

save('Wu-FTD_MIND_mytstat_DK308.mat', 'mytstat');
% plot surface map (Fig2a)
MINDmap = parcel_to_surface(filtered_tstat, 'aparcDK308_fsa5');
f = figure,
    plot_cortical(MINDmap, 'surface_name', 'fsa5', 'color_range', ...
                  [min(MINDmap) max(MINDmap)], 'cmap', 'RdBu_r')

%% MEAN MIND scatter plot (FigS1)
mean_MIND_degree_NC = MIND_degree_compared(find(group==1),:);
mean_MIND_degree_NC=mean(vertcat(mean_MIND_degree_NC),1);
mean_MIND_degree_FTD = MIND_degree_compared(find(group==2),:);
mean_MIND_degree_FTD=mean(vertcat(mean_MIND_degree_FTD),1);

figure
plot( mean_MIND_degree_NC',mytstat, 'o')
xlabel('Mean control MIND', 'FontSize', 14);
ylabel('Case-control t-value', 'FontSize', 14);
[R, p] = corrcoef(mean_MIND_degree_NC',mytstat);
text(0.5, 0.5, sprintf('R = %.4f, p = %.4f', R(1, 2), p(1, 2)), 'FontSize', 12)
grid on

[Spin_p_perm, null_dist] = perm_sphere_p(mean_MIND_degree_NC',mean_MIND_degree_FTD', perm_id, 'pearson');

%% the correlation between mean value of significant regions and clinical scales (Fig2e and FigS4)
SigregionMean_N = [];  
SigregionMean_P = [];
for n = 1:158
    
SigregionMean_N(n) = mean(MIND_degree_compared(n,find(mypval<0.05/308 & mytstat<0)));
SigregionMean_P(n) = mean(MIND_degree_compared(n,find(mypval<0.05/308 & mytstat>0)));
end
SigregionMean_N = SigregionMean_N';
SigregionMean_P = SigregionMean_P';
SigregionMean = SigregionMean';
tbl = table(subID,age, sex, education,TIV,group, FBI, FBI_apathy, FBI_disin,MMSE,CDR,SigregionMean_N,SigregionMean_P);
writetable(tbl, 'regionMean_FBI_WDNI.csv', 'Delimiter', ',');
% Calculate the correlation coefficient and p-value
tbl.sex = double(tbl.sex);
tbl = tbl(tbl.group == 2, :);
[r, p] = partialcorr(tbl.FBI_apathy, tbl{:,12}, tbl{:,2:5})

%% Reload and calculate the correlation between the discovery cohort and the validation cohort (FigS1)

clear
load('Wu-FTD_MIND_mytstat_DK308.mat');
load('NIFD_MIND_mytstat_DK68.mat');
WDNI_Tmap = mytstat;
NIFD_Tmap= mytstat;

[R, p] = corrcoef(WDNI_Tmap,NIFD_Tmap);
[Spin_p_perm, null_dist] = perm_sphere_p(WDNI_Tmap,NIFD_Tmap, perm_id, 'pearson');

%% PLS analysis
DK_gene = table2array(readtable('your allen pathway/Allen/abagen/allen_expression_DK308.csv'));
DK_gene(:, 1) = []; 
n_lh=152;
X=DK_gene(1:n_lh,:); 
Y=mytstat(1:n_lh,:); 
% z-score:
X=zscore(X);
Y=zscore(Y);

%perform full PLS and plot variance in Y explained by top 10 components

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
dim=10;
plot(1:dim,100*PCTVAR(2,1:dim),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on


%%% plot correlation of PLS component 1 with t-statistic (Fig3b)
figure
plot(XS(:,1),Y,'r.')
[R,p]=corrcoef(XS(:,1),Y) 
xlabel('PLS1 scores','FontSize',14);
ylabel('Case-control t-value(Zscore)','FontSize',14);
grid on

PLS1map = XS(:,1);
[Spin_p_perm, null_dist] = perm_sphere_p(Y,PLS1map, perm_id_lh, 'pearson');
Spin_p_perm  

PLSmap = parcel_to_surface(PLS1map, 'aparcDK308_fsa5');
f = figure,
    plot_cortical(PLSmap, 'surface_name', 'fsa5', 'color_range', ...
                  [-0.15 0.15], 'cmap', 'RdBu_r')

% permutation testing to assess significance of PLS result (FigS5)
rep=10000;
My_P = zeros(1,8);
My_R = zeros(1,8);
for dim=1:8
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=100*PCTVAR(2,1:dim);
Rsquared = temp(dim);
    for j=1:rep
      
        order=randperm(size(Y,1));
        Yp=Y(order,:);

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

        temp=100*PCTVAR(2,1:dim);
        Rsq(j) = temp(dim);
    end
dim
My_R(dim)=Rsquared
My_P(dim)=length(find(Rsq>=Rsquared))/rep
end
figure
plot(1:dim, My_P,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on

dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%% Bootstrap to get the gene list. This code has been modified based on https://github.com/SarahMorgan/Morphometric_Similarity_SZ
% this needs to be imported first
genes15631 = readtable('/media/dai/devC/.jmm_project/code/Allen/abagen/gene15631_list.csv');
genes=table2array(genes15631); 
geneindex=1:15631;

%number of bootstrap
bootnum=10000;

dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);


%align PLS components with desired direction for interpretability 
if R1(1,1)<0 
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);

%print out results
csvwrite('your pathway/MIND/outcome/PLS1_ROIscores_FTD_tmap_gene_pls.csv',XS(:,1));
csvwrite('your pathway/MIND/outcome/PLS2_ROIscores_FTD_tmap_gene_pls.csv',XS(:,2));

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; 
    Xr=X(myresample,:); 
    Yr=Y(myresample,:); 
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); 
      
    temp=stats.W(:,1);
    newW=temp(x1); 
    if corr(PLS1w,newW)<0 
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];
    
    temp=stats.W(:,2);
    newW=temp(x2); 
    if corr(PLS2w,newW)<0 
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; 
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);

%print out results
fid1 = fopen('your pathway/MIND/outcome/PLS1_geneWeights_FTD_tmap.csv','w')
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f, %.3e,%.3e\n', PLS1{i}, geneindex1(i), Z1(i), P_Z1(i));
end
fclose(fid1)

fid2 = fopen('your pathway/MIND/outcome/PLS2_geneWeights_FTD_tmap.csv','w')
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f, %.3e,%.3e\n', PLS2{i},geneindex2(i), Z2(i), P_Z2(i));
end
fclose(fid2)

% Then copy the results into sheets 8 and 9 of Wu-FTD_MIND_degree.xlsx

%% Identify overlapping genes between PLS1-weighted genes and FTD-related genes (Fig3e)
FTD_gene = readtable('your pathway/MIND/outcome/Genes_FTD_370.csv');

file_path = 'your pathway/MIND/outcome/Wu-FTD_MIND_degree.xlsx'; 
PLS1_P = readtable(file_path, 'Sheet', 8);
PLS1_N = readtable(file_path, 'Sheet', 9);

overlappingGenes_P = intersect(PLS1_P.Gene(:,1), FTD_gene.GeneSymbol(:,1));
overlappingGenes_N = intersect(PLS1_N.Gene(:,1), FTD_gene.GeneSymbol(:,1));

%% Examine the spatial correlation between a gene map and the MIND changes in FTD (Fig3f)

genes_of_interest = {'KANSL1'}; 

target_index = find(ismember(genes15631{:,1}, "KANSL1"));
target_gene_data = DK_gene(:, target_index);
Genemap = parcel_to_surface(target_gene_data, 'aparcDK308_fsa5');
f = figure,
    plot_cortical(Genemap, 'surface_name', 'fsa5', 'color_range', ...
                  [-max(Genemap) max(Genemap)], 'cmap', 'RdBu_r')


result_table = table();

for i = 1:length(overlappingGenes_P)
    target_index = find(strcmp(genes, overlappingGenes_P{i}));
    target_gene_data = DK_gene(1:152, target_index);
    gene_data = array2table(target_gene_data, 'VariableNames', {overlappingGenes_P{i}});
    result_table = [result_table, gene_data];
end

writetable(result_table, 'genes_data_lh_P.csv', 'Delimiter', ',');


GENE_name = {result_table.Properties.VariableNames{1:35}}';
r = [1:35]';
P_spin = [1:35;]';
result_corgene = table(GENE_name, r, P_spin, 'VariableNames',{'Neurosynth_name','r','P_spin'});

for i = 1:width(result_table)
    % load one map
    gene_map = result_table{:, i};
    % Spatial_corr
    r_gene = corr(gene_map, mytstat(1:152,:));
    % Spin permutation testing
    [Spin_p_perm, null_dist] = perm_sphere_p(gene_map,mytstat(1:152,:), perm_id_lh, 'pearson');
    
    result_corgene{i, 2} = r_gene;
    result_corgene{i, 3} = Spin_p_perm;
end
                         

%%  Neurosynth decode 

% load Neurosynth 123 cognition maps
% DK308
Neurosynth = readtable('your pathway/Annotation_brainmap/results_all_cognition123_annotations.csv');

Neurosynth_name = {Neurosynth.Properties.VariableNames{1:123}}';
r = [1:123]';
P_spin = [1:123;]';
result_Neurosynth = table(Neurosynth_name, r, P_spin, 'VariableNames',{'Neurosynth_name','r','P_spin'});

for i = 1:width(Neurosynth)
    % load one map
    Neurosynth_map = Neurosynth{:, i};
    % Spatial_corr
    r_Neurosynth = corr(Neurosynth_map, mytstat);
    % Spin permutation testing
    [Spin_p_perm, null_dist] = perm_sphere_p(Neurosynth_map,mytstat, perm_id, 'pearson');

    result_Neurosynth{i, 2} = r_Neurosynth;
    result_Neurosynth{i, 3} = Spin_p_perm;
end


mypath = 'your pathway/MIND/outcome/Decode/'; 
writetable(result_Neurosynth, [mypath,'cor_result_cognition123_WDNI_FTD_Tmap_annotations_DK308.csv']);



%% ENIGMA
% Stratify based on Yeo networks and Economo-Koskinas classes (Fig5a)
% fcn_ctd from https://github.com/netneurolab/hansen_genescognition
% fcn_ltd and fcn_psytd are modified gene lists based on fcn_ctd: Cortical_Layer_gene.csv and Psy6_classes_gens.csv
class_mean = economo_koskinas_spider(mytstat, 'parcellation', 'aparcDK308_fsa5','axis_range', [0.6 4.7])
class_mean = Yeo_spider(mytstat, 'parcellation', 'aparcDK308_fsa5','axis_range', [1 5])

%% Enrichment:Psychiatric disorders, cell types and cortical layers (Fig4c-e)
GeneListFull = readtable(file_path, 'Sheet', 7);  % Full PLS1 gene list
GeneListFull = string(GeneListFull{:, 1});
GeneListSub = readtable(file_path, 'Sheet', 9);  % 8: PLS1+ gene list;9: PLS1+ gene list
GeneListSub = string(GeneListSub{:, 1});

% cell;
[ctd, ctd_null, ctd_pvals, cellnames] = fcn_ctd(GeneListSub, GeneListFull, 10000);

% layers;
[ltd, ltd_null, ltd_pvals, cellnames] = fcn_ltd(GeneListSub, GeneListFull, 10000);

% six Psychiatric disorders genes
[psytd, psytd_null, psytd_pvals, cellnames] = fcn_psytd(GeneListSub, GeneListFull, 10000);


%% epicenter (Fig2c-d)
% mytstat is the FTD MIND network alterations map
% Identify cortical epicenter values (from functional connectivity)
% Load cortico-cortical functional connectivity data
[fc_ctx, fc_ctx_labels, ~, ~] = load_fc();

% Load cortico-cortical structural connectivity data
[sc_ctx, sc_ctx_labels, ~, ~] = load_sc();

fc_ctx_epi              = zeros(size(fc_ctx, 1), 1);
fc_ctx_epi_p            = zeros(size(fc_ctx, 1), 1);
for seed = 1:size(fc_ctx, 1)
    seed_conn           = fc_ctx(:, seed);
    r_tmp               = corrcoef(seed_conn, mytstat);
    fc_ctx_epi(seed)    = r_tmp(1, 2);
    fc_ctx_epi_p(seed)  = spin_test(seed_conn, mytstat, 'surface_name', 'fsa5', ...
                                    'parcellation_name', 'aparc', 'n_rot', 10000, 'type', 'pearson');
end

% Identify cortical epicenter values (from structural connectivity)
sc_ctx_epi              = zeros(size(sc_ctx, 1), 1);
sc_ctx_epi_p            = zeros(size(sc_ctx, 1), 1);
for seed = 1:size(sc_ctx, 1)
    seed_conn           = sc_ctx(:, seed);
    r_tmp               = corrcoef(seed_conn, mytstat);
    sc_ctx_epi(seed)    = r_tmp(1, 2);
    sc_ctx_epi_p(seed)  = spin_test(seed_conn, mytstat, 'surface_name', 'fsa5', ...
                                    'parcellation_name', 'aparc', 'n_rot', 10000, 'type', 'pearson');
end
    
    
% Selecting only regions with p < 0.1 (functional epicenters)
fc_ctx_epi_p_sig = zeros(length(fc_ctx_epi_p), 1);
fc_ctx_epi_p_sig(find(fc_ctx_epi_p < 0.05)) = fc_ctx_epi(fc_ctx_epi_p<0.05);
fc_ctx_epi_p_sig = -1*fc_ctx_epi_p_sig;
f = figure,
    plot_cortical(parcel_to_surface(fc_ctx_epi_p_sig, 'aparc_fsa5'), ...
                'color_range', [-0.5 0.5], 'cmap', 'GyRd_r') 
    
% Selecting only regions with p < 0.1 (structural epicenters)
sc_ctx_epi_p_sig = zeros(length(sc_ctx_epi_p), 1);
sc_ctx_epi_p_sig(find(sc_ctx_epi_p < 0.05)) = sc_ctx_epi(sc_ctx_epi_p<0.05);
sc_ctx_epi_p_sig = -1*sc_ctx_epi_p_sig;
f = figure,
    plot_cortical(parcel_to_surface(sc_ctx_epi_p_sig, 'aparc_fsa5'), ...
                'color_range', [-0.5 0.5], 'cmap', 'GyBu_r')


fc_epicenter = [fc_ctx_epi,fc_ctx_epi_p,fc_ctx_epi_p_sig];
sc_epicenter = [sc_ctx_epi,sc_ctx_epi_p,sc_ctx_epi_p_sig];

csvwrite('fc_epicenter_r_Pspin_sig005p.csv',fc_epicenter');
csvwrite('sc_epicenter_r_Pspin_sig005p.csv',sc_epicenter');
fc_ctx_trans_p_sig = csvread('trans_epicenter.csv');

% file from https://github.com/CNG-LAB/cngopen/tree/26235158c2a0f1847b6cf88b8146c23828251c65/transdiagnostic_gradients/Data/Generated_Data
load('your pathway/transdiagnostic_gradients_psychi/Data/Generated_Data/Epicenters.mat');

f = figure,
    plot_cortical(parcel_to_surface(fc_ctx_trans_p_sig, 'aparc_fsa5'), ...
                'color_range', [-1 1], 'cmap', 'GyRd')

fc_ctx_epi = csvread('fc_epicenter_r_Pspin_sig005p.csv');
fc_ctx_epi = fc_ctx_epi(1,:)';

sc_ctx_epi = csvread('sc_epicenter_r_Pspin_sig005p.csv');
sc_ctx_epi = sc_ctx_epi(1,:)';

col_epi = csvread('col_alteration_epicenter_FTD_psy_FC.csv');


temp_R = corr(fc_ctx_epi,fc_ctx_trans);
% Spin permutation testing
[P_spin, r_dist]   = spin_test(fc_ctx_epi,fc_ctx_trans, 'surface_name', 'fsa5', ...
                                   'parcellation_name', 'aparc', 'n_rot', 10000, ...
                                   'type', 'pearson');
