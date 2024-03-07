% Genetic decoding analysis
% set which of the AHBA donors should be included.
% written 2023 by Svenja Kuchenhoff, contributions by Mohamad Amin.
% input: sex difference maps + transcriptomic maps
% output: correlation between gene expr. maps and sex difference measures
% makes use of a previously loaded transcriptomic map in schaefer400
% parcellation from the BrainStat toolbox

loadold = 1;

if loadold == 0
    clear all
    close all
    loadold = 0;
end

female_donor = 0;
male_donor = 0;
all_donors = 1;
old_donors = 0;

% if bonferroni
bonf = 1;

saveall = 0;
count = 0;


% make a path to directory
addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';
scriptDir = fullfile(homeDir,'script');
addpath(genpath(scriptDir));
dataDir = fullfile(homeDir, 'data');
addpath(genpath(dataDir));
outDir = fullfile(homeDir,'output');
addpath(genpath(outDir));
figDir = fullfile(homeDir, 'figures/HCP/gendec');
addpath(figDir);

if loadold == 1
    load (fullfile(outDir, 'geneticdecoding_femdon.mat'));
    femresults = [Steroidreceptorresults;Steroidsynthesisresults];
    load (fullfile(outDir, 'geneticdecoding_maledon.mat'));
    maleresults = [Steroidreceptorresults;Steroidsynthesisresults];
    load (fullfile(outDir, 'geneticdecoding_allldons.mat'));
    allresults = [Steroidreceptorresults;Steroidsynthesisresults];

    idx_female_in_male = ismember(cellstr(femresults.Genes),cellstr(maleresults.Genes));
    idx_female_in_all = ismember(cellstr(femresults.Genes),cellstr(allresults.Genes));
    idx_all_in_male = ismember(cellstr(allresults.Genes),cellstr(maleresults.Genes));
    
    fem_male_corr_grad = corr(femresults.Gradient_corr(idx_female_in_male), maleresults.Gradient_corr, 'Type', 'Spearman');
    fem_male_corr_skew = corr(femresults.Skew_corr(idx_female_in_male), maleresults.Skew_corr, 'Type', 'Spearman');
    fem_male_corr_mean = corr(femresults.Mean_corr(idx_female_in_male), maleresults.Mean_corr, 'Type', 'Spearman');
    
    fem_all_corr_grad = corr(femresults.Gradient_corr(idx_female_in_all), allresults.Gradient_corr, 'Type', 'Spearman');
    fem_all_corr_skew = corr(femresults.Skew_corr(idx_female_in_all), allresults.Skew_corr, 'Type', 'Spearman');
    fem_all_corr_mean = corr(femresults.Mean_corr(idx_female_in_all), allresults.Mean_corr, 'Type', 'Spearman');
    
    male_all_corr_grad = corr(allresults.Gradient_corr(idx_all_in_male), maleresults.Gradient_corr, 'Type', 'Spearman');
    male_all_corr_skew = corr(allresults.Skew_corr(idx_all_in_male), maleresults.Skew_corr, 'Type', 'Spearman');
    male_all_corr_mean = corr(allresults.Mean_corr(idx_all_in_male), maleresults.Mean_corr, 'Type', 'Spearman');
    

    keyboard
end

        
%% Set up transcriptomic maps and load data 


load genes_schaefer400.mat
genelabels = genes_schaefer400.gene_names{1, 1}; 
genelabels = (string(genelabels))';

expression_schaefer400 = genes_schaefer400.expression{1, 1};

Gen_chrom = readtable('Genes_allbychrom.csv'); % doesnt work
% for some reason, I cannot read in the labels 'X' and 'Y' for chromosomes
Gen_sexrec = readtable('Genes_sexsteroidrec.csv');
Gen_sexsynth = readtable('Genes_sexsteroidsynth.csv');

keepX = find(Gen_chrom.chromosome == -1);
keepY = find(Gen_chrom.chromosome == -2);

Genes_of_interest = [Gen_sexrec.gene_symbol; "ESR2"; "ESRRA"; "ESRRG"; "ESRRB"; "GREB1"; "PGR"; "SHBG"];
% subselection of those genes I have maps for
[labels_steroidrec, index_sexrec, keep_steroidrec] = intersect(Genes_of_interest, genelabels);
transm_steroidrec = expression_schaefer400(:,keep_steroidrec);


Genes_synth_of_interest = [Gen_sexsynth.gene_symbol; "PIBF1"];
% subselection of those genes I have maps for
[labels_steroidsynth, index_sexsynth, keep_steroidsynth] = intersect(Genes_synth_of_interest, genelabels);
transm_steroidsynth = expression_schaefer400(:,keep_steroidsynth);


[labels_Xchrom, index_genchrom, keep_Xchroms] = intersect(Gen_chrom.gene_symbol(keepX), genelabels);
transm_Xchroms = expression_schaefer400(:,keep_Xchroms);

[labels_Ychrom, index_genchrom, keep_Ychroms] = intersect(Gen_chrom.gene_symbol(keepY), genelabels);
transm_Ychroms = expression_schaefer400(:,keep_Ychroms);

%keyboard

% add AHBA data seggregated by sex
all_genes = readtable('/ahba_data_by_sex/exp_df_all.csv');
gene_names_all = all_genes.Properties.VariableNames;
female_genes = readtable('/ahba_data_by_sex/exp_df_F.csv');
gene_names_female = female_genes.Properties.VariableNames;
male_genes = readtable('/ahba_data_by_sex/exp_df_M_corr_intensity.csv');
gene_names_male = male_genes.Properties.VariableNames;

[labels_steroidrec_fem, index_sexrec_fem, keep_steroidrec_fem] = intersect(Genes_of_interest, gene_names_female);
[labels_steroidsynth_fem, index_sexrec_fem, keep_steroidsynth_fem] = intersect(Genes_synth_of_interest, gene_names_female);

[labels_steroidrec_male, index_sexrec_male, keep_steroidrec_male] = intersect(Genes_of_interest, gene_names_male)
[labels_steroidsynth_male, index_sexrec_male, keep_steroidsynth_male] = intersect(Genes_synth_of_interest, gene_names_male);

[labels_steroidrec_all, index_sexrec_all, keep_steroidrec_all] = intersect(Genes_of_interest, gene_names_all);
[labels_steroidsynth_all, index_sexrec_all, keep_steroidsynth_all] = intersect(Genes_synth_of_interest, gene_names_all);


if female_donor == 1
    transm_steroidsynth = table2array(female_genes(:,keep_steroidsynth_fem));
    transm_steroidrec = table2array(female_genes(:,keep_steroidrec_fem));
    labels_steroidrec = labels_steroidrec_fem;
    labels_steroidsynth = labels_steroidsynth_fem;
elseif male_donor == 1
    transm_steroidsynth = table2array(male_genes(:,keep_steroidsynth_male));
    transm_steroidrec = table2array(male_genes(:,keep_steroidrec_male));
    labels_steroidrec = labels_steroidrec_male;
    labels_steroidsynth = labels_steroidsynth_male;
elseif all_donors == 1
    transm_steroidsynth = table2array(all_genes(:,keep_steroidsynth_all));
    transm_steroidrec = table2array(all_genes(:,keep_steroidrec_all));
    labels_steroidrec = labels_steroidrec_all;
    labels_steroidsynth = labels_steroidsynth_all;
end


%% for spatial specificity, compute baseline of 'brain genes'

if old_donors == 1
    imagesc(expression_schaefer400(1:200,:));
    all_brain_genes = (expression_schaefer400(1:200,:))';
elseif male_donor == 1
        male_array = table2array(male_genes);
        imagesc(male_array(1:200,2:end));
        all_brain_genes = (male_array(1:200,2:end))';
elseif female_donor == 1
        female_array = table2array(female_genes);
        imagesc(female_array(1:200,2:end));
        all_brain_genes = (female_array(1:200,2:end))';
elseif all_donors == 1
        all_array = table2array(all_genes);
        imagesc(all_array(1:200,2:end));
        all_brain_genes = (all_array(1:200,2:end))';
end


imagesc(all_brain_genes);
% first exclude nans, otherwise it doesnt work
% there are a few parcels for which I just don't have a value.
exclude_cols = find(~isnan(all_brain_genes(1,:)));
to_exclude = find(isnan(all_brain_genes(1,:)));
all_brain_genes = all_brain_genes(:,exclude_cols);
imagesc(all_brain_genes);

% run PCA of all genes.
% Rows of X correspond to observations (15631 genes) and columns correspond to
% variables (190 parcels in the left hemisphere)
[coeff, score, latent, tsquared, explained, mu] = pca(all_brain_genes);
% coeff > each columns = principal component in descending order of variance.
% latent = principal component variance -> first = .87
% explained = % total variance explained bz each PC
PC(exclude_cols) = coeff(:,1);
PC(to_exclude) = NaN;
PC(200:400) = NaN;
% now, correlation of PC brain genes with effect maps


%% Compute general overlap with sex-hormone related genes
load HCP_T1wT2w_sexdiffsmaps.mat

per_mom(1000,1) = load('gradient_perm.mat')
per_mom(1000,2) = load('profile_mean_perm.mat')
per_mom(1000,3) = load('skew_perm.mat')

% and addtional GLM.
repo = zeros(1000,3);
for measure = 1:3
    if measure == 1                
    namemoment = 'gradient';
    else if measure == 2
        namemoment = 'profile_mean';
        else if measure == 3
            namemoment = 'profile_skewness';
            end
        end
    end
    X = [transm_steroidsynth, transm_steroidrec];

    y = results.tvals(:,measure); %
    labels = [labels_steroidsynth; labels_steroidrec];
    % Convert to a table
    T = array2table(X, 'VariableNames', labels);
    % Add the response variable to the table
    T.ResponseVar = y;
    % Start with the response variable part of the formula
    formula = 'ResponseVar ~ ';
    % Add each predictor to the formula
    for i = 1:length(labels)
        formula = [formula labels{i}];
        if i < length(labels)
            formula = [formula ' + ']; % Add + between predictors
        end
    end
    mdl = fitglm(T, formula);
    trueF(:,measure) = mdl.devianceTest{:,'FStat'}(2);
    % Display the model summary
    namemoment(measure)
    for j = 1:1000
        j
        y = x_perm(:,j); % Your 400x1 response vector t_schaefer_400(:, 2); %
        labels = [labels_steroidsynth; labels_steroidrec];
        % Convert to a table
        T = array2table(X, 'VariableNames', labels);
        % Add the response variable to the table
        T.ResponseVar = y;
        % Start with the response variable part of the formula
        formula = 'ResponseVar ~ ';
        % Add each predictor to the formula
        for i = 1:length(labels)
            formula = [formula labels{i}];
            if i < length(labels)
                formula = [formula ' + ']; % Add + between predictors
            end
        end
        mdl = fitglm(T, formula);
        repo(j,measure) = mdl.devianceTest{:,'FStat'}(2);
    end
    %disp(mdl);
    %h = fdr_bh(mdl.Coefficients.pValue,0.05)
    sum(repo(:,measure) > trueF(:,measure))./1000
end

%% Genes of interest spatial correlations.

load HCP_T1wT2w_sexdiffsmaps.mat
for measure = 1:3      
    if measure == 1                
        namemoment = 'gradient';
        else if measure == 2
            namemoment = 'profile_mean';
            else if measure == 3
                namemoment = 'profile_skewness';
                end
            end
    end
    
    % this has been calculated in sex_diffs_3measures.mat
    t_schaefer_400(:, measure) = results.tvals(:,measure);

    %t_schaefer_400(:, measure) = results.Cohensd(:,measure);
    transm_steroidsynth(to_exclude,:) = NaN;
    transm_steroidsynth(200:400,:) = NaN;

    transm_steroidrec(to_exclude,:) = NaN;
    transm_steroidrec(200:400,:) = NaN;

    no_all_tests = size(transm_steroidrec, 2) + size(transm_steroidsynth, 2);

    % nan out these values plus the left hemisphere.
    t_schaefer_400(200:400,:) = NaN;
    t_schaefer_400(to_exclude) = NaN;

    % baseline correlation
    [rho,pval] = corr(t_schaefer_400(:, measure), PC', 'Type', 'Spearman', 'rows', 'complete');
        baselinecorr(1, measure) = rho;
        baselinecorr(2, measure) = pval;
      
    %compute correlations with Cohen's-d map, but only LH
    % for steroidreceptor transcriptomic maps
    for i = 1 : size(transm_steroidrec, 2)
        [rho,pval] = corr(t_schaefer_400(:, measure), transm_steroidrec(:,i), 'Type', 'Spearman', 'rows', 'complete');
        steroidreccorr(i,1, measure) = rho
        steroidreccorr(i,2, measure) = pval;
        % plot null distributions, correlation, and spin-p vals  
        if pval < 0.05 
            % do spin testing additionally (careful! This is whole-brain!)
            [p_spin, r_dist] = spin_test(t_schaefer_400(:,measure), transm_steroidrec(:,i), 'surface_name',...
            'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
            'type', 'Spearman');
            f = figure
            his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                               'facealpha', 1, 'linewidth', 0.5);
            l = line(rho, 5,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
            xlabel(['Null correlations' newline (sprintf('(sex diffs between %s and transcriptomic map %d)', namemoment, i))])
            legend(l,['{\it r}=' num2str(round(rho, 2)) newline ...
                              '{\it p}=' num2str(round(p_spin,3 ))])
            count = count + 1
            FigName{(count)} = sprintf('/HCP_Tw1T2w_sterrec_gene_spin_%s%d.png', namemoment, i);
            % store results
            steroid.recspin(i,1, measure) = num2cell(pval);
            steroidreccorr(i,3,measure) = p_spin;
            steroid.recspin = {p_spin, r_dist};
        end
    end

    % for steroidynthesis transcriptomic maps
    for i = 1 : size(transm_steroidsynth, 2)
        [rho,pval] = corr(t_schaefer_400(:, measure), transm_steroidsynth(:,i), 'Type', 'Spearman', 'rows', 'complete');
        steroidsynthcorr(i,1,measure) = rho;
        steroidsynthcorr(i,2, measure) = pval;
        % plot null distributions, correlation, and spin-p vals  
        if pval < 0.05 
           % do spin testing additionally (careful! This is whole-brain!)
            [p_spin, r_dist] = spin_test(t_schaefer_400(:,measure), transm_steroidsynth(:,i), 'surface_name',...
            'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
            'type', 'Spearman');
            f = figure
            his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                               'facealpha', 1, 'linewidth', 0.5);
            l = line(rho, 5,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
            xlabel(['Null correlations' newline (sprintf('(sex diffs between %s and transcriptomic map %d)', namemoment, i))])
            legend(l,['{\it r}=' num2str(round(rho, 2)) newline ...
                              '{\it p}=' num2str(round(p_spin,3 ))])
            count = count + 1
            FigName{(count)} = sprintf('//HCP_Tw1T2w_stersynth_gene_spin_%s%d.png', namemoment, i);
            % store results
            steroid.synthspin(i,1, measure) = num2cell(pval);
            steroidsynthcorr(i,3,measure) = p_spin;
            steroid.synthspin = {p_spin, r_dist};
        end
    end

    % for genes on the Y chromosome
    for i = 1 : size(transm_Ychroms, 2)
        [rho,pval] = corr(t_schaefer_400(1:200, measure), transm_Ychroms(1:200,i), 'Type', 'Spearman', 'rows', 'complete');
        chromsYcorr(i,1, measure) = rho;
        chromsYcorr(i,2, measure) = pval;
    end
    % for genes on the X chromosome
    for i = 1 : size(transm_Xchroms, 2)
        [rho,pval] = corr(t_schaefer_400(1:200, measure), transm_Xchroms(1:200,i), 'Type', 'Spearman', 'rows', 'complete');
        chromsXcorr(i,1, measure) = rho;
        chromsXcorr(i,2, measure) = pval;
    end

end

% save results
Genes = labels_steroidrec;
Gradient_corr = steroidreccorr(:,1,1);
Gradient_p = steroidreccorr(:,2,1);
Gradient_pspin = steroidreccorr(:,3,1);
Mean_corr = steroidreccorr(:,1,2);
Mean_p = steroidreccorr(:,2,2);
Mean_pspin = steroidreccorr(:,3,2);
Skew_corr = steroidreccorr(:,1,3);
Skew_p = steroidreccorr(:,2,3);
Skew_pspin = steroidreccorr(:,3,3);

Steroidreceptorresults = table(Genes, Gradient_corr, Gradient_p,Gradient_pspin, Mean_corr, Mean_p, Mean_pspin, Skew_corr, Skew_p, Skew_pspin);

Genes = labels_steroidsynth;
Gradient_corr = steroidsynthcorr(:,1,1);
Gradient_p = steroidsynthcorr(:,2,1);
Gradient_pspin = steroidsynthcorr(:,3,1);
Mean_corr = steroidsynthcorr(:,1,2);
Mean_p = steroidsynthcorr(:,2,2);
Mean_pspin = steroidsynthcorr(:,3,2);
Skew_corr = steroidsynthcorr(:,1,3);
Skew_p = steroidsynthcorr(:,2,3);
Skew_pspin = steroidsynthcorr(:,3,3);

Steroidsynthesisresults = table(Genes, Gradient_corr, Gradient_p,Gradient_pspin, Mean_corr, Mean_p, Mean_pspin, Skew_corr, Skew_p, Skew_pspin);

if old_donors == 1
    save(fullfile(outDir, 'geneticdecoding.mat'), 'Steroidreceptorresults', 'Steroidsynthesisresults', 'steroidreccorr', 'steroidsynthcorr', 'chromsYcorr', 'chromsXcorr', 'baselinecorr', 'coeff', 'PC');
elseif female_donor == 1
    save(fullfile(outDir, 'geneticdecoding_femdon.mat'), 'Steroidreceptorresults', 'Steroidsynthesisresults', 'steroidreccorr', 'steroidsynthcorr', 'chromsYcorr', 'chromsXcorr', 'baselinecorr', 'coeff', 'PC');
elseif male_donor == 1
    save(fullfile(outDir, 'geneticdecoding_maledon.mat'), 'Steroidreceptorresults', 'Steroidsynthesisresults', 'steroidreccorr', 'steroidsynthcorr', 'chromsYcorr', 'chromsXcorr', 'baselinecorr', 'coeff', 'PC');
elseif all_donors == 1
    save(fullfile(outDir, 'geneticdecoding_allldons.mat'), 'Steroidreceptorresults', 'Steroidsynthesisresults', 'steroidreccorr', 'steroidsynthcorr', 'chromsYcorr', 'chromsXcorr', 'baselinecorr', 'coeff', 'PC');
end




%% Plotting 

clear C yl xl pvals xpAll Ctmp Ctmp2 pvalstmp pvalstmp2 myLabeltmp myLabeltmp2

Ctmp(:,3) = Steroidsynthesisresults.Gradient_corr;
Ctmp(:,1) = Steroidsynthesisresults.Mean_corr;
Ctmp(:,2) = Steroidsynthesisresults.Skew_corr;
if bonf == 1
    pvalstmp(:,3) = Steroidsynthesisresults.Gradient_p;
    pvalstmp(:,1) = Steroidsynthesisresults.Mean_p;
    pvalstmp(:,2) = Steroidsynthesisresults.Skew_p;
elseif bonf == 0
    pvalstmp(:,3) = Steroidsynthesisresults.Gradient_pspin;
    pvalstmp(:,1) = Steroidsynthesisresults.Mean_pspin;
    pvalstmp(:,2) = Steroidsynthesisresults.Skew_pspin;
end
myLabeltmp = Steroidsynthesisresults.Genes;

Ctmp2(:,3) = Steroidreceptorresults.Gradient_corr;
Ctmp2(:,1) = Steroidreceptorresults.Mean_corr;
Ctmp2(:,2) = Steroidreceptorresults.Skew_corr;
if bonf == 1
    pvalstmp2(:,3) = Steroidreceptorresults.Gradient_p;
    pvalstmp2(:,1) = Steroidreceptorresults.Mean_p;
    pvalstmp2(:,2) = Steroidreceptorresults.Skew_p;
elseif bonf == 0
    pvalstmp2(:,3) = Steroidreceptorresults.Gradient_pspin;
    pvalstmp2(:,1) = Steroidreceptorresults.Mean_pspin;
    pvalstmp2(:,2) = Steroidreceptorresults.Skew_pspin;
end
myLabeltmp2 = Steroidreceptorresults.Genes; 

baseline_r(1) = baselinecorr(1,1);
baseline_r(2) = baselinecorr(1,2);
baseline_r(3) = baselinecorr(1,3);
baseline_p(1) = baselinecorr(2,1);
baseline_p(2) = baselinecorr(2,2);
baseline_p(3) = baselinecorr(2,3);

C = [Ctmp; Ctmp2; baseline_r];
pvals = [pvalstmp; pvalstmp2; baseline_p];
myLabel = [myLabeltmp; myLabeltmp2; 'baseline'];
no_all_tests = size(myLabel,1);

measureLabel = {'Profile Mean', 'Profile Skewness', 'Gradient'};

% Set [min,max] value of C to scale colors
clrLim = [-0.5,0.5];
% load('CorrColormap.mat') % Uncomment for custom CorrColormap
% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];


% Compute center of each circle
% This assumes the x and y values were not entered in imagesc()
x = 1 : 1 : size(C,1); % x edges
y = 1 : 1 : size(C,2); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(C==0)=nan; % eliminate cordinates for zero correlations

xpAll = xAll;
pvals = pvals';

if bonf == 1
    xpAll(pvals>0.05/no_all_tests)=nan; % eliminate cordinates for pspin bigger 0.05 
elseif bonf == 0
    xpAll(pvals==0)=nan; % eliminate cordinates for non-sig pspin vals
    xpAll(pvals>0.05)=nan; % eliminate cordinates for pspin bigger 0.05 
end

% Set color of each rectangle
% Set color scale
cmap = colormap(flipud(cbrewer('div','RdBu',11)));

% Create figure
fh = figure();
ax = axes(fh);
hold(ax,'on')
colormap(flipud(cbrewer('div','RdBu',11)))
tickvalues = 1:3;
x = zeros(1,3);
text(x, tickvalues, measureLabel, 'HorizontalAlignment', 'right','fontsize', 16, 'fontname', 'Calibri');
tickvalues = 1:length(C);
x = zeros(size(tickvalues));
x(:) = 4
text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',45, 'fontsize', 16, 'fontname', 'Calibri');


C = C';
%     scatter(xAll(:), yAll(:), 400.*abs(C(:)), C(:), 'filled', 'MarkerFaceAlpha', 1);
scatter(xAll(:), yAll(:), 800.*abs(C(:)), C(:), 'filled', 'MarkerFaceAlpha', 1);
scatter(xpAll(:), yAll(:), 800.*abs(C(:)), C(:), 'filled', 'MarkerFaceAlpha', 1,...
    'MarkerEdgeColor', [0 0 0], 'LineWidth', 4);

% enclose markers in a grid
n = size(C,1);
xl = xAll - 0.5;
%xadd = 26.5*ones(size(C,1),1);
xadd = (size(C,2) + 0.5)*ones(size(C,1),1);

xl = [xl,xadd];
xl = [xl;xl(1,:)];

yl = yAll + 0.5;
yadd = 0.5*ones(size(C,2),1)';
yl = [yadd;yl];
yl = [yl,yl(:,1)];

line(xl, yl, 'color', 'k') % horizontal lines
line(xl', yl', 'color', 'k') % vertical lines

axis(ax,'equal')
axis(ax,'tight')
set(ax,'YDir','Reverse', 'fontsize', 16, 'fontname', 'Calibri')
colorbar()
caxis(clrLim);
axis off
    
% black figure:
% Create figure
% fh = figure();
% ax = axes(fh);
% hold(ax,'on')
% colormap(flipud(cbrewer('div','RdBu',11)))
% tickvalues = 1:3;
% x = zeros(1,3);
% text(x, tickvalues, measureLabel, 'HorizontalAlignment', 'right','fontsize', 16, 'fontname', 'Calibri', 'color', 'w');
% tickvalues = 1:length(C);
% x = zeros(size(tickvalues));
% x(:) = 4
% text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',45, 'fontsize', 16, 'fontname', 'Calibri', 'color', 'w');
% 
% 
% set(gcf,'color','k'); % Set figure background to black
% 
% C = C';
% %     scatter(xAll(:), yAll(:), 400.*abs(C(:)), C(:), 'filled', 'MarkerFaceAlpha', 1);
% scatter(xAll(:), yAll(:), 800.*abs(C(:)), C(:), 'filled', 'MarkerFaceAlpha', 1);
% scatter(xpAll(:), yAll(:), 800.*abs(C(:)), C(:), 'filled', 'MarkerFaceAlpha', 1,...
%     'MarkerEdgeColor', [1 1 1], 'LineWidth', 4);
% 
% % enclose markers in a grid
% n = size(C,1);
% xl = xAll - 0.5;
% xadd = 26.5*ones(size(C,1),1);
% 
% xl = [xl,xadd];
% xl = [xl;xl(1,:)];
% 
% yl = yAll + 0.5;
% yadd = 0.5*ones(size(C,2),1)';
% yl = [yadd;yl];
% yl = [yl,yl(:,1)];
% 
% line(xl, yl, 'color', 'w') % horizontal lines
% line(xl', yl', 'color', 'w') % vertical lines
% 
% axis(ax,'equal')
% axis(ax,'tight')
% set(ax,'YDir','Reverse', 'fontsize', 16, 'fontname', 'Calibri', 'color', 'w')
% colorbar('color', 'w')
% caxis(clrLim);
% axis off

%% Save all results

disp(' ...saving results')
if saveall == 1
    save(fullfile(outDir, '/genetic_decoding_sexdiffs.mat'));
    %save((fullfile(outDir, sprintf('%s_sexdiffsmaps.mat', dataname))), 'T1T2moments', 'results', 'covariates', 'keep');    
    FigList = findobj(allchild(0), 'flat','Type','figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        savefig(FigHandle, fullfile(figDir, string(FigName(iFig))));
    end
end

disp(' ...done!:)')


