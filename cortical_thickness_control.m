% Cortical Thickness analysis
% Supplement 4

close all
clear all

saveall = 1;
count = 0;
takesinglemeas = 1; %no repeated measures bc it doesnt work
dataname = 'HCP_T1wT2w';
examine = 0;

% make a path to directory
addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';
scriptDir = fullfile(homeDir,'script');
addpath(genpath(scriptDir));
dataDir = fullfile(homeDir, 'data');
addpath(genpath(dataDir));
outDir = fullfile(homeDir,'output');
addpath(genpath(outDir));
figDir = fullfile(homeDir, 'figures/HCP/main');
addpath(figDir);

%%
load('HCP_T1wT2w_sexdiffsmaps.mat');
load('female_gradients.mat'); %1206 subjects
SN = female_gradients.SN;
D1 = female_gradients.D1;                             
unres = female_gradients.unres;                        
load('HCP400_CT.mat')
MPCnx1 = female_gradients.MPCnx1; %basis for gradient   
keep = find(squeeze(mean(MPCnx1(:,1,1:400),3))>0);
genderkeep = cellstr(unres.Gender(keep));
schaefer_400 = (fetch_parcellation('fsaverage5', 'schaefer', 400))';
 
genderkeep = cellstr(unres.Gender(keep));
agekeep = D1.Age_in_Yrs(keep);
icvkeep = unres.FS_IntraCranial_Vol(keep);
G = term(genderkeep);
M      =    term(agekeep) + term(genderkeep) + term(icvkeep);
slm = SurfStatLinMod(HCP400_CT(keep,:),M)
slm = SurfStatT(slm,G.F-G.M)

ct_sex = slm.t;

% FDR correction females - males
p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
% > if p > .95, very 
h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
h= h+hn; % complete mask

% compute effect sizes
Cohensd_ct = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
CohensdFDR = (Cohensd_ct.*h)'; % mask only those that were FDR corrected 

heri_ct = zeros(1,20484);

heri_ct = zeros(1,20484);
for i = 1:400 %200 parcels per hemisphere
    heri_ct(find(schaefer_400==i)) = Cohensd_ct(i);
end
f = figure,
BoSurfStatViewData(heri_ct,SN,'')
colormap(flipud(cbrewer('div','RdBu',11)))
SurfStatColLim([-1 1])

heri_ct = zeros(1,20484);
for i = 1:400 %200 parcels per hemisphere
    heri_ct(find(schaefer_400==i)) = CohensdFDR(i);
end
f = figure,
BoSurfStatViewData(heri_ct,SN,'')
colormap(flipud(cbrewer('div','RdBu',11)))
SurfStatColLim([-1 1])


% finally, check if cortical thickness relates to any of the 3 measures.
for m = 2
    if m == 1                
    namemoment = 'gradient'
    else if m == 2
        namemoment = 'profile mean'
        else if m == 3
            namemoment = 'profile skewness'
        end
    end
    end
      f=figure,
      scatter(Cohensd_ct,results.Cohensd(:,m),'k','filled'),lsline, set (gca, 'fontsize', 40, 'fontname', 'Calibri')
      [rho,pval] = corr(results.Cohensd(:,m), Cohensd_ct', 'Type', 'Spearman', 'rows', 'complete');
    
       [p_spin, r_dist] = spin_test(results.Cohensd(:,m)', Cohensd_ct, 'surface_name',...
            'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
            'type', 'spearman');
    
        f = figure,
        his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
             'facealpha', 1, 'linewidth', 0.5);
        l = line(rho, 3,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
        set (gca, 'fontsize', 40, 'fontname', 'Calibri')
        xlabel(['Null correlations' newline (sprintf('(%s)', namemoment))])
        xlim([-0.5, 0.5]);
        legend(l,['{\it r}=' num2str(round(rho, 2)) newline ...
                          '{\it p}=' num2str(round(p_spin,3 ))], 'Location', 'northwest')
end

