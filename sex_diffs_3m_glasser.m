% This script generates the main effects for the sex differences in
% microstructure paper, but for the glasser parcellation.

% input:subject + variable lists, layer profiles, individual gradients, 
%       MPCs/ mean MPC, HCP data
% output: sex differences in the 3 microstructural measures, figures
% functions:exportfigbo, BoSurfStatViewData, cbrewer, Surfstat, brainstat
% Supplement 1

loadold = 0;
if loadold == 0
    clear all
    close all
    loadold = 0;
end


saveall = 1;
count = 0;
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



%% loading data
if loadold == 0
    dataname = 'HCP_glasser';

    load('female_gradients.mat'); %1206 subjects
    D1 = female_gradients.D1;                                  
    unres = female_gradients.unres;                                        
    SN = female_gradients.SN;
    
    load('MPC_layers_retest2.mat')
    layers = MPC_layers_retest;

    % load glasser parcellation
    load('parcels360.mat')
    
    % gradient for mean matrices
    load('MPC_all_retest2.mat')
    MPCnx1 = MPC_all_retest(:, [2:181, 183:362], [2:181, 183:362]);
    
    MPC_schaffer = female_gradients.MPCnx1; %basis for gradient                   
    %make keep to keep those with complete data
    keep_schaffer = find(squeeze(mean(MPC_schaffer(:,1,1:400),3))>0);

    % CAREFUL! Run create_gradients_HCP.m before
    load('HCP_glasser_gradient.mat')
    c1_tx = cell2mat(gradient.grad1persubj);
    
    % first compute the measure before getting keep.
    T1T2moments(:,:,1) = c1_tx; % first moment is gradient
    
    %clear T1T2moments
    for i = 1:size(c1_tx,1)
        % subjects x parcels x moments
        T1T2moments(i,:,2) = mean(squeeze(layers(i,2:11,[2:181, 183:362]))); % second moment is mean
        T1T2moments(i,:,3) = skewness(squeeze(layers(i,2:11,[2:181, 183:362]))); % third moment is skewness
        %T1T2moments(i,:,3) = rescale(T1T2moments(i,:,3));
        % calculate mean of 10 MPC layers for each subject (for cov)
        meant1t2(i,:) = mean(squeeze(layers(i,2:11,[2:181, 183:362])));
    end

    %make keep to keep those with complete data
    keep_grad = find(mean(T1T2moments(:,:,1),2) ~= 0 & ~any(isnan(T1T2moments(:,:,1)), 2));
    keep_mean = find(mean(T1T2moments(:,:,2),2) ~= 0 & ~any(isnan(T1T2moments(:,:,2)), 2));
    keep_skew = find(mean(T1T2moments(:,:,3),2) ~= 0 & ~any(isnan(T1T2moments(:,:,3)), 2));
    %keep = intersect(keep_mean, keep_skew)
    keep = intersect(keep_schaffer, intersect(keep_grad, intersect(keep_mean, keep_skew)));
    
    T1T2moments = T1T2moments(keep,:,:);
    meant1t2= meant1t2(keep,:);

    % I unfortunately don't have all euler values. Thus, interpolate
    load('eulerHCPs1200.mat')
    euler_compl = fillmissing(eulerHCP.euler_HCP, 'movmean', 100);
    euler_no = euler_compl(keep)';
   

    %layers = female_gradients.MPC_layers1(keep,:,:); %T1T2w profiles
    ICV = unres.FS_IntraCranial_Vol(keep,:);
    ICV = ICV/1000000; 
    sex = cellstr(unres.Gender(keep,:));
    agekeep = D1.Age_in_Yrs(keep);
    keepmale = find(cell2mat(sex) == 'M');
    keepfemale = find(cell2mat(sex) == 'F');

    %make mean matrixes
    %meanMPC = squeeze(mean(MPCnx1(keep,:,:)));
    meanMPC = squeeze(mean(MPCnx1(keep,:,:)));
    meanMPC(eye(size(meanMPC))==1) = 0;



else if loadold == 1        
    % alternatively: load complete mat from running this script last time. 
    load('main_effects_sexdiffs.mat');
    keyboard
    end
end


%% Plot method.

f = figure,
BoSurfStatViewData(parcels360x,SN,'')
colormap((cbrewer('qual','Set3',11)))
count = count + 1
FigName{(count)} = 'glasser_parcellation.png';


f=figure,
%colormap(bone)
colormap(cbrewer('seq','Greens',11))
imagesc(meanMPC, [0,1])
%imagesc(meanMPC, [-0.5,3])
colorbar
set (gca, 'fontsize', 35, 'fontname', 'Calibri')
exportfigbo(f,[figDir sprintf('%s_meanMPC.png', dataname)],'png', 10)
count = count + 1
FigName{(count)} = sprintf('%s_meanMPC.fig', dataname);


% plot a sorted MPC for one subject
[B, I] = sort(mean(c1_tx));
for i = 1:size(c1_tx,2)
    x = I(i);
    for j = 1:size(c1_tx,2)
        y = I(j);
        sorted_M(i,j) = MPCnx1(size(MPCnx1, 1),x,y);
    end
end

f = figure 
%colormap(bone)
colormap(cbrewer('seq','Greens',11))
%imagesc(sorted_M, [-0.5,3])
imagesc(sorted_M, [0,1])
colorbar
set (gca, 'fontsize', 35, 'fontname', 'Calibri')
exportfigbo(f,[figDir, sprintf('%s_orderedMPC.png', dataname)],'png', 10)
count = count + 1
FigName{(count)} = sprintf('%s_orderedMPC.fig', dataname);


% check order: the first row should be white, the last one should be black.
f=figure,
imagesc(squeeze(mean(mean(layers,3))'))
%colormap((cbrewer('seq','Greys',99)))
colormap(cbrewer('seq','Greens',11))
colorbar
exportfigbo(f,[figDir, sprintf('%s_avg_intensity_acrosscortex.png', dataname)],'png', 10) 
count = count + 1
FigName{(count)} = sprintf('%s_avg_intensity_acrosscortex.fig', dataname);


%% ANALYSIS: Gradient, Mean Layer Profile + Skewness Layer Profile

% set up all covariates
meant1t2persub = mean(meant1t2,2); % mean over the whole cortex
meant1t2cov = term(meant1t2persub);  
icvterm = term(ICV);
sexterm = term(sex);
ageterm = term(agekeep);
eulerterm = term(euler_no);
%icvterm = term(unres.FS_IntraCranial_Vol(keep));

% save the covs 
covariates.meant1t2 = meant1t2persub;
covariates.icv = ICV;
covariates.sex = sex;
covariates.age = agekeep;
covariates.euler = euler_no;


%% Sex difference Analysis for all 3 measures.
for m = 1:size(T1T2moments, 3)
    %for m = 3
    if m == 1
        namemoment = 'gradient';
    else if m == 2
            namemoment = 'profile_mean';
        else if m == 3
                namemoment = 'profile_skewness';
            end
        end
    end

    to_surface = zeros(20484,1);
    % 2:181 + 183:362
    for i = 1:180 % to exclude medial wall
      to_surface(parcels360x==i+1) = mean(T1T2moments(:,i,m));
    end
    for i = 1:180 % to exclude medial wall
      to_surface(parcels360x==i+182) = mean(T1T2moments(:,i+180,m));
    end
    
    f = figure,
    BoSurfStatViewData(to_surface,SN,'')
    %colormap((cbrewer('seq','Greys',11)))
    colormap((cbrewer('seq','Greens',11)))
    clim = [prctile(mean(T1T2moments(:,:,m)),5) prctile(mean(T1T2moments(:,:,m)),95)]; % set colour limits
    SurfStatColLim(clim)
    exportfigbo(f,[figDir, sprintf('/%s_groupavg_%s', dataname, namemoment)],'png', 10)
    count = count + 1
    FigName{(count)} = sprintf('%s_groupavg_%s.fig', dataname, namemoment);


    % visualize how the measure looks like per layer.
    % show one high example for each measure, and one low example:
    % create group avarage of mean moment for all parcels,
    % then plot the profile layer-wise for the parcel of smallest and
    % biggest moment

    minlayers = 1;
    maxlayers = 3;

    minmean = find(mean(T1T2moments(:,:,m))==min(mean(T1T2moments(:,:,m))));
    maxmean = find(mean(T1T2moments(:,:,m))==max(mean(T1T2moments(:,:,m))));
    f=figure,
    hold on,
    subplot(1,4,1)
    %1. subplot shows layer-wise profile of parcel with smallest moment
    %imagesc(mean(squeeze(layers(:,:,minmean)))',[1.2 2.8])
    imagesc(mean(squeeze(layers(:,:,minmean)))')
    %colormap((cbrewer('seq','Greys',99)))
    colormap(cbrewer('seq','Greens',11))    
    set (gca, 'fontsize', 35, 'fontname', 'Calibri')
    caxis([minlayers maxlayers])
    yticks([])

    colorbar
    subplot(1,4,2)
    %2. subplot shows T1T2w intensity per layer
    plot(mean(squeeze(layers(:,:,minmean))),1:12, 'k', 'Linewidth', 4)
    set (gca, 'ydir', 'reverse', 'fontsize', 35, 'fontname', 'Calibri')
    axis([minlayers maxlayers 1 12])
    yticks([])

    subplot(1,4,3)
    %3. subplot shows layer-wise profile of parcel with biggest moment 
    %imagesc(mean(squeeze(layers(:,:,maxmean)))',[1.2 2.8])
    imagesc(mean(squeeze(layers(:,:,maxmean)))')
    colorbar
    set (gca, 'fontsize', 35, 'fontname', 'Calibri')
    caxis([minlayers maxlayers])
    yticks([])

    subplot(1,4,4)
    %4. subplot shows T1T2w intensity per layer
    plot(mean(squeeze(layers(:,:,maxmean))),1:12, 'k', 'Linewidth', 4)
    set (gca, 'ydir', 'reverse', 'fontsize', 35, 'fontname', 'Calibri')
    axis([minlayers maxlayers 1 12])
    yticks([])

    exportfigbo(f,[figDir, sprintf('/%s_method_%s', dataname, namemoment)],'png', 10)
    count = count + 1
    FigName{(count)} = sprintf('%s_method_%s.fig', dataname, namemoment);

        
    % set up linear model with covariates of sex, age and ICV 
    M = 1 + sexterm + ageterm + icvterm + eulerterm;        
    slm = SurfStatLinMod(T1T2moments(:,:,m),M);
    slm = SurfStatT(slm,sexterm.F-sexterm.M); %contrast is female - male

    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected
    
    % plot sex difference EFFECT SIZE FDR corrected
    for i = 1:180 % to exclude medial wall
      to_surface(parcels360x==i+1) =  Cohensd(i)*h(i);
    end
    for i = 1:180 % to exclude medial wall
      to_surface(parcels360x==i+182) =  Cohensd(i+180)*h(i+180);
    end

    
    f = figure,
    BoSurfStatViewData(to_surface,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    %clim = [prctile(mean(slm.t),5) prctile(mean(slm.t),95)]; % set colour limits
    SurfStatColLim([-1 1])
    %SurfStatColLim(clim)
    exportfigbo(f,[figDir, sprintf('/%s_fem-mal_FDRCohensD_%s', dataname, namemoment)],'png', 10)
    count = count + 1
    FigName{(count)} = sprintf('%s_fem-mal_FDRCohensD_%s.fig', dataname, namemoment);
    
    % identify mean of all positive and negative effects
    bigwomen = find(CohensdFDR>0);
    bigwomeneffect = mean(CohensdFDR(bigwomen)); %0.2603
    bigmen = find(CohensdFDR<0);
    bigmeneffect = mean(CohensdFDR(bigmen)); %-0.3154
    
    % store results
    results.tvals(:,m) = (slm.t)';
    results.FDR(:,m) = (slm.t.*h)';
    results.namemoment{m} = namemoment;
    results.Cohensd(:,m) = Cohensd';
    results.CohensdFDR(:,m) = CohensdFDR;
    results.bigwomeneffect(:,m) = bigwomeneffect;
    results.bigmeneffect(:,m) = bigmeneffect;

    descriptives.mean(1,m) = mean(mean(T1T2moments(keepfemale,:,m)));
    descriptives.mean(2,m) = mean(mean(T1T2moments(keepmale,:,m)));
    descriptives.std(1,m,:) = mean(std(T1T2moments(keepfemale,:,m)));
    descriptives.std(2,m,:) = mean(std(T1T2moments(keepmale,:,m)));
    descriptives.posmean(1,m) = mean(CohensdFDR(bigwomen));
    descriptives.negmean(2,m) = mean(CohensdFDR(bigmen));
    descriptives.posstd(1,m) = std(CohensdFDR(bigwomen));
    descriptives.negstd(2,m) = std(CohensdFDR(bigmen));
    
end


%% Add analysis: correlation between Mean and Skewness

% correlation between skewsness sex-diffs and mean sex-diffs: r = -.4120
mean_sex_diffs = results.Cohensd(:,2);
skew_sex_diffs = results.Cohensd(:,3);
corr_mean_skew_sexdiffs = corr(mean_sex_diffs, skew_sex_diffs)

% parcel-wise correlation between skew and mean, all males
mean_males = T1T2moments(keepmale, :, 2);
skew_males = T1T2moments(keepmale, :, 3);
corr_males_mean_skew_avg = corrcoef(mean_males, skew_males)
corr_males_mean_skew_mtx = corr(mean_males, skew_males);
corr_males_mean_skew = diag(corr_males_mean_skew_mtx);

for i = 1:180 % to exclude medial wall
  to_surface(parcels360x==i+1) =  corr_males_mean_skew(i);
end
for i = 1:180 % to exclude medial wall
  to_surface(parcels360x==i+182) =  corr_males_mean_skew(i+180);
end


f = figure,
BoSurfStatViewData(to_surface,SN,'')
colormap(flipud(cbrewer('div','RdBu',11)))
%clim = [prctile(mean(slm.t),5) prctile(mean(slm.t),95)]; % set colour limits
SurfStatColLim([-1 1])
%SurfStatColLim(clim)
exportfigbo(f,[figDir, sprintf('/%s_corr_males_mean_skew%s')],'png', 10)


% parcel-wise correlation between skew and mean, all females
mean_females = T1T2moments(keepfemale, :, 2);
skew_females = T1T2moments(keepfemale, :, 3);
corr_females_mean_skew_avg = corrcoef(mean_females, skew_females)
corr_females_mean_skew_mtx = corr(mean_females, skew_females);
corr_females_mean_skew = diag(corr_females_mean_skew_mtx);


for i = 1:180 % to exclude medial wall
  to_surface(parcels360x==i+1) =  corr_females_mean_skew(i);
end
for i = 1:180 % to exclude medial wall
  to_surface(parcels360x==i+182) =  corr_females_mean_skew(i+180);
end


f = figure,
BoSurfStatViewData(to_surface,SN,'')
colormap(flipud(cbrewer('div','RdBu',11)))
%clim = [prctile(mean(slm.t),5) prctile(mean(slm.t),95)]; % set colour limits
SurfStatColLim([-1 1])
%SurfStatColLim(clim)
exportfigbo(f,[figDir, sprintf('/%s_corr_females_mean_skew%s')],'png', 10)


% parcel-wise correlation between skew and mean, all 
mean_all = T1T2moments(:, :, 2);
skew_all = T1T2moments(:, :, 3);
corr_all_mean_skew_avg = corrcoef(mean_all, skew_all)
corr_all_mean_skew_mtx = corr(mean_all, skew_all);
corr_all_mean_skew = diag(corr_all_mean_skew_mtx);


for i = 1:180 % to exclude medial wall
  to_surface(parcels360x==i+1) =  corr_all_mean_skew(i);
end
for i = 1:180 % to exclude medial wall
  to_surface(parcels360x==i+182) =  corr_all_mean_skew(i+180);
end


f = figure,
BoSurfStatViewData(to_surface,SN,'')
colormap(flipud(cbrewer('div','RdBu',11)))
%clim = [prctile(mean(slm.t),5) prctile(mean(slm.t),95)]; % set colour limits
SurfStatColLim([-1 1])
%SurfStatColLim(clim)
exportfigbo(f,[figDir, sprintf('/%s_corr_all_mean_skew%s')],'png', 10)



%% Save all results


if saveall == 1
    disp(' ...saving results')
    save(fullfile(outDir, sprintf('%s_main_effects_sexdiffs_glasser.mat', dataname)));
    save((fullfile(outDir, sprintf('%s_sexdiffsmaps_glasser.mat', dataname))), 'T1T2moments', 'results', 'covariates', 'keep');    
    FigList = findobj(allchild(0), 'flat','Type','figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        savefig(FigHandle, fullfile(figDir, string(FigName(iFig))));
    end
end

disp(' ...done!:)')
