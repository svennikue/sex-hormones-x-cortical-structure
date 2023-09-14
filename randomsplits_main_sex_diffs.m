% this script loads the results from sex_diffs_3measures.m and creates random
% splits, weighted for men and women. 
% it does so 1000 times and outputs the correlations.

clear all

% make a path to directory
addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';
scriptDir = fullfile(homeDir,'script');
addpath(genpath(scriptDir));
dataDir = fullfile(homeDir, 'data');
addpath(genpath(dataDir));
outDir = fullfile(homeDir,'output');
addpath(genpath(outDir));
figDir = fullfile(homeDir, 'figures/HCP');
addpath(figDir);


%% loading data
% load('female_gradients_out.mat');
load HCP_T1wT2w_sexdiffsmaps.mat
results_HCP= results.tvals;

sexterm = term(covariates.sex);
sex = covariates.sex;
ageterm = term(covariates.age);
age = covariates.age;
icvterm = term(covariates.icv);
icv = covariates.icv;

keepmale = find(cell2mat(covariates.sex) == 'M');
keepfemale = find(cell2mat(covariates.sex) == 'F');

halfsplitsGrad = [];
halfsplitsMean = [];
halfsplitsSkew = [];
takeout20gradc = [];
takeout20meanc = [];
takeout20skewc = [];

% original t-values
    M = 1 + sexterm + ageterm + icvterm;
    slm_grad =  SurfStatLinMod(T1T2moments(:,:,1),M);
    slm_grad = SurfStatT(slm_grad,sexterm.F-sexterm.M);
    valsgrad = (slm_grad.t)';
    
    slm_mean =  SurfStatLinMod(T1T2moments(:,:,2),M);
    slm_mean = SurfStatT(slm_mean,sexterm.F-sexterm.M);
    valsmean = (slm_mean.t)';

    slm_skew =  SurfStatLinMod(T1T2moments(:,:,3),M);
    slm_skew = SurfStatT(slm_skew,sexterm.F-sexterm.M);
    valsskew = (slm_skew.t)';

for splits = 1:1000
    % do random split 
    indicesmale = randperm(size(keepmale,1));
    indicesfemale = randperm(size(keepfemale,1));
    
    set1male = keepmale(indicesmale(1:249));
    set2male = keepmale(indicesmale(250:499));
    
    set1female = keepfemale(indicesfemale(1:375));
    set2female = keepfemale(indicesfemale(298:594));
    
    set1 = [set1male; set1female];
    set2 = [set2male; set2female];
   
    % set up linear model with covariates for mean T1T2w layers
    sexterm1 = term(sex(set1));
    ageterm1 = term(age(set1));
    icvterm1 = term(icv(set1));
       
    sexterm2 = term(sex(set2));
    ageterm2 = term(age(set2));
    icvterm2 = term(icv(set2));

    % set up linear model: controlling for age, meanT1T2, ICV
    % females - males

    M1 = 1 + sexterm1 + ageterm1 + icvterm1;
    slmgrad1 =  SurfStatLinMod(T1T2moments(set1,:,1),M1);
    slmgrad1 = SurfStatT(slmgrad1,sexterm1.F-sexterm1.M); %contrast is female - male

    slmmean1 =  SurfStatLinMod((T1T2moments(set1,:,2)),M1);
    slmmean1 = SurfStatT(slmmean1,sexterm1.F-sexterm1.M); %contrast is female - male
    
    slmskew1 =  SurfStatLinMod((T1T2moments(set1,:,3)),M1);
    slmskew1 = SurfStatT(slmskew1,sexterm1.F-sexterm1.M); %contrast is female - male
    
    M2 = 1 + sexterm2 + ageterm2 + icvterm2;
    slmgrad2 = SurfStatLinMod((T1T2moments(set2,:,1)),M2);
    slmgrad2 = SurfStatT(slmgrad2,sexterm2.F-sexterm2.M); %contrast is female - male
    
    slmmean2 = SurfStatLinMod((T1T2moments(set2,:,2)),M2);
    slmmean2 = SurfStatT(slmmean2,sexterm2.F-sexterm2.M); %contrast is female - male
    
    slmskew2 = SurfStatLinMod((T1T2moments(set2,:,3)),M2);
    slmskew2 = SurfStatT(slmskew2,sexterm2.F-sexterm2.M); %contrast is female - male
    
    valsgrad1 = (slmgrad1.t)';
    valsgrad2 = (slmgrad2.t)';
    
    valsmean1 = (slmmean1.t)';
    valsmean2 = (slmmean2.t)';
    
    valsskew1 = (slmskew1.t)';
    valsskew2 = (slmskew2.t)';
    
    halfsplitrepg = corr(valsgrad1, valsgrad2);   
    halfsplitsGrad = [halfsplitsGrad, halfsplitrepg];
    
    halfsplitrepm = corr(valsmean1, valsmean2);   
    halfsplitsMean = [halfsplitsMean, halfsplitrepm];
    
    halfsplitreps = corr(valsskew1, valsskew2);   
    halfsplitsSkew = [halfsplitsSkew, halfsplitreps];

end

%% plot the result.
r1 = 1 - rand(500,1)*0.1;
r2 = 1 + rand(500,1)*0.1;
r = [r1;r2];

colors = colormap(flipud(cbrewer('seq','Greens',11)));

f=figure,
hold on
subplot(1,3,1)
hold on 
scatter(r, halfsplitsMean, 'MarkerEdgeColor', [0.4039         0    0.1216], 'Marker', '.');
% mean = 0.8622, range = 0.4871, std = 0.0563
boxplot(halfsplitsMean, 'widths', 0.5, 'Colors', [0.4039         0    0.1216], 'Symbol', 'w');
ylabel(['Spatial Correlation'])
xlabel(['Profile Mean'])
set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])

subplot(1,3,2)
hold on
scatter(r, halfsplitsSkew, 'MarkerEdgeColor', [0.0196    0.1882    0.3804], 'Marker', '.');
boxplot(halfsplitsSkew, 'widths', 0.5, 'Colors', [0.0196    0.1882    0.3804], 'Symbol', 'w');
ylabel(['Spatial Correlation'])
xlabel(['Profile Skewness'])
set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])

subplot(1,3,3)
hold on 
scatter(r, halfsplitsGrad, 'MarkerEdgeColor', [0    0.2667    0.1059], 'Marker', '.');
% mean = 0.5506, range = 0.3588, std = 0.0408
boxplot(halfsplitsGrad, 'widths', 0.5, 'Colors', [0    0.2667    0.1059], 'Symbol', 'w');
ylabel(['Spatial Correlation'])
xlabel(['Gradient'])
set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])


%% Save all results
disp(' ...saving results')

save(fullfile(outDir, 'randomsplits_HCP.mat'));
FigList = findobj(allchild(0), 'flat','Type','figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = sprintf('Figure%02d.fig', iFig);
    savefig(FigHandle, fullfile(figDir, FigName));
end


disp(' ...done!:)')

