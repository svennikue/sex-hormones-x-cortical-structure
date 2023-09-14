% random splits hippocampus

% this script loads the results from hippocampus_r/lh_sexdiffs.m and creates random
% splits, weighted for men and women. 
% it does so 1000 times and outputs the correlations.

% careful: one has to adjust the script in case the take out 20% is wished
% for. In the default, only the first figure ist usable.

clear all

load_old = 0;

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

hippoDir_RH = fullfile(dataDir, 'hipp_unfold_T1T2w/right_hippo');
hippoDir_LH = fullfile(dataDir, 'hipp_unfold_T1T2w');

addpath(genpath(hippoDir_RH));
addpath(genpath(hippoDir_LH));

%% loading data
files_RH = dir(fullfile(hippoDir_RH, '*.txt'));
files_LH = dir(fullfile(hippoDir_LH, '*.txt'));
hippo_LH = []; 
hippo_RH = []; 

for i=1:length(files_LH)
    filename = files_LH(i).name;
    hippo_LH = [hippo_LH, importdata(filename)]; % 7262 hippocampal vertices, 867 subjects
    filename = files_RH(i).name;
    hippo_RH = [hippo_RH, importdata(filename)]; % 7262 hippocampal vertices, 867 subjects
end

subjects = importdata(fullfile(dataDir, 'subjectListS1200_hipp_unfold.txt')); %867 subjects
load subjects_HCP.mat; % subjects, sex, unres.FS_IntraCranial_Vol 
[C, keep, IB] = intersect(unres.Subject, subjects); %unres is subjectlist
subjectlist = unres(keep,:); 
sex = cellstr(subjectlist.Gender);
icv =  subjectlist.FS_IntraCranial_Vol;

load subjects_HCP_ageBMI.mat; % BMI, age, handedness, menstrual stuff,...%D1
[C, keep, IB] = intersect(D1.Subject, subjects); %D1 are subjectvariables
subjectvars = D1(keep,:);
age = subjectvars.Age_in_Yrs;

keepmale = find(cell2mat(sex) == 'M');
keepfemale = find(cell2mat(sex) == 'F');


%%
% set up models
halfsplitsRH = [];
halfsplitsLH = [];

for splits = 1:1000
    % do random split 
    indicesmale = randperm(size(keepmale,1));
    indicesfemale = randperm(size(keepfemale,1));
    
    set1male = keepmale(indicesmale(1:183));
    set2male = keepmale(indicesmale(184:367));
    
    set1female = keepfemale(indicesfemale(1:250));
    set2female = keepfemale(indicesfemale(251:500));
    
    set1 = [set1male; set1female];
    set2 = [set2male; set2female];
   
    % set up linear model with covariates for mean T1T2w layers
    % include meanT1T2 per subject as covariate
    sexterm1 = term(sex(set1));
    ageterm1 = term(age(set1));
    icvterm1 = term(icv(set1));
       
    sexterm2 = term(sex(set2));
    ageterm2 = term(age(set2));
    icvterm2 = term(icv(set2));

    % set up linear model: controlling for age, meanT1T2, ICV
    % females - males
    M1 = 1 + sexterm1 + ageterm1 + icvterm1;

    slm_hipp_rh1 = SurfStatLinMod(hippo_RH(:,set1)',M1);
    slm_hipp_rh1 = SurfStatT(slm_hipp_rh1,sexterm1.F-sexterm1.M); %contrast is female - male
    
    slm_hipp_lh1 = SurfStatLinMod(hippo_LH(:,set1)',M1);
    slm_hipp_lh1 = SurfStatT(slm_hipp_lh1,sexterm1.F-sexterm1.M); %contrast is female - male
    
    M2 = 1 + sexterm2 + ageterm2 + icvterm2;
    slm_hipp_rh2 = SurfStatLinMod(hippo_RH(:,set2)',M2);
    slm_hipp_rh2 = SurfStatT(slm_hipp_rh2,sexterm2.F-sexterm2.M); %contrast is female - male
    slm_hipp_lh2 = SurfStatLinMod(hippo_LH(:,set2)',M2);
    slm_hipp_lh2 = SurfStatT(slm_hipp_lh2,sexterm2.F-sexterm2.M); %contrast is female - male
    
    hippo_lh1 = (slm_hipp_lh1.t);
    hippo_rh1 = (slm_hipp_rh1.t);
    hippo_lh2 = (slm_hipp_lh2.t);
    hippo_rh2 = (slm_hipp_rh2.t);
    
    haltsplitreg_lh = corr(hippo_lh1', hippo_lh2');
    haltsplitreg_rh = corr(hippo_rh1', hippo_rh2');

    halfsplitsRH = [halfsplitsRH,haltsplitreg_rh];
    halfsplitsLH = [halfsplitsLH,haltsplitreg_lh];
      
end


% to show the percentiles
CI_lh = prctile(halfsplitsLH, [5,95])
mean_lh = mean(halfsplitsLH)

CI_rh = prctile(halfsplitsRH, [5,95])
mean_rh = mean(halfsplitsRH)

%% plot results
r1 = 1 - rand(500,1)*0.1;
r2 = 1 + rand(500,1)*0.1;
r = [r1;r2];

colors = colormap(flipud(cbrewer('seq','Greens',11)));

f=figure,
hold on
subplot(1,2,1)
hold on 
scatter(r, halfsplitsLH, 'MarkerEdgeColor', [0.4039         0    0.1216], 'Marker', '.');
% mean = 0.8622, range = 0.4871, std = 0.0563
boxplot(halfsplitsLH, 'widths', 0.5, 'Colors', [0.4039         0    0.1216], 'Symbol', 'w');
ylabel(['Spatial Correlation'])
xlabel(['Hippocampus, LH'])
set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])

subplot(1,2,2)
hold on
scatter(r, halfsplitsRH, 'MarkerEdgeColor', [0.0196    0.1882    0.3804], 'Marker', '.');
boxplot(halfsplitsRH, 'widths', 0.5, 'Colors', [0.0196    0.1882    0.3804], 'Symbol', 'w');
ylabel(['Spatial Correlation'])
xlabel(['Hippocampus, RH'])
set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])



%% Save all results
disp(' ...saving results')

save(fullfile(outDir, 'randomsplits_hippo.mat'));
FigList = findobj(allchild(0), 'flat','Type','figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = sprintf('Figure%02d.fig', iFig);
    savefig(FigHandle, fullfile(figDir, FigName));
end


disp(' ...done!:)')
