% This script analysis if microstructural moments changes as a function of
% hormones. 
% This script analysis if microstructural moments changes as a function of
% hormones. 
% it takes the 5 models included in the manuscript and then runs a
% reliabilty anlyses with it: does the effect depend on something in the
% male sample?
% this scripts requires to have run sex_diffs_3measures.m before
% last edited 2024, written by Svenja Kuchenhoff.
% % Script for Supplement figure 8


loadold = 0;

if loadold == 0
    clear all
    close all
    loadold = 0;
end

saveall = 1;
plotfigs = 1;


% make a path to directory
addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';
scriptDir = fullfile(homeDir,'script');
addpath(genpath(scriptDir));
dataDir = fullfile(homeDir, 'data');
addpath(genpath(dataDir));
outDir = fullfile(homeDir,'output');
addpath(genpath(outDir));
figDir = fullfile(homeDir, 'figures/hormonal_figs');
addpath(figDir);

count = 0;


%% load data
% this scripts requires to have run sex_diffs_3measures.m before
load sexdiffsmaps.mat

% load schaefer400 parcellation
% careful! Now all parcels400 + 1 is schaefer_400
schaefer_400 = (fetch_parcellation('fsaverage5', 'schaefer', 400))';

% load sofies structure
load('female_gradients.mat'); %1206 subjects
D1 = female_gradients.D1;                             
unres = female_gradients.unres;                                            
SN = female_gradients.SN;
sex = cellstr(unres.Gender(keep,:));
agekeep = D1.Age_in_Yrs(keep);
keepmale = find(cell2mat(sex) == 'M'); %499 male
keepfemale = find(cell2mat(sex) == 'F'); %594 female

% create keep vectors: 
OC = D1.Menstrual_UsingBirthControl(keep); 
OC_estr_keep = find(D1.Menstrual_UsingBirthControl == 3);
OC_prog_keep = find(D1.Menstrual_UsingBirthControl == 4);

% all subjects for which OC is recorded
keep_OCrecorded = find(~isnan(OC)); % for 593 women, OC is recorded
% exclude all OC taking women from whole sample
excl_OC = find(OC ~= 1); % all men and all women that are NC
% include only women that take OC
keep_takeOC = find(OC == 1); %170 take OC
% exclude all naturally cycling women from whole sample
excl_natcyc = find(OC ~= 0); % all men and all women that take OC
% include only women that are naturally cycling
keep_NC = find(OC == 0); % n = 423

% prepare OC covariate  
OCcell = num2cell(OC(keep_OCrecorded));
OCcell(OC(keep_OCrecorded) == 1) = {'yes'};
OCcell(OC(keep_OCrecorded) == 0) = {'no'};

% prepare menstrual cycle measure
% D1(:,40) is Menstrual Cycle Length: 2= between 25-34 days
% D1(:,41) is Days since last Menstruation: all that are available
menses = D1.Menstrual_DaysSinceLast(keep);
subjectvars = D1(keep,:);


% exclude all that take birth control, that dont have a regular cycle,
% and that are more than 28 days since menstruation.

% women with a regular cycle shorter than 29 days = 284
% high estr = 184
% low estr = 100
% high prog = 113
% low prog = 171
% OC = 170
% women with regular cycles but not OC taking
keep_menses = mintersect(find(~isnan(menses)), find(menses<29), find(menses>-1), ...
    find(subjectvars.Menstrual_RegCycles==1), find(subjectvars.Menstrual_UsingBirthControl==0));

% overall 431 women with regular cycles or OC taking
keep_menses_andOC = mintersect(find(~isnan(menses)), find(menses<29), find(menses>-1), ...
    find(subjectvars.Menstrual_RegCycles==1), find(~isnan(OC)));

% prepare OC covariate
OCcov = OC(keep_menses_andOC);    

dayssincemenses = menses(keep_menses);
progconce = num2cell(dayssincemenses);
progconce(dayssincemenses < 15) = {'low'};
progconce(dayssincemenses > 14) = {'high'};

% prepare OC vs prog women covariate
progcov = menses(keep_menses_andOC);    
OC_vs_prog = num2cell(progcov);
OC_vs_prog(progcov < 15) = {'low'};
OC_vs_prog(progcov > 14) = {'high'};
OC_vs_prog(OCcov == 1) = {'OC'};


keeplowprog = mintersect(find(~isnan(menses)), find(menses<15), ...
    find(menses>-1), find(subjectvars.Menstrual_RegCycles==1), ...
    find(subjectvars.Menstrual_UsingBirthControl==0));
keepmenandlowprog = [keepmale; keeplowprog];

keephighprog = mintersect(find(~isnan(menses)), find(menses<29), ...
    find(menses>14), find(subjectvars.Menstrual_RegCycles==1), ...
    find(subjectvars.Menstrual_UsingBirthControl==0));
keepmenandhighprog = [keepmale; keephighprog];

%Definition phases: ovulation = days 7-14, premenstrual = 24-28
% luteal = after day 14, follicular = before day 15

estrconce = num2cell(dayssincemenses);
estrconce(dayssincemenses < 7) = {'low'};
estrconce(dayssincemenses > 6) = {'high'};
estrconce(dayssincemenses > 23) = {'low'};

% prepare OC vs estr women covariate
estrcov = menses(keep_menses_andOC);    
OC_vs_estr = num2cell(estrcov);
OC_vs_estr(progcov < 7) = {'low'};
OC_vs_estr(progcov > 6) = {'high'};
OC_vs_estr(progcov > 23) = {'low'};
OC_vs_estr(OCcov == 1) = {'OC'};


keephighestr = mintersect(find(~isnan(menses)), find(menses<23), ...
    find(menses>6), find(subjectvars.Menstrual_RegCycles==1), ...
    find(subjectvars.Menstrual_UsingBirthControl==0));
keepmenandhighestr = [keepmale; keephighestr];

% all menstrual cycle but not the high estrogen values.
% keep menses - keephighestr

keeplowestr = ~ismember(keep_menses, keephighestr);
keeplowestr = keep_menses(keeplowestr);  
keepmenandlowestr = [keepmale; keeplowestr];



%% ANALYSIS

%% Analyses.
% set up all covariates
sex = covariates.sex;
age = covariates.age;
icv = covariates.icv;


for splits = 1:1000
    modelname{1} = 'Males_OC_fem';
    randomIdxOC = randperm(length(keepmale),length(keep_takeOC)*2);
    keep_subsample_one{1} = [keepmale(randomIdxOC(1:(length(keep_takeOC)))); keep_takeOC];
    keep_subsample_two{1} = [keepmale(randomIdxOC((length(keep_takeOC)+1:length(keep_takeOC)*2))); keep_takeOC];

    modelname{2} = 'Males_high_prog';
    randomIdxlowProg = randperm(length(keepmale),length(keeplowprog)*2);
    keep_subsample_one{2} = [keepmale(randomIdxlowProg(1:(length(keeplowprog)))); keeplowprog];
    keep_subsample_two{2} = [keepmale(randomIdxlowProg((length(keeplowprog)+1:length(keeplowprog)*2))); keeplowprog];

    modelname{3} = 'Males_low_prog';
    randomIdxhighProg = randperm(length(keepmale),length(keephighprog)*2);
    keep_subsample_one{3} = [keepmale(randomIdxhighProg(1:(length(keephighprog)))); keephighprog];
    keep_subsample_two{3} = [keepmale(randomIdxhighProg((length(keephighprog)+1:length(keephighprog)*2))); keephighprog];
    
    modelname{4} = 'Males_high_estr';
    randomIdxlowestr = randperm(length(keepmale),length(keeplowestr)*2);
    keep_subsample_one{4} = [keepmale(randomIdxlowestr(1:(length(keeplowestr)))); keeplowestr];
    keep_subsample_two{4} = [keepmale(randomIdxlowestr((length(keeplowestr)+1:length(keeplowestr)*2))); keeplowestr];
     
    modelname{5} = 'Males_low_estr';
    randomIdxhighestr = randperm(length(keepmale),length(keephighestr)*2);
    keep_subsample_one{5} = [keepmale(randomIdxhighestr(1:(length(keephighestr)))); keephighestr];
    keep_subsample_two{5} = [keepmale(randomIdxhighestr((length(keephighestr)+1:length(keephighestr)*2))); keephighestr];
    
    for modelnumber = 1:5
        name_contrast = modelname{modelnumber};
        clear keep_one keep_two

        keep_one = keep_subsample_one{modelnumber};
        ageterm_one = term(age(keep_one));
        icvterm_one = term(icv(keep_one));   
        groupcompterm_one = term(sex(keep_one));

        keep_two = keep_subsample_two{modelnumber};
        ageterm_two = term(age(keep_two));
        icvterm_two = term(icv(keep_two));   
        groupcompterm_two = term(sex(keep_two));
    
        for moment = 1:size(T1T2moments, 3) 
            % 1 = gradient; 2 = mean; 3 = skewness
            M_one = 1 + groupcompterm_one + ageterm_one + icvterm_one;
            slm_one = SurfStatLinMod(T1T2moments(keep_one,:,moment),M_one);
            slm_one = SurfStatT(slm_one,groupcompterm_one.F-groupcompterm_one.M); %contrast is female - male
            Cohensd_one = 2*slm_one.t / sqrt(slm_one.df);
                
            M_two = 1 + groupcompterm_two + ageterm_two + icvterm_two;
            slm_two = SurfStatLinMod(T1T2moments(keep_two,:,moment),M_two);
            slm_two = SurfStatT(slm_two,groupcompterm_two.F-groupcompterm_two.M); %contrast is female - male
            Cohensd_two = 2*slm_two.t / sqrt(slm_two.df);
           
            halfsplitD.(name_contrast)(moment,splits) = corr(Cohensd_one', Cohensd_two');  
        end
    end
end


%% Plotting.
for model = 1:5
    name_contrast = modelname{model};

    r1 = 1 - rand(500,1)*0.1;
    r2 = 1 + rand(500,1)*0.1;
    r = [r1;r2];
    
    colors = colormap(flipud(cbrewer('seq','Greens',11)));
    
    f=figure,
    hold on
    subplot(1,3,1)
    hold on 
    scatter(r, halfsplitD.(name_contrast)(2,:), 'MarkerEdgeColor', [0.4039         0    0.1216], 'Marker', '.');
    % mean = 0.8622, range = 0.4871, std = 0.0563
    boxplot(halfsplitD.(name_contrast)(2,:), 'widths', 0.5, 'Colors', [0.4039         0    0.1216], 'Symbol', 'w');
    ylabel(['Spatial Correlation'])
    xlabel(['Profile Mean'])
    set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])
    
    subplot(1,3,2)
    hold on
    scatter(r, halfsplitD.(name_contrast)(3,:), 'MarkerEdgeColor', [0.0196    0.1882    0.3804], 'Marker', '.');
    boxplot(halfsplitD.(name_contrast)(3,:), 'widths', 0.5, 'Colors', [0.0196    0.1882    0.3804], 'Symbol', 'w');
    ylabel(['Spatial Correlation'])
    xlabel(['Profile Skewness'])
    set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])
    
    subplot(1,3,3)
    hold on 
    scatter(r, halfsplitD.(name_contrast)(1,:), 'MarkerEdgeColor', [0    0.2667    0.1059], 'Marker', '.');
    % mean = 0.5506, range = 0.3588, std = 0.0408
    boxplot(halfsplitD.(name_contrast)(1,:), 'widths', 0.5, 'Colors', [0    0.2667    0.1059], 'Symbol', 'w');
    ylabel(['Spatial Correlation'])
    xlabel(['Gradient'])
    set(gca, 'FontSize', 25, 'YLim', [0.4, 1], 'fontname', 'Calibri', 'xtick', [])
    
    sgtitle(name_contrast)
end

