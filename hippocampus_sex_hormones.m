% toolboxes used
% brainspace, cbrewer, export_fig, RainCloudPlots, spider_plot, surfstat
% input:'subjectListS1200_hipp_unfold.txt', 'subj_HCP_Lina.txt'
% output:

% This script computes sex-differences and hormonal effects on the
% hippocampus
% Supplement 6

loadold = 1;
if loadold == 0
    clear all
    loadold = 0;
end

saveall = 1;
smallsample = 1; %more rigorous group-making

% make a path to directory
addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';
scriptDir = fullfile(homeDir,'script');
addpath(genpath(scriptDir));
dataDir = fullfile(homeDir, 'data');
addpath(genpath(dataDir));
outDir = fullfile(homeDir,'output/hippocampus');
addpath(genpath(outDir));
hippoDir = fullfile(dataDir, 'hipp_unfold_T1T2w');
addpath(genpath(hippoDir));

%% loading data
if loadold == 0
    hippo = zeros(2,867,7262); %7262 values for the original resolution
    for hemi_lr = 1:2
        files(:,1) = dir(fullfile(hippoDir, '*.txt'));
        files(:,2) = dir(fullfile(hippoDir, 'right_hippo', '*.txt'));       
        name_files = (fullfile(dataDir, 'subjectListS1200_hipp_unfold.txt'));
        subjects = importdata(name_files); %867 subjects
        for i=1:length(files)
            filename = files(i,hemi_lr).name;
            hippo(hemi_lr,i,:) = importdata(filename); 
        end
    end
     
    
    load subjects_HCP.mat; % subjects, sex, unres.FS_IntraCranial_Vol 
    [C, keep, IB] = intersect(unres.Subject, subjects); %unres is subjectlist
    subjectlist = unres(keep,:); 
    sex = cellstr(subjectlist.Gender);
      
    icv =  subjectlist.FS_IntraCranial_Vol;
    load subjects_HCP_ageBMI.mat; % BMI, age, handedness, menstrual stuff,...%D1
    subjectvars = D1(keep,:);
    age = subjectvars.Age_in_Yrs; 
    load Campbell_limbic_mask.mat; % limbic cortex mask
    load fsaverage5.mat; % fs5 struct to plot 
    
    keepmale = find(cell2mat(sex) == 'M');
    keepfemale = find(cell2mat(sex) == 'F');
    

    % create keep vectors: 
    OC = subjectvars.Menstrual_UsingBirthControl; 
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
    OCcov = OC(keep_OCrecorded);    
    OCcell = num2cell(OCcov);
    OCcell(OCcov == 1) = {'yes'};
    OCcell(OCcov == 0) = {'no'};

    % prepare menstrual cycle measure
    % D1(:,40) is Menstrual Cycle Length: 2= between 25-34 days
    % D1(:,41) is Days since last Menstruation: all that are available

    menses = subjectvars.Menstrual_DaysSinceLast;
    
    if smallsample == 1
        % exclude all that take birth control, that dont have a regular cycle,
        % and that are more than 28 days since menstruation.
        % women with a regular cycle shorter than 29 days = 284
        % high estr = 184
        % low estr = 100
        % high prog = 113
        % low prog = 171
        % OC = 170
        % overall 431 women with regular cycles or OC taking
        keep_menses = mintersect(find(~isnan(menses)), find(menses<29), find(menses>-1), ...
            find(subjectvars.Menstrual_RegCycles==1), find(subjectvars.Menstrual_UsingBirthControl==0));

        dayssincemenses = menses(keep_menses);
        progconce = num2cell(dayssincemenses);
        progconce(dayssincemenses < 15) = {'low'};
        progconce(dayssincemenses > 14) = {'high'};

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

        keephighestr = mintersect(find(~isnan(menses)), find(menses<23), ...
            find(menses>6), find(subjectvars.Menstrual_RegCycles==1), ...
            find(subjectvars.Menstrual_UsingBirthControl==0));
        keepmenandhighestr = [keepmale; keephighestr];

        % all menstrual cycle but not the high estrogen values.
        % keep menses - keephighestr

        keeplowestr = ~ismember(keep_menses, keephighestr);
        keeplowestr = keep_menses(keeplowestr);  
        keepmenandlowestr = [keepmale; keeplowestr];

    elseif smallsample == 0
            % 120 -> liberal criterium to include more women 
            % 198 high estr
            % 84 low estr
            % 108 high prog
            % 174 low prog
            % 143 OC
            % 367 male
            % 356 NC           
            %only 398 women menses<35
            
            keep_menses = mintersect(find(~isnan(menses)), find(menses<35), ...
                find(menses>-1), find(subjectvars.Menstrual_UsingBirthControl==0));
   
            dayssincemenses = menses(keep_menses);

            progconce = num2cell(dayssincemenses);
            progconce(dayssincemenses < 17) = {'low'}; 
            progconce(dayssincemenses > 16) = {'high'};

            keeplowprog = mintersect(find(~isnan(menses)), find(menses<17), ...
                find(menses>-1), find(subjectvars.Menstrual_UsingBirthControl==0));
            keepmenandlowprog = [keepmale; keeplowprog];

            keephighprog = mintersect(find(~isnan(menses)), find(menses<35), ...
                find(menses>16), find(subjectvars.Menstrual_UsingBirthControl==0));
            keepmenandhighprog = [keepmale; keephighprog];

            estrconce = num2cell(dayssincemenses);
            estrconce(dayssincemenses < 7) = {'low'};
            estrconce(dayssincemenses > 6) = {'high'};
            estrconce(dayssincemenses > 23) = {'low'};

            keephighestr = mintersect(find(~isnan(menses)), find(menses<26), ...
                find(menses>5), find(subjectvars.Menstrual_UsingBirthControl==0));
            keepmenandhighestr = [keepmale; keephighestr];

            % all menstrual cycle but not the high estrogen values.
            % keep menses - keephighestr
            
            keeplowestr = ~ismember(keep_menses, keephighestr);
            keeplowestr = keep_menses(keeplowestr);  
            keepmenandlowestr = [keepmale; keeplowestr];
    end
    
elseif loadold == 1
    load hippo_lr_sex_horms_out.mat;
    keyboard
end


%% ANALYSIS 1: Hippocampal sex differences

for hipposexdiffs = 1 

    sexcov = term(cellstr(sex));
    agecov = term(age);
    ICVcov = term(icv);

    for hemi_lr = 1:2 
        M = 1 + sexcov + agecov + ICVcov;
        slm_hipp = SurfStatLinMod(squeeze(hippo(hemi_lr,:,:)),M);
        slm_hipp = SurfStatT(slm_hipp,sexcov.F-sexcov.M); %contrast is female - male
    
        % to save the hippocampal file
        hipposex(hemi_lr,:) = slm_hipp.t';
        
        % FDR correction
        p = 1-tcdf(slm_hipp.t,slm_hipp.df);
        h = fdr_bh(p,0.025);
        hn = fdr_bh(1-p,0.025);
        h= h+hn;
        fdr_hippo_tval(hemi_lr,:) = (slm_hipp.t.*h)';

        slm_hipp.df 
        Cohensd(hemi_lr,:) = (2*slm_hipp.t / sqrt(slm_hipp.df))';
        CohensdFDR(hemi_lr,:) = (squeeze(Cohensd(hemi_lr,:)).*h);
        
        if hemi_lr == 1
            save(fullfile(outDir, 'hipposex_lh.txt'), 'hipposex', '-ascii');
            save(fullfile(outDir, 'hipposexFDR_lh.txt'), 'fdr_hippo_tval', '-ascii');
            Cohensd_lh = Cohensd(hemi_lr,:);
            save(fullfile(outDir, 'hipposexCohensd_lh.txt'), 'Cohensd_lh', '-ascii');
            Cohensd_FDR_lh = CohensdFDR(hemi_lr,:);
            save(fullfile(outDir, 'hipposexFDRCohensd_lh.txt'), 'Cohensd_FDR_lh', '-ascii');
        elseif hemi_lr == 2   
            save(fullfile(outDir, 'hipposex_rh.txt'), 'hipposex', '-ascii');
            save(fullfile(outDir, 'hipposexFDR_rh.txt'), 'fdr_hippo_tval', '-ascii');
            Cohensd_rh = Cohensd(hemi_lr,:);
            save(fullfile(outDir, 'hipposexCohensd_rh.txt'), 'Cohensd_rh', '-ascii');
            Cohensd_FDR_rh = CohensdFDR(hemi_lr,:);
            save(fullfile(outDir, 'hipposexFDRCohensd_rh.txt'), 'Cohensd_FDR_rh', '-ascii');
        end
    
        % identify mean of all positive and negative effects
         bigwomen = find(squeeze(CohensdFDR(:,hemi_lr))<0); %-.3159 lh, 0.2554 lh
         bigwomeneffect = mean(CohensdFDR(bigwomen,hemi_lr)); %-.3315 rh; 0.2559 rh
         bigwomeneffect(hemi_lr,:) = mean(squeeze(CohensdFDR(hemi_lr,bigwomen(hemi_lr,:)))); %0.25
       
         % to include both hemispheres in one 
        CohensdFDR_all = [CohensdFDR(hemi_lr,:), CohensdFDR(hemi_lr,:)];
        bigwomen = find(CohensdFDR_all>0);
        bigwomeneffect = mean(CohensdFDR_all(bigwomen));
        bigmen = find(CohensdFDR_all<0);
        bigmeneffect = mean(CohensdFDR_all(bigmen));
    
        results.bigwomeneffect(:) = bigwomeneffect;
        results.bigmeneffect(:) = bigmeneffect;
        
        descriptives.mean(1) = mean(mean(CohensdFDR_all(keepfemale)));
        descriptives.mean(2) = mean(mean(CohensdFDR_all(keepmale)));
        descriptives.std(1,:) = mean(std(CohensdFDR_all(keepfemale)));
        descriptives.std(2,:) = mean(std(CohensdFDR_all(keepmale)));
        descriptives.posmean(1) = mean(CohensdFDR_all(bigwomen));
        descriptives.negmean(1) = mean(CohensdFDR_all(bigmen));
        descriptives.posstd(1) = std(CohensdFDR_all(bigwomen));
        descriptives.negstd(1) = std(CohensdFDR_all(bigmen));
    end
end

%keyboard
%% ANALYSIS 2: HORMONAL ANALYSIS

% how many models am I testing?
amount_models = 9;

% define inclusion per model so I can loop through them
keep_subsample{1} = keep_OCrecorded; % 1. to compare OC vs no OC women
modelname{1} = 'women-OCvsNC';
keep_subsample{2} = excl_OC; % 2. to compare men vs. NC women
modelname{2} = 'MenvsNC-women';
keep_subsample{3} = excl_natcyc; % 3. to compare men vs. OC women
modelname{3} = 'MenvsOC-women';
keep_subsample{4} = keep_menses; % keep only those NC women with regular cycle >> estr
modelname{4} = 'women-highvslow-estrogen';
keep_subsample{5} = keep_menses; % keep only those NC women with regular cycle >> prog
modelname{5} = 'women-highvslow-progesterone';
keep_subsample{6} = keepmenandhighprog; % keep men high prog women with regular cycle
modelname{6} = 'Menvshigh-prog';
keep_subsample{7} = keepmenandlowprog; % keep men and low prog women with regular cycle
modelname{7} = 'Menvslow-prog';
keep_subsample{8} = keepmenandhighestr; % keep men high estr women with regular cycle
modelname{8} = 'Menvshigh-estr';
keep_subsample{9} = keepmenandlowestr; % keep men and low estr women with regular cycle
modelname{9} = 'Menvslow-estr';

%for modelnumber = 8
for hemi_lr = 1:2 
    for modelnumber = 1:amount_models
        clear keep meant1t2part 
        keep = keep_subsample{modelnumber};
        
        agepart = age(keep);
        ageterm = term(agepart);
        
        icvpart = icv(keep);
        icvterm = term(icvpart);
        
        name_contrast = modelname{modelnumber};
        
        % in one model, it has to be sexterm, in the other OCterm.
        if modelnumber == 1
            groupcompterm = term(OCcell);
        elseif (modelnumber == 2) || (modelnumber == 3) || (modelnumber == 6) || (modelnumber == 7) || (modelnumber == 8) || (modelnumber == 9)
                % 2 is men vs. NC women, 3 is men vs. OC women, 
                % 6 is men vs. high prog, 7 is men vs. low prog
                % 8 is men vs. high estr, 9 is men vs. low estr
            sexpart = sex(keep);
            groupcompterm = term(sexpart);
        elseif modelnumber == 4
                    groupcompterm = term(estrconce);
        elseif modelnumber == 5 
                        groupcompterm = term(progconce);
        end
       
        M = 1 + groupcompterm + ageterm + icvterm;
        if modelnumber == 1  
            slm = SurfStatLinMod(squeeze(hippo(hemi_lr,keep, :)),M);
            slm = SurfStatT(slm,groupcompterm.no-groupcompterm.yes); %contrast is noOC - OCtakers
        elseif modelnumber == (modelnumber == 2) || (modelnumber == 3) || (modelnumber == 6) || (modelnumber == 7) || (modelnumber == 8) || (modelnumber == 9)
            % 2 is men vs. NC women, 3 is men vs. OC women,
            % 6 is men vs. high prog, 7 is men vs. low prog, 
            % 8 is men vs. high estr, 9 is men vs. low estr
            slm = SurfStatLinMod(squeeze(hippo(hemi_lr,keep, :)),M);
            slm = SurfStatT(slm,groupcompterm.F-groupcompterm.M); %contrast is female - male
        elseif modelnumber == 4 || modelnumber == 5 % 4 is estrogen, 5 is progesterone
                    slm = SurfStatLinMod(squeeze(hippo(hemi_lr,keep, :)),M);;
                    slm = SurfStatT(slm,groupcompterm.high-groupcompterm.low); %contrast is high - low
        end
    
        % Cohen's D
        % Cohen's D = 2t / sqrt(df)
        Cohensd = (2*slm.t / sqrt(slm.df))';
    
        % FDR correction
        p = 1-tcdf(slm.t,slm.df);
        h = fdr_bh(p,0.025);
        hn = fdr_bh(1-p,0.025);
        h= h+hn;
    
        % store results and set measure name for figures
        tvals_hippo{modelnumber, hemi_lr} = (slm.t)';
        tvals = (slm.t)';
        cohensd_hippo{modelnumber, hemi_lr} = Cohensd;
        tFDR_hippo{modelnumber, hemi_lr} = (slm.t .* h)';
        tFDR = (slm.t .* h)';
        dFDR_hippo{modelnumber, hemi_lr} = (Cohensd.*h');
        dFDR = (Cohensd.*h');

        % save the hippocampal file
        if hemi_lr == 1  
            save(fullfile(outDir, sprintf('lh_hippo_tvals_%s.txt', modelname{modelnumber})), 'tvals', '-ascii');
            save(fullfile(outDir, sprintf('lh_hippo_tFDR_%s.txt', modelname{modelnumber})), 'tFDR', '-ascii');
            save(fullfile(outDir, sprintf('lh_hippo_Cohensd_%s.txt', modelname{modelnumber})), 'Cohensd', '-ascii');
            save(fullfile(outDir, sprintf('lh_hippo_dFDR_%s.txt', modelname{modelnumber})), 'dFDR', '-ascii');
        elseif hemi_lr == 2
            save(fullfile(outDir, sprintf('rh_hippo_tvals_%s.txt', modelname{modelnumber})), 'tvals', '-ascii');
            save(fullfile(outDir, sprintf('rh_hippo_tFDR_%s.txt', modelname{modelnumber})), 'tFDR', '-ascii');
            save(fullfile(outDir, sprintf('rh_hippo_Cohensd_%s.txt', modelname{modelnumber})), 'Cohensd', '-ascii');
            save(fullfile(outDir, sprintf('rh_hippo_dFDR_%s.txt', modelname{modelnumber})), 'dFDR', '-ascii');
        end
    end
end

for modelnumber = 1:amount_models
    % then safe some descriptives
    % then include both hemispheres in one file, per model
    CohensdFDR_all_H(modelnumber,:) = [dFDR_hippo{modelnumber, 1}; dFDR_hippo{modelnumber, 2}];
    bigwomen_H = find(CohensdFDR_all_H(modelnumber,:)>0);
    bigwomeneffect_H(modelnumber) = mean(CohensdFDR_all_H(modelnumber,bigwomen_H));
    bigmen_H = find(CohensdFDR_all_H(modelnumber,:)<0);
    bigmeneffect_H(modelnumber) = mean(CohensdFDR_all_H(modelnumber,bigmen_H));
    
    resultsH.bigwomeneffect(modelnumber,:) = bigwomeneffect_H(modelnumber);
    resultsH.bigmeneffect(modelnumber,:) = bigmeneffect_H(modelnumber);
    
    descriptivesH.posmean(modelnumber) = mean(CohensdFDR_all_H(modelnumber,bigwomen_H));
    descriptivesH.negmean(modelnumber) = mean(CohensdFDR_all_H(modelnumber,bigmen_H));
    descriptivesH.posstd(modelnumber) = std(CohensdFDR_all_H(modelnumber,bigwomen_H));
    descriptivesH.negstd(modelnumber) = std(CohensdFDR_all_H(modelnumber,bigmen_H));
end


%% ANOVA part: post-hoc testing which results were different.
for i = 1 :size((CohensdFDR_all_H), 1)
    cohensd_mean_mat(:, i) = CohensdFDR_all_H(i,:)';
end
cohensd_mean_mat(:, 10) = CohensdFDR_all';
modelname{10} = 'all_sex_diffs'
 

[p, t, stats_d_mean] = anova1(cohensd_mean_mat)
[all_comps_mean_cohensd,mean_mean,h,gnames] = multcompare(stats_d_mean)

tbl = array2table(mean_mean, "RowNames", modelname, "VariableNames", ["Mean", "Standard Error"])



%% Save all results
disp(' ...saving results')
if saveall == 1
    save(fullfile(outDir, 'hippo_lr_sex_horms_out.mat'));
end

disp(' ...done!:)')