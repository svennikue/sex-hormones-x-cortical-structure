% This script analysis if microstructural moments changes as a function of
% hormones. 
% last edited 2023, written by Svenja Kuchenhoff.
% firstly, it compares men to women taking oral contraceptives (OC) vs.
% women that are naturally cycling (NC); and secondly, it compares men to
% women in different phases of their menstrual cycle. 
% this scripts requires to have run sex_diffs_3measures.m before

loadold = 0;

if loadold == 0
    clear all
    close all
    loadold = 0;
end

saveall = 1;
plotfigs = 1;
smallsample = 1; % always take the smaller one (final analysis)


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

%% loading data
if loadold == 0
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
    
    if smallsample == 1
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

    else if smallsample == 0
            % with a slightly different sample
            % not sure if these numbers are correct [22.05.23]
            % women with a regular cycle = 334 
            % high estr = 237
            % low estr = 97
            % high prog = 129
            % low prog = 205
            % OC = 170
            
            keep_menses = mintersect(find(~isnan(menses)), find(menses<35), ...
                find(menses>-1), find(subjectvars.Menstrual_UsingBirthControl==0));
            
            % overall 401 women with regular cycles or OC taking
            keep_menses_andOC = mintersect(find(~isnan(menses)), find(menses<29), find(menses>-1), ...
            find(subjectvars.Menstrual_RegCycles==1), find(~isnan(OC)));
            
            % prepare OC covariate
            OCcov = OC(keep_menses_andOC);    

            dayssincemenses = menses(keep_menses);

            progconce = num2cell(dayssincemenses);
            progconce(dayssincemenses < 17) = {'low'};
            progconce(dayssincemenses > 16) = {'high'};

            % prepare OC vs prog women covariate
            progcov = menses(keep_menses_andOC);    
            OC_vs_prog = num2cell(progcov);
            OC_vs_prog(progcov < 17) = {'low'};
            OC_vs_prog(progcov < 17) = {'low'};
            OC_vs_prog(OCcov == 1) = {'OC'};


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

            % prepare OC vs estr women covariate
            estrcov = menses(keep_menses_andOC);    
            OC_vs_estr = num2cell(estrcov);
            OC_vs_estr(estrcov < 7) = {'low'};
            OC_vs_estr(estrcov > 6) = {'high'};
            OC_vs_estr(estrcov > 23) = {'low'};
            OC_vs_estr(OCcov == 1) = {'OC'};

            keephighestr = mintersect(find(~isnan(menses)), find(menses<26), ...
                find(menses>5), find(subjectvars.Menstrual_UsingBirthControl==0));
            keepmenandhighestr = [keepmale; keephighestr];

            % all menstrual cycle but not the high estrogen values.
            % keep menses - keephighestr
            
            keeplowestr = ~ismember(keep_menses, keephighestr);
            keeplowestr = keep_menses(keeplowestr);  
            keepmenandlowestr = [keepmale; keeplowestr];
        end
    end
      
else if loadold == 1        
    % alternatively: load complete mat from running this script last time. 
    load('hormonal_effects_moments.mat');
    keyboard
    end
end


%% Analyses.
% set up all covariates
meant1t2 = covariates.meant1t2;
sex = covariates.sex;
age = covariates.age;
icv = covariates.icv;


% define inclusion per model so I can loop through them
keep_subsample{1} = keep_OCrecorded; % 1. to compare OC vs no OC women
modelname{1} = 'fem: OC vs NC';
keep_subsample{2} = excl_OC; % 2. to compare men vs. NC women
modelname{2} = 'Men vs NC women';
keep_subsample{3} = excl_natcyc; % 3. to compare men vs. OC women
modelname{3} = 'Men vs OC women';
keep_subsample{4} = keep_menses; % keep only those NC women with regular cycle >> estr
modelname{4} = 'fem: h v l estr';
keep_subsample{5} = keep_menses; % keep only those NC women with regular cycle >> prog
modelname{5} = 'fem: h v l prog';
keep_subsample{6} = keepmenandhighprog; % keep men high prog women with regular cycle
modelname{6} = 'Men vs high prog';
keep_subsample{7} = keepmenandlowprog; % keep men and low prog women with regular cycle
modelname{7} = 'Men vs low prog';
keep_subsample{8} = keepmenandhighestr; % keep men high estr women with regular cycle
modelname{8} = 'Men vs high estr';
keep_subsample{9} = keepmenandlowestr; % keep men and low estr women with regular cycle
modelname{9} = 'Men vs low estr';

% additional female models.
modelname{10} = 'fem: OC vs high estr';
keep_subsample{10} =  keep_menses_andOC; % 10. to compare OC vs high estr NC
modelname{11} = 'fem: OC vs high prog';
keep_subsample{11} =  keep_menses_andOC; % 11. to compare OC vs high prog NC
modelname{12} = 'fem: OC vs low estr';
keep_subsample{12} =  keep_menses_andOC; % 12. to compare OC vs low estr NC
modelname{13} = 'fem: OC vs low prog';
keep_subsample{13} =  keep_menses_andOC; % 13. to compare OC vs low prog NC


% for modelnumber = 9
for modelnumber = 1:9
%for modelnumber = 3
    clear keep meant1t2part 
    keep = keep_subsample{modelnumber};
    
    meant1t2part = meant1t2(keep);
    meant1t2cov = term(meant1t2part);
    
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
    elseif (modelnumber == 10) || (modelnumber == 12)
        groupcompterm = term(OC_vs_estr); 
    elseif (modelnumber == 11) || (modelnumber == 13)
        groupcompterm = term(OC_vs_prog);
    end
       
    % start with the analysis for all 3 measures.
    %for moment = 1 % or only choose one
    for moment = 1:size(T1T2moments, 3) 
        % ADJUST MODEL WITHOUT MEANT1T2!!
        M = 1 + groupcompterm + ageterm + icvterm;
         %M = 1 + meant1t2cov + groupcompterm + ageterm + icvterm;
        if modelnumber == 1  
            slm = SurfStatLinMod(T1T2moments(keep,:,moment),M);
            slm = SurfStatT(slm,groupcompterm.no-groupcompterm.yes); %contrast is noOC - OCtakers
        elseif modelnumber == (modelnumber == 2) || (modelnumber == 3) || (modelnumber == 6) || (modelnumber == 7) || (modelnumber == 8) || (modelnumber == 9)
            % 2 is men vs. NC women, 3 is men vs. OC women,
            % 6 is men vs. high prog, 7 is men vs. low prog, 
            % 8 is men vs. high estr, 9 is men vs. low estr
            slm = SurfStatLinMod(T1T2moments(keep,:,moment),M);
            slm = SurfStatT(slm,groupcompterm.F-groupcompterm.M); %contrast is female - male
        elseif modelnumber == 4 || modelnumber == 5 % 4 is estrogen, 5 is progesterone
            slm = SurfStatLinMod(T1T2moments(keep,:,moment),M);
            slm = SurfStatT(slm,groupcompterm.high-groupcompterm.low); %contrast is high - low
        elseif (modelnumber == 10) || (modelnumber == 11)
            slm = SurfStatLinMod(T1T2moments(keep_menses_andOC, :, moment), M);
            slm = SurfStatT(slm, groupcompterm.high - groupcompterm.OC); %contrast NC high - OC
        elseif (modelnumber == 12) || (modelnumber == 13)
            slm = SurfStatLinMod(T1T2moments(keep_menses_andOC, :, moment), M);
            slm = SurfStatT(slm, groupcompterm.low - groupcompterm.OC); %contrast NC low - OC 
        end
        
        
        % Cohen's D
        % Cohen's D = 2t / sqrt(df)
        Cohensd = 2*slm.t / sqrt(slm.df);
        
        % FDR correction
        p = 1-tcdf(slm.t,slm.df);
        h = fdr_bh(p,0.025);
        hn = fdr_bh(1-p,0.025);
        h= h+hn;
        
        % store results and set measure name for figures
        if moment == 1
            tvals_gradient{modelnumber} = slm.t;
            cohensd_gradient{modelnumber} = Cohensd;
            tFDR_gradient{modelnumber} = slm.t .* h;
            namemoment = 'Gradient';
        elseif moment == 2
            tvals_profilemean{modelnumber} = slm.t;
            cohensd_profilemean{modelnumber} = Cohensd;
            tFDR_profilemean{modelnumber} = slm.t .* h;
            namemoment = 'profile Mean';
        elseif moment == 3
            tvals_profileskew{modelnumber} = slm.t;
            cohensd_profileskew{modelnumber} = Cohensd;
            tFDR_profileskew{modelnumber} = slm.t .* h;
            namemoment = 'profile Skewness';
        end
        
             
        if plotfigs == 1    
%             % plot sex difference t values
%             vertices = zeros(20484,1);
%             for i = 1:400 %200 parcels per hemisphere
%                 vertices(find(schaefer_400==i)) = slm.t(i);
%             end   
%             f = figure,
% 
%             BoSurfStatViewData(vertices,SN,'')
%             if (modelnumber == 1) ||  (modelnumber == 4) || (modelnumber == 5)
%                 colormap(flipud(cbrewer('div','PuOr',11)))
%             else 
%                 colormap(flipud(cbrewer('div','RdBu',11)))
%             end
%             %clim = [prctile(mean(slm.t),5) prctile(mean(slm.t),95)]; % set colour limits
%             SurfStatColLim([-10 10])
%             %SurfStatColLim(clim)
%             title(sprintf('Tvals for contrast %s, for %s', name_contrast, namemoment));
%             exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_hormonaleff_tvals_%s_contrast_%s', namemoment, name_contrast)],'png', 10)
%             count = count + 1
%             FigName{(count)} = sprintf('HCP_Tw1T2w_hormonaleff_tvals_%s_constrast_%s.fig', namemoment, name_contrast);
% 
% 
%             %plot tval only of those parcels that survive FDR correction
%             vertices = zeros(20484,1);
%             for i = 1:400 %200 parcels per hemisphere
%                 vertices(find(schaefer_400==i)) = slm.t(i)*h(i);
%             end
% 
%             f = figure,
%             BoSurfStatViewData(vertices,SN,'')
%             if (modelnumber == 1) ||  (modelnumber == 4) || (modelnumber == 5)
%                 colormap(flipud(cbrewer('div','PuOr',11)))
%             else 
%                 colormap(flipud(cbrewer('div','RdBu',11)))
%             end  
%             SurfStatColLim([-10 10])
%             title(sprintf('FDR corr tvals for contrast %s, for %s', name_contrast, namemoment));
%             exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_hormonaleff_FDR_%s_%s', namemoment, name_contrast)],'png', 10)
%             count = count + 1
%             FigName{(count)} = sprintf('HCP_Tw1T2w_hormonaleff_FDR_%s_%s.fig', namemoment, name_contrast);
%         end
        
            % plot EFFECT SIZE FDR corrected 

            vertices = zeros(20484,1);
            for i = 1:400 %200 parcels per hemisphere
                CohensdFDR(i) = Cohensd(i)*h(i);
                vertices(find(schaefer_400==i)) = Cohensd(i)*h(i);
                %vertices(find(schaefer_400==i)) = Cohensd(i);
            end
    
            f = figure,
            BoSurfStatViewData(vertices,SN,'')
            if (modelnumber == 1) ||  (modelnumber == 4) || (modelnumber == 5) || (modelnumber == 10) || (modelnumber == 11) || (modelnumber == 12) || (modelnumber == 13)
                colormap(flipud(cbrewer('div','PuOr',11))) %for comparisons between females
            else 
                colormap(flipud(cbrewer('div','RdBu',11))) % for females vs males contrasts
            end
            %clim = [prctile(mean(slm.t),5) prctile(mean(slm.t),95)]; % set colour limits
            SurfStatColLim([-1 1])
            %SurfStatColLim(clim)
            title(sprintf('FDR corr. Cohens d, %s, for %s', name_contrast, namemoment));
            exportfigbo(f,[figDir, sprintf('/HCP_hormonaleff_FDRCohensD_%s%s', namemoment, name_contrast)],'png', 10)
            count = count + 1
            FigName{(count)} = sprintf('HCP_hormonaleff_FDRCohensD_%s%s.tiff', namemoment, name_contrast);
    
            imagewd = getframe(gcf); 
            imwrite(imagewd.cdata, fullfile(figDir, sprintf('HCP_hormonaleff_FDRCohensD_%s%s.tiff', namemoment, name_contrast)));
            imagewd = getframe(gcf); 
            imwrite(imagewd.cdata, fullfile(figDir, sprintf('HCP_hormonaleff_FDRCohensD_%s%s.png', namemoment, name_contrast)));


            % identify mean of all positive and negative effects
            bigwomen = find(CohensdFDR>0);
            bigwomeneffect(modelnumber) = mean(CohensdFDR(bigwomen)); %
            minwomeneffect(modelnumber) = min(CohensdFDR); %
            bigmen = find(CohensdFDR<0);
            bigmeneffect(modelnumber) = mean(CohensdFDR(bigmen)) %
            maxmeneffect(modelnumber) = max(CohensdFDR) %
        end
    end
end

%% Visualize effects.
% plot effect sizes per parcel for each comparison
% very stupidly, this is hard-coded
% so always adjust input data number and Titles
value = 7;
% value 1+2 = skewness tval&FDR, 3+4 = mean tval&FDR, 5+6 = gradient tval&FDR
% 7 is gradient Cohens d, 8 is mean Cohens d, 9 is skewness Cohens d


% also load 'normal' sex differences
load HCP_T1wT2w_sexdiffsmaps.mat

    % the numbers given to d{x} are the models you want to plot
    % 2 is men vs. NC women, 3 is men vs. OC women, 
    % 6 is men vs. high prog, 7 is men vs. low prog
    % 8 is men vs. high estr, 9 is men vs. low estr
if value == 1 
    namecontrast = 'Skewness';
    d{1} = tvals_profileskew{1,3}; % OC women
    d{2} = tvals_profileskew{1,9}; % low estr
    d{3} = tvals_profileskew{1,7}; % low prog
    d{4} = tvals_profileskew{1,6}; % high prog
    d{5} = tvals_profileskew{1,8}; % high estr 
    d{6} = results.Cohensd(:,3)';
        else if value == 2
                namecontrast = 'tFDR across cortex, Skewness';                  
                d{1} = tFDR_profileskew{1,3}; % OC women
                d{2} = tFDR_profileskew{1,9}; % low estr
                d{3} = tFDR_profileskew{1,7}; % low prog
                d{4} = tFDR_profileskew{1,6}; % high prog
                d{5} = tFDR_profileskew{1,8}; % high estr 
                d{6} = results.Cohensd(:,3)';
            else if value == 3
                    namecontrast = 'Mean';
                    d{1} = tvals_profilemean{1,3}; % OC women
                    d{2} = tvals_profilemean{1,9}; % low estr
                    d{3} = tvals_profilemean{1,7}; % low prog
                    d{4} = tvals_profilemean{1,6}; % high prog
                    d{5} = tvals_profilemean{1,8}; % high estr 
                    d{6} = results.Cohensd(:,2)';
                else if value == 4
                        namecontrast = 'tFDR across cortex, Mean';
                        d{1} = tFDR_profilemean{1,3}; % OC women
                        d{2} = tFDR_profilemean{1,9}; % low estr
                        d{3} = tFDR_profilemean{1,7}; % low prog
                        d{4} = tFDR_profilemean{1,6}; % high prog
                        d{5} = tFDR_profilemean{1,8}; % high estr 
                        d{6} = results.Cohensd(:,2)';
                    else if value == 5
                            namecontrast = 'Gradient'; % TVAL
                            d{1} = tvals_gradient{1,3}; % OC women
                            d{2} = tvals_gradient{1,9}; % low estr
                            d{3} = tvals_gradient{1,7}; % low prog
                            d{4} = tvals_gradient{1,6}; % high prog
                            d{5} = tvals_gradient{1,8}; % high estr 
                            d{6} = results.Cohensd(:,1)'; % sex diff gradient
                        else if value == 6
                                namecontrast = 'tFDR across cortex, Gradient'; %tval FDR
                                d{1} = tFDR_gradient{1,3}; % OC women
                                d{2} = tFDR_gradient{1,9}; % low estr
                                d{3} = tFDR_gradient{1,7}; % low prog
                                d{4} = tFDR_gradient{1,6}; % high prog
                                d{5} = tFDR_gradient{1,8}; % high estr 
                                d{6} = results.Cohensd(:,1)'; % sex diff gradient
                            else if value == 7
                                    namecontrast = 'Gradient'; %Cohnes D
                                    d{1} = cohensd_gradient{1,3}; % OC women
                                    d{2} = cohensd_gradient{1,9}; % low estr
                                    d{3} = cohensd_gradient{1,7}; % low prog
                                    d{4} = cohensd_gradient{1,6}; % high prog
                                    d{5} = cohensd_gradient{1,8}; % high estr 
                                    d{6} = results.Cohensd(:,1)';
                                else if value == 8
                                        namecontrast = 'Profile Mean';
                                        d{1} = cohensd_profilemean{1,3}; % OC women
                                        d{2} = cohensd_profilemean{1,9}; % low estr
                                        d{3} = cohensd_profilemean{1,7}; % low prog
                                        d{4} = cohensd_profilemean{1,6}; % high prog
                                        d{5} = cohensd_profilemean{1,8}; % high estr 
                                        d{6} = results.Cohensd(:,2)';
                                    else if value == 9 
                                            namecontrast = 'Profile Skewness';
                                            d{1} = cohensd_profileskew{1,3}; % OC women
                                            d{2} = cohensd_profileskew{1,9}; % low estr
                                            d{3} = cohensd_profileskew{1,7}; % low prog
                                            d{4} = cohensd_profileskew{1,6}; % high prog
                                            d{5} = cohensd_profileskew{1,8}; % high estr 
                                            d{6} = results.Cohensd(:,3)';
                                        end
                                    end
                                end
                            end
                        end         
                    end
                end              
            end
end

% also compute the difference between the groups 
% these tables also show the means per group
% all_comps matrix contains:
    % Group, Control Group, Lower Limit, Difference, Upper Limit, Pval

modelname{size(modelname,2)+1} = 'men v all women';


for i = 1 :size((cohensd_profilemean), 2)
    cohensd_mean_mat(:, i) = cohensd_profilemean{i};
end
cohensd_mean_mat(:, size(modelname,2)) = results.Cohensd(:,2)';


[p, t, stats_d_mean] = anova1(cohensd_mean_mat)
[all_comps_mean_cohensd,mean_mean,h,gnames] = multcompare(stats_d_mean)

tbl = array2table(mean_mean, "RowNames", modelname, "VariableNames", ["Mean", "Standard Error"])

for i = 1 :size((cohensd_profileskew), 2)
    cohensd_skew_mat(:, i) = cohensd_profileskew{i};
end

cohensd_skew_mat(:, size(modelname,2)) = results.Cohensd(:,3)';

[p, t, stats_d_skew] = anova1(cohensd_skew_mat)
[all_comps_skew_cohensd,mean_skew,h,gnames] = multcompare(stats_d_skew)


tbl = array2table(mean_skew, "RowNames", modelname, "VariableNames", ["Mean", "Standard Error"])



for i = 1 :size((cohensd_gradient), 2)
    cohensd_gradient_mat(:, i) = cohensd_gradient{i};
end
cohensd_gradient_mat(:, size(modelname,2)) = results.Cohensd(:,1)';


[p, t, stats_d_gradient] = anova1(cohensd_gradient_mat)
[all_comps_gradient_cohensd,mean_grad,h,gnames] = multcompare(stats_d_gradient)

tbl = array2table(mean_grad, "RowNames", modelname, "VariableNames", ["Mean", "Standard Error"])

%% additional scatter plot for revision


% modelname{1} = 'women: OC vs NC';
% modelname{2} = 'Men vs NC women';
% modelname{3} = 'Men vs OC women';
% modelname{4} = 'women: high vs low estrogen';
% modelname{5} = 'women: high vs low progesterone';
% modelname{6} = 'Men vs high prog';
% modelname{7} = 'Men vs low prog';
% modelname{8} = 'Men vs high estr';
% modelname{9} = 'Men vs low estr';

% for skew
% for i = 1:9
%     datatest(i,:) = cohensd_profileskew{i};
% end
% % select only those that I want
% data_wanted_models(1,:) = results.Cohensd(:,3)'; % main effect
% data_wanted_models(2,:) = datatest(3,:); % men vs OC
% for i = 6:9
%     data_wanted_models(i-3,:) = datatest(i,:);
% end

% for gradient
for i = 1:9
    datatest(i,:) = cohensd_gradient{i};
end
%select only those that I want
data_wanted_models(1,:) = results.Cohensd(:,1)'; % main effect
data_wanted_models(2,:) = datatest(3,:); % men vs OC
for i = 6:9
    data_wanted_models(i-3,:) = datatest(i,:);
end


%for mean
% for i = 1:9
%   datatest(i,:) = cohensd_profilemean{i};
% end 
% 
% % select only those that I want
% data_wanted_models(1,:) = results.Cohensd(:,2)'; % main effect
% data_wanted_models(2,:) = datatest(3,:); % men vs OC
% for i = 6:9
%     data_wanted_models(i-3,:) = datatest(i,:);
% end


[vertices, label, colortablel] = ...
  fs_read_annotation(['tpl-fsaverage_hemi-L_desc-types.annot']);
parcel_left = label;
label_left = label;
for i = 1:size(colortablel.table, 1)
  mycode = colortablel.table(i,5);
  parcel_left(find(parcel_left == mycode)) = i;
end
 
[vertices, label, colortabler] = ...
  fs_read_annotation(['tpl-fsaverage_hemi-R_desc-types.annot']);
parcel_right = label;
label_right = label;
for i = 1:size(colortabler.table, 1)
  mycode = colortabler.table(i,5);
  parcel_right(find(parcel_right == mycode)) = i;
end
 
types_garcia = [parcel_left(1:10242); parcel_right(1:10242)];
types_garcia = types_garcia';
 
       
% Generate sample data
data = data_wanted_models';
 
% Define colors for each data point
load Cortical_types.mat;
colors = types400;
 
% Create a plot matrix with scatter and histograms
figure;
 
model_names = {'All Females - Males','OC F - M',...
'High prog F - M','Low prog F - M','High Estr F - M','Low Estr F - M'};
 
% Scatter plots
for i = 1:size(data,2)
  for j = 1:size(data,2)
    subplot(size(data,2), size(data,2), (i - 1) * size(data,2) + j);
    if i == j
      % Diagonal - Histogram
      % Diagonal - Density Plot
      [f, xi] = ksdensity(data(:, i));
      plot(xi, f, 'LineWidth', 2, 'Color', 'k');
      title([model_names{i}]);
      xlim([-0.5 0.5])
      ylim([0 4])
      xline(0)
    else
      % Scatter plot
      h = scatter(data(:, i), data(:, j), 10, colors, 'filled'),lsline;
      colormap(colortabler.table(:,1:3)./256)
      xlim([-0.5 0.5])
      ylim([-0.5 0.5])
      yline(0)
      %title(['Variable ' num2str(i) ' vs Variable ' num2str(j)]);
      if i < size(data,2)
        xticks([]); % Remove x-axis labels for non-bottom subplots
      end
      if j > 1
        yticks([]); % Remove y-axis labels for non-leftmost subplots
      end
    end
  end
end





%% plot the results in the cloud-plots

% create mean and variance
means = cellfun(@mean, d);
variances = cellfun(@std, d);

[cb] = cbrewer('div', 'RdBu', 11); % define colours
fig_position = [200 200 600 400]; % coordinates for figures

% read into cell array of the appropriate dimensions
% make figure
data = d(1:6)';

if value == 7 | value == 8 | value== 9
    f  = figure('Position', fig_position);
    h   = rm_raincloud(data, cb(1,:));
    set(gca, 'XLim', [-0.8, 0.8], 'YTickLabel',...
        {'All women', 'High estrogen','High progesterone','Low progesterone','Low estrogen','OC'},...
        'fontsize', 20, 'fontname', 'Calibri');
    title([sprintf('%s', namecontrast)]);
    xlabel("Cohen's d per parcel");
else   
    f  = figure('Position', fig_position);
    h   = rm_raincloud(data, cb(1,:));
    set(gca, 'XLim', [-8, 8], 'YTickLabel',...
        {'All women', 'High estrogen','High progesterone','Low progesterone','Low estrogen','OC'},...
        'fontsize', 20, 'fontname', 'Calibri');
    title([sprintf('%s', namecontrast)]);
    subtitle('in men compared to women grouped by hormonal status');
    xlabel('t-values per parcel');
end

xline(0, '--');

% change one subset to new colour and alter dot size
% h.l(i,j) is the handle for the line connecting h.m(i,j) and h.m(i+1,j)

h.p{1, 1}.FaceColor         = cb(5,:);
h.s{1, 1}.MarkerFaceColor   = cb(5,:);
h.m(1, 1).MarkerFaceColor   = cb(5,:);

h.p{2, 1}.FaceColor         = cb(4,:);
h.s{2, 1}.MarkerFaceColor   = cb(4,:);
h.m(2, 1).MarkerFaceColor   = cb(4,:);

h.p{3, 1}.FaceColor         = cb(3,:);
h.s{3, 1}.MarkerFaceColor   = cb(3,:);
h.m(3, 1).MarkerFaceColor   = cb(3,:);

h.p{4, 1}.FaceColor         = cb(2,:);
h.s{4, 1}.MarkerFaceColor   = cb(2,:);
h.m(4, 1).MarkerFaceColor   = cb(2,:);

h.p{5, 1}.FaceColor         = cb(1,:);
h.s{5, 1}.MarkerFaceColor   = cb(1,:);
h.m(5, 1).MarkerFaceColor   = cb(1,:);

h.p{6, 1}.FaceColor         = cb(10,:);
h.s{6, 1}.MarkerFaceColor   = cb(10,:);
h.m(6, 1).MarkerFaceColor   = cb(10,:);

%% Save all results

if saveall == 1
    disp(' ...saving results')
    save(fullfile(outDir, 'hormonal_effects_moments.mat'));
    %save((fullfile(outDir, 'hormonalstuff.mat'), *addhere whatever I want to save*));    
    FigList = findobj(allchild(0), 'flat','Type','figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        savefig(FigHandle, fullfile(figDir, string(FigName(iFig))));
    end
end

disp(' ...done!:)')