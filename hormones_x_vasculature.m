% This script analysis is the changes due to hormones relate to 
% vasculature.
% input: vasculature maps, sex difference results, hormonal results.
% output: spatial correlations, spin-corrected pvalues and spin-rho
% distributions of correlation between effect maps and vasculature maps
% Supplement 11

loadold = 1;
saveall = 0;

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
    
    
if loadold == 0
    clear all
    close all
    loadold = 0;
    saveall = 1;
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

    load 'hormonal_effects_moments.mat'

    % the important variables are 
    % tvals_profilemean, tvals_profileskew, tvals_gradient
    % model 3 = OC women, model 6 = high prog, model 7 = low prog
    % model 8 = high estr, model 9 = low estr

    % load vascular maps
    D = fullfile(dataDir, 'vascular_stuff/');

    map = tvals_profileskew{8};
    atlas_a = readmatrix([D 'mean_Ved_ToF_ThreshtoMNIlin_400.csv']);
    atlas_v = readmatrix([D 'mean_Ved_swi_ThreshtoMNIlin_400.csv']);
    atlas = {atlas_a,atlas_v};
    modelnumber = 9;
    load HCP_T1wT2w_sexdiffsmaps.mat
    mainresults = results;

else if loadold == 1
        load vascular_x_hormones.mat
        keyboard
    end
end


%% run analysis and plot
% close all

for a = 1:2
    vascular_map = atlas{a};
    
    if a == 1
        nameatlas = 'arteries'
    else if a == 2
            nameatlas = 'veins'
        end
    end
    
    
    % plot atlas
    vertices = zeros(20484,1);
    for i = 1:400 %200 parcels per hemisphere
        vertices(find(schaefer_400==i)) = vascular_map(i);
    end      
    f = figure,
    BoSurfStatViewData(vertices,SN,'')
    colormap(cbrewer('seq','Reds',11))
    title(sprintf('Atlas of cerebral %s', nameatlas));
    exportfigbo(f,[figDir, sprintf('Atlas_of_cerebral_%s', nameatlas)],'fig', 10)
    count = count + 1
    FigName{(count)} = sprintf('Atlas_of_cerebral%s.png', nameatlas);
    
    for moment = 1:3
        if moment == 1 
            map_dvals = cohensd_gradient;
            namemoment = 'gradient'
        else if moment == 2
                map_dvals = cohensd_profilemean;
                namemoment = 'profile mean'
            else if moment == 3
                    map_dvals = cohensd_profileskew;
                    namemoment = 'profile skewness'
                end
            end
        end
       
        %first for main effects.
        map = mainresults.Cohensd(:,moment);
        [rho, pval] = corr(map,vascular_map, 'type', 'spearman');
%            if abs(rho) > .10
                %spin test 
                [p_spin, r_dist] = spin_test(map', vascular_map, 'surface_name',...
                'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
                'type', 'spearman');

                % plot null distributions, correlation, and spin-p vals 
                f = figure
                his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                                   'facealpha', 1, 'linewidth', 0.5);
                l = line(rho, 5,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
                xlabel(['Null correlations' newline (sprintf('(sex diffs between %s and vascular map %s)', namemoment, nameatlas))])
                legend(l,['{\it r}=' num2str(round(rho, 2)) newline ...
                                  '{\it p}=' num2str(round(p_spin,3 ))])
                count = count + 1
                FigName{(count)} = sprintf('/HCP_Tw1T2w_sterrec_gene_spin_%s%d.png', namemoment, i);

         % save values
         sexdiffresults.corr_atlas(a, moment, :) = rho; % save corr per model and measure
         sexdiffresults.pspin_atlas(a, moment, :) = p_spin; % save pval per model and measure
         sexdiffresults.pval(a,moment,:) = pval; % save the initial pval
    
                
        % then continue with the hormonal subgroups
        for currmodel = 1:modelnumber
            map = map_dvals{currmodel};
            % correlate stat map and vasculature
            [rho, pval] = corr(map',vascular_map, 'type', 'spearman')
            %if abs(rho) > .10
            %spin test 
            [p_spin, r_dist] = spin_test(map', vascular_map, 'surface_name',...
            'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
            'type', 'spearman');

            % plot null distributions, correlation, and spin-p vals    
            f = figure
            his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                            'facealpha', 1, 'linewidth', 0.5);
            l = line(rho, 1.7,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
            xlabel(['Null correlations' newline (sprintf('(%s and %s for model %d)', nameatlas, namemoment, currmodel))])
            xlim([-0.55, 0.55]);
            legend(l,['{\it r}=' num2str(round(rho, 2)) newline ...
                              '{\it p}=' num2str(round(p_spin,3 ))], 'Location', 'northwest')
            set (gca, 'fontsize', 20, 'fontname', 'Calibri');
            exportfigbo(f,[figDir, sprintf('/Vasc_%s_spin_%s_model%d.fig', nameatlas, namemoment, currmodel)],'png', 10)
            count = count + 1;
            FigName{(count)} = sprintf('/Vasc_%s_spin_%s_model%d.png', nameatlas, namemoment, currmodel);

            % store results
            hormresults.corr_atlas(a, moment, currmodel, :) = rho; % save corr per model and measure
            hormresults.pval_atlas(a, moment, currmodel, :) = p_spin; % save pval per model and measure
            hormresults.rdist_atlas(a, moment, currmodel, :) = rho;
            hormresults.pval_uncorr(a, moment, currmodel, :) = pval;
            %end
        end
    end
end

%% FDR Correction a la Benjamini Hochberg
pValues_horm = reshape(hormresults.pval_atlas, 1, []);
h_horm = fdr_bh(pValues_hrom,0.05); % if h is 1, then significant.
H_horm = reshape(h_horm, 2, 3, 9);

h_sexd = fdr_bh(sexdiffresults.pspin_atlas,0.05); % if h is 1, then significant.

% none of the relevant models (3, 6,7,8,9) or main effects survive FDR 

%% Save all results

disp(' ...saving results')
if saveall == 1
    save(fullfile(outDir, 'vascular_x_hormones.mat'));
    %save((fullfile(outDir, 'hormonalstuff.mat'), *addhere whatever I want to save*));    
    FigList = findobj(allchild(0), 'flat','Type','figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        savefig(FigHandle, fullfile(figDir, string(FigName(iFig))));
    end
end

disp(' ...done!:)')

