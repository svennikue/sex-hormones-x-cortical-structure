% This script analysis if microstructural moments changes as a function of
% cortical types. 
% firstly, it correlates sex difference tval maps with different cortical
% types. Second, it plots tvals (or d-values) per cortical type. Third, it compares sex 
% differences with BigBrain microstructural gradients

clear all
close all

plotfigs = 1;
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
figDir = fullfile(homeDir, 'figures/cortical_type_figs');
addpath(figDir);

count = 0;


%% Analyses.

load HCP_T1wT2w_sexdiffsmaps.mat
load Cortical_types.mat; %parcellation Paquola/Economo
parc400 = types400; % Parcellation = Cortical Types
cl = colortabler.table(2:7,1:3)./256; % colorsets

% prepare BigBrain gradients
profiles = read_histology_profile('template', 'fsaverage5');
schaefer_400 = fetch_parcellation('fsaverage5', 'schaefer', 400);
mpc = compute_mpc(profiles, schaefer_400);
gm_bb = compute_histology_gradients(mpc);
[surf_left, surf_right] = fetch_template_surface('fsaverage5');

% plot the 2 BigBrain gradients 
vertexwise_gradients = parcel2full(gm_bb.gradients{1}(:,1:2), schaefer_400);
vertexwise_gradients(isnan(vertexwise_gradients)) = inf; 
set(gca, 'fontname', 'Calibri');
obj = plot_hemispheres(vertexwise_gradients, ...
    {surf_left, surf_right}, ...
    'labeltext', {'BB Gradient 1', 'BB Gradient 2'});

colormap((cbrewer('seq','YlGn',11)))
%obj.colormaps([parula; .7 .7 .7])

%close all 


for measure = 1:3      
    if measure == 1                
        namemoment = 'gradient'
        else if measure == 2
            namemoment = 'profile mean'
            else if measure == 3
                namemoment = 'profile skewness'
                end
            end
    end
    
    % this has been calculated in sex_diffs_3measures.m
    d_vals(:, measure) = results.Cohensd(:,measure)';
    absolute_dvals(:, measure) = abs(d_vals(:, measure));
    
    
    % PT 1: Correlation between cortical types and tvals.
     %correlation between cortical types and tvals (difference male-female)
     
    % correlation between cortical types and tvals (difference male-female)
    % need to exclude type 1 here: cortical wall
    % exclude = find(parc400 ~= 1);
%   [rho,pval] = corr(t_schaefer_400(exclude, measure), parc400(exclude)', 'Type', 'Spearman');

    
    % instead of deleting exclude, put cortical wall to nan
    exclude = find(parc400 == 1);
    parc400(exclude) == NaN;
    d_vals(exclude, measure) = NaN;

    % check one thing! What if I also exlclude koniocortex for mean - will
    % there be a stronger correlation?
    % ok nevermind, the correlation would be the wrong way around anyways
%     exclude_konio = find(parc400 == 2);
%     parc400(exclude_konio) == NaN;
%     d_vals(exclude_konio,measure) = NaN;
%     [rho,pval] = corr(d_vals(:, measure), parc400', 'Type', 'Spearman', 'rows', 'complete');
% 

    
    % also correlate mean group map with cortical types
    % THIS IS NOT WHAT I WANT TO LOOK AT!
    % [rho,pval] = corr(mean(T1T2moments(:,exclude,measure))', parc400(exclude)', 'Type', 'Spearman')
    [rho,pval] = corr(d_vals(:, measure), parc400', 'Type', 'Spearman', 'rows', 'complete');

    % spin permutation testing
    [p_spin, r_dist] = spin_test(d_vals(:, measure), parc400', 'surface_name',...
        'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
        'type', 'Spearman');

    % plot null distributions, correlation, and spin-p vals    
    f = figure,
    his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
        'facecolor', cl(1,:), 'facealpha', 1, 'linewidth', 0.5);
    l = line(rho, 4.5,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
    set (gca, 'fontsize', 50, 'fontname', 'Calibri')
    xlabel(['Null correlations' newline (sprintf('(%s)', namemoment))])
    xlim([-0.5, 0.5]);
    legend(l,['{\it r}=' num2str(round(rho, 2)) newline ...
                      '{\it p}=' num2str(round(p_spin,3 ))], 'Location', 'northwest')
                  
    count = count + 1;
    FigName{(count)} = sprintf('/HCP_Tw1T2w_corticalType_Corr_Tval_spin_%s.png', namemoment);
    
    %store results
    results.corrs(measure,:) = rho;
    results.pvals(measure,:) = pval;
    results.spin(measure,:) = p_spin;
    
    % PT 2: PLOT TVALS PER CORTICAL TYPE
    %types400 = types400';
    
    % Rainbowplot: females vs. males, first gradient 
% Rainbowplot types for moment x sexdiff
%    to_brain_m = t_schaefer_400(:, measure);
    to_brain_m = d_vals(:,measure);
    fig_position = [200 200 600 400]; % size of the figures
    d = [];
    change_p.cw = to_brain_m(find(types400==1));
    change_p.ko = to_brain_m(find(types400==2));
    change_p.eu3 = to_brain_m(find(types400==3));
    change_p.eu2 = to_brain_m(find(types400==4));
    change_p.eu1 = to_brain_m(find(types400==5));
    change_p.dys = to_brain_m(find(types400==6));
    change_p.ag = to_brain_m(find(types400==7));
    d{1} = change_p.cw';
    d{2} = change_p.ko';
    d{3} = change_p.eu3';
    d{4} = change_p.eu2';
    d{5} = change_p.eu1';
    d{6} = change_p.dys';
    d{7} = change_p.ag';
    means = cellfun(@mean, d);
    variances = cellfun(@std, d);
    f = figure('Position', fig_position);        
    %h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
    %    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
    h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 0.5, 'dot_dodge_amount', 0.5, 'box_col_match', 0);
    h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 1, 'dot_dodge_amount', 1, 'box_col_match', 0);
    h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 1.5, 'dot_dodge_amount', 1.5, 'box_col_match', 0);
    h5 = raincloud_plot(d{5}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 2, 'dot_dodge_amount', 2, 'box_col_match', 0);
    h6 = raincloud_plot(d{6}, 'box_on', 1, 'color', cl(5,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 2.5, 'dot_dodge_amount', 2.5, 'box_col_match', 0);
    h7 = raincloud_plot(d{7}, 'box_on', 1, 'color', cl(6,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', 3, 'dot_dodge_amount', 3, 'box_col_match', 0);
    
    set(gca,'XLim', [-1 1], 'color', 'none', 'fontsize', 26, 'fontname', 'Calibri', 'ytick', []);
    title([sprintf('%s', namemoment)]);
    xlabel(['Cohens d']);
    box off
    exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_fem-mal_rainclouds_%s', namemoment)],'png', 10)
    count = count + 1;
    FigName{(count)} = sprintf('/HCP_Tw1T2w_corticalType_sexdiff_raincloud_%s.png', namemoment);
    
end

keyboard

%     h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
%     h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cl(2,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
%     h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cl(3,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
%     h5 = raincloud_plot(d{5}, 'box_on', 1, 'color', cl(4,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
%     h6 = raincloud_plot(d{6}, 'box_on', 1, 'color', cl(5,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', 0.95, 'dot_dodge_amount', 0.95, 'box_col_match', 0);
%     h7 = raincloud_plot(d{7}, 'box_on', 1, 'color', cl(6,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', 1.15, 'dot_dodge_amount', 1.15, 'box_col_match', 0);
%% PT 3: BigBrain Gradient Decoding
    % correlation with BigBrain gradients
    % Histological decoding.
for measure = 1:3      
    if measure == 1                
        namemoment = 'gradient'
        else if measure == 2
            namemoment = 'profile mean'
            else if measure == 3
                namemoment = 'profile skewness'
                end
            end
    end
    % this has been calculated in sex_diffs_3measures.m
%     t_schaefer_400(:, measure) = results.tvals(measure,:)';
%     t_schaefer_400 = abs(t_schaefer_400);
%     
% 
%     % BigBrain Gradient 1
%     t_bb_corr1 = corr(t_schaefer_400(:, measure), gm_bb.gradients{1}(:,1));
    
% this is just the correlation with the mean microstructural measure maps!
    t_bb_corr1 = corr(mean(T1T2moments(:,:,measure))', gm_bb.gradients{1}(:,1))
    t_bb_corr2 = corr(mean(T1T2moments(:,:,measure))', gm_bb.gradients{1}(:,2))
% this would be the correlation with the sex differences:
d_vals(:, measure) = results.Cohensd(:,measure)';
[rho,pval] = corr(d_vals(:, measure), gm_bb.gradients{1}(:,1)', 'Type', 'Spearman', 'rows', 'complete');
[rho,pval] = corr(d_vals(:, measure), gm_bb.gradients{2}(:,2)', 'Type', 'Spearman', 'rows', 'complete');

    
%     % spin test 
%     [p_spin, r_dist] = spin_test(t_schaefer_400(:, measure), gm_bb.gradients{1}(:,1), 'surface_name',...
%     'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
%     'type', 'pearson');
% 
%     % plot null distributions, correlation, and spin-p vals    
%     f = figure
%     his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
%                        'facecolor', cl(5,:), 'facealpha', 1, 'linewidth', 0.5);
%     l = line(t_bb_corr1, 1.7,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
%     xlabel(['Null correlations' newline (sprintf('(%s and BigBrain Gradient1)', namemoment))])
%     xlim([-0.55, 0.55]);
%     legend(l,['{\it r}=' num2str(round(t_bb_corr1, 2)) newline ...
%                       '{\it p}=' num2str(round(p_spin,3 ))], 'Location', 'northwest')
%     set (gca, 'fontsize', 20, 'fontname', 'Calibri');
%     exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_sexdiff_BB1_spin_%s.png', namemoment)],'png', 10)
%     count = count + 1;
%     FigName{(count)} = sprintf('/HCP_Tw1T2w_sexdiff_BB1_spin_%s.png', namemoment);
%     
%     % store results
%     results.corrBB1(measure,:) = t_bb_corr1; % save BB corr per measure
%     results.corrBB1spin = {p_spin, r_dist};
%     
%     
%     % plot with scatter plot
%     bestfit_coefficients1 = polyfit(t_schaefer_400(:,measure), gm_bb.gradients{1}(:,1), 1);
%     figure();
%     scatter(t_schaefer_400(:,measure), gm_bb.gradients{1}(:,1), 200, 'k', '.');
%     hold on
%     plot (t_schaefer_400(:,measure), polyval(bestfit_coefficients1, t_schaefer_400(:,measure)), 'k', 'LineWidth', 3);
%     xlabel(sprintf('t-values (%s)', namemoment));
%     ylabel('BigBrain Gradient 1')
%     set(gca, 'FontSize', 20, 'fontname', 'Calibri', 'YLim', [-1, 1], 'YTick', [-1, 1], 'Xtick', [0, 6], 'XLim', [0, 6])
%     text(0.2, 0.8, sprintf('r = %0.2f', t_bb_corr1), 'FontSize', 16);
%     exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_sexdiff_BB1_scatter_%s.png', namemoment)],'png', 10)
%     count = count + 1;
%     FigName{(count)} = sprintf('/HCP_Tw1T2w_sexdiff_BB1_scatter_%s.png', namemoment);
%     
% 
%     % BigBrain Gradient 2
%     t_bb_corr2 = corr(t_schaefer_400(:,measure), gm_bb.gradients{1}(:,2));
%     %corr_bbgrad2(:,measure) = t_bb_corr2; % save BB corr per measure    
%     % spin test 
%     [p_spin, r_dist] = spin_test(t_schaefer_400(:,measure), gm_bb.gradients{1}(:,2), 'surface_name',...
%     'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
%     'type', 'pearson');
%     % plot null distributions, correlation, and spin-p vals    
%     f = figure
%     his = histogram(r_dist, 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
%                        'facecolor', cl(5,:), 'facealpha', 1, 'linewidth', 0.5);
%     l = line(t_bb_corr2, 1.7 ,'Color', 'k','LineStyle', '-', 'Marker', 'o', 'MarkerFaceColor', 'k');
%     xlabel(['Null correlations' newline (sprintf('(%s and BigBrain Gradient2)', namemoment))])
%     xlim([-0.55, 0.55]);
%     legend(l,['{\it r}=' num2str(round(t_bb_corr2, 2)) newline ...
%                       '{\it p}=' num2str(round(p_spin,3 ))])
%     set (gca, 'fontsize', 20, 'fontname', 'Calibri');
%     exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_sexdiff_BB2_spin_%s.png', namemoment)],'png', 10)
%     count = count + 1;
%     FigName{(count)} = sprintf('/HCP_Tw1T2w_sexdiff_BB2_spin_%s.png', namemoment);
%     
%     
%     % store results
%     results.corrBB2(measure,:) = t_bb_corr2; % save BB corr per measure
%     results.corrBB2spin = {p_spin, r_dist};
%     
%     % plot with scatter plot        
%     bestfit_coefficients2 = polyfit(t_schaefer_400(:,measure), gm_bb.gradients{1}(:,2), 1);
%     figure();
%     scatter(t_schaefer_400(:,measure), gm_bb.gradients{1}(:,2), 200, 'k', '.');
%     hold on
%     plot (t_schaefer_400(:,measure), polyval(bestfit_coefficients2, t_schaefer_400(:,measure)), 'k', 'LineWidth', 3);
%     xlabel(sprintf('t-values (%s)', namemoment));
%     ylabel('BigBrain Gradient 2')
%     set(gca, 'FontSize', 20, 'fontname', 'Calibri', 'YLim', [-1, 1], 'YTick', [-1, 1], 'Xtick', [0 6], 'XLim', [0 6])
%     text(0.2, 0.8, sprintf('r = %0.2f', t_bb_corr2), 'FontSize', 16);
%     exportfigbo(f,[figDir, sprintf('/HCP_Tw1T2w_sexdiff_BB2_scatter_%s.png', namemoment)],'png', 10)
%     count = count + 1;
%     FigName{(count)} = sprintf('/HCP_Tw1T2w_sexdiff_BB2_scatter_%s.png', namemoment);
%     

end


%% Save all results

disp(' ...saving results')
if saveall == 1
    save(fullfile(outDir, 'cortical_types_sexdiff_results.mat'));
    %save((fullfile(outDir, 'hormonalstuff.mat'), *addhere whatever I want to save*));    
    FigList = findobj(allchild(0), 'flat','Type','figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        savefig(FigHandle, fullfile(figDir, string(FigName(iFig))));
    end
end

disp(' ...done!:)')