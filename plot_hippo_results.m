%%%%%% THIS SCRIPT ALLOWS TO PLOT ALL HIPPOCAMPAL FIGURES
%% SETTINGS:
% HEMISPHERES, SUBFIELDS, GROUP-AVE MAPS, SEXDIFFERENCES, HORMONAL
% CONTRASTS, WOMEN ONLY HORMONES
clear all

rh = 1;
lh = 0;


name = 'all_sexdiffs';

subfields = 0;
groupmean = 0;
sexdiff = 1;
low_estr = 0;
low_prog = 0;
high_estr = 0;
high_prog = 0;
OC = 0;

OCNC = 0;
highlowprog = 0;
highlowestr = 0;
FDR_corr = 0;


addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';

figDir = fullfile(homeDir, 'figures/hormonal_figs');
addpath(figDir);

%%
% Left hemisphere

if lh == 1
    if OCNC == 1
       data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_women-OCvsNC.txt');
       if FDR_corr == 0
           data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_Cohensd_women-OCvsNC.txt');
       end
    end

    if highlowprog == 1
       data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_women-highvslow-progesterone.txt');
       if FDR_corr == 0
           data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_Cohensd_women-highvslow-progesterone.txt');
       end
    end

    if highlowestr == 1
       data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_women-highvslow-estrogen.txt');
       if FDR_corr == 0
           data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_Cohensd_women-highvslow-estrogen.txt');
       end
    end    

    if sexdiff == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/hipposexFDRCohensd_lh.txt');
        %load /Users/skuech/Documents/my_projects/female_gradients/output/results_hippocampus_lh_sexdiffs.mat;
        % data = (results.CohensdFDR)';
    end

    if groupmean == 1
        groupmean_file = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_L_space-T1w_den-0p5mm_T1wT2w.shape.gii');
        data = groupmean_file.cdata;
    end

    if low_estr == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_Menvslow-estr.txt');
    end

    if low_prog == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_Menvslow-prog.txt');
    end

    if high_estr == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_Menvshigh-estr.txt');
    end


    if high_prog == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_Menvshigh-prog.txt');
    end

    if OC == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/lh_hippo_dFDR_MenvsOC-women.txt');
    end


       
    % JORDAN 
    % FOLDED
    
    window = true;
    smooth = 0;
    g = [];
    g.XLim = [0 1];
    figure; hold on;
    
    % addpath('/data/p_02542/dekrak/scripts/gifti-master')
    addpath(genpath('/Users/skuech/Documents/toolboxes'));
    % subfield labels
    % subfs = gifti('/data/p_02542/dekrak/data/groupAve/average_R_space-T1w_den-0p5mm_label-hipp_subfields.label.gii');
    subfs = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_L_space-T1w_den-0p5mm_label-hipp_subfields.label.gii');
    
  
    % rotate surface to corobl
    rot = eye(4);
    
    % surface left hipp 
    slh = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii');
    v = [slh.vertices, zeros(length(slh.vertices),1)];
    v = rot'*v';
    slh.vertices = v(1:3,:)';
    tuplh = mean(slh.vertices(:,2));
    slh.vertices(:,2) = slh.vertices(:,2) -tuplh; % align middle
    
    if subfields == 1
        lh = plot_gifti(slh,subfs.cdata);
        colormap(subfs.labels.rgba(2:end,1:3));
        clim([min(subfs.cdata), max(subfs.cdata)]);

    elseif sexdiff == 1 || low_prog == 1|| low_estr == 1 || high_estr == 1 || OC == 1 || high_prog == 1
        % if with data
        lh = plot_gifti(slh,data);
        colormap(flipud(cbrewer('div','RdBu',11)));
        clim([-1,1]);

    elseif OCNC == 1 ||highlowprog == 1 || highlowestr == 1
        lh = plot_gifti(slh,data);
        colormap(flipud(cbrewer('div','PuOr',11)));
        clim([-1,1]);

    elseif groupmean == 1
        lh = plot_gifti(slh,data);
        colormap(cbrewer('seq','Greens',11));
        clim([(prctile(data,5)) (prctile(data,95))]);
    end
    
    
    
    
    % draw borders
    slh = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii');
    v = [slh.vertices, zeros(length(slh.vertices),1)];
    v = rot'*v';
    slh.vertices = v(1:3,:)';
    tuplh = mean(slh.vertices(:,2));
    slh.vertices(:,2) = slh.vertices(:,2) -tuplh; % align middle
    p = plot_giftibordersSofie(slh, subfs.cdata);
    %[p, A] = plot_giftiborders(slh, subfs.cdata);
    
    
    % UNFOLDED 
    subfs = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_L_space-T1w_den-0p5mm_label-hipp_subfields.label.gii');
    sluh = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_L_space-unfolded_den-0p5mm_label-hipp_midthickness.surf.gii');
    sluh.vertices(:,[1 2 3]) = sluh.vertices(:,[2 1 3]);
    sluh.vertices(:,1) = -sluh.vertices(:,1); % flip left
    g = gca;
    sluh.vertices(:,1) = sluh.vertices(:,1) + (g.XLim(1)-max(sluh.vertices(:,1))) -1; % translate to left edge
    sluh.vertices(:,2) = sluh.vertices(:,2) -mean(sluh.vertices(:,2)); % align middle
    
    % if with data
    if subfields == 1
        luh = plot_gifti(sluh, subfs.cdata);
        colormap(subfs.labels.rgba(2:end,1:3));
        clim([min(subfs.cdata), max(subfs.cdata)]);
    elseif sexdiff == 1 || low_prog == 1|| low_estr == 1 || high_estr == 1 || OC == 1 || high_prog == 1
        luh = plot_gifti(sluh, data);
        colormap(flipud(cbrewer('div','RdBu',11)));
        clim([-1,1])
    elseif OCNC == 1 ||highlowprog == 1 || highlowestr == 1
        lh = plot_gifti(sluh,data);
        colormap(flipud(cbrewer('div','PuOr',11)));
        clim([-1,1]);
    elseif groupmean == 1
        lh = plot_gifti(sluh,data);
        colormap(cbrewer('seq','Greens',11));
        clim([(prctile(data,5)) (prctile(data,95))]);
    end
    
    
    
    % UNFOLDED BORDERS
    p = plot_giftiborders(sluh, subfs.cdata);
    
    colorbar;
    light;
    
    figDir = fullfile(homeDir, 'figures/hormonal_figs');
    imagewd = getframe(gcf); 
    imwrite(imagewd.cdata, fullfile(figDir, sprintf('/lh_%s.tiff', name)));
    imagewd = getframe(gcf); 
    imwrite(imagewd.cdata, fullfile(figDir, sprintf('/lh_%s.png', name)));

    close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RIGHT HEMISPHERE

if rh == 1
    if OCNC == 1
       data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_women-OCvsNC.txt');
    end

    if highlowprog == 1
       data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_women-highvslow-progesterone.txt');
    end

    if highlowestr == 1
       data = load ('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_women-highvslow-estrogen.txt');
    end  


    if sexdiff == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/hipposexFDRCohensd_rh.txt');
        % data = (results.CohensdFDR)';
    end

    if groupmean == 1
        groupmean_file = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_R_space-T1w_den-0p5mm_T1wT2w.shape.gii');
        data = groupmean_file.cdata;
    end

    if low_estr == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_Menvslow-estr.txt');
    end

    if low_prog == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_Menvslow-prog.txt');
    end

    if high_estr == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_Menvshigh-estr.txt');
    end

        
    if high_estr == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_Menvshigh-estr.txt');
    end


    if high_prog == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_Menvshigh-prog.txt');
    end

    if OC == 1
        data = load('/Users/skuech/Documents/my_projects/female_gradients/output/hippocampus/rh_hippo_dFDR_MenvsOC-women.txt');
    end

    
    % JORDAN 
    % FOLDED
    
    window = true;
    smooth = 0;
    g = [];
    g.XLim = [0 1];
    figure; hold on;
    
    % addpath('/data/p_02542/dekrak/scripts/gifti-master')
    addpath(genpath('/Users/skuech/Documents/toolboxes'));
    % subfield labels
    % subfs = gifti('/data/p_02542/dekrak/data/groupAve/average_R_space-T1w_den-0p5mm_label-hipp_subfields.label.gii');
    subfs = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_R_space-T1w_den-0p5mm_label-hipp_subfields.label.gii');
    
    
    % rotate surface to corobl
    rot = eye(4);
    
    % surface left hipp 
    slh = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii');
    v = [slh.vertices, zeros(length(slh.vertices),1)];
    v = rot'*v';
    slh.vertices = v(1:3,:)';
    tuplh = mean(slh.vertices(:,2));
    slh.vertices(:,2) = slh.vertices(:,2) -tuplh; % align middle
    
    if subfields == 1
        lh = plot_gifti(slh,subfs.cdata);
        colormap(subfs.labels.rgba(2:end,1:3))
        clim([min(subfs.cdata), max(subfs.cdata)])
    elseif sexdiff == 1 || low_prog == 1|| low_estr == 1 || high_estr == 1 || OC == 1 || high_prog == 1
        % if with data
        lh = plot_gifti(slh,data);
        colormap(flipud(cbrewer('div','RdBu',11)))
        clim([-1,1])

    elseif OCNC == 1 ||highlowprog == 1 || highlowestr == 1
        lh = plot_gifti(slh,data);
        colormap(flipud(cbrewer('div','PuOr',11)));
        clim([-1,1]);

    elseif groupmean == 1
        lh = plot_gifti(slh,data);
        colormap(cbrewer('seq','Greens',11));
        clim([(prctile(data,5)) (prctile(data,95))]);
    end
    
    
    
    
    % draw borders
    slh = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii');
    v = [slh.vertices, zeros(length(slh.vertices),1)];
    v = rot'*v';
    slh.vertices = v(1:3,:)';
    tuplh = mean(slh.vertices(:,2));
    slh.vertices(:,2) = slh.vertices(:,2) -tuplh; % align middle
    p = plot_giftibordersSofie(slh, subfs.cdata);
    %[p, A] = plot_giftiborders(slh, subfs.cdata);
    
    
    % UNFOLDED 
    subfs = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_R_space-T1w_den-0p5mm_label-hipp_subfields.label.gii');
    sluh = gifti('/Users/skuech/Documents/my_projects/female_gradients/script/Hippunfold_downsample/Sofie_plotting_hippo/average_R_space-unfolded_den-0p5mm_label-hipp_midthickness.surf.gii');
    sluh.vertices(:,[1 2 3]) = sluh.vertices(:,[2 1 3]);
    % sluh.vertices(:,1) = -sluh.vertices(:,1); % flip left
    g = gca;
    sluh.vertices(:,1) = sluh.vertices(:,1) + (g.XLim(1)-max(sluh.vertices(:,1))) -1; % translate to left edge
    sluh.vertices(:,2) = sluh.vertices(:,2) -mean(sluh.vertices(:,2)); % align middle
    
    if subfields == 1
        colormap(subfs.labels.rgba(2:end,1:3));
        clim([min(subfs.cdata), max(subfs.cdata)]);
        luh = plot_gifti(sluh, subfs.cdata);
    elseif sexdiff == 1 || low_prog == 1|| low_estr == 1 || high_estr == 1 || OC == 1 || high_prog == 1
        % if with data
        luh = plot_gifti(sluh, data);
        colormap(flipud(cbrewer('div','RdBu',11)))
        clim([-1,1])
    elseif OCNC == 1 ||highlowprog == 1 || highlowestr == 1
        lh = plot_gifti(sluh,data);
        colormap(flipud(cbrewer('div','PuOr',11)));
        clim([-1,1]);
    elseif groupmean == 1
        lh = plot_gifti(sluh,data);
        colormap(cbrewer('seq','Greens',11));
        clim([(prctile(data,5)) (prctile(data,95))]);
    end
    
    
    
    % UNFOLDED BORDERS
    p = plot_giftiborders(sluh, subfs.cdata);
    
    colorbar;
    light;

    figDir = fullfile(homeDir, 'figures/hormonal_figs');
    imagewd = getframe(gcf); 
    imwrite(imagewd.cdata, fullfile(figDir, sprintf('/rh_%s.tiff', name)));
    imwrite(imagewd.cdata, fullfile(figDir, sprintf('/rh_%s.png', name)));
end
