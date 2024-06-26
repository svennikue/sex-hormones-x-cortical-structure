% Gradient computation for HCP data
% requires: female_gradients.mat (= 1206 HCP dataset in T1T2w profiles)
% Functions: schaefer parcellation toolbox, exportfigbo, brainspace toolbox
% First step for main analysis: needed for Figure 1,2

count = 0;
parcellation_glasser = true;

% make a path to directory
addpath(genpath('/Users/skuech/Documents/toolboxes'));
homeDir = '/Users/skuech/Documents/my_projects/female_gradients/';
scriptDir = fullfile(homeDir,'script');
addpath(genpath(scriptDir));
dataDir = fullfile(homeDir, 'data');
addpath(genpath(dataDir));
outDir = fullfile(homeDir,'output');
addpath(genpath(outDir));
figDir = fullfile(homeDir, 'figures');
addpath(figDir);


load('female_gradients.mat'); %1206 subjects
SN = female_gradients.SN;

if parcellation_glasser == true
    load('MPC_all_retest2.mat')
    MPCnx1 = MPC_all_retest(:, [2:181, 183:362], [2:181, 183:362]);
    dataname = 'HCP_glasser';
    load("parcels360.mat")
    glasser = (fetch_parcellation('fsaverage5', 'glasser', 360))'; 
else
    dataname = 'HCP_T1wT2w';
    MPCnx1 = female_gradients.MPCnx1; %basis for gradient
    % load schaefer400 parcellation
    % careful! Now all parcels400 + 1 is schaefer_400
    schaefer_400 = (fetch_parcellation('fsaverage5', 'schaefer', 400))';
end


%make keep to keep those with complete data
keep = find(squeeze(mean(MPCnx1(:,1,1:size(MPCnx1,2)),3))>0);

%make mean matrixes
meanMPC = squeeze(mean(MPCnx1(keep,:,:)));
meanMPC(eye(size(meanMPC))==1) = 0;

% plot mean
f = figure 
colormap(bone)
imagesc(meanMPC, [-0.5,3])
colorbar
set (gca, 'fontsize', 20, 'fontname', 'Calibri')
exportfigbo(f,[figDir, sprintf('%mean_MPC.png', dataname)],'png', 10)
count = count + 1
FigName{(count)} = sprintf('%mean_MPC.fig', dataname);

%% ANALYSIS 1.1.1: SEX COMPARISON, MPC GRADIENT computation  
% Calculate gradient 
% alignment of MPC to microstructural mean MPC data
c1_tx =  zeros(size(MPCnx1, 1),size(MPCnx1,2)); %1206 subjects x 400 parcels
for i = 1:size(MPCnx1, 1) %for each subject
    i
    try
        gm = GradientMaps('kernel','na','approach','dm','align','pa'); % initializes gradient computation
        gm = gm.fit({meanMPC,squeeze(MPCnx1(i,:,:))}); 
        % computes gradients of all provided matrices (all n) and aligns
        % the gradient of the current person to the mean MPC gradient
        % input are microstructuralcovariance matrices

        fc2mpc_t  = gm.aligned{2}(:,1); %first gradient of person i, procrustes 
        % aligned to mean MPC matrix. gm.aligned{1} would be the gradient of 
        % mean MPC matrix (first element provided in gm.fit)
        
        % check how the results are changed if I don't rescale.
        % c1_tx(i,:) = -fc2mpc_t;
        c1_tx(i,:) = rescale(-fc2mpc_t); % rescale aligned individual gradient 0-1
    catch
        % for those subjects that create errors, stays at 0
    end
end

if parcellation_glasser == false
    % plot a sorted MPC for one subject
    [B, I] = sort(mean(c1_tx));
    for i = 1:size(MPCnx1,2)
        x = I(i);
        for j = 1:size(MPCnx1,2)
            y = I(j);
            sorted_M(i,j) = MPCnx1(size(MPCnx1, 1),x,y);
        end
    end
    
    f = figure 
    colormap(bone)
    imagesc(sorted_M, [-0.5,3])
    colorbar
    set (gca, 'fontsize', 20, 'fontname', 'Calibri')
    exportfigbo(f,[figDir, sprintf('%s_orderedMPC.png', dataname)],'png', 10)
    count = count + 1
    FigName{(count)} = sprintf('%s_orderedMPC.fig', dataname);
    
    
    %mean MPC gradient
    for i = 1:size(MPCnx1,2) %200 parcels per hemisphere
        vertices(find(schaefer_400==i)) = mean(c1_tx(keep,i));
    end
    
    
    % plot the mean MPC gradient for all subjects 
    f = figure,
    BoSurfStatViewData(vertices,SN,'')
    colormap((cbrewer('seq','Greens',11)))
    SurfStatColLim([0 1])
    exportfigbo(f,[figDir, sprintf('%s_MPC_G1.png', dataname)],'png', 10)
    count = count + 1
    FigName{(count)} = sprintf('%s_MPC_G1.fig', dataname);
end

if parcellation_glasser == true
    to_surface = zeros(20484,1);
    % 2:181 + 183:362
    for i = 1:180 % to exclude medial wall
      to_surface(parcels360x==i+1) = mean(c1_tx(keep,i));
    end
    for i = 1:180 % to exclude medial wall
      to_surface(parcels360x==i+182) = mean(c1_tx(keep,i+180));
    end

    f = figure,
    BoSurfStatViewData(to_surface,SN,'')
    colormap((cbrewer('seq','Greens',11)))
    SurfStatColLim([0 1])
    exportfigbo(f,[figDir, sprintf('%s_MPC_G1.png', dataname)],'png', 10)
    FigName{(count)} = sprintf('%s_MPC_G1.fig', dataname);
end


gradient.grad1persubj = {c1_tx};
gradient.gm = {gm};




%% Save all results
disp(' ...saving results')

save(fullfile(outDir, sprintf('%s_gradient.mat', dataname)), 'gradient');
FigList = findobj(allchild(0), 'flat','Type','figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    savefig(FigHandle, fullfile(figDir, string(FigName(iFig))));
end


disp(' ...done!:)')
