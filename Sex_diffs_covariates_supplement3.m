% Supplement 3
% Written by Sofie Valk, 2024

%spiderplot
cd('/Users/sofievalk/Downloads/github_repo/')
% toolboxes used
% brainspace, surfstat, colorbrewer
addpath(genpath('/Users/sofievalk/Documents/GitHub/BrainSpace/matlab/'));

% make a path to directory
mkdir('/Users/sofievalk/Dropbox/my_projects/female_gradients/')
RPATH = '/Users/sofievalk/Dropbox/my_projects/female_gradients/'

cd(RPATH)
load('female_gradients.mat')
%make keep to keep those with complete data
keep_mpc = find(squeeze(mean(female_gradients.MPCnx1(:,1,1:400),3))>0);
keep = keep_mpc;

%make mean matrixes
meanMPC = squeeze(mean(female_gradients.MPCnx1(keep,:,:)));
meanMPC(eye(size(meanMPC))==1) = 0;

 f=figure,
 imagesc(squeeze(mean(mean(MPC_layers1(keep,:,:)),3))',[1.1 2.5])
 colormap((cbrewer('seq','Greys',99)))
 colorbar
 exportfigbo(f,[RPATH 'meant1t2.method.png'],'png', 10)

     
     
     
% alignment of MPC to microstructural mean MPC data similar to function
c1_tx =  zeros(1206,400);
for i = 1:1206
    i
    try
        gm = GradientMaps('kernel','na','approach','dm','align','pa');
        gm = gm.fit({meanMPC,squeeze(female_gradients.MPCnx1(i,:,:))});
        
        fc2mpc_t  = gm.aligned{2}(:,1);
        
        c1_tx(i,:) = rescale(-fc2mpc_t);
    catch
    end
end

for gradient_analysis = 1
 
    genderkeep  = cellstr(unres.Gender(keep));
    agekeep     = D1.Age_in_Yrs(keep);
    familyHCPk  = cellstr(D1.Family_ID(keep));
    twinHCPk    = cellstr(D1.ZygosityGT(keep));
    familyHCP   = term(familyHCPk);
    twinHCP     = term(twinHCPk);
    icvkeep     = unres.FS_IntraCranial_Vol(keep);
    G           = term(genderkeep);
    M           =       term(agekeep) + term(genderkeep) + term(icvkeep) + ( 1 + twinHCP )*( 1 + random(familyHCP) ) + I;

    slm         = SurfStatLinMod(c1_tx(keep,:),M);
    slm         = SurfStatT(slm,G.F-G.M);
    
    % FDR correction females - males
    p   = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h   = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn  = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h   = h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1))    = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'mpcG1.gender.Family_controlled.t.png'],'png', 10)
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1))    = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'mpcG1.gender.Family_controlled.FDR.png'],'png', 10)
   
end
   
for skewness_analysis = 1
    
    for i =1:1206
        skew(i,:) = skewness(squeeze(female_gradients.MPC_layers1(i,2:11,:)));
    end
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = mean(skew(keep,i));
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = mean(skew(keep,i+200));
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap((cbrewer('seq','Greens',99)))
  
    
    
    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    familyHCPk = cellstr(D1.Family_ID(keep));
    familyHCP = term(familyHCPk);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M           =       term(agekeep) + term(genderkeep) + term(icvkeep);
    slm =  SurfStatLinMod(skew(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'skewt1t2.gender.Family_controlled.t.png'],'png', 10)
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'skewt1t2.gender.Family_controlled.FDR.png'],'png', 10)
  
    
end

for meant1wt2w  =  1
          
    for i =1:1206
        meant1t2(i,:) = mean(squeeze(female_gradients.MPC_layers1(i,2:11,:)));
    end
    
    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    familyHCPk = cellstr(D1.Family_ID(keep));
    familyHCP = term(familyHCPk);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M           =       term(agekeep) + term(genderkeep) + term(icvkeep) + ( 1 + twinHCP )*( 1 + random(familyHCP) ) + I;
    slm =  SurfStatLinMod(meant1t2(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'meant1t2.gender.Family_controlled.t.png'],'png', 10)
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'meant1t2.gender.Family_controlled.FDR.png'],'png', 10)
    
end


load('/Users/sofievalk/Dropbox/my_projects/Juelich/2019.StrucCOV/matlab.mat', 'HCP400_CT')

meant1wt2w_ctx = zeros(1206,400);
for j = 1
    measure_used = meant1t2;
    for i = 1:400
        seed = HCP400_CT(keep,i);
        M    = 1 + term(seed);
        slm  = SurfStatLinMod(measure_used(keep,i),M);
        meant1wt2w_ctx(keep,i) = measure_used(keep,i) - slm.X*slm.coef;
    end
end

for meant1wt2wCTX  =  1

    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    familyHCPk = cellstr(D1.Family_ID(keep));
    familyHCP = term(familyHCPk);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M = 1 + term(agekeep) + term(genderkeep) + term(icvkeep);
    slm =  SurfStatLinMod(meant1wt2w_ctx(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'meant1t2.gender.CTX_controlled.t.png'],'png', 10)
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'meant1t2.gender.CTX_controlled.FDR.png'],'png', 10)
    
end

skew_ctx = zeros(1206,400);
for j = 1
    measure_used = skew;
    for i = 1:400
        seed = HCP400_CT(keep,i);
        M    = 1 + term(seed);
        slm  = SurfStatLinMod(measure_used(keep,i),M);
        skew_ctx(keep,i) = measure_used(keep,i) - slm.X*slm.coef;
    end
end

for skewCTX  =  1

    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    familyHCPk = cellstr(D1.Family_ID(keep));
    familyHCP = term(familyHCPk);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M = 1 + term(agekeep) + term(genderkeep) + term(icvkeep);
    slm =  SurfStatLinMod(skew_ctx(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'skew.gender.CTX_controlled.t.png'],'png', 10)
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'skew.CTX_controlled.FDR.png'],'png', 10)
    
end


cttx_ctx = zeros(1206,400);
for j = 1
    measure_used = c1_tx;
    for i = 1:400
        seed = HCP400_CT(keep,i);
        M    = 1 + term(seed);
        slm  = SurfStatLinMod(measure_used(keep,i),M);
        cttx_ctx(keep,i) = measure_used(keep,i) - slm.X*slm.coef;
    end
end

for cttxCTX  =  1

    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    familyHCPk = cellstr(D1.Family_ID(keep));
    familyHCP = term(familyHCPk);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M = 1 + term(agekeep) + term(genderkeep) + term(icvkeep);
    slm =  SurfStatLinMod(cttx_ctx(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'mpc.gender.CTX_controlled.t.png'],'png', 10)
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'mpc.CTX_controlled.FDR.png'],'png', 10)
    
end

%% do cortical thickness analysis
load('/Users/sofievalk/Dropbox/my_projects/Juelich/2019.StrucCOV/matlab.mat', 'HCP400_CT')
load('/Users/sofievalk/Dropbox/my_projects/Juelich/2019.StrucCOV/matlab.mat', 'unres')
load('/Users/sofievalk/Dropbox/my_projects/Juelich/2019.StrucCOV/matlab.mat', 'D1')
load('/Users/sofievalk/Dropbox/my_projects/Juelich/2019.StrucCOV/matlab.mat', 'parcels400')


   for cortcal_thickness = 1     
    for i =1:1206
        meant1t2(i,:) = mean(squeeze(female_gradients.MPC_layers1(i,2:11,:)));
    end
    
    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M           =       term(agekeep) + term(genderkeep) + term(icvkeep);
    slm =  SurfStatLinMod(meant1t2(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    t1wt2wsex = slm.t;
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])

    genderkeep = cellstr(unres.Gender(keep));
    agekeep = D1.Age_in_Yrs(keep);
    icvkeep = unres.FS_IntraCranial_Vol(keep);
    G = term(genderkeep);
    M           =       term(agekeep) + term(genderkeep) + term(icvkeep);
    slm =  SurfStatLinMod(HCP400_CT(keep,:),M)
    slm = SurfStatT(slm,G.F-G.M)
    
    ct_sex = slm.t;
    
    % FDR correction females - males
    p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
    % > if p > .95, very 
    h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
    hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
    h= h+hn; % complete mask

    % compute effect sizes
    Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
    CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected 
    
    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = Cohensd(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = Cohensd(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'CTX_sex_difference.d.png'],'png', 10)

    heri_ct = zeros(1,20484);
    for i = 1:200
        heri_ct(:,find(parcels400==i+1)) = CohensdFDR(i);
    end
    for i = 1:200
        heri_ct(:,find(parcels400==i+1001)) = CohensdFDR(i+200);
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(flipud(cbrewer('div','RdBu',11)))
    SurfStatColLim([-1 1])
    exportfigbo(f,[RPATH 'CTX_sex_difference.FDR.png'],'png', 10)

    f=figure,
    scatter(ct_sex,t1wt2wsex,'k','filled'),lsline,
    exportfigbo(f,[RPATH 'CTX_T1wT2w_sex_difference.png'],'png', 10)

    
     [p_spin, r_dist] = spin_test(t1wt2wsex', ct_sex, 'surface_name',...
                'fsa5', 'parcellation_name', 'schaefer_400', 'n_rot', 1000, ...
                'type', 'spearman');

   end



%% data of svenja for the seaborn plot

load('/Users/sofievalk/Desktop/hormonal_effects_struct.mat')

% modelname{1} = 'women: OC vs NC';
% modelname{2} = 'Men vs NC women';
% modelname{3} = 'Men vs OC women';
% modelname{4} = 'women: high vs low estrogen';
% modelname{5} = 'women: high vs low progesterone';
% modelname{6} = 'Men vs high prog';
% modelname{7} = 'Men vs low prog';
% modelname{8} = 'Men vs high estr';
% modelname{9} = 'Men vs low estr';

for i = 1:9
    datatest(i,:) = cohensd_profilemean{i};
end

[S,AX,BigAx,H,HAx] = plotmatrix(datatest')
S(3).Color = ;

% Generate sample data
data = randn(100, 3);

% Define colors for each data point
colors = jet(100);

% Create a plot matrix with scatter
figure;

% Using the scatter function

subplot(2, 2, 2);
h = scatter(data(:, 1), data(:, 2), 50, colors, 'filled');
title('Variable 1 vs Variable 2');

subplot(2, 2, 3);
h = scatter(data(:, 2), data(:, 1), 50, colors, 'filled');
title('Variable 2 vs Variable 1');

% Customize the markers in the plot matrix
for i = 1:length(h)
    h(i).MarkerEdgeColor = 'k';  % Set marker edge color to black
end

% Add labels to the whole plot matrix
sgtitle('Scatter Plot Matrix with Varying Circle Colors');

% Add labels to individual scatter plots
for i = 1:4
    subplot(2, 2, i);
    xlabel(['Variable ' num2str(mod(i, 2) + 1)]);
    ylabel(['Variable ' num2str(floor((i - 1) / 2) + 1)]);
end

for i = 1:2
    subplot(2, 2, i.^2);
    histogram(data(:, i), 'Normalization', 'probability', 'FaceColor', 'b');
    title(['Histogram of Variable ' num2str(i)]);
end

[vertices, label, colortablel] = ...
    fs_read_annotation(['/Users/sofievalk/Desktop/cng/code_methods/tpl-fsaverage_hemi-L_desc-types.annot']);
parcel_left = label;
label_left = label;
for i = 1:size(colortablel.table, 1)
    mycode = colortablel.table(i,5);
    parcel_left(find(parcel_left == mycode)) = i;
end

[vertices, label, colortabler] = ...
    fs_read_annotation(['/Users/sofievalk/Desktop/cng/code_methods/tpl-fsaverage_hemi-R_desc-types.annot']);
parcel_right = label;
label_right = label;
for i = 1:size(colortabler.table, 1)
    mycode = colortabler.table(i,5);
    parcel_right(find(parcel_right == mycode)) = i;
end

types_garcia = [parcel_left(1:10242); parcel_right(1:10242)];
types_garcia = types_garcia';

            
% Generate sample data
data = datatest';

% Define colors for each data point
colors = types400;

% Create a plot matrix with scatter and histograms
figure;

model_names = {'women: OC vs NC','Men vs NC women','Men vs OC women','women: high vs low estrogen','women: high vs low progesterone',...
'Men vs high prog','Men vs low prog','Men vs high estr','Men vs low estr'};

% Scatter plots
for i = 1:9
    for j = 1:9
        subplot(9, 9, (i - 1) * 9 + j);
        if i == j
            % Diagonal - Histogram
            % Diagonal - Density Plot
            [f, xi] = ksdensity(data(:, i));
            plot(xi, f, 'LineWidth', 2, 'Color', 'k');
            title([model_names{i}]);
        else
            % Scatter plot
            h = scatter(data(:, i), data(:, j), 10, colors, 'filled'),lsline;
            colormap(colortabler.table(:,1:3)./256)
            %title(['Variable ' num2str(i) ' vs Variable ' num2str(j)]);
            if i < 9
                xticks([]); % Remove x-axis labels for non-bottom subplots
            end
            if j > 1
                yticks([]); % Remove y-axis labels for non-leftmost subplots
            end
        end
    end
end


%% glasser HCP 
% load the mask

load('/Users/sofievalk/Downloads/glasser_hcp_ind.mat')


 for i =1:1206
        skew_gl(i,:) = skewness(squeeze(glasser_hcp_ind.I(i,2:11,:)));
 end
   
 skew_gl(isnan(skew_gl)) = 0;
 
 
to_surface = zeros(20484,1);
for i = 1:360
    to_surface(glasser_fsa5==i) = nanmean(skew_gl(10,i));
end

f = figure,
BoSurfStatViewData(to_surface,SN,'')
colormap((cbrewer('seq','Greens',99)))
BoSurfStatColLim([1.5 2])


for i =1:1206
    skew_gl(i,:) = mean(squeeze(MPC_layers_retest(i,:,:)));
end

skew_gl(isnan(skew_gl)) = 0;

keep = intersect(keep_mpc, find(mean(skew_gl,2)~=0))
genderkeep = cellstr(unres.Gender(keep));
agekeep = D1.Age_in_Yrs(keep);
familyHCPk = cellstr(D1.Family_ID(keep));
familyHCP = term(familyHCPk);
icvkeep = unres.FS_IntraCranial_Vol(keep);
G = term(genderkeep);
M = 1 + term(agekeep) + term(genderkeep) + term(icvkeep);
slm =  SurfStatLinMod(skew_gl(keep,:),M)
slm = SurfStatT(slm,G.F-G.M)

% FDR correction females - males
p = 1-tcdf(slm.t,slm.df); % probability of a sample having a larger t-score than the t-score of the sample.
% > if p > .95, very
h = fdr_bh(p,0.025); % if h = 1, then pval is significant, first tail
hn = fdr_bh(1-p,0.025); % if h = 1, then pval is significant, other tail
h= h+hn; % complete mask

% compute effect sizes
Cohensd = 2*slm.t / sqrt(slm.df); %Cohen's d = 2t / Sqrt(df) (See Rosenthal (1994) and Howell (2013), p.649)
CohensdFDR = (Cohensd.*h)'; % mask only those that were FDR corrected
  
% 2:181 + 183:362
to_surface = zeros(20484,1);
for i = 1:180
    to_surface(parcels360x==i+1) = mean(skew_gl(keep,i));
end
for i = 1:180
    to_surface(parcels360x==i+182) = mean(skew_gl(keep,i+180));
end

f = figure,
BoSurfStatViewData(to_surface,SN,'')
colormap(flipud(cbrewer('div','RdBu',99)))
SurfStatColLim([-0.5 0.5])


for i = 1:12
    for j = 1:12
        for k = 1:360
            r = corr(squeeze(glasser_hcp_ind.I(keep,i,k)),squeeze(glasser_hcp_ind.I(keep,j,k)),'type','spearman');
            layer_r(i,j,k) = r;
        end
    end
end

f= figure,
subplot(6,2, 1);
imagesc(squeeze(layer_r(1,:,:)),[0.5 1])
subplot(6,2, 2);
imagesc(squeeze(layer_r(2,:,:)),[0.5 1])
subplot(6,2, 3);
imagesc(squeeze(layer_r(3,:,:)),[0.5 1])
subplot(6,2, 4);
imagesc(squeeze(layer_r(4,:,:)),[0.5 1])
subplot(6,2, 5);
imagesc(squeeze(layer_r(5,:,:)),[0.5 1])
subplot(6,2, 6);
imagesc(squeeze(layer_r(6,:,:)),[0.5 1])
subplot(6,2, 7);
imagesc(squeeze(layer_r(7,:,:)),[0.5 1])
subplot(6,2, 8);
imagesc(squeeze(layer_r(8,:,:)),[0.5 1])
subplot(6,2, 9);
imagesc(squeeze(layer_r(9,:,:)),[0.5 1])
subplot(6,2, 10);
imagesc(squeeze(layer_r(10,:,:)),[0.5 1])
subplot(6,2, 11);
imagesc(squeeze(layer_r(11,:,:)),[0.5 1])
subplot(6,2, 12);
imagesc(squeeze(layer_r(12,:,:)),[0.5 1])


for j = 1 : 12
    to_surface = zeros(20484,1);
    for i = 1:360
        to_surface(glasser_fsa5==i) = mean(squeeze(layer_r(j,:,i)));
    end
    
    f = figure,
    BoSurfStatViewData(to_surface,SN,'')
    colormap(flipud(cbrewer('div','RdBu',99)))
end


%% test svenja


skew_gl = []
for i =1:1206
    skew_gl(i,:) = mean(squeeze(MPC_layers_retest(i,:,:)));
end

skew_gl(isnan(skew_gl)) = 0;

keep = find(mean(skew_gl,2)~=0)


% 2:181 + 183:362
to_surface = zeros(20484,1);
for i = 1:180
    to_surface(parcels360x==i+1) = mean(skew_gl(keep,i+1));
end
for i = 1:180
    to_surface(parcels360x==i+182) = mean(skew_gl(keep,i+182));
end


f = figure,
BoSurfStatViewData(to_surface,SN,'')


