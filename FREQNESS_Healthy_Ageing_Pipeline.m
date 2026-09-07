% This script contains the pipeline for the 'Ageing reshapes the resting-state network landscape of the human brain' publication. 
% GED is carried out on the whole sample of participants.
% The processing of NIFTI files (spatial activation patterns) and
% eigenspectra (explained variance over frequenciez) is performed on
% 'Young' and 'Old' groups, separately. 
% Statistical testing compares the eigenspectra from the two groups.

% Nikita Kudriavtsev, Aarhus (DK), 06/09/2026

%% Perform GENERALISED EIGENDECOMPOSITION (GED)

% Setting up cluster (parallel computing)
% clusterconfig('scheduler', 'none');  % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); % If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 8); % slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')

% Choose whether to test whole sample of frequencies (set to 1), or just
% the target stimulation frequency (set to 0)
allfrex_flag = 1; 

% NOTE: with the current approach, we anchor the width of the filter to the
% first value based on Rosso et al., (2021a,2023a), and then separately
% compute the other widths as logarithmic functions of the frequency bins 

% Settings
S = struct();     % initialize input structure to submit to cluster
S.time = 304;     % in seconds (to match the original FREQNESS dataset)
S.targetfrex = 2.439; % stimulation frequency: hard-coded 1 / 0.410 s (range is 0.400 and 0.420 but 0.410 is the mode) 
S.fwhm  = .35;        % filter width
S.shr   = 0.01;       % shrinking proportion for regularization; value between 0 and 1 (until 0.001 it works for restoring the full rank of the matrix)
S.ncomps2keep = 100;  % how many components to keep (given the huge amount of components = N voxels)
S.comps2see   = [1 2 3]; % list components of interest, of which you want to save nifti image
S.path_script    = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia'; % path of the present script
S.path_source    = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1_TSA2021_SR250Hz/';%
S.path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/Aging/Aging_GED'; 
S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz/']; 
S.srate = 250; % 250 Hz is the standard after downsampling

% Set all frequencies to test, and associated filter widths
if allfrex_flag == 1 

    % Set with respect to stimulation frequency
    nfrex_above = 80; % tested for even numbers
    nfrex_below = 6;   % tested for even numbers
    % Compute all frequencies
    frex_above = linspace(S.targetfrex, S.targetfrex * (nfrex_above/2), nfrex_above);
    frex_below = linspace(S.targetfrex / 2, S.targetfrex / (2 * nfrex_below), nfrex_below);
    frex_all = [frex_below(end:-1:1), frex_above];
    % compute all filter widths
    fwhm_above = logspace(log10(S.fwhm), log10(S.fwhm * (nfrex_above)), nfrex_above);
    fwhm_below = logspace(log10(S.fwhm), log10(S.fwhm / (nfrex_below)), nfrex_below+1);
    fwhm_all   = [fwhm_below(end:-1:2) fwhm_above];

    % Visualize frewquencies and filters
    figure
    plot(frex_all, fwhm_all, 's-')
    xlabel('Center frequencies (Hz)')
    ylabel('Filter FWHM (Hz)')
    title('Filter width as a function of center frequency')

    % Run GED with parallel computing for a set of frequencies
    nfrex = length(frex_all);
    for frexi = 1:nfrex
        S.targetfrex  = frex_all(frexi); % overwrite fields of input structure
        S.fwhm        = fwhm_all(frexi);
        S.name_folder = ['GEDoutput_frex_' num2str(S.targetfrex) 'Hz/'];
        % Running with parallel computing
        jobid = job2cluster(@FREQNESS_AdvScience_GED_SingleSubjects,S); % always give a structure as input
        
        disp(['Submitting GED for frequency = ' num2str(S.targetfrex) ' Hz and width = ' num2str(S.fwhm) ' Hz on cluster'])
        
    end

else
    
    % Run GED with parallel computing for single stimulation frequency
    jobid = job2cluster(@FREQNESS_AdvScience_GED_SingleSubjects,S); % always give a structure as input

end

% repeat the same for the second dataset.

% The sections below use the same paths, workspace variables and external
% functions as the original analysis script. Run the sections as before.

%% Setting up age groups and collecting eigenvalues
% Run the original dataset selection and collection sections for each dataset.

%% setting up age groups to match the GED optput
data_switch = 1 %this allows switching between outputs by the previous section and combining them into one structure
if data_switch==1
    GED_dir = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/Aging/Aging_GEDGEDoutput_frex_';
    curr_freq = 0.813;
    curr_freq_dir = dir([GED_dir num2str(curr_freq) 'Hz/*.mat']);
    age_info = readtable('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/TSA2021_Nikita.xlsx'); % loading behavioural data
    age_info = age_info(:,1:4);
    IDs_formated = compose('%04d', age_info.Subject);
else
    GED_dir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Aging_GED/ResultsGEDoutput_frex_';
    curr_freq = 0.813;
    curr_freq_dir = dir([GED_dir num2str(curr_freq) 'Hz/*.mat']);
    age_info = readtable('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/Dement_Fold/Participant_categories_AuditMemDement.xlsx');
    age_info = age_info(:,1:4);
    IDs_formated = compose('%04d', age_info.ID);
end

prefix = 'SUBJ_';
suffix = '.mat';
age_info.ID_forInd = IDs_formated+suffix;
IDs_formated = prefix + IDs_formated + suffix;
IDs_present = string({curr_freq_dir.name}).';

logical_vector_relevan_part = ismember(IDs_formated,IDs_present);
age_info = age_info(logical_vector_relevan_part,:);
%logical indexing of age groups
if data_switch==1
    Y_Part = age_info.Age<27;
    O_Part = age_info.Age>58;
else
    Y_Part = ismember(age_info.Group,'A');
    O_Part = ismember(age_info.Group,{'B', 'B_MCI'});
end

idx_young = age_info.ID_forInd(Y_Part);
idx_old = age_info.ID_forInd(O_Part);
idx_all = {idx_young; idx_old};

%%
path_GEDoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/Aging/';

ngroups = length(idx_all);

% Suffixes for saving the output
suffix_group = {'_young'; '_old'};

% Select components of interest
comps2see = [1:10];
stimfrex = 2.439;

% Get all_frex from home path
list_GED = dir(path_GEDoutput);
list_GED = list_GED( [list_GED.isdir] );
list_GED = list_GED(endsWith({list_GED.name}, 'Hz')); % added since GED and PCA coexist in the same folder

% Initialize vectors for storing frequencies
nfrex = length(list_GED);
all_frex = zeros(nfrex,1);
% Initialize topEvals
topEvals_avg = zeros(length(comps2see),nfrex); % avg of top components across subjects, per frequency
topEvals_avg_group = zeros([size(topEvals_avg), ngroups]); % for assignment to group

for groupi = 1:ngroups
% Import and process top eigenvalues, per frequency
for frexi = 1:nfrex
    
    % Get folder name
    this_folder = [path_GEDoutput, list_GED(frexi).name '/'];
    % Retrieve frequency
    this_frex = regexp(this_folder, '([\d.]+)Hz', 'tokens');
    % Assign frequency to vector
    all_frex(frexi) = str2double(this_frex{1,1});
    
    % All subjects within frequency frexi
    fulllist = dir([this_folder '/SUBJ*.mat']);
    
    % Filter participants by age
    age_finder = {fulllist(:).name};
    age_filter = contains(age_finder, idx_all{groupi});
    list_subs = fulllist(age_filter);
    nsubs = length(list_subs);
    
    % Process evals from subjects
    for subi = 1:nsubs % over subjects
        
        % Get eigenvalues (already sorted and normalized)
        load([this_folder, list_subs(subi).name],'evals')
        nevals = length(evals);
        
        if subi == 1 && frexi == 1 && groupi == 1
            % Initialize matrix for all participants, once upon first iteration
            evals_all = zeros(nevals,40,nfrex); % for Rest and Listen conditions
            topEvals_all_group = zeros([size(evals_all), ngroups]);
        end
        
        evals_all(:,subi,frexi) = evals; % assign for computing average and for exporting

    end
    
    % Compute grand-averages over subjects (just for visualization)
    topEvals_avg(:,frexi) = squeeze( mean(evals_all(comps2see,:,frexi),2) );
    
    disp(['Processing evals of frex ' this_frex{1,1} 'Hz'])
    
end
% Re-sort vector of frequencies (strings follow alphabetic order)
[all_frex, idx_sort] = sort(all_frex, 'ascend');
topEvals_avg = topEvals_avg(:,idx_sort);

% Final assignment
topEvals_avg_group(:,:,groupi) = topEvals_avg; % only for inspection
topEvals_all_group(:,:,:,groupi) = evals_all;  % for actual testing, as contains all subjects per group

end

% Save output
save([landscape_output '/NetworkLandscape_Groups_TempSeq_Ages.mat'], 'topEvals_all_group', 'all_frex', 'idx_sort')

%% Statistics

clc
clear all
close all
%%
path_GEDoutput = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/'; % This loads output from previous block

% Load data
data_switch = 2
if data_switch==1
    load([path_GEDoutput 'NetworkLandscape_Groups_TempSeq_Ages.mat'])
    Y_subNum = 37;
    O_subNum = 40;
    Dat_num = 1;
else
    load([path_GEDoutput 'NetworkLandscape_Groups_2.mat'])
    Y_subNum = 48;
    O_subNum = 40;
    Dat_num = 2;
end

%size(topEvals_all_group)

% Reduce the eigenspectrum
last_frex = 30; % consider the eigenspectrum up to this frequency
[~,idx_last] = min(abs(all_frex-last_frex));
last_comp = 3; % consider up to this component
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:); % re-sort frequencies
topEvals_all_group = topEvals_all_group(1:last_comp,:,1:idx_last,:);

% Initialize vector to concatenate stats for FDR correction
all_Pevals = [];
all_Tevals = [];

% Test differences in the eigenspectra, across frequencies, between the two groups

% Select components of interest
comp2test = 1;
ncomps = length(comp2test);
nfrex = size(topEvals_all_group,3);

% Test eigenvalues
range2test = 1:idx_last; % based on reference PNAS paper (2016)
all_frex = all_frex(range2test);

% Initialize output
[Pevals,Tevals] = deal(zeros(ncomps,length(all_frex)));

% Run stats
for frexi = 1:length(all_frex)
    [~,p,~,stats] = ttest2( squeeze(topEvals_all_group(comp2test,1:Y_subNum,frexi,1)) , squeeze(topEvals_all_group(comp2test,1:O_subNum,frexi,2)));
    Pevals(ncomps,frexi) = p;
    Tevals(ncomps,frexi) = stats.tstat;
end

% Concatenate vector for FDR correction on all conditions
% NOTE: This is redundant when you test one component at a time, but if you
% want to loop over components, then you need to concatenate all P and T
% vals in a row vector for FDR correction
all_Pevals = [all_Pevals, Pevals];
all_Tevals = [all_Tevals, Tevals];
%FDR correction
thresh_FDR = zeros(size(all_Pevals,1),1); %FDR thresholds (independently for each GED component)
P_above_thresh = zeros(size(all_Pevals)); %frequencies (variance explained) which are significant after FDR correction
for ii = 1:size(all_Pevals,1) %over GED components 
    thresh_FDR(ii,1) = fdr(all_Pevals(ii,:)); %FDR correction (getting the threshold)
    P_above_thresh(ii,all_Pevals(ii,:)<=thresh_FDR(ii,1)) = 1; %getting the significant frequencies (independently for each component)
end

% Find significant frequencies after FDR
cat_all_frex = repmat(all_frex,[4 1]);
frex_postFDR = cat_all_frex(find(P_above_thresh));
% Display
for frexi = 1:length(frex_postFDR)
    disp(['The significant frequencies after FDR are: ' num2str(frex_postFDR(frexi)) ' Hz'])
end

% --- Extract group matrices: [subjects x freqs] ---
G1 = squeeze(topEvals_all_group(comp2test,1:Y_subNum,:,1));  % group 1 (37 subs)
G2 = squeeze(topEvals_all_group(comp2test,1:O_subNum,:,2));  % group 2 (all subs present)

% --- Group 1 stats (95% CI, NaN-robust) ---
m1  = nanmean(G1,1);
n1  = sum(~isnan(G1),1);
sd1 = nanstd(G1,0,1);
se1 = sd1 ./ max(sqrt(n1),1);
t1  = tinv(0.975, max(n1-1,1));     % two-sided 95% CI
ci1 = t1 .* se1;
lo1 = m1 - ci1; hi1 = m1 + ci1;

% --- Group 2 stats (95% CI, NaN-robust) ---
m2  = nanmean(G2,1);
n2  = sum(~isnan(G2),1);
sd2 = nanstd(G2,0,1);
se2 = sd2 ./ max(sqrt(n2),1);
t2  = tinv(0.975, max(n2-1,1));
ci2 = t2 .* se2;
lo2 = m2 - ci2; hi2 = m2 + ci2;

% --- Plot shaded bands + mean lines (blue/red like before) ---
xx = [all_frex, fliplr(all_frex)];
xx=xx.';
% force row vectors
% --- before calling fill ---
x   = all_frex(:)';   % 1×F
lo1 = lo1(:)';  hi1 = hi1(:)';
lo2 = lo2(:)';  hi2 = hi2(:)';

xx1 = [x fliplr(x)];          % 1×(2F)
yy1 = [lo1 fliplr(hi1)];      % 1×(2F)
xx2 = [x fliplr(x)];
yy2 = [lo2 fliplr(hi2)];

figure; hold on
fill(xx1, yy1, [0 0.447 0.741], 'FaceAlpha',0.38, 'EdgeColor','none');
fill(xx2, yy2, [0.85 0.325 0.098], 'FaceAlpha',0.38, 'EdgeColor','none');

% If you have the means, plot them (better than (lo+hi)/2):
plot(x, m1, 'Color',[0 0.447 0.741], 'LineWidth',4)
plot(x, m2, 'Color',[0.85 0.325 0.098], 'LineWidth',4)
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)
set(gca, 'XTick', all_frex(1:3:last_frex));
set(gca, 'XTickLabels', string(all_frex(1:3:last_frex)));
xtickangle(45)

xlabel('Frequency (Hz)', 'FontSize', 18)
ylabel('Explained variance (%)', 'FontSize', 18)
title(['Dataset#' num2str(Dat_num) ': Comp. ' num2str(comp2test) ' Eigenspectrum (group means \pm 95% CI)'], 'FontSize', 24)
legend({'Youth 95% CI', 'Elderly 95% CI','Youth mean','Elderly mean'}, 'Location','best'); box on

%% Computing center of mass

%% jobs to cluster to compute my info
addpath('//scratch7/MINDLAB2023_MEG-AuditMemDement/nikita');
%S.GED_dir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Aging_GED/';
S.GED_dir = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/Aging/';
dir_list_files = dir([S.GED_dir 'Aging*']);
S.curr_freq = 0;
S.outdir_coodr = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Comp_coord/TempSeq_MNI_cent_mass/';
S.comp_saved = 1; 
S.thresh_swtitch = 1;
for freq_sig = [1:51 61 70 79]%7:12%50%length(dir_list_files)
    S.curr_freq = dir_list_files(freq_sig).name;
    clusterconfig('scheduler', 'cluster'); %set the cluster
    % clusterconfig('scheduler', 'none'); %set locally
    clusterconfig('long_running', 0); %set the most typical cue (see our cluster guide in labook for details)
    clusterconfig('slot', 1); %set amount of memory; between 1 and 12 (each slot is 8gb of ram)
    jobid = job2cluster(@center_of_mass,S)
    disp([num2str(freq_sig) '/' num2str(length(dir_list_files))]);
end
%
addpath('//scratch7/MINDLAB2023_MEG-AuditMemDement/nikita');
S.GED_dir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Aging_GED/';
%S.GED_dir = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/Aging/';
dir_list_files = dir([S.GED_dir 'Results*']);
S.curr_freq = 0;
S.outdir_coodr = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Comp_coord/AudMem_MNI_cent_mass/';
S.comp_saved = 1; 
S.thresh_swtitch = 1;
for freq_sig = [1:51 61 70 79]%7:12%50%length(dir_list_files)
    S.curr_freq = dir_list_files(freq_sig).name;
    clusterconfig('scheduler', 'cluster'); %set the cluster
    % clusterconfig('scheduler', 'none'); %set locally
    clusterconfig('long_running', 0); %set the most typical cue (see our cluster guide in labook for details)
    clusterconfig('slot', 1); %set amount of memory; between 1 and 12 (each slot is 8gb of ram)
    jobid = job2cluster(@center_of_mass,S)
    disp([num2str(freq_sig) '/' num2str(length(dir_list_files))]);
end

%% Quadratic Renyi entropy - both datasets combined

%% BOTH DATASETS COMBINED ANALYSIS SECTION

path_GEDoutput = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/'; % This loads output from previous block
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita')
% Load data

landscape = load([path_GEDoutput 'NetworkLandscape_Groups_TempSeq_Ages.mat'])
topEvals_all_group = landscape.topEvals_all_group;
all_frex = landscape.all_frex;
idx_sort = landscape.idx_sort;
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:);
Y_subNum = 37;
O_subNum = 40;
Dat_num = 1;
comp_incl = 3559;
freq_retained = 80;
Y_FREQ_TS = [];
Y_topEvals_TS = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
O_topEvals_TS = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
Y_holder_TS = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
Y_holder_TS = permute(Y_holder_TS,[1 3 2]);
Y_FREQ_TS.evals = Y_holder_TS;
Y_FREQ_TS.frex = all_frex(1:freq_retained);
O_FREQ_TS = [];
O_holder_TS = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
O_holder_TS = permute(O_holder_TS,[1 3 2]);
O_FREQ_TS.evals = O_holder_TS;
O_FREQ_TS.frex = all_frex(1:freq_retained);

landscape = load([path_GEDoutput 'NetworkLandscape_Groups_2.mat'])
topEvals_all_group = landscape.topEvals_all_group;
all_frex = landscape.all_frex;
idx_sort = landscape.idx_sort;
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:);
Y_subNum = 48;
O_subNum = 40;
Dat_num = 2;
Y_FREQ_AMD = [];
Y_topEvals_AMD = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
O_topEvals_AMD = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
Y_holder_AMD = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
Y_holder_AMD = permute(Y_holder_AMD,[1 3 2]);
Y_FREQ_AMD.evals = Y_holder_AMD;
Y_FREQ_AMD.frex = all_frex(1:freq_retained);
O_FREQ_AMD = [];
O_holder_AMD = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
O_holder_AMD = permute(O_holder_AMD,[1 3 2]);
O_FREQ_AMD.evals = O_holder_AMD;
O_FREQ_AMD.frex = all_frex(1:freq_retained);

Y_FREQ_BOTH = [];
Y_holder_BOTH = cat(3,Y_holder_TS,Y_holder_AMD);
Y_FREQ_BOTH.evals = Y_holder_BOTH;
Y_FREQ_BOTH.frex = all_frex(1:freq_retained);

O_FREQ_BOTH = [];
O_holder_BOTH = cat(3,O_holder_TS,O_holder_AMD);
O_FREQ_BOTH.evals = O_holder_BOTH;
O_FREQ_BOTH.frex = all_frex(1:freq_retained);
% Test differences in the eigenspectra, across frequencies, between the two groups

[ED_Y, H2_Y] = FREQNESS_EntropyLandscape(Y_FREQ_BOTH);
[ED_O, H2_O] = FREQNESS_EntropyLandscape(O_FREQ_BOTH);
H2_Y(isinf(H2_Y)) = NaN;
H2_O(isinf(H2_O)) = NaN;

%%
all_frex = all_frex(1:freq_retained);
%
% --- Extract group matrices: [subjects x freqs] ---

G1 = H2_Y';  % group 1 (37 subs)
G2 = H2_O';  % group 2 (all subs present) 

% group_comparison
freq_stat = zeros(freq_retained,1);
for freq = 1:freq_retained
    [pval,tval] = ranksum(G1(:,freq),G2(:,freq));
    freq_stat(freq,1) = pval;
end

[sig_ent_freq,crit_p] = fdr_bh(freq_stat,0.01,'pdep','yes');

% --- Group 1 stats (95% CI, NaN-robust) ---
m1  = nanmean(G1,1);
n1  = sum(~isnan(G1),1);
sd1 = nanstd(G1,0,1);
se1 = sd1 ./ max(sqrt(n1),1);
t1  = tinv(0.975, max(n1-1,1));     % two-sided 95% CI
ci1 = t1 .* se1;
lo1 = m1 - ci1; hi1 = m1 + ci1;

% --- Group 2 stats (95% CI, NaN-robust) ---
m2  = nanmean(G2,1);
n2  = sum(~isnan(G2),1);
sd2 = nanstd(G2,0,1);
se2 = sd2 ./ max(sqrt(n2),1);
t2  = tinv(0.975, max(n2-1,1));
ci2 = t2 .* se2;
lo2 = m2 - ci2; hi2 = m2 + ci2;

% --- Plot shaded bands + mean lines (blue/red like before) ---
xx = [all_frex, fliplr(all_frex)];
xx=xx.';
% force row vectors
% --- before calling fill ---
x   = all_frex(:)';   % 1×F
lo1 = lo1(:)';  hi1 = hi1(:)';
lo2 = lo2(:)';  hi2 = hi2(:)';

xx1 = [x fliplr(x)];          % 1×(2F)
yy1 = [lo1 fliplr(hi1)];      % 1×(2F)
xx2 = [x fliplr(x)];
yy2 = [lo2 fliplr(hi2)];

figure; hold on
fill(xx1, yy1, [0 0.447 0.741], 'FaceAlpha',0.38, 'EdgeColor','none');
fill(xx2, yy2, [0.85 0.325 0.098], 'FaceAlpha',0.38, 'EdgeColor','none');

% If you have the means, plot them (better than (lo+hi)/2):
plot(x, m1, 'Color',[0 0.447 0.741], 'LineWidth',4)
plot(x, m2, 'Color',[0.85 0.325 0.098], 'LineWidth',4)
scatter(all_frex(find(sig_ent_freq>0)),sig_ent_freq(find(sig_ent_freq>0))*2.5, 'r', '*')

set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)

set(gca, 'XTick', all_frex(1:3:freq_retained));
set(gca, 'XTickLabels', string(all_frex(1:3:freq_retained)));
xtickangle(45)
ylim([2.5 4]);

xlabel('Frequency (Hz)', 'FontSize', 18)
ylabel('Quadratic Renyi entropy (H2)', 'FontSize', 18)
title(['Entropy over all eigenvalues'], 'FontSize', 24)  
legend({'Youth 95% CI', 'Elderly 95% CI','Youth mean','Elderly mean'}, 'Location','best'); box on

%% Group average activation patterns for the figures
% Kept here while Y_FREQ_BOTH is still a structure, before the numeric analyses.

%% combining cross frequencies all my info
freq_ind = [1:6 15 24 34 51:54 7:14 16:23 25:33 35:42];
num_freq = 40;
freq_ind_select = freq_ind(1:num_freq);
Dataset_switch = 1; % Change to upload connectivity matrices from a certain dataset 0 = TemoSeqAges, 1 = Dementia

%Y/O separation
participant_info = readtable('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/TSA2021_Nikita.xlsx'); % loading behavioural data
Valid_dat = participant_info(:,:); %defining a dataset with participants that have all the behavioural tests
Young = Valid_dat.Subject(Valid_dat.Age < 27);
Old = Valid_dat.Subject(Valid_dat.Age > 55);
Y_Selection = ismember(Valid_dat.Subject,Young);
O_Selection = ismember(Valid_dat.Subject,Old);
%Frequency indices

Y_Selection_TempSEQ = Y_Selection;
O_Selection_TempSEQ = O_Selection;

%%
%Y/O separation
participant_info = readtable('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/Descript_Dat/Participant_categories_AuditMemDement_ages.xlsx'); % loading behavioural data
% part 3,8,59,65 have no resting state GED
Val_part = ismember(participant_info.ID,[3 6 59 65]);
Valid_participant_info  = participant_info(~Val_part,:);
Y_Selection = ismember(Valid_participant_info.Group,{'A'});
%O_Selection = ismember(Valid_participant_info.Group,{'B'});
O_Selection = ismember(Valid_participant_info.Group,{'B','B_MCI'});
D_Selection = ismember(Valid_participant_info.Group,{'C'});
%%
Y_Selection_AudMem = Y_Selection;
O_Selection_AudMem = O_Selection;

Y_Selection_BOTH = [Y_Selection_TempSEQ;Y_Selection_AudMem];
O_Selection_BOTH = [O_Selection_TempSEQ;O_Selection_AudMem];

Y_Selection = Y_Selection_BOTH;
O_Selection = O_Selection_BOTH;

%% Creating group average activation patterns for the figures
path_GEDoutputGrAvg = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Age_GED_mapss';
frex_to_plot = [6,12,13,14,15,25];
comp2avrg = 1;
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); % getting the mask for creating the figure

SS = size(maskk.img);

for freqforfigi = frex_to_plot
    AudMemDir = dir(['/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Aging_GED/ResultsGEDoutput_frex_' , num2str(Y_FREQ_BOTH.frex(freqforfigi)) , 'Hz/*.mat']);
    AudMemAllActPatt = nan(length(AudMemDir),3559);
    for tempi = 1:length(AudMemDir)
        tmp = load([AudMemDir(tempi).folder,'/',AudMemDir(tempi).name],'GEDmap');
        tmp = tmp.GEDmap;
        tmp = squeeze(tmp(comp2avrg,:));
        tmp = abs(tmp);
        tmp = tmp/sum(tmp);
        AudMemAllActPatt(tempi,:) = tmp;
        disp(tempi)
    end

    TempSeqDir = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/Aging/Aging_GEDGEDoutput_frex_' , num2str(Y_FREQ_BOTH.frex(freqforfigi)) , 'Hz/*.mat']);
    TempSeqAllActPatt = nan(length(TempSeqDir),3559);
    for tempi = 1:length(TempSeqDir)
        tmp = load([TempSeqDir(tempi).folder,'/',TempSeqDir(tempi).name],'GEDmap');
        tmp = tmp.GEDmap;
        tmp = squeeze(tmp(comp2avrg,:));
        tmp = abs(tmp);
        tmp = tmp/sum(tmp);
        TempSeqAllActPatt(tempi,:) = tmp;
        disp(tempi)
    end
    
    ActPattAll = cat(1,TempSeqAllActPatt,AudMemAllActPatt);
    
    Y_ActPattAll = ActPattAll(Y_Selection,:);
    Y_ActPattAll = nanmean(Y_ActPattAll,1);
    Y_thresh = mean(Y_ActPattAll) + std(Y_ActPattAll);
    surv_Vox_Y = find(Y_ActPattAll>Y_thresh);
    dumimg = zeros(SS(1),SS(2),SS(3),1);
    for i = 1:length(surv_Vox_Y) % over brain sources
        sourci = surv_Vox_Y(i);
        dumm = find(maskk.img == sourci); % finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); % getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = Y_ActPattAll(sourci); % storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    save_nii(nii,[path_GEDoutputGrAvg '/YoungAverage_Comp_',num2str(comp2avrg),'frex_' num2str(Y_FREQ_BOTH.frex(freqforfigi)) '.nii']); % printing image
    
    O_ActPattAll = ActPattAll(O_Selection,:);
    O_ActPattAll = nanmean(O_ActPattAll,1);
    O_thresh = mean(O_ActPattAll) + std(O_ActPattAll);
    surv_Vox_O = find(O_ActPattAll>O_thresh);
    dumimg = zeros(SS(1),SS(2),SS(3),1);
    for i = 1:length(surv_Vox_O) % over brain sources
        sourci = surv_Vox_O(i);
        dumm = find(maskk.img == sourci); % finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); % getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = O_ActPattAll(sourci); % storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    save_nii(nii,[path_GEDoutputGrAvg '/OldAverage_Comp_',num2str(comp2avrg),'frex_' num2str(Y_FREQ_BOTH.frex(freqforfigi)) '.nii']); % printing image
    
end

%% Leading component eigenspectrum - both datasets combined

%%
path_GEDoutput = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/'; % This loads output from previous block
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita')
% Load data

landscape = load([path_GEDoutput 'NetworkLandscape_Groups_TempSeq_Ages.mat'])
topEvals_all_group = landscape.topEvals_all_group;
all_frex = landscape.all_frex;
idx_sort = landscape.idx_sort;
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:);
Y_subNum = 37;
O_subNum = 40;
Dat_num = 1;
comp_incl = 3559;
freq_retained = 80;
Y_FREQ_TS = [];
Y_topEvals_TS = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
O_topEvals_TS = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
Y_holder_TS = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
Y_holder_TS = permute(Y_holder_TS,[1 3 2]);
Y_FREQ_TS.evals = Y_holder_TS;
Y_FREQ_TS.frex = all_frex(1:freq_retained);
O_FREQ_TS = [];
O_holder_TS = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
O_holder_TS = permute(O_holder_TS,[1 3 2]);
O_FREQ_TS.evals = O_holder_TS;
O_FREQ_TS.frex = all_frex(1:freq_retained);

landscape = load([path_GEDoutput 'NetworkLandscape_Groups_2.mat'])
topEvals_all_group = landscape.topEvals_all_group;
all_frex = landscape.all_frex;
idx_sort = landscape.idx_sort;
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:);
Y_subNum = 48;
O_subNum = 40;
Dat_num = 2;
Y_FREQ_AMD = [];
Y_topEvals_AMD = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
O_topEvals_AMD = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
Y_holder_AMD = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
Y_holder_AMD = permute(Y_holder_AMD,[1 3 2]);
Y_FREQ_AMD.evals = Y_holder_AMD;
Y_FREQ_AMD.frex = all_frex(1:freq_retained);
O_FREQ_AMD = [];
O_holder_AMD = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
O_holder_AMD = permute(O_holder_AMD,[1 3 2]);
O_FREQ_AMD.evals = O_holder_AMD;
O_FREQ_AMD.frex = all_frex(1:freq_retained);
%%
Y_FREQ_BOTH = [];
Y_holder_BOTH = cat(3,Y_holder_TS,Y_holder_AMD);
Y_FREQ_BOTH.evals = Y_holder_BOTH;
Y_FREQ_BOTH.frex = all_frex(1:freq_retained);

O_FREQ_BOTH = [];
O_holder_BOTH = cat(3,O_holder_TS,O_holder_AMD);
O_FREQ_BOTH.evals = O_holder_BOTH;
O_FREQ_BOTH.frex = all_frex(1:freq_retained);
%%
% Test differences in the eigenspectra, across frequencies, between the two groups
Y_FREQ_BOTH = Y_FREQ_BOTH.evals;
O_FREQ_BOTH = O_FREQ_BOTH.evals;
%%
% Initialize vector to concatenate stats for FDR correction
all_Pevals = [];
all_Tevals = [];

% Select components of interest
comp2test = 1;
ncomps = length(comp2test);
nfrex = size(Y_FREQ_BOTH,2);

% Test eigenvalues
load([path_GEDoutput 'NetworkLandscape_Groups_TempSeq_Ages.mat'], 'all_frex')
range2test = 1:nfrex; % based on reference PNAS paper (2016)
all_frex = all_frex(range2test);

% Initialize output
[Pevals,Tevals] = deal(zeros(ncomps,length(all_frex)));

% Run stats
for frexi = 1:length(all_frex)
    [~,p,~,stats] = ttest2( squeeze(Y_FREQ_BOTH(comp2test,frexi,:)) , squeeze(O_FREQ_BOTH(comp2test,frexi,:)));
    Pevals(ncomps,frexi) = p;
    Tevals(ncomps,frexi) = stats.tstat;
end

% Concatenate vector for FDR correction on all conditions
% NOTE: Thid is redundant when you test one component at a time, but if you
% want to loop over components, then you need to concatenate all P and T
% vals in a row vector for FDR correction
all_Pevals = [all_Pevals, Pevals];
all_Tevals = [all_Tevals, Tevals];
%FDR correction
thresh_FDR = zeros(size(all_Pevals,1),1); %FDR thresholds (independently for each GED component)
P_above_thresh = zeros(size(all_Pevals)); %frequencies (variance explained) which are significant after FDR correction
for ii = 1:size(all_Pevals,1) %over GED components 
    thresh_FDR(ii,1) = fdr(all_Pevals(ii,:)); %FDR correction (getting the threshold)
    P_above_thresh(ii,all_Pevals(ii,:)<=thresh_FDR(ii,1)) = 1; %getting the significant frequencies (independently for each component)
end

% Find significant frequencies after FDR
cat_all_frex = repmat(all_frex,[4 1]);
frex_postFDR = cat_all_frex(find(P_above_thresh));
% Display
for frexi = 1:length(frex_postFDR)
    disp(['The significant frequencies after FDR are: ' num2str(frex_postFDR(frexi)) ' Hz'])
end

% --- Extract group matrices: [subjects x freqs] ---
G1 = squeeze(Y_FREQ_BOTH(comp2test,:,:))';  % group 1 (37 subs)
G2 = squeeze(O_FREQ_BOTH(comp2test,:,:))';  % group 2 (all subs present)

% --- Group 1 stats (95% CI, NaN-robust) ---
m1  = nanmean(G1,1);
n1  = sum(~isnan(G1),1);
sd1 = nanstd(G1,0,1);
se1 = sd1 ./ max(sqrt(n1),1);
t1  = tinv(0.975, max(n1-1,1));     % two-sided 95% CI
ci1 = t1 .* se1;
lo1 = m1 - ci1; hi1 = m1 + ci1;

% --- Group 2 stats (95% CI, NaN-robust) ---
m2  = nanmean(G2,1);
n2  = sum(~isnan(G2),1);
sd2 = nanstd(G2,0,1);
se2 = sd2 ./ max(sqrt(n2),1);
t2  = tinv(0.975, max(n2-1,1));
ci2 = t2 .* se2;
lo2 = m2 - ci2; hi2 = m2 + ci2;

% --- Plot shaded bands + mean lines (blue/red like before) ---
xx = [all_frex, fliplr(all_frex)];
xx=xx.';
% force row vectors
% --- before calling fill ---
x   = all_frex(:)';   % 1×F
lo1 = lo1(:)';  hi1 = hi1(:)';
lo2 = lo2(:)';  hi2 = hi2(:)';

xx1 = [x fliplr(x)];          % 1×(2F)
yy1 = [lo1 fliplr(hi1)];      % 1×(2F)
xx2 = [x fliplr(x)];
yy2 = [lo2 fliplr(hi2)];

figure; hold on
fill(xx1, yy1, [0 0.447 0.741], 'FaceAlpha',0.38, 'EdgeColor','none');
fill(xx2, yy2, [0.85 0.325 0.098], 'FaceAlpha',0.38, 'EdgeColor','none');

% If you have the means, plot them (better than (lo+hi)/2):
plot(x, m1, 'Color',[0 0.447 0.741], 'LineWidth',4)
plot(x, m2, 'Color',[0.85 0.325 0.098], 'LineWidth',4)
scatter(all_frex(find(P_above_thresh>0)),P_above_thresh(find(P_above_thresh>0)), 'r', '*')
xlim([0 freq_retained])
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 14)

xlabel('Frequency (Hz)', 'FontSize', 18)
ylabel('Explained variance (%)', 'FontSize', 18)
title(['Comp. ' num2str(comp2test) ' Eigenspectrum (group means \pm 95% CI)'], 'FontSize', 24)
legend({'Youth 95% CI', 'Elderly 95% CI','Youth mean','Elderly mean'}, 'Location','best'); box on
set(gca, 'XTick', all_frex(1:3:freq_retained));
set(gca, 'XTickLabels', string(all_frex(1:3:freq_retained)));
xtickangle(45)

%% Center of mass - both datasets combined

%% BOTH DATASETS COMBINED
%% combining cross frequencies all my info
freq_ind = [1:6 15 24 34 51:54 7:14 16:23 25:33 35:42];
num_freq = 40;
freq_ind_select = freq_ind(1:num_freq);
Dataset_switch = 1; % Change to upload connectivity matrices from a certain dataset 0 = TemoSeqAges, 1 = Dementia

num_part = 78;
coord_dir = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Comp_coord/TempSeq_MNI_cent_mass/Aging*');
all_freq_coord = zeros(num_freq,3,num_part);
for freq_looper = 1:num_freq %length(coord_dir)
    j = freq_ind_select(freq_looper);
    tmp = load([coord_dir(j).folder '/' coord_dir(j).name]);
    this_freq_coord = tmp.coord_part_info;
    all_freq_coord(freq_looper,:,:) = this_freq_coord;
end
%Y/O separation
participant_info = readtable('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/TSA2021_Nikita.xlsx'); % loading behavioural data
Valid_dat = participant_info(:,:); %defining a dataset with participants that have all the behavioural tests
Young = Valid_dat.Subject(Valid_dat.Age < 27);
Old = Valid_dat.Subject(Valid_dat.Age > 55);
Y_Selection = ismember(Valid_dat.Subject,Young);
O_Selection = ismember(Valid_dat.Subject,Old);
%Frequency indices
tmp = load('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/NetworkLandscape_Groups_TempSeq_Ages.mat', 'all_frex');
all_frex = tmp.all_frex;

all_freq_coord_TempSeq = all_freq_coord;
Y_Selection_TempSEQ = Y_Selection;
O_Selection_TempSEQ = O_Selection;

num_part = 100;
coord_dir = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/Comp_coord/AudMem_MNI_cent_mass/Result*');
all_freq_coord = zeros(num_freq,3,num_part);
for freq_looper = 1:num_freq %length(coord_dir)
    j = freq_ind_select(freq_looper);
    tmp = load([coord_dir(j).folder '/' coord_dir(j).name]);
    this_freq_coord = tmp.coord_part_info;
    all_freq_coord(freq_looper,:,:,:) = this_freq_coord;
end
%Y/O separation
participant_info = readtable('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/Descript_Dat/Participant_categories_AuditMemDement_ages.xlsx'); % loading behavioural data
% part 3,8,59,65 have no resting state GED
Val_part = ismember(participant_info.ID,[3 6 59 65]);
Valid_participant_info  = participant_info(~Val_part,:);
Y_Selection = ismember(Valid_participant_info.Group,{'A'});
O_Selection = ismember(Valid_participant_info.Group,{'B'});
%O_Selection = ismember(Valid_participant_info.Group,{'B','B_MCI'});
D_Selection = ismember(Valid_participant_info.Group,{'C'});

all_freq_coord_AudMem = all_freq_coord;
Y_Selection_AudMem = Y_Selection;
O_Selection_AudMem = O_Selection;

%Frequency indices
tmp = load('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/NetworkLandscape_Groups_2.mat', 'all_frex');
all_frex = tmp.all_frex;

all_freq_coord_BOTH = cat(3,all_freq_coord_TempSeq,all_freq_coord_AudMem);
Y_Selection_BOTH = [Y_Selection_TempSEQ;Y_Selection_AudMem];
O_Selection_BOTH = [O_Selection_TempSEQ;O_Selection_AudMem];

all_freq_coord = all_freq_coord_BOTH;
Y_Selection = Y_Selection_BOTH;
O_Selection = O_Selection_BOTH;
%%
%% Computing group difference and plotting CENTER OF MASS
stats_storage = zeros((num_freq),3);
for freq = 1:num_freq
    for dim = 1:3
        cur_array = squeeze(all_freq_coord(freq,dim,:));
        cur_array_Y = cur_array(Y_Selection);
        cur_array_O = cur_array(O_Selection);
        [pval,h,tval] = ranksum(cur_array_Y,cur_array_O);
        stats_storage((freq),dim) = pval;
        
    end
end

%
stats_storage_fdr = zeros(size(stats_storage));
%%
for dim = 1:3
        [h,crit_p] = fdr_bh(stats_storage(:,dim),0.05,'pdep','yes');% fdr correction
        stats_storage_fdr(:,dim) = h;
        crit_p
        all_frex(h)
end
%%
% plotting
% create per frequency difference
diff_table = zeros(num_freq,3);
for freq = 1:num_freq
    for dim = 1:3
        cur_array = squeeze(all_freq_coord(freq,dim,:));
        cur_array_Y_mean = mean(cur_array(Y_Selection));
        cur_array_O_mean = mean(cur_array(O_Selection));
        diff_table(freq,dim) = cur_array_Y_mean - cur_array_O_mean;
    end
end
tmp = load('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/NetworkLandscape_Groups_2.mat', 'all_frex');
all_frex = tmp.all_frex;
figure
for dim = 1:3
    subplot(3,1,dim)
    plot(diff_table(:,dim))
    hold on
    plot(stats_storage_fdr(:,dim)*2)
    hold off
    set(gca, 'XTick', 1:num_freq);
    set(gca, 'XTickLabels', string(all_frex(1:num_freq)));
    xtickangle(45)
end

%%
tmp = load('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/RQAstuff/MNI152_8mm_coord_dyi.mat')
coord_MEG_MNI = tmp.MNI8;

%% Center of mass old young plot
n_freq_plot = 40;
Cmap = parula(n_freq_plot);
figure; hold on
for dim = 1:3
    if dim == 1
        coord_lab = 'X axis';
    elseif dim == 2
        coord_lab = 'Y axis';
    else
        coord_lab = 'Z axis';
    end
    x_all_Y = [];
    x_all_O = [];
    y_all_Y = [];
    y_all_O = [];
    c_all = [];
    for freq = 1:num_freq
        cur_array = squeeze(all_freq_coord(freq,dim,:));
        cur_array_Y = cur_array(Y_Selection);
        cur_array_O = cur_array(O_Selection);
        
        coord_vec = cur_array(O_Selection);
        
        x_shuffle_magnitude = 0.1;
        y_shuffle_magnitude = 0.001;
        %adding jitters for visual readability
        xcoord_shuffle = x_shuffle_magnitude*randn(size(coord_vec));
        %coord_vec = coord_vec+xcoord_shuffle;
        ycoord_Y = freq + y_shuffle_magnitude*randn(size(cur_array_Y))+0.13;
        ycoord_O = freq + y_shuffle_magnitude*randn(size(cur_array_O))-0.13;
        
        %appending to one combined vector
        x_all_Y = [x_all_Y;cur_array_Y];
        x_all_O = [x_all_O;cur_array_O];
        y_all_Y = [y_all_Y;ycoord_Y];
        y_all_O = [y_all_O;ycoord_O];
        %third vector stores color code for each dot of the scatterplot (same dimentionality as my combined vectors)
        c_all = [c_all; repmat(Cmap(freq,:),numel(coord_vec),1)];

    end
    
    subplot(1,3,dim)
    scatter(x_all_Y, y_all_Y, 5, 'b', 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);
    hold on
    scatter(x_all_O, y_all_O, 5, 'r', 'filled','MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);
    stats_storage_fdr_thisdim = stats_storage_fdr(:,dim);
    scatter(stats_storage_fdr_thisdim(find(stats_storage_fdr_thisdim>0))*min(coord_MEG_MNI(:,dim)), find(stats_storage_fdr_thisdim>0), 'r', '*');
    hold off
    xlim([min(coord_MEG_MNI(:,dim)) max(coord_MEG_MNI(:,dim))]);
    
    ylim([min(y_all_Y)-0.2 max(y_all_Y)]);
    set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 14)
    set(gca, 'YTick', 1:n_freq_plot);
    set(gca, 'YTickLabels', '');
    xlabel(['Coordinate ' coord_lab], 'FontSize', 18)
    
    if dim == 1
        set(gca, 'YTickLabels', string(all_frex(1:num_freq)));
        ylabel('Frequency', 'FontSize', 18)
    end
end

%% and not the same for canonical frequency bands CENTER OF MASS

%split coord info between young and old
all_freq_coord_Y = squeeze(all_freq_coord(:,:,Y_Selection));
all_freq_coord_O = squeeze(all_freq_coord(:,:,O_Selection));

band_limits =[];

delta = 1:7;
theta = 8:11;
aplha = 12:15;
beta = 16:30;
gamma = 31:40;
band_limits{1} = delta;
band_limits{2} = theta;
band_limits{3} = aplha;
band_limits{4} = beta;
band_limits{5} = gamma;

band_pval_store = zeros(5,6);
for band_range = 1:5
    curr_range = band_limits{band_range};
    %young
    npart = size(all_freq_coord_Y,3);
    store_fit_Y = zeros(npart,3,2);
    for part = 1:npart
        for dim = 1:3
            cur_array = squeeze(all_freq_coord_Y(curr_range,dim,part));
            x = 1:numel(cur_array);
            x = x-mean(x);
            p = polyfit(x, cur_array',1);
            store_fit_Y(part,dim,1) = p(1);
            store_fit_Y(part,dim,2) = p(2);
        end
    end
    
    %old
    npart = size(all_freq_coord_O,3);
    store_fit_O = zeros(npart,3,2);
    for part = 1:npart
        for dim = 1:3
            cur_array = squeeze(all_freq_coord_O(curr_range,dim,part));
            x = 1:numel(cur_array);
            x = x-mean(x);
            p = polyfit(x, cur_array',1);
            store_fit_O(part,dim,1) = p(1);
            store_fit_O(part,dim,2) = p(2);
        end
    end
    slope_intersect_swticth = 2;
    [pval,~,~] = ranksum(store_fit_Y(:,1,slope_intersect_swticth),store_fit_O(:,1,slope_intersect_swticth));
    band_pval_store(band_range,1) = pval;
    band_pval_store(band_range,2) = median(squeeze(store_fit_Y(:,1,slope_intersect_swticth)))-median(squeeze(store_fit_O(:,1,slope_intersect_swticth)));
    [pval,~,~] = ranksum(store_fit_Y(:,2,slope_intersect_swticth),store_fit_O(:,2,slope_intersect_swticth));
    band_pval_store(band_range,3) = pval;
    band_pval_store(band_range,4) = median(squeeze(store_fit_Y(:,2,slope_intersect_swticth)))-median(squeeze(store_fit_O(:,2,slope_intersect_swticth)));
    [pval,~,~] = ranksum(store_fit_Y(:,3,slope_intersect_swticth),store_fit_O(:,3,slope_intersect_swticth));
    band_pval_store(band_range,5) = pval;
    band_pval_store(band_range,6) = median(squeeze(store_fit_Y(:,3,slope_intersect_swticth)))-median(squeeze(store_fit_O(:,3,slope_intersect_swticth)));
    
end

%% Individual alpha peaks and post-alpha decay - both datasets combined

%%
path_GEDoutput = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/'; % This loads output from previous block
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita')
% Load data

landscape = load([path_GEDoutput 'NetworkLandscape_Groups_TempSeq_Ages.mat'])
topEvals_all_group = landscape.topEvals_all_group;
all_frex = landscape.all_frex;
idx_sort = landscape.idx_sort;
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:);
Y_subNum = 37;
O_subNum = 40;
Dat_num = 1;
comp_incl = 3559;
freq_retained = 42;
Y_FREQ_TS = [];
Y_topEvals_TS = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
O_topEvals_TS = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
Y_holder_TS = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
Y_holder_TS = permute(Y_holder_TS,[1 3 2]);
Y_FREQ_TS.evals = Y_holder_TS;
Y_FREQ_TS.frex = all_frex(1:freq_retained);
O_FREQ_TS = [];
O_holder_TS = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
O_holder_TS = permute(O_holder_TS,[1 3 2]);
O_FREQ_TS.evals = O_holder_TS;
O_FREQ_TS.frex = all_frex(1:freq_retained);

landscape = load([path_GEDoutput 'NetworkLandscape_Groups_2.mat'])
topEvals_all_group = landscape.topEvals_all_group;
all_frex = landscape.all_frex;
idx_sort = landscape.idx_sort;
topEvals_all_group = topEvals_all_group(:,:,idx_sort,:);
Y_subNum = 48;
O_subNum = 40;
Dat_num = 2;
Y_FREQ_AMD = [];
Y_topEvals_AMD = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
O_topEvals_AMD = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
Y_holder_AMD = topEvals_all_group(1:comp_incl,1:Y_subNum,1:freq_retained,1);
Y_holder_AMD = permute(Y_holder_AMD,[1 3 2]);
Y_FREQ_AMD.evals = Y_holder_AMD;
Y_FREQ_AMD.frex = all_frex(1:freq_retained);
O_FREQ_AMD = [];
O_holder_AMD = topEvals_all_group(1:comp_incl,1:O_subNum,1:freq_retained,2);
O_holder_AMD = permute(O_holder_AMD,[1 3 2]);
O_FREQ_AMD.evals = O_holder_AMD;
O_FREQ_AMD.frex = all_frex(1:freq_retained);

Y_FREQ_BOTH = [];
Y_holder_BOTH = cat(3,Y_holder_TS,Y_holder_AMD);
Y_FREQ_BOTH.evals = Y_holder_BOTH;
Y_FREQ_BOTH.frex = all_frex(1:freq_retained);

O_FREQ_BOTH = [];
O_holder_BOTH = cat(3,O_holder_TS,O_holder_AMD);
O_FREQ_BOTH.evals = O_holder_BOTH;
O_FREQ_BOTH.frex = all_frex(1:freq_retained);

Y_FREQ_BOTH = Y_FREQ_BOTH.evals;
O_FREQ_BOTH = O_FREQ_BOTH.evals;
%% Alpha peaks and decay 
comp_of_int = 1;
search_window_start = 11;
search_window_end = 16;
%finding top values areounf the alpha peak
peak_val_Y  =  max(squeeze(Y_FREQ_BOTH(comp_of_int,search_window_start:search_window_end,:)));
peak_val_O  = max(squeeze(O_FREQ_BOTH(comp_of_int,search_window_start:search_window_end,:)));

a_peak_Y = zeros(size(peak_val_Y));
for i = 1:length(a_peak_Y)
    landsc = squeeze(Y_FREQ_BOTH(comp_of_int,search_window_start:search_window_end,i));
    a_peak_Y(i) = find(landsc == peak_val_Y(i))+search_window_start-1;
end

a_peak_O = zeros(size(peak_val_O));
for i = 1:length(a_peak_O)
    landsc = squeeze(O_FREQ_BOTH(comp_of_int,search_window_start:search_window_end,i));
    a_peak_O(i) = find(landsc == peak_val_O(i))+search_window_start-1;
end

figure
histogram(a_peak_Y)
hold on
histogram(a_peak_O)
hold off
set(gcf, 'Color', 'w')

addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/RQAstuff/Scripts')
%%
[pval,h,tval] = ranksum(a_peak_Y,a_peak_O)
%%
a_peak_Y_Freq = all_frex(a_peak_Y);
a_peak_O_Freq = all_frex(a_peak_O);
figure
histogram(a_peak_Y_Freq)
hold on
histogram(a_peak_O_Freq)
hold off
set(gcf, 'Color', 'w')

%%
pval_stor = zeros(6,10);
for length_fit = 5%4:13
    Gr = 1;
    if Gr == 1
        a_peak = a_peak_Y;
        nsubs = size(Y_FREQ_BOTH,3);
    else
        a_peak = a_peak_O;
        nsubs = size(O_FREQ_BOTH,3);
    end
    decayCoeff = nan(nsubs,1);
    goodFit.R2 = nan(nsubs,1);
    A_sub      = nan(nsubs,1);   % store amplitudes for optional individual plots
    
    for subi = 1:nsubs
        currr_alpha = a_peak(subi);
        idx_last = currr_alpha +length_fit;
        y = Y_FREQ_BOTH(comp_of_int,currr_alpha:idx_last,subi);
        
        % Keep only positive, finite values
        mask = y > 0 & isfinite(y);
        if sum(mask) < 3
            % Not enough data points to fit reliably
            continue
        end
        
        x_fit = currr_alpha:idx_last;
        x_sub = x_fit(mask);
        y_sub = y(mask);

        % Linearize: log(y) = log(A) - lambda * x
        logy = log(y_sub);
        X    = [x_sub(:) ones(numel(x_sub),1)];
        
        beta = X \ logy';          % [slope; intercept]
        lambda = -beta(1);        % decay coefficient
        A      = exp(beta(2));    % amplitude (not returned, but used for plotting)
        
        % Store decay coefficient
        decayCoeff(subi) = lambda;
        A_sub(subi)      = A;
        
        % Compute R² in the original (non-log) space
        y_hat = A * exp(-lambda * x_sub);
        sse   = sum((y_sub - y_hat).^2);
        sst   = sum((y_sub - mean(y_sub)).^2);
        if sst > 0
            goodFit.R2(subi) = 1 - sse/sst;
        else
            goodFit.R2(subi) = NaN;
        end
    end
    dec_Y = decayCoeff;
    med_dec_Y = median(dec_Y,'omitnan');
    %
    Gr = 2;
    if Gr == 1
        a_peak = a_peak_Y;
        nsubs = size(Y_FREQ_BOTH,3);
    else
        a_peak = a_peak_O;
        nsubs = size(O_FREQ_BOTH,3);
    end
    decayCoeff = nan(nsubs,1);
    goodFit.R2 = nan(nsubs,1);
    A_sub      = nan(nsubs,1);   % store amplitudes for optional individual plots
    
    for subi = 1:nsubs
        currr_alpha = a_peak(subi);
        idx_last = currr_alpha +length_fit;
        y = O_FREQ_BOTH(comp_of_int,currr_alpha:idx_last,subi);
        
        % Keep only positive, finite values
        mask = y > 0 & isfinite(y);
        if sum(mask) < 3
            % Not enough data points to fit reliably
            continue
        end
        
        x_fit = currr_alpha:idx_last;
        x_sub = x_fit(mask);
        y_sub = y(mask);
        
        % Linearize: log(y) = log(A) - lambda * x
        logy = log(y_sub);
        X    = [x_sub(:) ones(numel(x_sub),1)];
        
        beta = X \ logy';          % [slope; intercept]
        lambda = -beta(1);        % decay coefficient
        A      = exp(beta(2));    % amplitude (not returned, but used for plotting)
        
        % Store decay coefficient
        decayCoeff(subi) = lambda;
        A_sub(subi)      = A;
        
        % Compute R² in the original (non-log) space
        y_hat = A * exp(-lambda * x_sub);
        sse   = sum((y_sub - y_hat).^2);
        sst   = sum((y_sub - mean(y_sub)).^2);
        if sst > 0
            goodFit.R2(subi) = 1 - sse/sst;
        else
            goodFit.R2(subi) = NaN;
        end
    end
    dec_O = decayCoeff;
    med_dec_O = median(dec_O,'omitnan');
    %
    %figure
    %histogram(dec_O)
    % hold on
    %  histogram(dec_O)
    %  hold off
    pos_stor = length_fit - 3;
    [pval,h,tval] = ranksum(dec_Y,dec_O);
    pval_stor(1,pos_stor) = pval;
    pval_stor(2,pos_stor) = tval.zval;
    pval_stor(3,pos_stor) = tval.ranksum;
    [pval,h,tval] = ranksum(a_peak_Y,a_peak_O);
    pval_stor(4,pos_stor) = pval;
    pval_stor(5,pos_stor) = tval.zval;
    pval_stor(6,pos_stor) = tval.ranksum;
    
end
disp('DONE')
%%
figure
histogram(dec_Y)
hold on
histogram(dec_O)
hold off
set(gcf, 'Color', 'w')
