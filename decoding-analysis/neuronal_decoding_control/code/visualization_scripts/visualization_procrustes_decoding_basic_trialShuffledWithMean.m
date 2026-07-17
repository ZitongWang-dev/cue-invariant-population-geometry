%{
Filename: visualization_procrustes_decoding_basic_trialShuffledWithMean.m
Author: Zitong Wang
Date: 2026-07-16

Description:
    Visualizes the structure-destroying (correspondence-shuffled) Procrustes
    control produced by procrustes_decoding_basic_trialShuffledWithMean.m.
    Loads pre-computed results, combines complementary stimulus pairs, computes
    means/SEMs, and plots the key accuracy columns with error bars and a 2%
    chance reference line.

    Result columns (from procrustes_decoding_basic_trialShuffledWithMean.m):
      1 genAcc               - stim2 within-condition 10-fold CV accuracy (ceiling)
      2 acc_true_rot         - TRUE correspondence, rotation-only PT (vs 1:50)
      3 acc_shuffle_target   - shuffled correspondence, rotation-only, vs TARGET labels (1:50)
      4 acc_shuffle_identity - shuffled correspondence, rotation-only, vs SHUFFLED identity (pi)
      5 acc_notransform      - raw stim1 test into stim2 model, no transform (vs 1:50)
      6 acc_chance           - random-permutation chance baseline

    Reading the figure:
      acc_shuffle_target HIGH  -> the rotation can force an arbitrary correspondence
                                  onto the target configuration (overfitting signature).
      acc_shuffle_identity ~2% -> the trials' true pre-rotation identity is not preserved,
                                  i.e. the correspondence really was destroyed.
%}

clc; clear;
%% Configuration
monkey = 'KO';                 % 'FR' or 'KO'
areas  = {'V1','V2'};          % use {'V1'} to plot a single area

% NOTE: must match save_path in procrustes_decoding_basic_trialShuffledWithMean.m.
% If that leaf is renamed to '..._results' (sibling-pipeline convention), update here too.
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_basic_trialShuffledWithMean'];

% Column indices of the [N x 6] result matrices
col_self           = 1;
col_true_rot       = 2;
col_shuffle_target = 3;
col_shuffle_id     = 4;
col_notransform    = 5;
col_chance         = 6;

% Measures to plot (order = bar order within each stimulus pair)
meas_cols  = [col_self col_notransform col_true_rot col_shuffle_target col_shuffle_id];
meas_names = {'self-decoding', ...
              'non-PT-decoding', ...
              'PT-decoding(only rotation)', ...
              'shuffled-corr PT vs target labels', ...
              'shuffled-corr PT vs identity labels'};
meas_style = {'-','--','-',':','-.'};

% Complementary stimulus pairs, combined vertically
pair_defs  = {'acec','ecac'; 'ecex','exec'; 'acex','exac'};
pair_names = {'EC-AC','EC-EX','AC-EX'};

%% Load, combine complementary pairs, compute mean & SEM
n_area = numel(areas);
n_pair = size(pair_defs,1);
n_meas = numel(meas_cols);

acc_group   = zeros(n_pair, n_meas*n_area);
er_group    = zeros(n_pair, n_meas*n_area);
chance_mean = zeros(n_pair, n_area);   % empirical chance column, for reporting

for p = 1:n_pair
    for a = 1:n_area
        combined = [result_loader(monkey,areas{a},pair_defs{p,1},file_location) ; ...
                    result_loader(monkey,areas{a},pair_defs{p,2},file_location)];
        nRep = size(combined,1);
        m    = mean(combined,1);
        sem  = std(combined,0,1)/sqrt(nRep);
        for k = 1:n_meas
            bar_idx = (k-1)*n_area + a;               % measure-major, area-minor
            acc_group(p,bar_idx) = m(meas_cols(k));
            er_group(p,bar_idx)  = sem(meas_cols(k));
        end
        chance_mean(p,a) = m(col_chance);
    end
end

% Sanity report: empirical chance should sit at ~0.02
disp('Empirical chance column (rows = stimulus pairs, cols = areas):');
disp(array2table(chance_mean,'RowNames',pair_names,'VariableNames',areas));

%% Plot bar graph with error bars and 2% chance line
pairs = categorical(pair_names);
pairs = reordercats(pairs, pair_names);   % <-- keep this order

figure('Position',[100 100 1080 540]);
b = bar(pairs, acc_group);

% Colors by area, line styles by measure
area_cols = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
if n_area > size(area_cols,1), area_cols = lines(n_area); end

legend_labels = cell(1, n_meas*n_area);
for k = 1:n_meas
    for a = 1:n_area
        i = (k-1)*n_area + a;
        b(i).FaceColor = area_cols(a,:);
        b(i).LineStyle = meas_style{k};
        b(i).LineWidth = 1.5;
        legend_labels{i} = [areas{a} ' ' meas_names{k}];
    end
end

ylim([0 1]); yline(0.02,'--');
legend([legend_labels, {'theoretical chance level: 2%'}],'Location','eastoutside');
ylabel('Decoding accuracy'); xlabel('Stimulus pair');
title([monkey ' correspondence-shuffled control, Z-scored (rotation-only)']);

% Error bars
[ng, nb] = size(acc_group);
x = nan(nb, ng);
for i = 1:nb, x(i,:) = b(i).XEndPoints; end
hold on;
errorbar(x', acc_group, er_group,'r','linestyle','none','LineWidth',1,'HandleVisibility','off');
hold off;

% saveas(gcf, fullfile('..','..','results','figures', ...
%     ['procrustes_decoding_basic_trialShuffledWithMean_' monkey '.png']));

%%
function organized_result = result_loader(monkey,vp,stimpair,file_location)
decoding_result = load(strcat(pwd,'\',file_location,'\',monkey,'\',vp,'\',stimpair,'_','results','.mat'),[stimpair,'_','results']);
decoding_result = decoding_result.([stimpair,'_','results']);
organized_result = decoding_result{1};
end
