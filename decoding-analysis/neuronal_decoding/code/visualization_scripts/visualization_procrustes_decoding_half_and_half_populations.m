%{
Filename: visualization_procrustes_decoding_half_and_half_populations.m
Author: Zitong Wang
Date: 2025-07-10

Description:
    Visualizes half‑and‑half population Procrustes decoding. Loads results for
    AC, EC, EX stimuli (V1 & V2), computes mean/SEM, and plots decoding
    accuracies (rotation‑only PT plus controls) in a single bar graph.
%}

%% Works on: half‑half decoding with partial transformations
clc
clear
%% Parameters
monkey = 'KO';
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_half_and_half_populations_results']
%% Load result matrices
V1_ac_result = result_loader(monkey,'V1','ac',file_location);
V2_ac_result = result_loader(monkey,'V2','ac',file_location);

V1_ec_result = result_loader(monkey,'V1','ec',file_location);
V2_ec_result = result_loader(monkey,'V2','ec',file_location);

V1_ex_result = result_loader(monkey,'V1','ex',file_location);
V2_ex_result = result_loader(monkey,'V2','ex',file_location);

%% Aliases for clarity (no data changes)
V1_ac = V1_ac_result;  V2_ac = V2_ac_result;
V1_ec = V1_ec_result;  V2_ec = V2_ec_result;
V1_ex = V1_ex_result;  V2_ex = V2_ex_result;

total_repeat = length(V1_ac);
calc = @(mat)[mean(mat,1); std(mat,0,1)/sqrt(total_repeat)];
V1_ac_ms = calc(V1_ac); V2_ac_ms = calc(V2_ac);
V1_ex_ms = calc(V1_ex); V2_ex_ms = calc(V2_ex);
V1_ec_ms = calc(V1_ec); V2_ec_ms = calc(V2_ec);

partial_idx = 6; % rotation‑only column
acc_group = [V1_ac_ms(1,1) V2_ac_ms(1,1) V1_ac_ms(1,2) V2_ac_ms(1,2) V1_ac_ms(1,partial_idx) V2_ac_ms(1,partial_idx) V1_ac_ms(1,4) V2_ac_ms(1,4) V1_ac_ms(1,8) V2_ac_ms(1,8); 
             V1_ex_ms(1,1) V2_ex_ms(1,1) V1_ex_ms(1,2) V2_ex_ms(1,2) V1_ex_ms(1,partial_idx) V2_ex_ms(1,partial_idx) V1_ex_ms(1,4) V2_ex_ms(1,4) V1_ex_ms(1,8) V2_ex_ms(1,8); 
             V1_ec_ms(1,1) V2_ec_ms(1,1) V1_ec_ms(1,2) V2_ec_ms(1,2) V1_ec_ms(1,partial_idx) V2_ec_ms(1,partial_idx) V1_ec_ms(1,4) V2_ec_ms(1,4) V1_ec_ms(1,8) V2_ec_ms(1,8)];

er_group = [V1_ac_ms(2,1) V2_ac_ms(2,1) V1_ac_ms(2,2) V2_ac_ms(2,2) V1_ac_ms(2,partial_idx) V2_ac_ms(2,partial_idx) V1_ac_ms(2,4) V2_ac_ms(2,4) V1_ac_ms(2,8) V2_ac_ms(2,8); 
            V1_ex_ms(2,1) V2_ex_ms(2,1) V1_ex_ms(2,2) V2_ex_ms(2,2) V1_ex_ms(2,partial_idx) V2_ex_ms(2,partial_idx) V1_ex_ms(2,4) V2_ex_ms(2,4) V1_ex_ms(2,8) V2_ex_ms(2,8); 
            V1_ec_ms(2,1) V2_ec_ms(2,1) V1_ec_ms(2,2) V2_ec_ms(2,2) V1_ec_ms(2,partial_idx) V2_ec_ms(2,partial_idx) V1_ec_ms(2,4) V2_ec_ms(2,4) V1_ec_ms(2,8) V2_ec_ms(2,8)];

%% Plot bar graph
cats = categorical({'AC','EX','EC'});
cats = reordercats(cats,{'AC','EX','EC'});
figure('Position',[100 100 1080 540]);
b = bar(cats, acc_group);

cols = repmat([0 0.4470 0.7410; 0.8500 0.3250 0.0980],5,1);
for i = 1:10, b(i).FaceColor = cols(i,:); end
b(3).LineStyle = '--'; b(4).LineStyle = '--';
b(7).LineStyle = ':';  b(8).LineStyle = ':';
b(9).LineStyle = '-.';  b(10).LineStyle = '-.';
for i = 3:10, b(i).LineWidth = 1.5; end

ylim([0 0.6]); yline(0.02,'--');
legend({'V1 self','V2 self','V1 non‑PT','V2 non‑PT','V1 rot‑only','V2 rot‑only','V1 PT‑ctrl','V2 PT‑ctrl','V1 non‑PT‑ctrl','V2 non‑PT‑ctrl','2 % chance'},'Location','best');
ylabel('Decoding accuracy'); xlabel('Stimulus');
title([monkey ' half‑&‑half decoding (rotation‑only)']);

% Error bars
[ng, nb] = size(acc_group);
x = nan(nb, ng);
for i = 1:nb, x(i,:) = b(i).XEndPoints; end
hold on; errorbar(x', acc_group, er_group,'r','linestyle','none','LineWidth',1,'HandleVisibility','off'); hold off;

%% Helper
function out = result_loader(monkey,vp,stim,fileloc)
    tmp = load([pwd '\' fileloc '\' monkey '\' vp '\' stim '_results.mat'], [stim '_results']);
    out = tmp.([stim '_results']);
end
