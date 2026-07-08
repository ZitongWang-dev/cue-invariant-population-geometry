%{
Filename: visualization_procrustes_decoding_half_populations_incremental.m
Author: Zitong Wang
Date: 2025-07-10

Description:
    Visualizes incremental half‑and‑half Procrustes decoding. Loads V1 & V2
    accuracy matrices across neuron counts, then plots mean ± SEM curves for
    self‑decoding vs. rotation‑only PT decoding.

%}

%% Works on: incremental half‑half decoding with partial transformations
clc
clear

%% Parameters
monkey = 'FR';
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_half_and_half_populations_incremental_results']
% Hard‑coded neuron‑number table
neuron_num_coode = { [48]; [48 112];   % FR‑V1 / FR‑V2
                     [109]; [109 146]}; % KO‑V1 / KO‑V2
name2idx = struct('FRV1',1,'FRV2',2,'KOV1',3,'KOV2',4);

%% X‑axis (unique neuron counts shown)
make_x = @(lst) unique([(1:floor(max(lst)/20))*10, floor(max(lst)/2), floor(min(lst)/2)]);
V1_xaxis = make_x(neuron_num_coode{name2idx.([monkey 'V1'])});
V2_xaxis = make_x(neuron_num_coode{name2idx.([monkey 'V2'])});

%% Load incremental results (each cell → one neuron count)
V1_ac = result_loader(monkey,'V1','ac',file_location);
V1_ec = result_loader(monkey,'V1','ec',file_location);
V1_ex = result_loader(monkey,'V1','ex',file_location);
V2_ac = result_loader(monkey,'V2','ac',file_location);
V2_ec = result_loader(monkey,'V2','ec',file_location);
V2_ex = result_loader(monkey,'V2','ex',file_location);

% Combine: rows concatenated, cols incremental counts × 4 metrics
V1_all = [cell2mat(V1_ac); cell2mat(V1_ec); cell2mat(V1_ex)];
V2_all = [cell2mat(V2_ac); cell2mat(V2_ec); cell2mat(V2_ex)];

v1_step = size(V1_all,2)/4; v2_step = size(V2_all,2)/4; repN = size(V2_all,1);
V1_split = mat2cell(V1_all, repN, ones(1,4)*v1_step);
V2_split = mat2cell(V2_all, repN, ones(1,4)*v2_step);

%% Plot curves (self vs. PT‑rotation)
figure('Position',[100 100 1080 540]); hold on;
lw = 1.5;
errorbar(V1_xaxis, mean(V1_split{1},1), std(V1_split{1},0,1)/sqrt(repN), 'Color',[0 0.4470 0.7410],'LineWidth',lw);
errorbar(V1_xaxis, mean(V1_split{3},1), std(V1_split{3},0,1)/sqrt(repN), 'Color',[0 0.4470 0.7410],'LineStyle','--','LineWidth',lw);
errorbar(V2_xaxis, mean(V2_split{1},1), std(V2_split{1},0,1)/sqrt(repN), 'Color',[0.8500 0.3250 0.0980],'LineWidth',lw);
errorbar(V2_xaxis, mean(V2_split{3},1), std(V2_split{3},0,1)/sqrt(repN), 'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',lw);
hold off;
legend('V1 self','V1 PT‑rot','V2 self','V2 PT‑rot','Location','northwest');
title([monkey ' incremental half‑half decoding']);
xlabel('Number of neurons'); ylabel('Decoding accuracy'); ylim([0 0.4]); xlim([0 max(V2_xaxis)+5]);


%%
function organized_result = result_loader(monkey,vp,stimpair,file_location)
decoding_result = load(strcat(pwd,'\',file_location,'\',monkey,'\',vp,'\',stimpair,'_','results','.mat'),[stimpair,'_','results']); 
decoding_result = decoding_result.([stimpair,'_','results']);


neuron_length = length(decoding_result);
organized_result = cell(1,4); % 1 model acc; 2 no-trans decoding; 3 trans-decoding(only rotation); 4 PT-control
for i =1:neuron_length
   organized_result{1} = [organized_result{1} decoding_result{i}(:,1)];
   organized_result{2} = [organized_result{2} decoding_result{i}(:,2)];
   organized_result{3} = [organized_result{3} decoding_result{i}(:,6)]; % 6th column in the data
   organized_result{4} = [organized_result{4} decoding_result{i}(:,4)];
end

end