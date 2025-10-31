%{
Filename: visualization_procrustes_decoding_basic_incremental.m
Author: Zitong Wang
Date: 2025-07-10

Description:
    Visualizes incremental‑neuron Procrustes decoding results. Loads accuracy
    matrices for V1 and V2 across increasing neuron counts, plots mean ± SEM
    curves (self‑decoding vs. rotation‑only PT decoding), and extracts summary
    tables for quick reference.
%}
%% Works on: incremental decoding with partial transformations
clc
clear

%%
monkey = 'KO';
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_basic_incremental_results']
% Hard‑coded neuron lists per monkey/area
neuron_num_coode = { [48]; [48 112];        % FR‑V1 / FR‑V2
                     [109]; [109 146] };    % KO‑V1 / KO‑V2
name2idx = struct('FRV1',1,'FRV2',2,'KOV1',3,'KOV2',4);

%% Generate X‑axes (unique neuron counts actually present in results)
v1_neuron_num_list = neuron_num_coode{name2idx.([monkey 'V1'])};
V1_xaxis = [];
idx = 1;
while 10*2^(idx-1) < max(v1_neuron_num_list)
    V1_xaxis = [V1_xaxis 10*2^(idx-1)]; %#ok<AGROW>
    idx = idx+1;
end
V1_xaxis = unique([V1_xaxis v1_neuron_num_list]);

v2_neuron_num_list = neuron_num_coode{name2idx.([monkey 'V2'])};
V2_xaxis = [];
idx = 1;
while 10*2^(idx-1) < max(v2_neuron_num_list)
    V2_xaxis = [V2_xaxis 10*2^(idx-1)]; %#ok<AGROW>
    idx = idx+1;
end
V2_xaxis = unique([V2_xaxis v2_neuron_num_list]);

%% Load incremental‑results (each cell = one neuron count)
V1_acec = result_loader(monkey,'V1','acec',file_location);
V1_acex = result_loader(monkey,'V1','acex',file_location);
V1_ecac = result_loader(monkey,'V1','ecac',file_location);
V1_ecex = result_loader(monkey,'V1','ecex',file_location);
V1_exac = result_loader(monkey,'V1','exac',file_location);
V1_exec = result_loader(monkey,'V1','exec',file_location);

V2_acec = result_loader(monkey,'V2','acec',file_location);
V2_acex = result_loader(monkey,'V2','acex',file_location);
V2_ecac = result_loader(monkey,'V2','ecac',file_location);
V2_ecex = result_loader(monkey,'V2','ecex',file_location);
V2_exac = result_loader(monkey,'V2','exac',file_location);
V2_exec = result_loader(monkey,'V2','exec',file_location);

%% Combine matrices (rows concatenated, columns incremental neuron counts)
V1_combined = [cell2mat(V1_acec); cell2mat(V1_acex); cell2mat(V1_ecac); ...
               cell2mat(V1_ecex); cell2mat(V1_exac); cell2mat(V1_exec)];
V2_combined = [cell2mat(V2_acec); cell2mat(V2_acex); cell2mat(V2_ecac); ...
               cell2mat(V2_ecex); cell2mat(V2_exac); cell2mat(V2_exec)];

v1_incr = size(V1_combined,2)/4;   % 4 analysis types per neuron‑count
v2_incr = size(V2_combined,2)/4;
repN    = size(V2_combined,1);
V1_split = mat2cell(V1_combined, repN, ones(1,4)*v1_incr);
V2_split = mat2cell(V2_combined, repN, ones(1,4)*v2_incr);

%% Plot mean ± SEM curves (self vs PT rotation‑only)
figure('Position',[100 100 1080 540]); hold on;
lw = 1.5;
errorbar(V1_xaxis, mean(V1_split{1},1), std(V1_split{1},0,1)/sqrt(repN), ...
         'Color',[0 0.4470 0.7410],'LineWidth',lw);
errorbar(V1_xaxis, mean(V1_split{3},1), std(V1_split{3},0,1)/sqrt(repN), ...
         'Color',[0 0.4470 0.7410],'LineStyle','--','LineWidth',lw);
errorbar(V2_xaxis, mean(V2_split{1},1), std(V2_split{1},0,1)/sqrt(repN), ...
         'Color',[0.8500 0.3250 0.0980],'LineWidth',lw);
errorbar(V2_xaxis, mean(V2_split{3},1), std(V2_split{3},0,1)/sqrt(repN), ...
         'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',lw);
hold off;
legend('V1 self‑decoding','V1 PT‑decoding','V2 self‑decoding','V2 PT‑decoding','Location','northwest');
title([monkey ' averaged decoding accuracy (incremental neurons)']);
xlabel('Number of neurons'); ylabel('Decoding accuracy'); ylim([0 0.6]); xlim([0 max(V2_xaxis)+5]);

%% Extract mean & SEM arrays for quick reference
v1_self         = [mean(V1_split{1},1); std(V1_split{1},0,1)/sqrt(repN)];
v1_pt_rotation  = [mean(V1_split{3},1); std(V1_split{3},0,1)/sqrt(repN)];
v2_self         = [mean(V2_split{1},1); std(V2_split{1},0,1)/sqrt(repN)];
v2_pt_rotation  = [mean(V2_split{3},1); std(V2_split{3},0,1)/sqrt(repN)];

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