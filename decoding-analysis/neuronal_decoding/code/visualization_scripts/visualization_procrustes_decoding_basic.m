%{
Filename: visualization_procrustes_decoding_basic.m
Author: Zitong Wang
Date: 2025-07-10

Description:
    Visualizes Procrustes decoding results for V1 and V2.
    Loads pre‑computed results, combines complementary stimulus pairs, computes
    means/SEMs, and plots rotation‑only decoding accuracies with error bars and a
    2 % chance reference line.
%}

clc; clear;
%%
monkey = 'FR';
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_basic_results'];
%% Load result matrices
V1_acec_result = result_loader(monkey,'V1','acec',file_location);
V2_acec_result = result_loader(monkey,'V2','acec',file_location);

V1_ecac_result = result_loader(monkey,'V1','ecac',file_location);
V2_ecac_result = result_loader(monkey,'V2','ecac',file_location);

V1_acex_result = result_loader(monkey,'V1','acex',file_location);
V2_acex_result = result_loader(monkey,'V2','acex',file_location);

V1_exac_result = result_loader(monkey,'V1','exac',file_location);
V2_exac_result = result_loader(monkey,'V2','exac',file_location);

V1_ecex_result = result_loader(monkey,'V1','ecex',file_location);
V2_ecex_result = result_loader(monkey,'V2','ecex',file_location);

V1_exec_result = result_loader(monkey,'V1','exec',file_location);
V2_exec_result = result_loader(monkey,'V2','exec',file_location);

%% Combine complementary stimulus pairs vertically
V1_acec_combined = [V1_acec_result ; V1_ecac_result];
V2_acec_combined = [V2_acec_result ; V2_ecac_result];
V1_acex_combined = [V1_acex_result ; V1_exac_result];
V2_acex_combined = [V2_acex_result ; V2_exac_result];
V1_ecex_combined = [V1_ecex_result ; V1_exec_result];
V2_ecex_combined = [V2_ecex_result ; V2_exec_result];

%% Compute mean & SEM
nRep = length(V1_acec_combined);
calcStats = @(mat)[mean(mat,1); std(mat,0,1)/sqrt(nRep)];
V1_acec_ms = calcStats(V1_acec_combined);
V2_acec_ms = calcStats(V2_acec_combined);
V1_acex_ms = calcStats(V1_acex_combined);
V2_acex_ms = calcStats(V2_acex_combined);
V1_ecex_ms = calcStats(V1_ecex_combined);
V2_ecex_ms = calcStats(V2_ecex_combined);

partial_idx = 6; % rotation‑only transformation index

%% Assemble accuracy & error matrices for bar plot
acc_group = [V1_acec_ms(1,1) V2_acec_ms(1,1) V1_acec_ms(1,2) V2_acec_ms(1,2) V1_acec_ms(1,partial_idx) V2_acec_ms(1,partial_idx) V1_acec_ms(1,4) V2_acec_ms(1,4) V1_acec_ms(1,8) V2_acec_ms(1,8); 
             V1_ecex_ms(1,1) V2_ecex_ms(1,1) V1_ecex_ms(1,2) V2_ecex_ms(1,2) V1_ecex_ms(1,partial_idx) V2_ecex_ms(1,partial_idx) V1_ecex_ms(1,4) V2_ecex_ms(1,4) V1_ecex_ms(1,8) V2_ecex_ms(1,8); 
             V1_acex_ms(1,1) V2_acex_ms(1,1) V1_acex_ms(1,2) V2_acex_ms(1,2) V1_acex_ms(1,partial_idx) V2_acex_ms(1,partial_idx) V1_acex_ms(1,4) V2_acex_ms(1,4) V1_acex_ms(1,8) V2_acex_ms(1,8)];

er_group = [V1_acec_ms(2,1) V2_acec_ms(2,1) V1_acec_ms(2,2) V2_acec_ms(2,2) V1_acec_ms(2,partial_idx) V2_acec_ms(2,partial_idx) V1_acec_ms(2,4) V2_acec_ms(2,4) V1_acec_ms(2,8) V2_acec_ms(2,8); 
            V1_ecex_ms(2,1) V2_ecex_ms(2,1) V1_ecex_ms(2,2) V2_ecex_ms(2,2) V1_ecex_ms(2,partial_idx) V2_ecex_ms(2,partial_idx) V1_ecex_ms(2,4) V2_ecex_ms(2,4) V1_ecex_ms(2,8) V2_ecex_ms(2,8); 
            V1_acex_ms(2,1) V2_acex_ms(2,1) V1_acex_ms(2,2) V2_acex_ms(2,2) V1_acex_ms(2,partial_idx) V2_acex_ms(2,partial_idx) V1_acex_ms(2,4) V2_acex_ms(2,4) V1_acex_ms(2,8) V2_acex_ms(2,8)];

%% Plot bar graph with error bars and 2 % chance line
pairs = categorical({'EC-AC','EC-EX','AC-EX'});
pairs = reordercats(pairs, {'EC-AC','EC-EX','AC-EX'});   % <-- keep this order

figure('Position',[100 100 1080 540]);
b = bar(pairs, acc_group);

% Colors & styles
cols = repmat([0 0.4470 0.7410; 0.8500 0.3250 0.0980],5,1);
for i = 1:10, b(i).FaceColor = cols(i,:); end
b(3).LineStyle = '--'; b(4).LineStyle = '--';
b(7).LineStyle = ':';  b(8).LineStyle = ':';
b(9).LineStyle = '-.'; b(10).LineStyle = '-.';
for i = 3:10, b(i).LineWidth = 1.5; end

ylim([0 0.6]); yline(0.02,'--');
legend('V1 self-decoding','V2 self-decoding','V1 non-PT-decoding','V2 non-PT-decoding','V1 PT-decoding(only rotation)','V2 PT-decoding(only rotation)','V1 PT control','V2 PT control','V1 non-PT control','V2 non-PT control','theoretical chance level: 2%')
ylabel('Decoding accuracy'); xlabel('Stimulus pair');
title([monkey ' neuron‑shuffled, Z‑scored (rotation‑only)']);

% Error bars
[ng, nb] = size(acc_group);
x = nan(nb, ng);
for i = 1:nb, x(i,:) = b(i).XEndPoints; end
hold on; errorbar(x', acc_group, er_group,'r','linestyle','none','LineWidth',1,'HandleVisibility','off'); hold off;
%%
function organized_result = result_loader(monkey,vp,stimpair,file_location)
decoding_result = load(strcat(pwd,'\',file_location,'\',monkey,'\',vp,'\',stimpair,'_','results','.mat'),[stimpair,'_','results']); 
decoding_result = decoding_result.([stimpair,'_','results']);
organized_result = decoding_result{1};
end