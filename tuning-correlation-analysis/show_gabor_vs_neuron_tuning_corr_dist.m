%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gabor

monkey='gabor_pyramid_resps';
neurons={'complex'};
gabor_type = neurons{1};
for neuron_group=1:length(neurons)
    % group name to write in the csv file
    neuron=neurons{neuron_group};
    
    %tuning_data_path='tuning_all\';
    %tuning_data_path=[monkey '\no_jitter\3_scale_16_orientation\'];
    tuning_data_path=[monkey '\5x5_jitter_avg_area\3_scale_16_orientation\'];
    
    
    type_pairs={{'ec','ac'},{'ec','ex'},{'ac','ex'}};
    sigs=[];
    
    for cnt=0:2
        stim_type=type_pairs{cnt+1};

        load([tuning_data_path neuron '\' stim_type{1} '.mat'])
        tuning_1=center_resps';
        load([tuning_data_path neuron '\' stim_type{2} '.mat'])
        tuning_2=center_resps';
        load([tuning_data_path neuron '\bad_channels.mat'])
        [num_neurons, ~] = size(tuning_1);
        chosen_layer=1:num_neurons;
        
        % remove bad channels
        length(setdiff(chosen_layer,all_bad_channels))
        tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
        tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
        selected_channels=setdiff(chosen_layer,all_bad_channels);
        
%         % remove 1 outlier
%         tuning_geo_product=tuning_1.*tuning_2;
%         [max_vals,max_indices]=max(tuning_geo_product,[],2);
%         tmp_tuning_1=[];
%         tmp_tuning_2=[];
%         for neuron_idx=1:length(max_indices)
%             curr_max_idx = max_indices(neuron_idx);
%             tmp_row = tuning_1(neuron_idx,:);
%             tmp_row(curr_max_idx)=[];
%             tmp_tuning_1=[tmp_tuning_1;tmp_row];
%             tmp_row = tuning_2(neuron_idx,:);
%             tmp_row(curr_max_idx)=[];
%             tmp_tuning_2=[tmp_tuning_2;tmp_row];
%         end
%         tuning_1=tmp_tuning_1;
%         tuning_2=tmp_tuning_2;
%         % remove outlier end
        
        
        [tuning_corr,p_vals]=corr(tuning_1',tuning_2','Type','Pearson');
        
        ident_diag=diag(tuning_corr,0);
        x_edge=-1:0.05:1;
         
%         % plot histograms - all
%         subplot(3,3,1+cnt)
%         ident_diag=diag(tuning_corr,0);
%         x_edge=-1:0.05:1;
%         histogram(ident_diag,x_edge)
%         xlim([-1,1])
%         ylim([0,30])
%         xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)))), ' count=' num2str(length(selected_channels))])
%         ylabel('Counts')
%         title([ neuron ' Gabor ' stim_type{1} ' vs ' stim_type{2}])
        
        % significant correlations
        curr_sigs=diag(p_vals,0)<=0.05;
        
        % negative significant correlations
        curr_sigs_neg=(diag(p_vals,0)<=0.05)&(ident_diag<0);
        
        % positive significant correlations
        curr_sigs_pos=(diag(p_vals,0)<=0.05)&(ident_diag>0);
        % plot histograms - pos sig only
        subplot(3,3,1+cnt)
        histogram(ident_diag(curr_sigs_pos),x_edge)
        xlim([-1,1])
        ylim([0,30])
        xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)&curr_sigs_pos))) ' count=' num2str(length(curr_sigs_pos(curr_sigs_pos>0))) '/' num2str(length(selected_channels)) '(' num2str(length(curr_sigs_pos(curr_sigs_pos>0)) / length(selected_channels)) ')'])
        ylabel('Counts')
        title([monkey ' ' neuron ' ' stim_type{1} ' vs ' stim_type{2}])
        
        sigs_pos_channels=selected_channels(curr_sigs_pos);
        sigs=[sigs, curr_sigs_pos];
    end
    %set(gcf,'position',[0,0,1000,500])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corr between type

datas_ko=[
    %v1
    {'KO041015'}
    {'KO041515'}
    {'KO042115'}
    {'KO043015'}
    {'KO060315'}
    %v2
    {'KO031915'}
    {'KO041715'}
    {'KO042515'}
    {'KO042815'}
    {'KO050115'}
    {'KO060815'}
    {'KO060915'}
    ];

datas_fr=[
    %v1
    {'FR041014'}
    {'FR081414'}
    {'FR082114'}
    {'FR082814'}
    {'FR082914'}
    %v2
    {'FR041814'}
    {'FR060414'}
    {'FR060514'}
    {'FR072614'}
    {'FR073014'}
    {'FR090414'}
    {'FR090914'}
    ];

monkey='ko';
datas=datas_ko;
t_window='40-250';
neurons={'v1','v2'};
csv_col_names={'v1_all','v2_all'};

fraction_values=[];
count_values=[];
n_samples=[];
for neuron_group=1:length(neurons)
    % group name to write in the csv file
    neuron=neurons{neuron_group};
    name_to_write=csv_col_names{neuron_group};
    
    tuning_data_path='tuning_all\';
    
    type_pairs={{'ec','ac'},{'ec','ex'},{'ac','ex'}};
    sigs=[];
    
    for cnt=0:2
        stim_type=type_pairs{cnt+1};
        
        % load stimulus tuning data
        neurons_layer_d=[];
        neurons_layer_4=[];
        neurons_layer_s=[];
        neurons_layer_all=[];
        if strcmp(neuron,'v1')
            %v1
            for i=1:5
                data=datas{i};
                load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s1.mat']);
                if strcmp(monkey,'fr')
                    layer_4=l4;
                end
                neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24];
                neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
                neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24];
                load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s2.mat']);
                if strcmp(monkey,'fr')
                    layer_4=l4;
                end
                neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24+12];
                neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
                neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24+12];
                neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
            end
        else
            %v2
            for i=6:12
                data=datas{i};
                load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s1.mat']);
                if strcmp(monkey,'fr')
                    layer_4=l4;
                end
                neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24];
                neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
                neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24];
                load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s2.mat']);
                if strcmp(monkey,'fr')
                    layer_4=l4;
                end
                neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24+12];
                neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
                neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24+12];
                neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
            end
        end
        
        % select the neuron group
        % s: superficial
        % 4: layer_4
        % d: deep
        % all: all layers combined
        chosen_layer=neurons_layer_all;
        
        % read tuning data
        load([tuning_data_path t_window '\' monkey '_' stim_type{1} '_all_tuning.mat'])
        tuning_1=all_tuning;
        load([tuning_data_path t_window '\' monkey '_' stim_type{2} '_all_tuning.mat'])
        tuning_2=all_tuning;
        
        % remove bad channels
        length(setdiff(chosen_layer,all_bad_channels))
        tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
        tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
        selected_channels=setdiff(chosen_layer,all_bad_channels);
        
%         % remove 1 outlier
%         tuning_geo_product=tuning_1.*tuning_2;
%         [max_vals,max_indices]=max(tuning_geo_product,[],2);
%         tmp_tuning_1=[];
%         tmp_tuning_2=[];
%         for neuron_idx=1:length(max_indices)
%             curr_max_idx = max_indices(neuron_idx);
%             tmp_row = tuning_1(neuron_idx,:);
%             tmp_row(curr_max_idx)=[];
%             tmp_tuning_1=[tmp_tuning_1;tmp_row];
%             tmp_row = tuning_2(neuron_idx,:);
%             tmp_row(curr_max_idx)=[];
%             tmp_tuning_2=[tmp_tuning_2;tmp_row];
%         end
%         tuning_1=tmp_tuning_1;
%         tuning_2=tmp_tuning_2;
%         % remove outlier end
        
        
        [tuning_corr,p_vals]=corr(tuning_1',tuning_2','Type','Pearson');
        ident_diag=diag(tuning_corr,0);
        x_edge=-1:0.05:1;
        
%         % plot histograms - all
%         subplot(3,3,1+neuron_group*3+cnt)
%         ident_diag=diag(tuning_corr,0);
%         x_edge=-1:0.05:1;
%         histogram(ident_diag,x_edge)
%         xlim([-1,1])
%         ylim([0,30])
%         xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)))) ' count=' num2str(length(selected_channels))])
%         ylabel('Counts')
%         title([monkey ' ' neuron ' ' stim_type{1} ' vs ' stim_type{2}])
        
        % significant correlations
        curr_sigs=diag(p_vals,0)<=0.05;
        
        % negative significant correlations
        curr_sigs_neg=(diag(p_vals,0)<=0.05)&(ident_diag<0);
        
        % positive significant correlations
        curr_sigs_pos=(diag(p_vals,0)<=0.05)&(ident_diag>0);
        % plot histograms - pos sig only
        subplot(3,3,1+neuron_group*3+cnt)
        histogram(ident_diag(curr_sigs_pos),x_edge)
        xlim([-1,1])
        ylim([0,30])
        xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)&curr_sigs_pos))) ' count=' num2str(length(curr_sigs_pos(curr_sigs_pos>0))) '/' num2str(length(selected_channels)) '(' num2str(length(curr_sigs_pos(curr_sigs_pos>0)) / length(selected_channels)) ')' ])
        ylabel('Counts')
        title([monkey ' ' neuron ' ' stim_type{1} ' vs ' stim_type{2}])
        
        sigs_pos_channels=selected_channels(curr_sigs_pos);
        sigs=[sigs, curr_sigs_pos];
    end
    set(gcf,'position',[0,0,1920,962])
    print([gabor_type '_' monkey '_pos_sig.eps'],'-depsc','-painters','-loose',gcf)
end