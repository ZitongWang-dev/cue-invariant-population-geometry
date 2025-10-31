monkey='gabor_pyramid_resps';
neurons={'simple'};
csv_col_names=neurons;

fraction_values=[];
count_values=[];
n_samples=[];
for neuron_group=1:length(neurons)
    % group name to write in the csv file
    neuron=neurons{neuron_group};
    name_to_write=csv_col_names{neuron_group};
    
    %tuning_data_path='tuning_all\';
    tuning_data_path=[ monkey '\3x3_jitter_avg_area\3_scale_16_orientation\'];
    
    
    type_pairs={'ec','ac','ex'};
    
    for cnt=0:2
        stim_type=type_pairs{cnt+1};

        load([tuning_data_path 'odd' '\' stim_type '.mat'])
        tuning_1=center_resps';
        load([tuning_data_path 'even' '\' stim_type '.mat'])
        tuning_2=center_resps';
        [num_neurons, ~] = size(tuning_1);
        load([tuning_data_path 'odd' '\bad_channels.mat'])
        chosen_layer=1:num_neurons;
        
        % remove bad channels
        length(setdiff(chosen_layer,all_bad_channels))
        tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
        tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
        selected_channels=setdiff(chosen_layer,all_bad_channels);
        
        odd_pos = tuning_1';
        odd_pos(odd_pos<0) = 0;
        odd_neg = -1.*tuning_1';
        odd_neg(odd_neg<0) = 0;
        even_pos = tuning_2';
        even_pos(even_pos<0) = 0;
        even_neg = -1.*tuning_2';
        even_neg(even_neg<0) = 0;
        center_resps = [odd_pos,odd_neg,even_pos,even_neg];
        save(['gabor_pyramid_resps\3x3_jitter_avg_area\3_scale_16_orientation\' 'simple' '\' stim_type '.mat'],'center_resps')
    end
end