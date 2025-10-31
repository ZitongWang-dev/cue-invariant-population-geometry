% stimtype='ec';
% num_blocks=2;
% gabor_size_comp=21;
% gabor_size_orient=21;
% load(['descriptor_scores/complexity/stimulus_' stimtype '_complexity_axis_' num2str(num_blocks^2) '_blocks_even_' num2str(gabor_size_comp) '.mat'])
% 
% [~,cx_rankorder]=sort(complexity_scores,'descend');
% cx_idx=zeros(1,50);
% for i=1:50
%     cx_idx(cx_rankorder(i))=51-i;
% end
% 
% % datas=[
% %     %v1
% %     {'KO041015'}
% %     {'KO041515'}
% %     {'KO042115'}
% %     {'KO043015'}
% %     {'KO060315'}
% %     %v2
% %     {'KO031915'} 
% %     {'KO041715'}
% %     {'KO042515'}
% %     {'KO042815'}
% %     {'KO050115'}
% %     {'KO060815'}
% %     {'KO060915'}
% % ]
% % 
% % neurons_layer_d=[];
% % neurons_layer_4=[];
% % neurons_layer_s=[];
% % neurons_layer_all=[];
% % 
% % for i=1:5
% %     data=datas{i};
% %     load(['layers_ko\KO_' data(3:end) '_s1.mat']);
% %     neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24];
% %     neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
% %     neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24];
% %     load(['layers_ko\KO_' data(3:end) '_s2.mat']);
% %     neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24+12];
% %     neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
% %     neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24+12];
% %     neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
% % end
% % 
% % % for i=6:12
% % %     data=datas{i};
% % %     load(['layers_ko\KO_' data(3:end) '_s1.mat']);
% % %     neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24];
% % %     neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
% % %     neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24];
% % %     load(['layers_ko\KO_' data(3:end) '_s2.mat']);
% % %     neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24+12];
% % %     neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
% % %     neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24+12];
% % %     neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
% % % end
% % 
% % monkey='KO';
% % stim_type='ex';
% % setting='v1_all';
% % neuron_group='V1';
% % t_window='110-180';
% % combined_tuning=[];
% % load(['tuning_all\' t_window '\ko_' stim_type '_all_tuning.mat'])
% %  all_tuning=all_tuning(setdiff(neurons_layer_all,all_bad_channels),:);
% %  %all_tuning(all_bad_channels,:)=[];
% %  all_tuning=all_tuning';
% % %all_tuning(all_bad_channels,:)=[];
% % combined_tuning=[combined_tuning,all_tuning];
% 
% datas=[
%     %v1
%     {'FR041014'}
%     {'FR081414'}
%     {'FR082114'}
%     {'FR082814'}
%     {'FR082914'}
%     %v2
%     {'FR041814'}
%     {'FR060414'}
%     {'FR060514'}
%     {'FR072614'}
%     {'FR073014'}
%     {'FR090414'}
%     {'FR090914'}
% ]
% 
% neurons_layer_d=[];
% neurons_layer_4=[];
% neurons_layer_s=[];
% neurons_layer_all=[];
% % 
% % for i=1:5
% %     data=datas{i};
% %     load(['layers_fr\FR_' data(3:end) '_s1.mat']);
% %     neurons_layer_d=[neurons_layer_d (ld)+(i-1)*24];
% %     neurons_layer_4=[neurons_layer_4 (l4)+(i-1)*24];
% %     neurons_layer_s=[neurons_layer_s (ls)+(i-1)*24];
% %     load(['layers_fr\FR_' data(3:end) '_s2.mat']);
% %     neurons_layer_d=[neurons_layer_d (ld)+(i-1)*24+12];
% %     neurons_layer_4=[neurons_layer_4 (l4)+(i-1)*24+12];
% %     neurons_layer_s=[neurons_layer_s (ls)+(i-1)*24+12];
% %     neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
% % end
% 
% for i=6:12
%     data=datas{i};
%     load(['layers_fr\FR_' data(3:end) '_s1.mat']);
%     neurons_layer_s=[neurons_layer_s (ld)+(i-1)*24];
%     neurons_layer_4=[neurons_layer_4 (l4)+(i-1)*24];
%     neurons_layer_d=[neurons_layer_d (ls)+(i-1)*24];
%     load(['layers_fr\FR_' data(3:end) '_s2.mat']);
%     neurons_layer_s=[neurons_layer_s (ld)+(i-1)*24+12];
%     neurons_layer_4=[neurons_layer_4 (l4)+(i-1)*24+12];
%     neurons_layer_d=[neurons_layer_d (ls)+(i-1)*24+12];
%     neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
% end
% 
% monkey='FR';
% stim_type='ec';
% setting='v2_all';
% neuron_group='V2';
% t_window='40-110';
% %t_window='110-180';
% load(['tuning_all\' t_window '\fr_' stim_type '_all_tuning.mat'])
% all_tuning=all_tuning(setdiff(neurons_layer_all,all_bad_channels),:);

%layers={'conv1','pool1','conv2','pool2','conv3','conv4','conv5','pool5'};
%layers={'conv1_2','pool1','conv2_1','conv2_2','pool2','conv3_1','conv3_4','pool3','conv4_1','conv4_4','pool4','conv5_1'};
%layers={'max_pooling2d_1','activation_4_relu','activation_7_relu','activation_10_relu','activation_13_relu','activation_16_relu'};
stim_types={'ec','ac','ex'};
% nn_name='vgg';
% layers={'pool1','pool2','pool3','pool4','pool5'};
%nn_name='Alexnet';
%layers={'pool1','pool2','pool5'};

nn_name='gabor_pyramid_resps';
%layers={'odd_sq','even_sq','odd','even','complex'};
layers={'simple','complex'};

for layer_index=1:length(layers)
    nn_setting=layers{layer_index};
    all_bad_channels=[];
    n_bad_neurons=[];
    for stim_type_index=1:length(stim_types)
        stim_type=stim_types{stim_type_index};
        %load(['D:\15386\corr_mat\' nn_name '\no_jitter\3_scale_16_orientation\' nn_setting '\' stim_type '.mat'])

        load(['D:\15386\corr_mat\' nn_name '\3x3_jitter_avg_area\3_scale_16_orientation\' nn_setting '\' stim_type '.mat'])
        all_tuning=center_resps';
        
        %all_tuning(all_bad_channels,:)=[];
        all_tuning=all_tuning';
        %all_tuning(all_bad_channels,:)=[];
        combined_tuning=[];
        combined_tuning=[combined_tuning,all_tuning];
        
        
        [sorted_tuning,tuning_rankorder]=sort(all_tuning,'descend');
        min_tuning=min(sorted_tuning);
        max_tuning=max(sorted_tuning);
        normalized_tuning=(sorted_tuning-min_tuning)./(max_tuning-min_tuning);
        diff=(normalized_tuning(1:end-1,:)-normalized_tuning(2:end,:))./2;
        sparsity=1-(sum(normalized_tuning)-sum(diff))./50;
        
        n_bad_neuron=sum(isnan(sparsity));
        n_bad_neurons=[n_bad_neurons,n_bad_neuron];
        all_bad_channels=[all_bad_channels;isnan(sparsity)];
    end
    all_bad_channels=sum(all_bad_channels);
    all_bad_channels=find(all_bad_channels==3);
    save([ nn_name '\3x3_jitter_avg_area\3_scale_16_orientation\' nn_setting '\' 'bad_channels.mat'],'all_bad_channels')
end