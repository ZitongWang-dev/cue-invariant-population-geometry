%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corr between type

% monkey='Alexnet';
% %neurons={'conv1','pool1','conv2','pool2','conv3','conv4','conv5','pool5'}';
% %neurons_label={'conv1','pool1','conv2','pool2','conv3','conv4','conv5','pool5'}';
% n_neuron=96;
% neurons={'pool1','pool2','pool5'}';
% neurons_label={'pool1','pool2','pool5'}';

monkey='vgg';
% neurons={'conv1_2','pool1','conv2_1','conv2_2','pool2','conv3_1','conv3_4','pool3','conv4_1','conv4_4','pool4','conv5_1'}';
% neurons_label={'conv1\_2','pool1','conv2\_1','conv2\_2','pool2','conv3\_1','conv3\_4','pool3','conv4\_1','conv4\_4','pool4','conv5\_1'}';
n_neuron=64;
neurons={'pool1','pool2','pool3','pool4','pool5'}';
neurons_label={'pool1','pool2','pool3','pool4','pool5'}';

% monkey='resnet-bottlenecks';
% %monkey='resnet';
% neurons={'max_pooling2d_1','activation_4_relu','activation_7_relu','activation_10_relu','activation_13_relu','activation_16_relu'}';
% neurons_label={'max\_pooling2d\_1','activation\_4\_relu','activation\_7\_relu','activation\_10\_relu','activation\_13\_relu','activation\_16\_relu'}';
% n_neuron=63;


curr_stim_set=1:50;
mat_diag_idx=[1];
dec_idx=length(curr_stim_set):-1:2;
for i=1:length(curr_stim_set)-1
    mat_diag_idx=[mat_diag_idx,1+sum(dec_idx(1:i))];
end

corr_between_rsas_by_layer=[];
for layer=1:length(neurons)
    
    neuron=neurons{layer};
    tuning_data_path=[ monkey '\no_jitter\' neuron '\'];
    
    type_pairs={{'ec','ac'},{'ec','ex'},{'ac','ex'}};
    corr_between_rsas_all=[];
    for cnt=0:2
        stim_type=type_pairs{cnt+1};
        
        % read tuning data
        load([tuning_data_path stim_type{1} '.mat'])
        tuning_1=center_resps';
        load([tuning_data_path stim_type{2} '.mat'])
        tuning_2=center_resps';
        load([tuning_data_path 'bad_channels.mat'])
        
        chosen_layer=1:length(center_resps);
        
        % remove bad channels
        length(setdiff(chosen_layer,all_bad_channels))
        tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
        tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
        selected_channels=setdiff(chosen_layer,all_bad_channels);

        size_tuning=size(tuning_1);
        neuron_sampling=randperm(size_tuning(1));
        tuning_1=tuning_1(neuron_sampling(1:n_neuron),:);
        tuning_2=tuning_2(neuron_sampling(1:n_neuron),:);

        
        [tuning_corr,p_vals]=corr(tuning_1',tuning_2','Type','Pearson');
        ident_diag=diag(tuning_corr,0);       
        
        m1=(tuning_1'-mean(tuning_1'))';
        m2=(tuning_2'-mean(tuning_2'))';
        rsam=corr(m1,m2,'Type','Pearson');
        
        rsam1=corr(m1,m1,'Type','Pearson');
        half1=rsam1(tril(true(size(rsam1))));
        half1=half1(setdiff(1:length(half1),mat_diag_idx));
        rsam2=corr(m2,m2,'Type','Pearson');
        half2=rsam2(tril(true(size(rsam2))));
        half2=half2(setdiff(1:length(half2),mat_diag_idx));
        
        corr(half1,half2);
        corr_between_rsas_all=[corr_between_rsas_all,corr(half1,half2)];
    end
    
    corr_between_rsas_all;
    corr_between_rsas_by_layer=[corr_between_rsas_by_layer;corr_between_rsas_all];
    
end

%plot(1:length(neurons),corr_between_rsas_by_layer)
bar(corr_between_rsas_by_layer')
%legend({'ec-ac','ec-ex','ac-ex'});
legend(neurons_label);
%xlabel('Layer')
ylabel('Correlation')
ylim([0,1])
% xlim([0,length(neurons)+1])
% xticks(1:length(neurons))
xlim([0,4])
xticks(1:3)
%xticklabels(neurons_label)
xticklabels({'ec-ac','ec-ex','ac-ex'})
title([monkey ' corrrelation between RSAs'])
set(gcf,'position',[0,0,2000,1000])
print(['rsa\cnns\neuron_num_controlled\20230908\rsacorr_bar_' monkey '.png'],'-dpng','-painters','-loose',gcf)
print(['rsa\cnns\neuron_num_controlled\20230908\rsacorr_bar_' monkey '.eps'],'-depsc','-painters','-loose',gcf)

% subplot(1,2,1)
% hold on
% bar(corr_between_rsas_all)
% errorbar([1,2,3],chance_means,chance_stderrs)
% ylim([-0.1,1])
% title([monkey ' ' neuron ' ec-ac ec-ex ac-ex'])
% subplot(1,2,2)
% hold on
% bar(corr_between_rsas_sig_only)
% errorbar([1,2,3],chance_means_sig_only,chance_stderrs_sig_only)
% ylim([-0.1,1])
% title([monkey ' ' neuron ' ec-ac ec-ex ac-ex pos sig neurons only'])
% set(gcf,'position',[0,0,2000,1000])
% print(['rsa\decoding_' monkey '_' neuron '_' stim_type{1} '_' stim_type{2} '.png'],'-dpng','-painters','-loose',gcf)