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
neuron='v2';
t_window='40-250';
% group name to write in the csv file
name_to_write='v2_all';

tuning_data_path='tuning_all\';


type_pairs={{'ec','ac'},{'ec','ex'},{'ac','ex'}};%,{'ec','ec'},{'ac','ac'},{'ex','ex'}};
sigs=[];
corr_between_rsas_all=[];
corr_between_rsas_sig_only=[];
chance_means=[];
chance_stderrs=[];
chance_means_sig_only=[];
chance_stderrs_sig_only=[];
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
    
    [tuning_corr,p_vals]=corr(tuning_1',tuning_2','Type','Pearson');
    ident_diag=diag(tuning_corr,0);
    
    % significant correlations
    curr_sigs=diag(p_vals,0)<=0.05;
    
    % negative significant correlations
    curr_sigs_neg=(diag(p_vals,0)<=0.05)&(ident_diag<0);
    
    % positive significant correlations
    curr_sigs_pos=(diag(p_vals,0)<=0.05)&(ident_diag>0);
    
    subplot(1,2,1)
    m1=(tuning_1'-mean(tuning_1'))';
    m2=(tuning_2'-mean(tuning_2'))';
    rsam=corr(m1,m2,'Type','Pearson');
%     imagesc(rsam);
%     xlabel(['Stim(' stim_type{1} ')'])
%     ylabel(['Stim(' stim_type{2} ')'])
%     title('All Neurons')
    
    rsam1=corr(m1,m1,'Type','Pearson');
    half1=rsam1(tril(true(size(rsam1))));
    rsam2=corr(m2,m2,'Type','Pearson');
    half2=rsam2(tril(true(size(rsam2))));

    shuffled_corr=[];
    for shuffle_idx=1:100
        shuffled_idx_1=randperm(length(half1));
        shuffled_idx_2=randperm(length(half2));
        half1_shuffled=half1(shuffled_idx_1);
        half2_shuffled=half2(shuffled_idx_2);
        shuffled_corr=[shuffled_corr,corr(half1_shuffled,half2_shuffled)];
    end

    
    corr(half1,half2)
    corr_between_rsas_all=[corr_between_rsas_all,corr(half1,half2)];
    
    chance_means=[chance_means, mean(shuffled_corr)];
    chance_stderrs=[chance_stderrs, std(shuffled_corr)/sqrt(length(shuffled_corr))];
    
    subplot(1,2,2)
    m1=(tuning_1(curr_sigs_pos,:)'-mean(tuning_1(curr_sigs_pos,:)'))';
    m2=(tuning_2(curr_sigs_pos,:)'-mean(tuning_2(curr_sigs_pos,:)'))';
    rsam=corr(m1,m2,'Type','Pearson');
    %imagesc(corr(tuning_1(curr_sigs_pos,:),tuning_2(curr_sigs_pos,:),'Type','Pearson'));
    
    rsam1=corr(m1,m1,'Type','Pearson');
    half1=rsam1(tril(true(size(rsam1))));
    rsam2=corr(m2,m2,'Type','Pearson');
    half2=rsam2(tril(true(size(rsam2))));
    
    shuffled_corr=[];
    for shuffle_idx=1:100
        shuffled_idx_1=randperm(length(half1));
        shuffled_idx_2=randperm(length(half2));
        half1_shuffled=half1(shuffled_idx_1);
        half2_shuffled=half2(shuffled_idx_2);
        shuffled_corr=[shuffled_corr,corr(half1_shuffled,half2_shuffled)];
    end
    
    corr(half1,half2)
    corr_between_rsas_sig_only=[corr_between_rsas_sig_only,corr(half1,half2)];
    
    
    chance_means_sig_only=[chance_means_sig_only, mean(shuffled_corr)];
    chance_stderrs_sig_only=[chance_stderrs_sig_only, std(shuffled_corr)/sqrt(length(shuffled_corr))];
    
%     imagesc(rsam);
%     xlabel(['Stim(' stim_type{1} ')'])
%     ylabel(['Stim(' stim_type{2} ')'])
%     title('Positive Significant Neurons')
%     colorbar()
%     caxis([-1,1])
%     set(gcf,'position',[0,0,2000,1000])
    %print(['svm_decoding\rsa\decoding_' monkey '_' stim_type{1} '_' stim_type{2} '.png'],'-dpng','-painters','-loose',gcf)
end

corr_between_rsas_all
corr_between_rsas_sig_only

subplot(1,2,1)
hold on
bar(corr_between_rsas_all)
errorbar([1,2,3],chance_means,chance_stderrs)
ylim([-0.1,1])
title([monkey ' ' neuron ' ec-ac ec-ex ac-ex'])
subplot(1,2,2)
hold on
bar(corr_between_rsas_sig_only)
errorbar([1,2,3],chance_means_sig_only,chance_stderrs_sig_only)
ylim([-0.1,1])
title([monkey ' ' neuron ' ec-ac ec-ex ac-ex pos sig neurons only'])
set(gcf,'position',[0,0,2000,1000])
print(['rsa\decoding_' monkey '_' neuron '_' stim_type{1} '_' stim_type{2} '.png'],'-dpng','-painters','-loose',gcf)