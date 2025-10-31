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
t_window='40-110';
neurons={'v1','v2'};
csv_col_names={'v1_all','v2_all'};

% monkey='vgg';
% neurons={'pool1','pool2','pool3','pool4','pool5'};
% csv_col_names=neurons;

% monkey='Alexnet';
% neurons={'pool1','pool2','pool5'};
% csv_col_names=neurons;

monkey='gabor_pyramid_resps';
neurons={'complex'};
csv_col_names=neurons;

out_file_path='.\';
f=fopen([out_file_path  monkey '_pools.csv'],'w');

out_str=['type,ec-ac,ec-ex,ac-ex,ec-ac-ec-ex,ec-ac-ac-ex,ec-ex-ac-ex,ec-ac-ex,all'];
fprintf(f,'%s\n',out_str);

fraction_values=[];
count_values=[];
n_samples=[];
for neuron_group=1:length(neurons)
    % group name to write in the csv file
    neuron=neurons{neuron_group};
    name_to_write=csv_col_names{neuron_group};
    
    %tuning_data_path='tuning_all\';
    tuning_data_path=[ monkey '\no_jitter\3_scale_16_orientation\'];
    
    
    type_pairs={{'ec','ac'},{'ec','ex'},{'ac','ex'}};
    sigs=[];
    
    for cnt=0:2
        stim_type=type_pairs{cnt+1};
        
%         % load stimulus tuning data
%         neurons_layer_d=[];
%         neurons_layer_4=[];
%         neurons_layer_s=[];
%         neurons_layer_all=[];
%         if strcmp(neuron,'v1')
%             %v1
%             for i=1:5
%                 data=datas{i};
%                 load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s1.mat']);
%                 if strcmp(monkey,'fr')
%                     layer_4=l4;
%                 end
%                 neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24];
%                 neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
%                 neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24];
%                 load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s2.mat']);
%                 if strcmp(monkey,'fr')
%                     layer_4=l4;
%                 end
%                 neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24+12];
%                 neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
%                 neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24+12];
%                 neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
%             end
%         else
%             %v2
%             for i=6:12
%                 data=datas{i};
%                 load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s1.mat']);
%                 if strcmp(monkey,'fr')
%                     layer_4=l4;
%                 end
%                 neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24];
%                 neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
%                 neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24];
%                 load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s2.mat']);
%                 if strcmp(monkey,'fr')
%                     layer_4=l4;
%                 end
%                 neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24+12];
%                 neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
%                 neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24+12];
%                 neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
%             end
%         end
%         
%         % select the neuron group
%         % s: superficial
%         % 4: layer_4
%         % d: deep
%         % all: all layers combined
%         chosen_layer=neurons_layer_all;
%         
%         
%         % read tuning data
%         load([tuning_data_path t_window '\' monkey '_' stim_type{1} '_all_tuning.mat'])
%         tuning_1=all_tuning;
%         load([tuning_data_path t_window '\' monkey '_' stim_type{2} '_all_tuning.mat'])
%         tuning_2=all_tuning;

        load([tuning_data_path neuron '\' stim_type{1} '.mat'])
        tuning_1=center_resps;
        load([tuning_data_path neuron '\' stim_type{2} '.mat'])
        tuning_2=center_resps;
        load([tuning_data_path neuron '\bad_channels.mat'])
        chosen_layer=1:length(center_resps);
        
        % remove bad channels
        length(setdiff(chosen_layer,all_bad_channels))
        tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
        tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
        selected_channels=setdiff(chosen_layer,all_bad_channels);
        
        shuffled_tuning_1=tuning_1';
        [m,n]=size(shuffled_tuning_1);
        [~,idx]=sort(rand(m,n));
        shuffled_tuning_1=shuffled_tuning_1(sub2ind([m, n], idx, ones(m,1)*(1:n)));

        shuffled_tuning_2=tuning_2';
        [m,n]=size(shuffled_tuning_2);
        [~,idx]=sort(rand(m,n));
        shuffled_tuning_2=shuffled_tuning_2(sub2ind([m, n], idx, ones(m,1)*(1:n)));
        
        tuning_1=shuffled_tuning_1';
        tuning_2=shuffled_tuning_2';

        
        [tuning_corr,p_vals]=corr(tuning_1',tuning_2','Type','Pearson');
        
        % plot histograms
        subplot(3,3,1+cnt)
        ident_diag=diag(tuning_corr,0);
        x_edge=-1:0.05:1;
        histogram(ident_diag,x_edge)
        xlim([-1,1])
        ylim([0,30])
        xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag))))])
        ylabel('Counts')
        title([monkey ' ' neuron ' neuron response correlation ' stim_type{1} ' vs ' stim_type{2}])
        
        % significant correlations
        curr_sigs=diag(p_vals,0)<=0.05;
        
        % negative significant correlations
        curr_sigs_neg=(diag(p_vals,0)<=0.05)&(ident_diag<0);
        
        subplot(3,3,4+cnt)
        histogram(ident_diag(curr_sigs_neg),x_edge)
        xlim([-1,1])
        ylim([0,30])
        xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)&curr_sigs_neg)))])
        ylabel('Counts')
        title(['neg sig ' stim_type{1} ' vs ' stim_type{2}])
        
        % positive significant correlations
        curr_sigs_pos=(diag(p_vals,0)<=0.05)&(ident_diag>0);
        
        subplot(3,3,7+cnt)
        histogram(ident_diag(curr_sigs_pos),x_edge)
        xlim([-1,1])
        ylim([0,30])
        xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)&curr_sigs_pos)))])
        ylabel('Counts')
        title(['pos sig ' stim_type{1} ' vs ' stim_type{2}])
        
        sigs_pos_channels=selected_channels(curr_sigs_pos);
        %save(['svm_decoding\significant_pairs\' upper(monkey) '_' upper(neuron) '_' upper(stim_type{1}) 'pred' upper(stim_type{2}) '.mat'],'sigs_pos_channels')
        %save(['svm_decoding\significant_pairs\' upper(monkey) '_' upper(neuron) '_' upper(stim_type{2}) 'pred' upper(stim_type{1}) '.mat'],'sigs_pos_channels')
        sigs=[sigs, curr_sigs_pos];
    end
    set(gcf,'position',[0,0,1000,500])
    cnts=[sum(sigs), sum(sigs(:,1)==1 & sigs(:,2)==1), sum(sigs(:,1)==1 & sigs(:,3)==1), sum(sigs(:,2)==1 & sigs(:,3)==1), sum(sigs(:,1)==1 & sigs(:,2)==1 & sigs(:,3)==1)]
    fracs=cnts./length(sigs)
    
    %%%%%%%%%%%%%%%%%%
    % write fractions into csv files
    %%%%%%%%%%%%%%%%%%
    
    out_str=[name_to_write];
    for i=1:length(cnts)
        out_str=[out_str ',' num2str(fracs(i))];
    end
    out_str=[out_str ',' num2str(1)];
    fprintf(f,'%s\n',out_str);
    
    out_str=[name_to_write '_n'];
    for i=1:length(cnts)
        out_str=[out_str ',' num2str(cnts(i))];
    end
    out_str=[out_str ',' num2str(length(sigs))];
    fprintf(f,'%s\n',out_str);
    
    out_str=[name_to_write '_SSE'];
    for i=1:length(cnts)
        out_str=[out_str ',' num2str(sqrt(fracs(i)*(1-fracs(i))/length(sigs)))];
    end
    out_str=[out_str];
    fprintf(f,'%s\n',out_str);

    out_str=[name_to_write '_Wilson_score_interval_lower'];
    sig_level_z_val=1.96;
    for i=1:length(cnts)
        out_str=[out_str ',' num2str((cnts(i)+sig_level_z_val^2/2)/(length(sigs)+sig_level_z_val^2)-sig_level_z_val/(length(sigs)+sig_level_z_val^2)*sqrt(cnts(i)*(length(sigs)-cnts(i))/length(sigs)+sig_level_z_val^2/4))];
    end
    out_str=[out_str];
    fprintf(f,'%s\n',out_str);    

    out_str=[name_to_write '_Wilson_score_interval_higher'];
    for i=1:length(cnts)
        out_str=[out_str ',' num2str((cnts(i)+sig_level_z_val^2/2)/(length(sigs)+sig_level_z_val^2)+sig_level_z_val/(length(sigs)+sig_level_z_val^2)*sqrt(cnts(i)*(length(sigs)-cnts(i))/length(sigs)+sig_level_z_val^2/4))];
    end
    out_str=[out_str];
    fprintf(f,'%s\n',out_str); 
    
    fraction_values=[fraction_values;fracs];
    count_values=[count_values;cnts];
    n_samples=[n_samples;length(sigs)];
end

% for col=1:length(fraction_values)
%     p_est=(count_values(1,col)+count_values(2,col))/(n_samples(1)+n_samples(2));
%     z_values=(fraction_values(1,:)-fraction_values(2,:))./sqrt(p_est*(1-p_est)*(1/n_samples(1)+1/n_samples(2)));
%     p_val_stat_test=min(normcdf(z_values),normcdf(-z_values));
% end
% 
% out_str='z_values (v1 vs v2)';
% for i=1:length(cnts)
%     out_str=[out_str ',' num2str(z_values(i))];
% end
% out_str=[out_str];
% fprintf(f,'%s\n',out_str);
% 
% out_str='p_values (v1 vs v2)';
% for i=1:length(cnts)
%     out_str=[out_str ',' num2str(p_val_stat_test(i))];
% end
% out_str=[out_str];
% fprintf(f,'%s\n',out_str);

fclose(f);