net=alexnet;%vgg19; %resnet50 %alexnet
%layers={'conv1','pool1','conv2','pool2','conv3','conv4','conv5','pool5'};
layers={'pool1','pool2','pool5'};
%layers={'pool1','pool2','pool3','pool4','pool5'};
%stim_sizes=[11,19,51,67,99,131,163,195];
stim_sizes=[19,67,195];
%stim_sizes=[6,16,52,124,224];
for l=1:length(layers)
    curr_layer=layers{l};
    stim_dir='stimuli\';
    % new_size=51;
    new_size=stim_sizes(l);
    stim_dir=['resized_stim/227/' num2str(new_size) '/'];
    
    nstimuli=51;
    stimlist=cell(1,nstimuli-1);
    stim_cnt=0;
    
    sizes = [7, 11, 14, 21, 35, 42, 49, 56, 60];
    %sizes = [35, 42, 49, 56];
    orientations = 0:10:170;
    gabors=cell(1,length(orientations)*length(sizes));
    
    stimtype='ac';
    if strcmp(stimtype,'ec')
        for idx=1:nstimuli-1
            stim_cnt=stim_cnt+1;
            stimlist{stim_cnt}=[stimtype '_' num2str(idx,'%03d') '_l.png' ];
        end
    else
        for idx=1:nstimuli-1
            stim_cnt=stim_cnt+1;
            stimlist{stim_cnt}=[stimtype '_' num2str(idx,'%03d') '_sl.png' ];
        end
    end
    
    %mean_resps=[];
    center_resps=[];
    % random position should be here, not below for no jittering
    pos1=randi(2)-1;
    pos2=randi(2)-1;
    for stim_idx=1:nstimuli-1
        curr_stim=stimlist{stim_idx};
        im=imread([stim_dir curr_stim]);
        %     if strcmp(stimtype,'ec')
        %         extended_im=ones(292).*128;
        %         extended_im(81:292-80,81:292-80)=im;
        %         im=extended_im;
        %     end
        im_rgb=repmat(im,[1 1 3]);
        resps=activations(net,im_rgb,curr_layer,'OutputAs','channels');
        
        % the random jittering is questionable: should jitter on the image
        % not the neuron
        %pos1=randi(2)-1;
        %pos2=randi(2)-1;
        resp_size=size(resps);
        if mod(resp_size(1),2)==0
            %even input
            curr_center_resps=reshape(resps(resp_size(1)/2+pos1,resp_size(1)/2+pos2,:),1,resp_size(3));
        else
            %odd input
            curr_center_resps=reshape(resps((resp_size(1)+1)/2,(resp_size(1)+1)/2,:),1,resp_size(3));
        end
        center_resps=[center_resps;curr_center_resps];
        
        %     curr_mean_resps=reshape(mean(mean(resps)),1,resp_size(3));
        %     mean_resps=[mean_resps;curr_mean_resps];
    end
    save(['corr_mat\Alexnet\no_jitter\' curr_layer '\' stimtype '.mat'],'center_resps')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % corr between type
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
% % ];
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
% ];
% 
% sigs=[]
% for cnt=0:2
%     type_pairs={{'ec','ac'},{'ec','ex'},{'ac','ex'}};
%     stim_type=type_pairs{cnt+1};
%     
%     monkey='fr';
%     neuron='v1';
%     t_window='110-250';
%     
%     neurons_layer_d=[];
%     neurons_layer_4=[];
%     neurons_layer_s=[];
%     neurons_layer_all=[];
%     if strcmp(neuron,'v1')
%         %v1
%         for i=1:5
%             data=datas{i};
%             %load(['layers_ko\KO_' data(3:end) '_s1.mat']);
%             load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s1.mat']);
%             if strcmp(monkey,'fr')
%                 layer_4=l4;
%             end
%             neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24];
%             neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
%             neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24];
%             %load(['layers_ko\KO_' data(3:end) '_s2.mat']);
%             load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s2.mat']);
%             %layer_4=l4;
%             if strcmp(monkey,'fr')
%                 layer_4=l4;
%             end
%             neurons_layer_d=[neurons_layer_d (1:min(layer_4)-1)+(i-1)*24+12];
%             neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
%             neurons_layer_s=[neurons_layer_s (max(layer_4)+1:12)+(i-1)*24+12];
%             neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
%         end
% %         tuning_1=tuning_1(setdiff(neurons_layer_d,all_bad_channels),:);
% %         tuning_2=tuning_2(setdiff(neurons_layer_d,all_bad_channels),:);
%     else
%         %v2
%         for i=6:12
%             data=datas{i};
%             %load(['layers_ko\KO_' data(3:end) '_s1.mat']);
%             load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s1.mat']);
%             %layer_4=l4;
%             if strcmp(monkey,'fr')
%                 layer_4=l4;
%             end
%             neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24];
%             neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24];
%             neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24];
%             %load(['layers_ko\KO_' data(3:end) '_s2.mat']);
%             load(['layers_' monkey '\' upper(monkey) '_' data(3:end) '_s2.mat']);
%             %layer_4=l4;
%             if strcmp(monkey,'fr')
%                 layer_4=l4;
%             end
%             neurons_layer_s=[neurons_layer_s (1:min(layer_4)-1)+(i-1)*24+12];
%             neurons_layer_4=[neurons_layer_4 (layer_4)+(i-1)*24+12];
%             neurons_layer_d=[neurons_layer_d (max(layer_4)+1:12)+(i-1)*24+12];
%             neurons_layer_all=[neurons_layer_all (i-1)*24+1:i*24];
%         end
% %         tuning_1=tuning_1(setdiff(neurons_layer_d,all_bad_channels),:);
% %         tuning_2=tuning_2(setdiff(neurons_layer_d,all_bad_channels),:);
%     end
%         
%     chosen_layer=neurons_layer_all;
%     %name_to_write='v1_late';
%     name_to_write='512';
%     
%     %gabor_setting='even_sq';
% %     tuning_1=0;
% %     load(['D:\15386\corr_mat\tuning_all\40-110\' monkey '_' stim_type{1} '_all_tuning.mat'])
% %     tuning_1=tuning_1+all_tuning/2;
% %     load(['D:\15386\corr_mat\tuning_all\110-180\' monkey '_' stim_type{1} '_all_tuning.mat'])
% %     tuning_1=tuning_1+all_tuning/2;
% %     tuning_2=0;
% %     load(['D:\15386\corr_mat\tuning_all\40-110\' monkey '_' stim_type{2} '_all_tuning.mat'])
% %     tuning_2=tuning_2+all_tuning/2;
% %     load(['D:\15386\corr_mat\tuning_all\110-180\' monkey '_' stim_type{2} '_all_tuning.mat'])
% %     tuning_2=tuning_2+all_tuning/2;
% 
%     load(['D:\15386\corr_mat\tuning_all\' t_window '\' monkey '_' stim_type{1} '_all_tuning.mat'])
%     tuning_1=all_tuning;
%     load(['D:\15386\corr_mat\tuning_all\' t_window '\' monkey '_' stim_type{2} '_all_tuning.mat'])
%     tuning_2=all_tuning;
%     
% %     load(['D:\15386\corr_mat\gabor_pyramid_resps\odd\' stim_type{1} '.mat'])
% %     tuning_1=double(center_resps)';
% %     load(['D:\15386\corr_mat\gabor_pyramid_resps\odd\' stim_type{2} '.mat'])
% %     tuning_2=double(center_resps)';
% 
% %     if strcmp(neuron,'v1')
% %         %v1
% %         tuning_1=tuning_1(setdiff(1:120,all_bad_channels),:);
% %         tuning_2=tuning_2(setdiff(1:120,all_bad_channels),:);
% %     else
% %         %v2
% %         tuning_1=tuning_1(setdiff(121:288,all_bad_channels),:);
% %         tuning_2=tuning_2(setdiff(121:288,all_bad_channels),:);
% %     end
% 
% %     load(['D:\15386\corr_mat\tuning_all\' t_window '\' monkey '_' stim_type{1} '_all_tuning.mat'])
% %     tuning_1=all_tuning;
% %     load(['D:\15386\corr_mat\tuning_all\' t_window '\' monkey '_' stim_type{2} '_all_tuning.mat'])
% %     tuning_2=all_tuning;
% %     load(['D:\15386\corr_mat\tuning_all\' t_window '\' monkey '_' stim_type{3} '_all_tuning.mat'])
% %     tuning_3=all_tuning;
% 
%     length(setdiff(chosen_layer,all_bad_channels))
%     tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
%     tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
%     
%     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for NN
% %     load(['Alexnet\' name_to_write '\' stim_type{1} '.mat']);
% %     tuning_1=center_resps';
% %     load(['Alexnet\' name_to_write '\' stim_type{2} '.mat']);
% %     tuning_2=center_resps';
%     
%     load(['D:\15386\corr_mat\Sparse_coding_filters\responses\' name_to_write '_' stim_type{1} '.mat'])
%     tuning_1=center_resps';
%     load(['D:\15386\corr_mat\Sparse_coding_filters\responses\' name_to_write '_' stim_type{2} '.mat'])
%     tuning_2=center_resps';
%     
%     subplot(3,3,1+cnt)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %shuffling
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % perms_1=[];
%     % for i=1:length(tuning_1)*10
%     %     perms_1=[perms_1;randperm(50)];
%     % end
%     %
%     % perms_2=[];
%     % for i=1:length(tuning_2)*10
%     %     perms_2=[perms_2;randperm(50)];
%     % end
%     %
%     % tuning_1_5t=[tuning_1;tuning_1;tuning_1;tuning_1;tuning_1];
%     % tuning_2_5t=[tuning_2;tuning_2;tuning_2;tuning_2;tuning_2];
%     %
%     % tuning_1_10t=[tuning_1_5t;tuning_1_5t];
%     % tuning_2_10t=[tuning_2_5t;tuning_2_5t];
%     %
%     % [tuning_corr_raw,p_vals_raw]=corr(tuning_1',tuning_2');
%     %
%     % tuning_1=tuning_1_10t(perms_1);
%     % tuning_2=tuning_2_10t(perms_2);
%     % % tuning_1=tuning_1(perms_1);
%     % % tuning_2=tuning_2(perms_2);
%     
%     
%     %[tuning_corr,p_vals]=corr(tuning_1',tuning_2');
%     [tuning_corr,p_vals]=corr(tuning_1',tuning_2','Type','Pearson');
%     
%     % ident_diag_raw=diag(tuning_corr_raw,0);
%     % ident_diag=diag(tuning_corr,0);
%     
%     ident_diag=diag(tuning_corr,0);
%     x_edge=-1:0.05:1;
%     histogram(ident_diag,x_edge)
%     xlim([-1,1])
%     ylim([0,30])
%     %ylim([0,64])
%     xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag))))])
%     ylabel('Counts')
%     title([monkey ' ' neuron ' neuron response correlation ' stim_type{1} ' vs ' stim_type{2}])
%     %title(['gabor 52 ' stim_type{1} ' vs ' stim_type{2}])
%     % % imagesc(tuning_corr)
%     % % colorbar()
%     % % caxis([-1,1])
%     % % title([monkey ' ' stim_type{1} ' vs ' stim_type{2} ' ' t_window])
%     %set(gcf,'position',[0,0,1000,1000])
%     
%     curr_sigs=diag(p_vals,0)<=0.05;
%     
%     curr_sigs_neg=(diag(p_vals,0)<=0.05)&(ident_diag<0);
%     subplot(3,3,4+cnt)
%     histogram(ident_diag(curr_sigs_neg),x_edge)
%     xlim([-1,1])
%     ylim([0,30])
%     xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)&curr_sigs_neg)))])
%     ylabel('Counts')
%     title(['neg sig ' stim_type{1} ' vs ' stim_type{2}])
%     
%     
%     curr_sigs_pos=(diag(p_vals,0)<=0.05)&(ident_diag>0);
%     subplot(3,3,7+cnt)
%     histogram(ident_diag(curr_sigs_pos),x_edge)
%     xlim([-1,1])
%     ylim([0,30])
%     xlabel(['Tuning Correlation' newline 'mean=' num2str(mean(ident_diag(~isnan(ident_diag)&curr_sigs_pos)))])
%     ylabel('Counts')
%     title(['pos sig ' stim_type{1} ' vs ' stim_type{2}])
%     
%     
%     sigs=[sigs, curr_sigs_pos];
% end
% set(gcf,'position',[0,0,1000,500])
% cnts=[sum(sigs), sum(sigs(:,1)==1 & sigs(:,2)==1), sum(sigs(:,1)==1 & sigs(:,3)==1), sum(sigs(:,2)==1 & sigs(:,3)==1), sum(sigs(:,1)==1 & sigs(:,2)==1 & sigs(:,3)==1)]
% fracs=cnts./length(sigs)
% 
% % f=fopen(['csvs\' 'SC' '\' name_to_write '.csv'],'w');
% % 
% % out_str=['type,ec-ac,ec-ex,ac-ex,ec-ac-ec-ex,ec-ac-ac-ex,ec-ex-ac-ex,ec-ac-ex,all'];
% % %fprintf(f,'%s\n',out_str);
% % 
% % out_str=[name_to_write];
% % for i=1:length(cnts)
% %     out_str=[out_str ',' num2str(fracs(i))];
% % end
% % out_str=[out_str ',' num2str(1)];
% % fprintf(f,'%s\n',out_str);
% % 
% % out_str=[name_to_write '_n'];
% % for i=1:length(cnts)
% %     out_str=[out_str ',' num2str(cnts(i))];
% % end
% % out_str=[out_str ',' num2str(length(sigs))];
% % fprintf(f,'%s\n',out_str);
% % 
% % fclose(f);