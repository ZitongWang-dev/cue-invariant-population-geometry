% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % corr between type
%
% datas_ko_v1={
%     %v1
% { {{'KO041015' '2'}} {{'KO041015' '3'}} {{'KO041015' '4'}} }
% { {{'KO041515' '2'}} {{'KO041515' '4'}} {{'KO041515' '3'}} }
% { {{'KO042115' '4'}} {{'KO042115' '2'}} {{'KO042115' '3'}} }
% { {{'KO043015' '4'}} {{'KO043015' '2'}} {{'KO043015' '3'}} }
% { {{'KO060315' '4'}} {{'KO060315' '3'}} {{'KO060315' '2'}} }
%
% };
% datas_ko_v2={
%     %v2
% %{ {{'KO031915' '2'}} {{'KO031915' '3'}} {{'KO031915' '4'}} }
% { {{'KO041715' '4'}} {{'KO041715' '3'}} {{'KO041715' '2'}} }
% { {{'KO042515' '4'}} {{'KO042515' '2'}} {{'KO042515' '3'}} }
% { {{'KO042815' '4'}} {{'KO042815' '3'}} {{'KO042815' '2'}} }
% { {{'KO050115' '4'}} {{'KO050115' '3'}} {{'KO050115' '2'}} }
% { {{'KO060815' '4'}} {{'KO060815' '3'}} {{'KO060815' '2'}} }
% { {{'KO060915' '4'}} {{'KO060915' '2'}} {{'KO060915' '3'}} }
%
% };
%
% datas_fr_v1={
%     %v1
% %{ {{'FR041014' '2'}} {{'FR041014' '4'}} {{'FR041014' '3'}}}
% %{ {{'FR041014' '2'}} {{'FR041014' '4'}} {{'FR041014' '3'}}}
% { {{'FR081414' '4'}} {{'FR081414' '2'}} {{'FR081414' '3'}}}
% { {{'FR082114' '4'}} {{'FR082114' '3'}} {{'FR082114' '2'}}}
% { {{'FR082814' '4'}} {{'FR082814' '3'}} {{'FR082814' '2'}}}
% { {{'FR082914' '4'}} {{'FR082914' '2'}} {{'FR082914' '3'}}}
% };
%
% datas_fr_v2={
%     %v2
% { {{'FR041814' '3'}} {{'FR041814' '4'}} {{'FR041814' '2'}}}
% { {{'FR060414' '4'}} {{'FR060414' '2'}} {{'FR060414' '3'}}}
% { {{'FR060514' '3'}} {{'FR060514' '4'}} {{'FR060514' '2'}}}
% { {{'FR072614' '2'}} {{'FR072614' '3'}} {{'FR072614' '4'}}}
% { {{'FR073014' '3'}} {{'FR073014' '4'}} {{'FR073014' '2'}}}
% %{ {{'FR090414' '4'}} {{'FR090414' '3'}} {{'FR090414' '2'}}}
% %{ {{'FR090414' '4'}} {{'FR090414' '3'}} {{'FR090414' '2'}}}
% { {{'FR090914' '2'}} {{'FR090914' '3'}} {{'FR090914' '4'}}}
% };
%
% n_stim=50;
% n_trial=10;
%
% monkey=upper('ko');
% brain_region=upper('v1');
% type_dict={'EC','AC','EX'};
% datas=datas_fr_v1;
% psth_path='psths\';
% session_trials_ec=[];
% session_trials_ac=[];
% session_trials_ex=[];
% n_good_channels=0;
% for i=1:length(datas)
%     load([psth_path datas{i}{1}{1}{1} '_' datas{i}{1}{1}{2} '.mat'])
%     n_good_channel=length(setdiff(1:24,bad_channel));
%     selected_channels=psths(1:50,:,setdiff(1:24,bad_channel));
%     for j=1:n_good_channel
%         session_trials_ec=[session_trials_ec, { selected_channels(:,:,j)}];
%     end
%
%     load([psth_path datas{i}{2}{1}{1} '_' datas{i}{2}{1}{2} '.mat'])
%     n_good_channel=length(setdiff(1:24,bad_channel));
%     selected_channels=psths(1:50,:,setdiff(1:24,bad_channel));
%     for j=1:n_good_channel
%         session_trials_ac=[session_trials_ac, { selected_channels(:,:,j)}];
%     end
%
%     load([psth_path datas{i}{3}{1}{1} '_' datas{i}{3}{1}{2} '.mat'])
%     n_good_channel=length(setdiff(1:24,bad_channel));
%     selected_channels=psths(1:50,:,setdiff(1:24,bad_channel));
%     for j=1:n_good_channel
%         session_trials_ex=[session_trials_ex, { selected_channels(:,:,j)}];
%     end
%
%     n_good_channels=n_good_channels+n_good_channel;
% end
%
% session_trials={session_trials_ec,session_trials_ac,session_trials_ex};
% tunings={};
%
% for i=1:3
%     curr_session=session_trials{i};
%     curr_tuning=zeros(n_good_channels,n_stim,n_trial);
%     for j=1:n_good_channels
%         curr_cell=curr_session{j};
%         for k=1:n_trial
%             for l=1:n_stim
%                 curr_tuning(j,l,k)=mean(curr_cell{l,k}(340:550));
%             end
%         end
%     end
%     tunings{i}=curr_tuning;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datas_ko=[
    %v1
    {'KO041015'}
    {'KO041515'}
    {'KO042115'}
    {'KO043015'}
    {'KO060315'}
    %v2
    %{'KO031915'}
    {'KO041715'}
    {'KO042515'}
    {'KO042815'}
    {'KO050115'}
    {'KO060815'}
    {'KO060915'}
    ];

datas_fr=[
    %v1
    %{'FR041014'}
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
    %{'FR090414'}
    {'FR090914'}
    ];

monkey='ko';
layer_name='test';
if strcmp(monkey,'ko')
    datas=datas_ko;
else
    datas=datas_fr;
end

neuron='v1';

% load stimulus tuning data
neurons_layer_d=[];
neurons_layer_4=[];
neurons_layer_s=[];
neurons_layer_all=[];
if strcmp(neuron,'v1')
    %v1
    if strcmp(monkey,'ko')
        i_max=5;
    else
        i_max=4;
    end
    for i=1:i_max
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
    for i=i_max+1:i_max+6
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
chosen_layer=neurons_layer_s;

load(['tuning_all\40-250\' monkey '_ec_all_tuning.mat'])
layer_all=setdiff(neurons_layer_all,all_bad_channels);
chosen_layer=setdiff(chosen_layer,all_bad_channels);
chosen_layer_idx_v1=[];
for it=chosen_layer
    chosen_layer_idx_v1=[chosen_layer_idx_v1,find(chosen_layer==it)];
end

neuron='v2';
% load stimulus tuning data
neurons_layer_d=[];
neurons_layer_4=[];
neurons_layer_s=[];
neurons_layer_all=[];
if strcmp(neuron,'v1')
    %v1
    if strcmp(monkey,'ko')
        i_max=5;
    else
        i_max=4;
    end
    for i=1:i_max
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
    for i=i_max+1:i_max+6
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
chosen_layer=neurons_layer_s;

load(['tuning_all\40-250\' monkey '_ec_all_tuning.mat'])
layer_all=setdiff(neurons_layer_all,all_bad_channels);
chosen_layer=setdiff(chosen_layer,all_bad_channels);
chosen_layer_idx_v2=[];
for it=chosen_layer
    chosen_layer_idx_v2=[chosen_layer_idx_v2,find(chosen_layer==it)];
end

n_neuron=min(length(chosen_layer_idx_v1),length(chosen_layer_idx_v2));

means=[];
stds=[];


for neuron_setting={'v1','v2'}
    neuron=neuron_setting{1};
    load(['tuning_all\' monkey '_' neuron '_bytrial_40-250.mat'])
    if strcmp(neuron,'v1')
        chosen_layer_idx=chosen_layer_idx_v1;
    else
        chosen_layer_idx=chosen_layer_idx_v2;       
    end
    
    half_rsam_corrs=[];
    upper_bounds=[];
    
    n_trial=10;
    parfor dup=1:100
        
        curr_randperm=randperm(n_trial);
        i_trials_1=curr_randperm(1:5);
        i_trials_2=curr_randperm(6:10);
        
        n_neuron=size(tunings{1});
        n_neuron=n_neuron(1);
        neuron_sampling=randperm(n_neuron);
        neuron_sampling=neuron_sampling(1:n_neuron);
        
        tuning_ec_1=mean(tunings{1}(neuron_sampling,:,i_trials_1),3);
        tuning_ec_2=mean(tunings{1}(neuron_sampling,:,i_trials_2),3);
        tuning_ec_1=tuning_ec_1(chosen_layer_idx,:);
        tuning_ec_2=tuning_ec_2(chosen_layer_idx,:);
        
        tuning_ac_1=mean(tunings{2}(neuron_sampling,:,i_trials_1),3);
        tuning_ac_2=mean(tunings{2}(neuron_sampling,:,i_trials_2),3);
        tuning_ac_1=tuning_ac_1(chosen_layer_idx,:);
        tuning_ac_2=tuning_ac_2(chosen_layer_idx,:);
        
        tuning_ex_1=mean(tunings{3}(neuron_sampling,:,i_trials_1),3);
        tuning_ex_2=mean(tunings{3}(neuron_sampling,:,i_trials_2),3);
        tuning_ex_1=tuning_ex_1(chosen_layer_idx,:);
        tuning_ex_2=tuning_ex_2(chosen_layer_idx,:);
        
        m_ec_1=(tuning_ec_1'-mean(tuning_ec_1'))';
        m_ac_1=(tuning_ac_1'-mean(tuning_ac_1'))';
        m_ex_1=(tuning_ex_1'-mean(tuning_ex_1'))';
        
        m_ec_2=(tuning_ec_2'-mean(tuning_ec_2'))';
        m_ac_2=(tuning_ac_2'-mean(tuning_ac_2'))';
        m_ex_2=(tuning_ex_2'-mean(tuning_ex_2'))';
        
        rsam_ec_1=corr(m_ec_1,m_ec_1,'Type','Pearson');
        rsam_ac_1=corr(m_ac_1,m_ac_1,'Type','Pearson');
        rsam_ex_1=corr(m_ex_1,m_ex_1,'Type','Pearson');
        
        rsam_ec_2=corr(m_ec_2,m_ec_2,'Type','Pearson');
        rsam_ac_2=corr(m_ac_2,m_ac_2,'Type','Pearson');
        rsam_ex_2=corr(m_ex_2,m_ex_2,'Type','Pearson');
        
        half_ec_1=rsam_ec_1(tril(true(size(rsam_ec_1))));
        half_ac_1=rsam_ac_1(tril(true(size(rsam_ac_1))));
        half_ex_1=rsam_ex_1(tril(true(size(rsam_ex_1))));
        
        half_ec_2=rsam_ec_2(tril(true(size(rsam_ec_2))));
        half_ac_2=rsam_ac_2(tril(true(size(rsam_ac_2))));
        half_ex_2=rsam_ex_2(tril(true(size(rsam_ex_2))));
        
        half_rsam_corrs=[half_rsam_corrs;[corr(half_ec_1,half_ac_1),corr(half_ec_1,half_ex_1),corr(half_ac_1,half_ex_1)]];
        upper_bounds=[upper_bounds;[corr(half_ec_1,half_ec_2),corr(half_ac_1,half_ac_2),corr(half_ex_1,half_ex_2)]];
        
    end
    
    mean_corrs=mean(half_rsam_corrs);
    mean_ub=mean(upper_bounds);
    mean_corrs=mean_corrs./[mean(mean_ub([1,2])),mean(mean_ub([1,3])),mean(mean_ub([2,3]))];
    std_corrs=std(half_rsam_corrs);
    
    means=[means;mean_corrs];
    stds=[stds;std_corrs];
    [h1,p1]=ttest2(half_rsam_corrs(:,1),half_rsam_corrs(:,2));
    [h2,p2]=ttest2(half_rsam_corrs(:,1),half_rsam_corrs(:,3));
    [h3,p3]=ttest2(half_rsam_corrs(:,2),half_rsam_corrs(:,3));
    neuron
    [h1,h2,h3]
    
    if strcmp(neuron,'v1')
        dist1=half_rsam_corrs;
    end
    
    if strcmp(neuron,'v2')
        dist2=half_rsam_corrs;
    end
end

'v1 vs v2'
[h1,p1]=ttest2(dist1(:,1),dist2(:,1));
[h2,p2]=ttest2(dist1(:,2),dist2(:,2));
[h3,p3]=ttest2(dist1(:,3),dist2(:,3));
[h1,h2,h3]

means=means';
stds=stds';

hb=bar(means);
hold on
for pos = 1:size(means,2)
    xpos=hb(pos).XData+hb(pos).XOffset;
    errorbar(xpos,means(:,pos),stds(:,pos),'LineStyle','none','Color','black')
end
ylim([0,1])
xlim([0,4])
xticks([1,2,3])
xticklabels({'ec-ac','ec-ex','ac-ex'})
title([monkey ' v1&v2 ' layer_name ', #neuron=' num2str(n_neuron)])
legend({'v1','v2'})


set(gcf,'position',[0,0,2000,1000])

print(['rsa_by_trial\ratio_to_ub\neuron_num_controlled\decoding_' monkey '_' layer_name '.png'],'-dpng','-painters','-loose',gcf)