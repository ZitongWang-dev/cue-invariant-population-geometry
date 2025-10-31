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

monkey='ko';
neuron='v1';
load(['tuning_all\' monkey '_' neuron '_bytrial_40-250.mat'])

half_rsam_corrs=[];
upper_bounds=[];
half_rsam_corrs_pos_only=[];
upper_bounds_pos_only=[];

n_trial=10;
parfor dup=1:100
    
    curr_randperm=randperm(n_trial);
    i_trials_1=curr_randperm(1:5);
    i_trials_2=curr_randperm(6:10);
    
    n_neuron=size(tunings{1});
    n_neuron=n_neuron(1);
    neuron_sampling=randperm(n_neuron);
    neuron_sampling=neuron_sampling(1:100);
    
    tuning_ec_1=mean(tunings{1}(neuron_sampling,:,i_trials_1),3);
    tuning_ec_2=mean(tunings{1}(neuron_sampling,:,i_trials_2),3);

    tuning_ac_1=mean(tunings{2}(neuron_sampling,:,i_trials_1),3);
    tuning_ac_2=mean(tunings{2}(neuron_sampling,:,i_trials_2),3);
    
    tuning_ex_1=mean(tunings{3}(neuron_sampling,:,i_trials_1),3);
    tuning_ex_2=mean(tunings{3}(neuron_sampling,:,i_trials_2),3);
        
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
    
    % pos only
    %%%%%%%%%%%%%%%%%%%
    [tuning_corr_1,p_vals_1]=corr(tuning_ec_1',tuning_ac_1','Type','Pearson');
    ident_diag_1=diag(tuning_corr_1,0);
    [tuning_corr_2,p_vals_2]=corr(tuning_ec_1',tuning_ex_1','Type','Pearson');
    ident_diag_2=diag(tuning_corr_2,0);
    [tuning_corr_3,p_vals_3]=corr(tuning_ac_1',tuning_ex_1','Type','Pearson');
    ident_diag_3=diag(tuning_corr_3,0);
        
    % positive significant correlations
    curr_sigs_pos_1=(diag(p_vals_1,0)<=0.05)&(ident_diag_1>0);
    curr_sigs_pos_2=(diag(p_vals_2,0)<=0.05)&(ident_diag_2>0);
    curr_sigs_pos_3=(diag(p_vals_3,0)<=0.05)&(ident_diag_3>0);
    
    m_ec_1_ps=(tuning_ec_1(curr_sigs_pos_1,:)'-mean(tuning_ec_1(curr_sigs_pos_1,:)'))';
    m_ac_1_ps=(tuning_ac_1(curr_sigs_pos_2,:)'-mean(tuning_ac_1(curr_sigs_pos_2,:)'))';
    m_ex_1_ps=(tuning_ex_1(curr_sigs_pos_3,:)'-mean(tuning_ex_1(curr_sigs_pos_3,:)'))';

    m_ec_2_ps=(tuning_ec_2(curr_sigs_pos_1,:)'-mean(tuning_ec_2(curr_sigs_pos_1,:)'))';
    m_ac_2_ps=(tuning_ac_2(curr_sigs_pos_2,:)'-mean(tuning_ac_2(curr_sigs_pos_2,:)'))';
    m_ex_2_ps=(tuning_ex_2(curr_sigs_pos_3,:)'-mean(tuning_ex_2(curr_sigs_pos_3,:)'))';
    
    rsam_ec_1_ps=corr(m_ec_1_ps,m_ec_1_ps,'Type','Pearson');
    rsam_ac_1_ps=corr(m_ac_1_ps,m_ac_1_ps,'Type','Pearson');
    rsam_ex_1_ps=corr(m_ex_1_ps,m_ex_1_ps,'Type','Pearson');

    rsam_ec_2_ps=corr(m_ec_2_ps,m_ec_2_ps,'Type','Pearson');
    rsam_ac_2_ps=corr(m_ac_2_ps,m_ac_2_ps,'Type','Pearson');
    rsam_ex_2_ps=corr(m_ex_2_ps,m_ex_2_ps,'Type','Pearson');
    
    half_ec_1_ps=rsam_ec_1_ps(tril(true(size(rsam_ec_1_ps))));
    half_ac_1_ps=rsam_ac_1_ps(tril(true(size(rsam_ac_1_ps))));
    half_ex_1_ps=rsam_ex_1_ps(tril(true(size(rsam_ex_1_ps))));
    
    half_ec_2_ps=rsam_ec_2_ps(tril(true(size(rsam_ec_2_ps))));
    half_ac_2_ps=rsam_ac_2_ps(tril(true(size(rsam_ac_2_ps))));
    half_ex_2_ps=rsam_ex_2_ps(tril(true(size(rsam_ex_2_ps))));
    
    half_rsam_corrs_pos_only=[half_rsam_corrs_pos_only;[corr(half_ec_1_ps,half_ac_1_ps),corr(half_ec_1_ps,half_ex_1_ps),corr(half_ac_1_ps,half_ex_1_ps)]];
    upper_bounds_pos_only=[upper_bounds_pos_only;[corr(half_ec_1_ps,half_ec_2_ps),corr(half_ac_1_ps,half_ac_2_ps),corr(half_ex_1_ps,half_ex_2_ps)]];
end

%% plot raw
% 
% mean_corrs=mean(half_rsam_corrs);
% stderr_corrs=std(half_rsam_corrs)/sqrt(length(half_rsam_corrs));
% 
% mean_ub=mean(upper_bounds);
% stderr_ub=std(upper_bounds)/sqrt(length(upper_bounds));
% 
% subplot(1,2,1)
% hold on
% bar(mean_corrs)
% errorbar([1,2,3],mean_corrs,stderr_corrs)
% line([-1,5],[mean_ub(1),mean_ub(1)],'Color','red')
% line([-1,5],[mean_ub(1)+stderr_ub(1),mean_ub(1)+stderr_ub(1)],'Color','red','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(1)-stderr_ub(1),mean_ub(1)-stderr_ub(1)],'Color','red','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(2),mean_ub(2)],'Color','green')
% line([-1,5],[mean_ub(2)+stderr_ub(2),mean_ub(2)+stderr_ub(2)],'Color','green','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(2)-stderr_ub(2),mean_ub(2)-stderr_ub(2)],'Color','green','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(3),mean_ub(3)],'Color','blue')
% line([-1,5],[mean_ub(3)+stderr_ub(3),mean_ub(3)+stderr_ub(3)],'Color','blue','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(3)-stderr_ub(3),mean_ub(3)-stderr_ub(3)],'Color','blue','LineStyle',':','HandleVisibility','off')
% ylim([0,1])
% xlim([0,4])
% xticks([1,2,3])
% xticklabels({'ec-ac','ec-ex','ac-ex'})
% title([monkey ' ' neuron ' ec-ac ec-ex ac-ex'])
% legend({'mean corr','stderr','ec ub','ac ub','ex ub'})
% 
% 
% mean_corrs=mean(half_rsam_corrs_pos_only);
% stderr_corrs=std(half_rsam_corrs_pos_only)/sqrt(length(half_rsam_corrs_pos_only));
% 
% mean_ub=mean(upper_bounds_pos_only);
% stderr_ub=std(upper_bounds_pos_only)/sqrt(length(upper_bounds_pos_only));
% 
% subplot(1,2,2)
% hold on
% bar(mean_corrs)
% errorbar([1,2,3],mean_corrs,stderr_corrs)
% line([-1,5],[mean_ub(1),mean_ub(1)],'Color','red')
% line([-1,5],[mean_ub(1)+stderr_ub(1),mean_ub(1)+stderr_ub(1)],'Color','red','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(1)-stderr_ub(1),mean_ub(1)-stderr_ub(1)],'Color','red','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(2),mean_ub(2)],'Color','green')
% line([-1,5],[mean_ub(2)+stderr_ub(2),mean_ub(2)+stderr_ub(2)],'Color','green','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(2)-stderr_ub(2),mean_ub(2)-stderr_ub(2)],'Color','green','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(3),mean_ub(3)],'Color','blue')
% line([-1,5],[mean_ub(3)+stderr_ub(3),mean_ub(3)+stderr_ub(3)],'Color','blue','LineStyle',':','HandleVisibility','off')
% line([-1,5],[mean_ub(3)-stderr_ub(3),mean_ub(3)-stderr_ub(3)],'Color','blue','LineStyle',':','HandleVisibility','off')
% ylim([0,1])
% xlim([0,4])
% xticks([1,2,3])
% xticklabels({'ec-ac','ec-ex','ac-ex'})
% title([monkey ' ' neuron ' ec-ac ec-ex ac-ex pos sig neurons only'])
% legend({'mean corr','stderr','ec ub','ac ub','ex ub'})
% set(gcf,'position',[0,0,2000,1000])
% print(['D:\15386\corr_mat\rsa_by_trial\ratio_to_ub\decoding_' monkey '_' neuron '.png'],'-dpng','-painters','-loose',gcf)

%% plot ratio

mean_corrs=mean(half_rsam_corrs);
mean_ub=mean(upper_bounds);
mean_corrs=mean_corrs./[mean(mean_ub([1,2])),mean(mean_ub([1,3])),mean(mean_ub([2,3]))];
stderr_corrs=std(half_rsam_corrs)/sqrt(length(half_rsam_corrs));


subplot(1,2,1)
hold on
bar(mean_corrs)
errorbar([1,2,3],mean_corrs,stderr_corrs)
ylim([0,1])
xlim([0,4])
xticks([1,2,3])
xticklabels({'ec-ac','ec-ex','ac-ex'})
title([monkey ' ' neuron ' ec-ac ec-ex ac-ex'])
legend({'mean corr','stderr'})


mean_corrs=mean(half_rsam_corrs_pos_only);
mean_ub=mean(upper_bounds_pos_only);
mean_corrs=mean_corrs./[mean(mean_ub([1,2])),mean(mean_ub([1,3])),mean(mean_ub([2,3]))];

stderr_corrs=std(half_rsam_corrs_pos_only)/sqrt(length(half_rsam_corrs_pos_only));


subplot(1,2,2)
hold on
bar(mean_corrs)
errorbar([1,2,3],mean_corrs,stderr_corrs)
ylim([0,1])
xlim([0,4])
xticks([1,2,3])
xticklabels({'ec-ac','ec-ex','ac-ex'})
title([monkey ' ' neuron ' ec-ac ec-ex ac-ex pos sig neurons only'])
legend({'mean corr','stderr'})
set(gcf,'position',[0,0,2000,1000])
print(['rsa_by_trial\ratio_to_ub\neuron_num_controlled\decoding_' monkey '_' neuron '.png'],'-dpng','-painters','-loose',gcf)