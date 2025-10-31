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
n_neuron=100; %40

stim_all_idx=1:50;
stim_complex_idx=[8,10,25,29,37];
stim_simple_idx=[1,2,3,4,5,6,7,9,12,13,14,15,16,17,20,21,23,24,40,46];
stim_middle_idx=[11,18,19,22,26,27,28,30,31,32,33,34,35,36,38,39,41,42,43,44,45,47,48,49,50];

stim_set_name='all';
curr_stim_set=stim_all_idx;


mat_diag_idx=[1];
dec_idx=length(curr_stim_set):-1:2;
for i=1:length(curr_stim_set)-1
    mat_diag_idx=[mat_diag_idx,1+sum(dec_idx(1:i))];
end
%mat_diag_idx=[];


means=[];
stds=[];
lower_cis=[];
upper_cis=[];
final_ubs=[];
for neuron_setting={'v1','v2'}
    neuron=neuron_setting{1};
    load(['tuning_all\' monkey '_' neuron '_bytrial_40-250.mat'])
    
    half_rsam_corrs=[];
    upper_bounds=[];
    half_rsam_corrs_pos_only=[];
    upper_bounds_pos_only=[];
    
    n_trial=10;
    parfor dup=1:1000
        
        curr_randperm=randperm(n_trial);
        i_trials_1=curr_randperm(1:5);
        i_trials_2=curr_randperm(6:10);
        
        n_neuron=size(tunings{1});
        n_neuron=n_neuron(1);
        neuron_sampling=randperm(n_neuron);
        neuron_sampling=neuron_sampling(1:n_neuron);
        
        tuning_ec_1=mean(tunings{1}(neuron_sampling,:,i_trials_1),3);
        tuning_ec_2=mean(tunings{1}(neuron_sampling,:,i_trials_2),3);
        tuning_ec_1=tuning_ec_1(:,curr_stim_set)
        tuning_ec_2=tuning_ec_2(:,curr_stim_set);
        
        tuning_ac_1=mean(tunings{2}(neuron_sampling,:,i_trials_1),3);
        tuning_ac_2=mean(tunings{2}(neuron_sampling,:,i_trials_2),3);
        tuning_ac_1=tuning_ac_1(:,curr_stim_set);
        tuning_ac_2=tuning_ac_2(:,curr_stim_set);
        
        tuning_ex_1=mean(tunings{3}(neuron_sampling,:,i_trials_1),3);
        tuning_ex_2=mean(tunings{3}(neuron_sampling,:,i_trials_2),3);
        tuning_ex_1=tuning_ex_1(:,curr_stim_set)
        tuning_ex_2=tuning_ex_2(:,curr_stim_set);
        
        
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
        
        half_ec_1=half_ec_1(setdiff(1:length(half_ec_1),mat_diag_idx));
        half_ac_1=half_ac_1(setdiff(1:length(half_ac_1),mat_diag_idx));
        half_ex_1=half_ex_1(setdiff(1:length(half_ex_1),mat_diag_idx));
        
        half_ec_2=rsam_ec_2(tril(true(size(rsam_ec_2))));
        half_ac_2=rsam_ac_2(tril(true(size(rsam_ac_2))));
        half_ex_2=rsam_ex_2(tril(true(size(rsam_ex_2))));

        half_ec_2=half_ec_2(setdiff(1:length(half_ec_2),mat_diag_idx));
        half_ac_2=half_ac_2(setdiff(1:length(half_ac_2),mat_diag_idx));
        half_ex_2=half_ex_2(setdiff(1:length(half_ex_2),mat_diag_idx));
        
        half_rsam_corrs=[half_rsam_corrs;[corr(half_ec_1,half_ac_1),corr(half_ec_1,half_ex_1),corr(half_ac_1,half_ex_1)]];
        upper_bounds=[upper_bounds;[corr(half_ec_1,half_ec_2),corr(half_ac_1,half_ac_2),corr(half_ex_1,half_ex_2)]];
        
%         % pos only
%         %%%%%%%%%%%%%%%%%%%
%         [tuning_corr_1,p_vals_1]=corr(tuning_ec_1',tuning_ac_1','Type','Pearson');
%         ident_diag_1=diag(tuning_corr_1,0);
%         [tuning_corr_2,p_vals_2]=corr(tuning_ec_1',tuning_ex_1','Type','Pearson');
%         ident_diag_2=diag(tuning_corr_2,0);
%         [tuning_corr_3,p_vals_3]=corr(tuning_ac_1',tuning_ex_1','Type','Pearson');
%         ident_diag_3=diag(tuning_corr_3,0);
%         
%         % positive significant correlations
%         curr_sigs_pos_1=(diag(p_vals_1,0)<=0.05)&(ident_diag_1>0);
%         curr_sigs_pos_2=(diag(p_vals_2,0)<=0.05)&(ident_diag_2>0);
%         curr_sigs_pos_3=(diag(p_vals_3,0)<=0.05)&(ident_diag_3>0);
%         
%         m_ec_1_ps=(tuning_ec_1(curr_sigs_pos_1,:)'-mean(tuning_ec_1(curr_sigs_pos_1,:)'))';
%         m_ac_1_ps=(tuning_ac_1(curr_sigs_pos_2,:)'-mean(tuning_ac_1(curr_sigs_pos_2,:)'))';
%         m_ex_1_ps=(tuning_ex_1(curr_sigs_pos_3,:)'-mean(tuning_ex_1(curr_sigs_pos_3,:)'))';
%         
%         m_ec_2_ps=(tuning_ec_2(curr_sigs_pos_1,:)'-mean(tuning_ec_2(curr_sigs_pos_1,:)'))';
%         m_ac_2_ps=(tuning_ac_2(curr_sigs_pos_2,:)'-mean(tuning_ac_2(curr_sigs_pos_2,:)'))';
%         m_ex_2_ps=(tuning_ex_2(curr_sigs_pos_3,:)'-mean(tuning_ex_2(curr_sigs_pos_3,:)'))';
%         
%         rsam_ec_1_ps=corr(m_ec_1_ps,m_ec_1_ps,'Type','Pearson');
%         rsam_ac_1_ps=corr(m_ac_1_ps,m_ac_1_ps,'Type','Pearson');
%         rsam_ex_1_ps=corr(m_ex_1_ps,m_ex_1_ps,'Type','Pearson');
%         
%         rsam_ec_2_ps=corr(m_ec_2_ps,m_ec_2_ps,'Type','Pearson');
%         rsam_ac_2_ps=corr(m_ac_2_ps,m_ac_2_ps,'Type','Pearson');
%         rsam_ex_2_ps=corr(m_ex_2_ps,m_ex_2_ps,'Type','Pearson');
%         
%         half_ec_1_ps=rsam_ec_1_ps(tril(true(size(rsam_ec_1_ps))));
%         half_ac_1_ps=rsam_ac_1_ps(tril(true(size(rsam_ac_1_ps))));
%         half_ex_1_ps=rsam_ex_1_ps(tril(true(size(rsam_ex_1_ps))));
%         
%         half_ec_2_ps=rsam_ec_2_ps(tril(true(size(rsam_ec_2_ps))));
%         half_ac_2_ps=rsam_ac_2_ps(tril(true(size(rsam_ac_2_ps))));
%         half_ex_2_ps=rsam_ex_2_ps(tril(true(size(rsam_ex_2_ps))));
%         
%         half_rsam_corrs_pos_only=[half_rsam_corrs_pos_only;[corr(half_ec_1_ps,half_ac_1_ps),corr(half_ec_1_ps,half_ex_1_ps),corr(half_ac_1_ps,half_ex_1_ps)]];
%         upper_bounds_pos_only=[upper_bounds_pos_only;[corr(half_ec_1_ps,half_ec_2_ps),corr(half_ac_1_ps,half_ac_2_ps),corr(half_ex_1_ps,half_ex_2_ps)]];
    end
    
    mean_corrs=mean(half_rsam_corrs);
    mean_ub=mean(upper_bounds);
    %mean_corrs=mean_corrs./[mean(mean_ub([1,2])),mean(mean_ub([1,3])),mean(mean_ub([2,3]))];
    final_ub=[mean(mean_ub([1,2])),mean(mean_ub([1,3])),mean(mean_ub([2,3]))];
    
    %ratio_dists=sort(half_rsam_corrs./[mean(mean_ub([1,2])),mean(mean_ub([1,3])),mean(mean_ub([2,3]))],1);
    ratio_dists=sort(half_rsam_corrs);
    lower_ci=ratio_dists(25,:);
    upper_ci=ratio_dists(975,:);
    
    std_corrs=std(half_rsam_corrs);
    %mean_corrs=mean_ub
    
    means=[means;mean_corrs];
    stds=[stds;std_corrs];
    lower_cis=[lower_cis;lower_ci];
    upper_cis=[upper_cis;upper_ci];
    final_ubs=[final_ubs;final_ub];

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
[h1,p11]=ttest2(dist1(:,1),dist2(:,1));
[h2,p12]=ttest2(dist1(:,2),dist2(:,2));
[h3,p13]=ttest2(dist1(:,3),dist2(:,3));
p1s=[p11;p12;p13]

p21=ranksum(dist1(:,1),dist2(:,1));
p22=ranksum(dist1(:,2),dist2(:,2));
p23=ranksum(dist1(:,3),dist2(:,3));
p2s=[p21;p22;p23]


means=means';
stds=stds';
lower_cis=lower_cis';
upper_cis=upper_cis';
lower_cis=means-lower_cis;
upper_cis=upper_cis-means;
final_ubs=final_ubs';

%means=final_ubs';

%hb=bar(means);
bars_to_plot=[final_ubs, means];
hb=bar(bars_to_plot);

set(hb(1,3),'LineStyle','--','FaceColor','#0072BD')
set(hb(1,4),'LineStyle','--','FaceColor','#D95319')


hold on
for pos = [3,4]%1:size(means,2)
    xpos=hb(pos).XData+hb(pos).XOffset;
    %errorbar(xpos,means(:,pos),stds(:,pos),'LineStyle','none','Color','black')
    %errorbar(xpos,means(:,pos),lower_cis(:,pos),upper_cis(:,pos),'LineStyle','none','Color','black')
    errorbar(xpos,bars_to_plot(:,pos),lower_cis(:,pos-2),upper_cis(:,pos-2),'LineStyle','none','Color','black')
end
ylim([-0.1,1])
xlim([0,4])
xticks([1,2,3])
xticklabels({['ec-ac'],'ec-ex','ac-ex'})
title([monkey ' ' stim_set_name ' v1&v2 ec-ac ec-ex ac-ex RSA corr. upper bounds'])
%title([monkey ' ' stim_set_name ' v1&v2 ec-ac ec-ex ac-ex'])
%legend({'v1','v2'})
legend({'v1 upper bound','v2 upper bound', 'v1', 'v2', '95% CI'})

cilb=means-lower_cis
ciub=means+upper_cis

set(gcf,'position',[0,0,2000,1000])
%print(['rsa_by_trial\ratio_to_ub\neuron_num_controlled\decoding_' monkey '.png'],'-depsc','-painters','-loose',gcf)
print(['rsa_' monkey ' ' stim_set_name '_raw_n_ub.png'],'-dpng','-painters','-loose',gcf)
print(['rsa_' monkey ' ' stim_set_name '_raw_n_ub.eps'],'-depsc','-painters','-loose',gcf)