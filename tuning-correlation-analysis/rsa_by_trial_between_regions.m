monkey='ko';
means=[];
stds=[];
n_neuron=40;


neuron=neuron_setting{1};
load(['tuning_all\' monkey '_v1_bytrial_40-250.mat'])
tunings_v1=tunings;
load(['tuning_all\' monkey '_v2_bytrial_40-250.mat'])
tunings_v2=tunings;

half_rsam_corrs=[];
upper_bounds=[];
half_rsam_corrs_pos_only=[];
upper_bounds_pos_only=[];

n_trial=10;
parfor dup=1:100
    
    curr_randperm=randperm(n_trial);
    i_trials_1=curr_randperm(1:5);
    i_trials_2=curr_randperm(6:10);
    
    n_neuron=size(tunings_v1{1});
    n_neuron=n_neuron(1);
    neuron_sampling=randperm(n_neuron);
    neuron_sampling=neuron_sampling(1:n_neuron);
    
    tuning_ec_1=mean(tunings_v1{1}(neuron_sampling,:,i_trials_1),3);
    tuning_ec_2=mean(tunings_v2{1}(neuron_sampling,:,i_trials_1),3);
    
    tuning_ac_1=mean(tunings_v1{2}(neuron_sampling,:,i_trials_1),3);
    tuning_ac_2=mean(tunings_v2{2}(neuron_sampling,:,i_trials_1),3);
    
    tuning_ex_1=mean(tunings_v1{3}(neuron_sampling,:,i_trials_1),3);
    tuning_ex_2=mean(tunings_v2{3}(neuron_sampling,:,i_trials_1),3);
    
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
    
    half_rsam_corrs=[half_rsam_corrs;[corr(half_ec_1,half_ec_2),corr(half_ac_1,half_ac_2),corr(half_ex_1,half_ex_2)]];
end

mean_corrs=mean(half_rsam_corrs);
std_corrs=std(half_rsam_corrs);

means=[means;mean_corrs];
stds=[stds;std_corrs];
[h1,p1]=ttest(half_rsam_corrs(:,1),half_rsam_corrs(:,2));
[h2,p2]=ttest(half_rsam_corrs(:,1),half_rsam_corrs(:,3));
[h3,p3]=ttest(half_rsam_corrs(:,2),half_rsam_corrs(:,3));
neuron
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
xticklabels({'ec','ac','ex'})
title([monkey ' v1 vs. v2'])
legend({'mean','std'})


set(gcf,'position',[0,0,2000,1000])
print(['rsa_by_trial\ratio_to_ub\neuron_num_controlled\between_regions\decoding_' monkey '.png'],'-dpng','-painters','-loose',gcf)