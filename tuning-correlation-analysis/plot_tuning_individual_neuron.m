%stim_path='stimuli\';
stim_path='stimuli\';
nstimuli=51;

stim_type='ec';
stimlist_ec=cell(1,nstimuli-1);
stim_cnt=0;
if strcmp(stim_type,'ec')
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist_ec{stim_cnt}=[stim_type '_' num2str(idx,'%03d') '_l.png' ];
    end
else
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist_ec{stim_cnt}=[stim_type '_' num2str(idx,'%03d') '_sl.png' ];
    end
end

stim_type='ex';
stimlist_ex=cell(1,nstimuli-1);
stim_cnt=0;
if strcmp(stim_type,'ec')
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist_ex{stim_cnt}=[stim_type '_' num2str(idx,'%03d') '_l.png' ];
    end
else
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist_ex{stim_cnt}=[stim_type '_' num2str(idx,'%03d') '_sl.png' ];
    end
end

stim_type='ac';
stimlist_ac=cell(1,nstimuli-1);
stim_cnt=0;
if strcmp(stim_type,'ec')
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist_ac{stim_cnt}=[stim_type '_' num2str(idx,'%03d') '_l.png' ];
    end
else
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist_ac{stim_cnt}=[stim_type '_' num2str(idx,'%03d') '_sl.png' ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corr between types
datas=[
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


monkey='ko';
neuron='v1';
t_window='40-250';

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

chosen_layer=neurons_layer_all;
name_to_write=[monkey '_' neuron];

stim_type={'ec','ac','ex'};

% load(['tuning_all\' t_window '\' monkey '_' stim_type{1} '_all_tuning.mat'])
% tuning_1=all_tuning;
% load(['tuning_all\' t_window '\' monkey '_' stim_type{2} '_all_tuning.mat'])
% tuning_2=all_tuning;
% load(['tuning_all\' t_window '\' monkey '_' stim_type{3} '_all_tuning.mat'])
% tuning_3=all_tuning;
% 
% 
% length(setdiff(chosen_layer,all_bad_channels))
% tuning_1=tuning_1(setdiff(chosen_layer,all_bad_channels),:);
% tuning_2=tuning_2(setdiff(chosen_layer,all_bad_channels),:);
% tuning_3=tuning_3(setdiff(chosen_layer,all_bad_channels),:);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for gabors
gabor_type='simple';
load(['gabor_pyramid_resps\5x5_jitter_avg_area\3_scale_16_orientation\' gabor_type '\' stim_type{1} '.mat'])
tuning_1=center_resps';
load(['gabor_pyramid_resps\5x5_jitter_avg_area\3_scale_16_orientation\' gabor_type '\' stim_type{2} '.mat'])
tuning_2=center_resps';
load(['gabor_pyramid_resps\5x5_jitter_avg_area\3_scale_16_orientation\' gabor_type '\' stim_type{3} '.mat'])
tuning_3=center_resps';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for gabors
% gabor_type='64';
% load(['Sparse_coding_filters\responses\' gabor_type '_' stim_type{1} '.mat'])
% tuning_1=center_resps';
% load(['Sparse_coding_filters\responses\' gabor_type '_' stim_type{2} '.mat'])
% tuning_2=center_resps';
% load(['Sparse_coding_filters\responses\' gabor_type '_' stim_type{3} '.mat'])
% tuning_3=center_resps';

[tuning_corr_1_2,p_vals_1_2]=corr(tuning_1',tuning_2','Type','Pearson');
[tuning_corr_1_3,p_vals_1_3]=corr(tuning_1',tuning_3','Type','Pearson');
[tuning_corr_2_3,p_vals_2_3]=corr(tuning_2',tuning_3','Type','Pearson');

[sorted_tuning_1,order_1]=sort(tuning_1,2,'descend');
[sorted_tuning_2,order_2]=sort(tuning_2,2,'descend');
[sorted_tuning_3,order_3]=sort(tuning_3,2,'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% remove top stim
% tuning_1(:,[order_1(1),order_2(1),order_3(1)])=[];
% tuning_2(:,[order_1(1),order_2(1),order_3(1)])=[];
% tuning_3(:,[order_1(1),order_2(1),order_3(1)])=[];

size_tuning=size(tuning_1);
n_neuron=size_tuning(1);
%n_neuron=length(tuning_1);

corrs_pair_1=diag(tuning_corr_1_2,0);
corrs_pair_2=diag(tuning_corr_1_3,0);
corrs_pair_3=diag(tuning_corr_2_3,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Histogram
% x_edge=-1:0.05:1;
% y_max=20;
% curr_fig=figure();
% subplot(1,3,1)
% histogram(corrs_pair_1,x_edge);
% title(["EC vs AC" ['Average=' num2str(mean(corrs_pair_1)) newline 'Std error=' num2str(std(corrs_pair_1)/sqrt(n_neuron))]])
% xlabel("Pearson Correlation")
% ylabel("Count")
% ylim([0,y_max])
% subplot(1,3,2)
% histogram(corrs_pair_2,x_edge);
% title(["EC vs EX" ['Average=' num2str(mean(corrs_pair_2)) newline 'Std error=' num2str(std(corrs_pair_2)/sqrt(n_neuron))]])
% xlabel("Pearson Correlation")
% ylabel(['Count' newline  ])
% ylim([0,y_max])
% subplot(1,3,3)
% histogram(corrs_pair_3,x_edge);
% title(["AC vs EX" ['Average=' num2str(mean(corrs_pair_3)) newline 'Std error=' num2str(std(corrs_pair_3)/sqrt(n_neuron))]])
% xlabel("Pearson Correlation")
% ylabel("Count")
% ylim([0,y_max])
% 
% set(gcf,'position',[0,0,1800,600])
% saveas(curr_fig,['tuning_curves_20220118\tuning_' monkey '_' neuron ],'png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Tuning & Scatter
sz=200;
for i=[163,12,177]%[1,38,45]%1:n_neuron
    curr_fig=figure();
    subplot(2,3,1)
    mapFRstimuli('ex','stimuli',0,stimlist_ec(order_1(i,:)),sorted_tuning_1(i,:),10,0,sorted_tuning_1(i,1),stim_path,1);
    xticks([]);
    yticks([]);
    title(upper(stim_type{1}))
    subplot(2,3,2)
    mapFRstimuli('ac','stimuli',0,stimlist_ac(order_2(i,:)),sorted_tuning_2(i,:),10,0,sorted_tuning_2(i,1),stim_path,1);
    xticks([]);
    yticks([]);
    title(upper(stim_type{2}))
    subplot(2,3,3)
    mapFRstimuli('ex','stimuli',0,stimlist_ex(order_3(i,:)),sorted_tuning_3(i,:),10,0,sorted_tuning_3(i,1),stim_path,1);
    xticks([]);
    yticks([]);
    title(upper(stim_type{3}))
    
    subplot(2,3,4)
    hold on
    scatter((tuning_1(i,:)-min(tuning_1(i,:)))./(max(tuning_1(i,:))-min(tuning_1(i,:))),(tuning_2(i,:)-min(tuning_2(i,:)))./(max(tuning_2(i,:))-min(tuning_2(i,:))),sz,'.');
    p=polyfit((tuning_1(i,:)-min(tuning_1(i,:)))./(max(tuning_1(i,:))-min(tuning_1(i,:))),(tuning_2(i,:)-min(tuning_2(i,:)))./(max(tuning_2(i,:))-min(tuning_2(i,:))),1);
    fit_x = linspace(0,1);
    fit_y = polyval(p,fit_x);
    plot(fit_x,fit_y);
    axis equal
    xlabel(upper(stim_type{1}))
    ylabel(upper(stim_type{2}))
    xlim([0,1])
    ylim([0,1])
    %title([stim_type{1} ' vs ' stim_type{2} newline 'corr=' num2str(tuning_corr_1_2(i,i)) newline 'is\_significant=' num2str(p_vals_1_2(i,i)<0.05)])
    title([upper(stim_type{1}) ' vs ' upper(stim_type{2}) newline 'corr=' num2str(tuning_corr_1_2(i,i)) newline 'p val=' num2str(p_vals_1_2(i,i))])
    
    subplot(2,3,5)
    hold on
    scatter((tuning_1(i,:)-min(tuning_1(i,:)))./(max(tuning_1(i,:))-min(tuning_1(i,:))),(tuning_3(i,:)-min(tuning_3(i,:)))./(max(tuning_3(i,:))-min(tuning_3(i,:))),sz,'.');
    p=polyfit((tuning_1(i,:)-min(tuning_1(i,:)))./(max(tuning_1(i,:))-min(tuning_1(i,:))),(tuning_3(i,:)-min(tuning_3(i,:)))./(max(tuning_3(i,:))-min(tuning_3(i,:))),1);
    fit_x = linspace(0,1);
    fit_y = polyval(p,fit_x);
    plot(fit_x,fit_y);
    axis equal
    xlabel(upper(stim_type{1}))
    ylabel(upper(stim_type{3}))
    xlim([0,1])
    ylim([0,1])
    %title([stim_type{1} ' vs ' stim_type{3} newline 'corr=' num2str(tuning_corr_1_3(i,i)) newline 'is\_significant=' num2str(p_vals_1_3(i,i)<0.05)])
    title([upper(stim_type{1}) ' vs ' upper(stim_type{3}) newline 'corr=' num2str(tuning_corr_1_3(i,i)) newline 'p val=' num2str(p_vals_1_3(i,i))])
    
    subplot(2,3,6)
    hold on
    scatter((tuning_2(i,:)-min(tuning_2(i,:)))./(max(tuning_2(i,:))-min(tuning_2(i,:))),(tuning_3(i,:)-min(tuning_3(i,:)))./(max(tuning_3(i,:))-min(tuning_3(i,:))),sz,'.');
    p=polyfit((tuning_2(i,:)-min(tuning_2(i,:)))./(max(tuning_2(i,:))-min(tuning_2(i,:))),(tuning_3(i,:)-min(tuning_3(i,:)))./(max(tuning_3(i,:))-min(tuning_3(i,:))),1);
    fit_x = linspace(0,1);
    fit_y = polyval(p,fit_x);
    plot(fit_x,fit_y);
    axis equal
    xlabel(upper(stim_type{2}))
    ylabel(upper(stim_type{3}))
    xlim([0,1])
    ylim([0,1])
    %title([stim_type{2} ' vs ' stim_type{3} newline 'corr=' num2str(tuning_corr_2_3(i,i)) newline 'is\_significant=' num2str(p_vals_2_3(i,i)<0.05)])
    title([upper(stim_type{2}) ' vs ' upper(stim_type{3}) newline 'corr=' num2str(tuning_corr_2_3(i,i)) newline 'p val=' num2str(p_vals_2_3(i,i))])
    
    set(gcf,'position',[0,0,2000,700])
    %saveas(curr_fig,['tuning_curves_20220118\tuning_' monkey '_' neuron '_' num2str(i)],'png')
    %print(['tuning_curves_20220118\sc' '_' gabor_type '_' num2str(i) '.png'],'-dpng','-painters','-loose',gcf)
    print([ gabor_type '_' num2str(i) '.eps'],'-depsc','-painters','-loose',gcf)
    close(curr_fig)
    
%     curr_fig=figure();
%     plot(sorted_tuning_1(i,:))
%     hold on
%     plot(sorted_tuning_2(i,:))
%     plot(sorted_tuning_3(i,:))
%     xlabel('Rank Order')
%     ylabel('Average Firing Rate')
%     legend({'EC','AC','EX'})
%     print(['tuning_curves_20220118\tuning_curve' '_' neuron '_' num2str(i) '.png'],'-dpng','-painters','-loose',gcf)
%     close(curr_fig)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% for tuning curves
% plot(sorted_tuning_1(i,:))
% hold on
% plot(sorted_tuning_2(i,:))
% plot(sorted_tuning_3(i,:))
% xlabel('Rank Order')
% ylabel('Average Firing Rate')
% legend({'EC','AC','EX'})
% print(['tuning_curves_20220118\even_' monkey '_' neuron '_' num2str(i) '.eps'],'-depsc','-painters','-loose',gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot stim patches in the orginal order
% mapFRstimuli('ex','stimuli',0,stimlist_ex,sorted_tuning_1(i,:),10,0,sorted_tuning_1(i,1),'stimuli\',1);
% xticks([]);
% yticks([]);
% set(gcf,'position',[0,0,1000,500])
% print(['tuning_curves_20220118\EX_original_order.eps'],'-depsc','-painters','-loose',gcf)
