%stim_dir='stimuli_small\pngs\';
stim_dir='stimuli_small\pngs\';

nstimuli=51;
stimlist=cell(1,nstimuli-1);
stim_cnt=0;

% sizes = [7, 18, 19, 52, 151];
%sizes = [7, 11, 14, 21, 35, 42, 49, 56, 60];

%sizes = [20,30,40,50,60];
%orientations = 0:10:179;

%sizes=[20,40,60];
sizes=[20];
orientations = 0:12:180;

gabors=cell(1,length(orientations)*length(sizes));

stimtype='ac';
gabor_type='complex';
if strcmp(stimtype,'ec')
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist{stim_cnt}=[stimtype '_' num2str(idx,'%03d') '_l0.png' ];
    end
else
    for idx=1:nstimuli-1
        stim_cnt=stim_cnt+1;
        stimlist{stim_cnt}=[stimtype '_' num2str(idx,'%03d') '_l0.png' ];
    end
end

%new_size=52;

idx=0;
for i=1:length(sizes)
    for j=1:length(orientations)
        idx=idx+1;
        gabors{idx}=imread(sprintf('gabors/%d_%d_odd.png',sizes(i),orientations(j)));
%         idx=idx+1;
%         gabors{idx}=imread(sprintf('gabors/%d_%d_even.png',sizes(i),orientations(j)));
    end
end

gabors_2=cell(1,length(orientations)*length(sizes));
idx=0;
for i=1:length(sizes)
    for j=1:length(orientations)
%         idx=idx+1;
%         gabors{idx}=imread(sprintf('gabors/%d_%d_odd.png',sizes(i),orientations(j)));
        idx=idx+1;
        gabors_2{idx}=imread(sprintf('gabors/%d_%d_even.png',sizes(i),orientations(j)));
    end
end

center_resps=[];
for stim_idx=1:nstimuli-1
    curr_stim=stimlist{stim_idx};
    im=imread([stim_dir curr_stim]);
    if ~strcmp(stimtype,'ec')
        im=im(41:end-40,41:end-40);
    end
%     tmp_size=35;
%     im=imresize(im,[tmp_size,tmp_size]);
%     new_im=uint8(ones(227,227)).*128;
%     new_im(114-round(tmp_size/2):113+tmp_size-round(tmp_size/2),114-round(tmp_size/2):113+tmp_size-round(tmp_size/2))=im;
%     imshow(new_im)
%     im=double(new_im);
    im=double(im);
    
    curr_stim_resps=[];
    for idx=1:length(orientations)*length(sizes)
%         curr_filter=double(gabors_2{idx});
%         gabor_size=length(curr_filter);
%         curr_filter=curr_filter./255.*2-1;
%         coef=conv2(im,curr_filter);
%         coef=coef(gabor_size:end-gabor_size+1,gabor_size:end-gabor_size+1);
%         coef_size=length(coef);
%         %coef=coef.^2;
%         %coef=abs(coef);
%         
%         % random jittering
%         %pos1=randi(5)-3;
%         %pos2=randi(5)-3;
%         
%         %curr_stim_resps=[curr_stim_resps coef(round(coef_size/2)+pos1,round(coef_size/2)+pos2)];
%         
%         % average/max of random jittering
%         pos1=-2;
%         pos2=2;
%         % no jittering
%         jitter_range=0;
%         
%         curr_stim_resps=[curr_stim_resps max(max(coef(round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range,round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range)))];
        
        curr_filter_1=double(gabors{idx});
        gabor_size=length(curr_filter_1);
        curr_filter_1=curr_filter_1./255.*2-1;
        % L2 nornalization
        curr_filter_1=(curr_filter_1)./(sum(sum(curr_filter_1.^2)));
        coef=conv2(im,curr_filter_1);
        coef=coef(gabor_size:end-gabor_size+1,gabor_size:end-gabor_size+1);
        coef_size=length(coef);
        coef_sq=coef.^2;
        
        
        curr_filter_2=double(gabors_2{idx});
        gabor_size=length(curr_filter_2);
        curr_filter_2=curr_filter_2./255.*2-1;
        % L2 nornalization
        curr_filter_2=(curr_filter_2)./(sum(sum(curr_filter_2.^2)));
        coef_2=conv2(im,curr_filter_2);
        coef_2=coef_2(gabor_size:end-gabor_size+1,gabor_size:end-gabor_size+1);
        coef_size=length(coef_2);
        coef_2_sq=coef_2.^2;
        %coef=sqrt(coef+coef_2);
        if strcmp(gabor_type, 'odd')
            coef = coef;
        elseif strcmp(gabor_type, 'even')
            coef = coef_2;
        elseif strcmp(gabor_type, 'odd_sq')
            coef = coef_sq;
        elseif strcmp(gabor_type, 'even_sq')
            coef = coef_2_sq;
        elseif strcmp(gabor_type, 'complex')
            coef = coef_sq+coef_2_sq;
        end
        
        %pos1=randi(5)-3;
        %pos2=randi(5)-3;
        %curr_stim_resps=[curr_stim_resps coef(round(coef_size/2)+pos1,round(coef_size/2)+pos2)];
        
        %jitter_range=0;
        %curr_stim_resps=[curr_stim_resps max(max(coef(round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range,round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range)))];
        
        %no jitter
        %curr_stim_resps=[curr_stim_resps max(max(coef(round(coef_size/2),round(coef_size/2))))];
        
        %with jitter, 10 avg
%         n_rep = 10;
%         avg_resp = 0;
%         for rep=1:n_rep
%             pos1=randi(5)-3;
%             pos2=randi(5)-3;
%             avg_resp = avg_resp+coef(round(coef_size/2)+pos1,round(coef_size/2)+pos2)/n_rep;
%         end
% 
%         curr_stim_resps=[curr_stim_resps avg_resp];
        
        jitter_range = 1;
        %with jitter, average area
        %curr_stim_resps=[curr_stim_resps mean(mean(coef(round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range,round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range)))];
        
        %trials
        %center_range_resp_mat = coef(round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range,round(coef_size/2)-jitter_range:round(coef_size/2)+jitter_range);
        %curr_stim_resps=[curr_stim_resps, center_range_resp_mat(:)];

        center_range_resp_mat = coef(round(coef_size/2)-jitter_range:round(coef_size/2),round(coef_size/2)-jitter_range:round(coef_size/2));
        curr_stim_resps=[curr_stim_resps, center_range_resp_mat(:)];

    end
    %avg
    %center_resps=[center_resps;curr_stim_resps];
    %trials
    center_resps = cat(3,center_resps,curr_stim_resps);
end
%save(['gabor_pyramid_resps\no_jitter\3_scale_16_orientation\' gabor_type '\' stimtype '.mat'],'center_resps')
%save(['gabor_pyramid_resps\2x2_jitter_trials\3_scale_16_orientation\' gabor_type '\' stimtype '.mat'],'center_resps')