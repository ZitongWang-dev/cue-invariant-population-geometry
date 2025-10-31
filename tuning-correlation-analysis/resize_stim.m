%stim_dir='stimuli\';
stim_dir='stimuli\';

nstimuli=51;
stimlist=cell(1,nstimuli-1);
stim_cnt=0;

%sizes = [7, 28,52];

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
% for idx=1:nstimuli-1
%     stim_cnt=stim_cnt+1;
%     stimlist{stim_cnt}=[stimtype '_' num2str(idx,'%03d') '_l0.png' ];
% end

%stim_sizes=[11,51,99,131,163,195];
stim_sizes=[224];
for s=1:length(stim_sizes)
    new_size=stim_sizes(s);
    for stim_idx=1:nstimuli-1
        curr_stim=stimlist{stim_idx};
        im=imread([stim_dir curr_stim]);
        
        if ~strcmp(stimtype,'ec')
            im=im(81:end-80,81:end-80);
        end
        im=imresize(im,[new_size,new_size]);
        
        new_im=im;    %newsize=224
        %new_im=uint8(ones(224,224)).*128;
        %new_im(112-round(new_size/2):111+new_size-round(new_size/2),112-round(new_size/2):111+new_size-round(new_size/2))=im;
%             new_im=uint8(ones(227,227)).*128;
%             new_im(114-round(new_size/2):113+new_size-round(new_size/2),114-round(new_size/2):113+new_size-round(new_size/2))=im;
        imshow(new_im)
        imwrite(mat2gray(new_im,[0,255]), ['resized_stim/224/' num2str(new_size) '/' curr_stim]);
        
        %     im_rgb=repmat(im,[1 1 3]);
        %     resps=activations(net,im_rgb,'conv1_2','OutputAs','channels');
        %     resp_size=size(resps);
        %     curr_mean_resps=reshape(mean(mean(resps)),1,resp_size(3));
        %     mean_resps=[mean_resps;curr_mean_resps];
    end
end