function [hdlfig, im_item_size]=mapFRstimuli(type,display_type,level,stimlist,meanFRlist,nbcol,border,maxFR,root_path,subplot_use);
%function [im_item_size]=mapFRstimuli(type,display_type,level,stimlist,meanFRlist,nbcol,border,maxFR,root_path,subplot_use);

%[im_item_size]=mapFRstimuli(type,display,stimlist,meanFRlist,nbcol,maxFR,root_path);
%
%   plot FR map of list of stimuli
%
%
%   type: type of stimuli
%   display_type: what to display: map, edge, stimuli
%   leve: level of the concepts -1 (156 stimuli) 0 (50 stimuli) >0
%   (sublevels)
%   ....
%   nbcol: number of stimuli per row (number of columns)
%   subplot_use: if use in a subplot of map then circhsift individual maps
%
% Corentin CNBC 03/13/13

%nb_stim=length(stimlist);
nb_stim=size(meanFRlist,2);%to match number in meanFRlist if different from the total number of stimuli

if level==-1,
    switch type
        case 'ec' %edge concept with 157 stimuli
            stim_path=[root_path];
            
            im_res=[];
            border=0;%0 %to display only stimuli;%border witdth
            
            for i = 1:nb_stim
                item_name=stimlist{i};
                im_item=imread(strcat(stim_path , item_name));
                
                im_item_size=size(im_item,1)+border;
                meanFR=meanFRlist(i);
                
                switch display_type
                    case 'edge'
                        %to display mean FR
                        %adjust color with meanFR
                        ind1=find(im_item~=max(max(im_item)));ind2=find(im_item==max(max(im_item)));
                        im_item(ind1)=meanFR;
                        im_item(ind2)=0;
                        
                        im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=im_item;
                        range=[0 maxFR];cmap='jet';
                        
                    case 'map'
                        %to display only meanFR like appearance stimuli
                        c_mask=meanFR*ones(im_item_size-border,im_item_size-border);
                        im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=c_mask;
                        range=[0 maxFR];cmap='jet';
                        
                    case 'stimuli'
                        %to display only stimuli
                        im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=im_item;
                        range=[0 255];cmap='gray';
                end
            end
            
            %hdlfig=figure;
            imagesc(im_res,range);colormap(cmap);axis square;
            
            
        case 'ap' %appearance concept with 157 stimuli
            stim_path=[root_path];
            im_res=[];im_res_item=[];b=0;%border witdth
            
            for i = 1:nb_stim
                item_name=stimlist{i};
                im_item=double(imread(strcat(stim_path , item_name)));
                
                im_item_size=size(im_item,1)+border;
                meanFR=meanFRlist(i);
                
                switch display_type
                    case 'map'
                        %to display mean FR
                        c_mask=meanFR*ones(im_item_size-border,im_item_size-border);
                        im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=c_mask;
                        range=[0 maxFR];cmap='jet';
                        
                    case 'stimuli'
                        %to display only stimuli
                        im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=im_item;
                        range=[0 255];cmap='gray';
                end
            end
            
            %hdlfig=figure;
            imagesc(im_res,range);colormap(cmap);axis square;                        
    end
    
    
else
    
    
    %concepts for diferent level with 51 stimuli
    if level==0,
        stim_path=[root_path ];
    else
        stim_path=[root_path ];
    end    
    im_res=[];
    border=0;%0;%border witdth
    
    for i = 1:nb_stim
        item_name=stimlist{i};
        im_item=double(imread(strcat(stim_path , item_name)));
        
        im_item_size=size(im_item,1)+border;
        meanFR=meanFRlist(i);
 
        switch display_type
            
            case 'map'
                %to display mean FR
                c_mask=meanFR*ones(im_item_size-border,im_item_size-border);
                im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=c_mask;
                range=[0 maxFR];cmap='jet';
                
            case 'edge'
                %to display stimuli (edge)
                %adjust color with meanFR
                ind1=find(im_item~=max(max(im_item)));ind2=find(im_item==max(max(im_item)));
                im_item(ind1)=meanFR;
                im_item(ind2)=0;
                im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=im_item;
                range=[0 maxFR];cmap='jet';
                
            case 'stimuli'
                %to display only stimuli
                if strcmp(type,'ec')
                    border=5;im_item_size=size(im_item,1)+border;
                    if i==1,im_res=128*ones(im_item_size+im_item_size*floor((nb_stim)/nbcol),im_item_size*nbcol);end;
                end
               
                im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=im_item;
                range=[0 255];cmap='gray';
                
            case 'stimuli_label'
                %to display only stimuli with  acolor coded label (ex: level number)
                im_item_size_l=im_item_size+border_label;
                im_item_size_c=im_item_size;
                
                if strcmp(type,'ec')
                    border=5;im_item_size=size(im_item,1)+border;
                    if i==1,im_res=128*ones(im_item_size_l+im_item_size_l*floor((nb_stim)/nbcol),im_item_size*nbcol);end;
                end
               
                im_res(1+im_item_size*floor((i-1)/nbcol):im_item_size+im_item_size*floor((i-1)/nbcol)-border, 1+im_item_size*(mod(i-1,nbcol)):im_item_size+im_item_size*(mod(i-1,nbcol))-border)=im_item;
                range=[0 255];cmap='gray';
                
        end
    end
    %original plot setting
    %%%%%%%%%%%%%%%
    %hdlfig=figure;
    if subplot_use
        imagesc(circshift(im_res,2),range);colormap(cmap);%axis square;
    else
        imagesc(im_res,range);colormap(cmap);%axis square;
    end
    
    %return figure handler
    %%%%%%%%%%%%%%%%
    %hdlfig=figure;
    %if subplot_use
    %    imagesc(circshift(im_res,2),range);colormap(cmap);axis square;
    %else
    %    imagesc(im_res,range);colormap(cmap);axis square;
    %end
    
end