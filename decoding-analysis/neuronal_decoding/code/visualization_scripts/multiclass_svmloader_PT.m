function [data,label] = multiclass_svmloader_PT(array,timewindow)
% array: 4d array
% data: row-each data points, column-neurons, each block-spike count in
% certain time range

[neuron_num,condition_num] = size(array);
% change this repeat extraction in later version, since repeats can vary in
% each condition
repeat = 10;
%create data&label
data = zeros(condition_num*repeat,neuron_num);
label = zeros(condition_num*repeat,1);



% Ten labels corresponding to ten images

conditon2label = [1:condition_num];
%trick: Since the noise_pattern_num is 10, every ten consecutive conditon_idx have the same noise pattern, different labels (label for images)


% iterate through conditons

for con = 1:condition_num
    con_temp = cell2mat(array(:,con)); %take a column of condition in the 4d array, dim:[neuronxrepeat,1500 ms]
    window = con_temp(:,timewindow(1):timewindow(2));
    spike_count = sum(window,2); %data with repeat scattered in this vector
    data_temp = reshape(spike_count,repeat,neuron_num); %dim [repeat, neuron]
    data(con*repeat-9:con*repeat,:) = data_temp;
    label(con*repeat-9:con*repeat,:) = conditon2label(con);
end

end