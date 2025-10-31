function [data,label] = spike_loader(array,timewindow,repeat)

% FR-V1-ac%ex, some neurons have few trials, 5 trials instead of 10
% neurons with fewer trials were removed here
% sequence: ac - ec -ex
%            1    2   3

[neuron_num,condition_num] = size(array);
% change this repeat extraction in later version, since repeats can vary in
% each condition

%create data&label
data = zeros(neuron_num,condition_num*repeat);




% Ten labels corresponding to ten images
labels = [1:condition_num*repeat];
label = fix((labels+(repeat-1))/repeat);




% iterate through conditons

for con = 1:condition_num
    con_temp = cell2mat(array(:,con)); %take a column of condition in the 4d array, dim:[neuronxrepeat,1500 ms]
    window = con_temp(:,timewindow(1):timewindow(2));
    spike_count = sum(window,2); %data with repeat scattered in this vector
    data_temp = reshape(spike_count,repeat,neuron_num); %dim [repeat, neuron]
    data_temp = data_temp'; %dim [neuron, repeat]
    data(:,con*repeat-(repeat-1):con*repeat) = data_temp;

end


end
