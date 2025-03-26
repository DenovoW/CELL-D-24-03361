clear all
close all

stim_time = [123,124] %frame名称后缀的序号,0..end e.g. 139,140,278,279,402,403,617,618
stim_time = stim_time + 1
event_time = [123] %e.g. 139,278,401,617
roi = [7,10,11,14,18,42,45,46,47,57,58,62,63,64,66,69]
cd 'I:\2p imaging\20241119_nG6s_aCOPN5\analyse\slice2\DPCPX_light2_0001'
data = csvread('I:\2p imaging\20241119_nG6s_aCOPN5\analyse\slice2\DPCPX_light2_0001\Results_analyse.csv')

% [m,n] = size(data);
a = zeros(length(data),1);
a(event_time,1) = 1;
d = data(:,roi); %取roi的那几列
a(stim_time,:) = []; d(stim_time,:) = []; %创建打标的前三列a；a和原始数据均删除光干扰frame

one = find(a(:,1)==1);
one_base_10 = one - 10; one_base_5 = one - 5; %ctrl time is -10 ~ -5
one_after_15 = one + 15; %test time is 0-15

m = size(d,2);
n = size(one,1);
for j = 1:m
    for i = 1:n
        F_ctrl(j,i) = mean(d(one_base_10(i):one_base_5(i),j));
        evt_time = one(i):one_after_15(i);
        deltaF(i,:,j) = (d(evt_time,j) - F_ctrl(j,i))';
    end
    dF(:,:,j) = mean(deltaF(:,:,j),1);
    max_df(j,:) = max(dF(:,:,j));
end

d = [a,d];
result = [roi',max_df];
csvwrite('result_deltaF.csv',d)
csvwrite('roi_maxDeltaF.csv',result)