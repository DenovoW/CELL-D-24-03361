'
/clear all
close all
% filename, start, t(the end time for analysis)
% baseline =
% mean(d_sort(1:100,:))这里的100要看具体的文件更改，短文件改短一点,hm3d一般就取到100，copn可能取到20左右吧
% check plot of cell type first
% then change the celltype manually, delete the '%' before csvwrite
% delete the '%' before savefig,saveas

filename = 'I:\2p imaging\20220524_nG6s_aHm3d\analyse\slice2 CNO 5uM'
cd (filename);
results_directory= filename; 
figures_directory=fullfile(results_directory,'nonresponsive_CellType');
if exist(figures_directory,'dir')~=7
    mkdir(figures_directory);
end
figures_visibility='on'; % either 'on' or 'off' (in any case figures are saved)
data = csvread('Results_analyse.csv');
roi_response = csvread('roi_SerialNum_activated.csv');
[m,n] = size(data);
roi_total = 1:n;
roi = setdiff(roi_total,roi_response);

[m,n] = size(roi);
t = 250;
start = 0;
d = data(1:t,roi);%取roi的那几列
d_sort = sort(d);
% baseline = mean(d)
baseline = mean(d_sort(1:100,:));
ma = max(d);
sd = std(d);
value = [sd;baseline;(baseline + sd);ma;roi];
    
for i = 1:n
    if value(1,i) < 30 
        value(6,i) = 0; %silent
    else value(6,i) = 1;
    end
end

pos2 = find(value(6,:) == 0);
silent = value(5,pos2)
non_silent = setdiff(roi,silent);
[~, idx] = ismember(non_silent, value(5,:));
[m,n] = size(idx);
for i = 1:n
    [pks,locs] = findpeaks(d(:,idx(i)),'MinPeakHeight',value(3,idx(i)),'MinPeakDistance',10);
    findpeaks(d(:,idx(i)),'MinPeakHeight',value(3,idx(i)),'MinPeakDistance',10);
    c{1,i} = pks;
    c{2,i} = [start;locs;t];
    c{3,i} = diff(c{2,i});
    c{4,i} = std(c{3,i});
    c{5,i} = non_silent(i); %c = {[pks];[start;locs;t];intervals;sd of intervals;serial num}
    if c{4,i} < 10
        c{6,i} = 2;
    else c{6,i} = 1; %c = {[pks];[start;locs;t];intervals;sd of intervals;serial num；1/2}
    end
    hold on
    plot(1:t,repmat(value(3,idx(i)),1,t))
    xlabel(num2str(c{5,i}))
    savefig(fullfile(figures_directory,num2str(c{5,i})))
    saveas(gcf,fullfile(figures_directory,num2str(c{5,i})),'png')
    figure
end

 % 获取第6行值为 2 的列的索引
sixth_row = cell2mat(c(6,:)); % 将第6行转换为矩阵
indices = find(sixth_row == 2);

% 输出对应第五行的值
regular = cell2mat(c(5,indices))
irregular = setdiff(non_silent,regular)
indices2 = ismember(value(5,:), regular);
value(6, indices2) = 2; %silent 0, irregular 1, regular 2

csvwrite('Serial_num_silent.csv',silent)
csvwrite('Serial_num_regular.csv',regular)
csvwrite('Serial_num_irregular.csv',irregular)
csvwrite('CellType.csv',[value(5,:)',value(6,:)'])

