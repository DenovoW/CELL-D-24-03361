
clearvars;
close all;

%things that must be changed:
% cd 'file path'
% title('A')
% results_directory, figures_directory,csvwrite
% baseline, drug_start/end

% things that can be changed:
% frame, f
% pks threshold, 0.05 as default

baseline = 720; % frame number of baseline
drug_start = baseline + 1; % the first frame of drug arrived
drug_end = baseline + 360; % 1.5 min after drug arrived

% % % % set results path
results_directory = '/share/home/huhailan/wangjunying/result/20231221/slice1_again';     
figures_directory = fullfile(results_directory,'NE_TTX');
if exist(figures_directory,'dir') ~= 7
    mkdir(figures_directory);
end
figures_visibility = 'off'; % either 'on' or 'off' (in any case figures are saved)

% % % % % create roi，128*128 in 2048*2048
% % no GUI to load data
cd '/share/home/huhailan/wangjunying/20231221/slice1/NE_5uM_TTX';
imgDir = dir('*.tif'); %read all the tif in direction
for i = 1 : numel(imgDir)
%     numsort(i) = str2double(imgDir(i).name(52:end-4))
    nam_str = strsplit(imgDir(i).name,{'_','.'})
    numsort(i) = str2num(nam_str{end-1})
end
[~,index] = sort(numsort)
newimgDir = imgDir(index)

init = [2048 2048]
final = [16 16]
for k = 1 : length(newimgDir)
    pic = imread([newimgDir(k).name])
    size(pic)
    c1 = zeros(final(1),1)
    c1 = c1+init(1)/final(1)
    c2 = zeros(final(2),1)
    c2 = c2+init(2)/final(2)
    c_pic = mat2cell(pic,c1,c2)
    rc_pic = cellfun(@(x) reshape(x, 1,init(1)*init(2)/final(1)/final(2)).', c_pic, 'UniformOutput', false)
    mc_pic = cellfun(@mean,rc_pic,'UniformOutput',false)
    final_pic = cell2mat(mc_pic)

    final_pic_double = im2double(final_pic)% uint16->double for calculate or it will overflow
    a_double(:,:,k) = final_pic_double
end

% % use GUI to load data
% cd '/Users/wangjunying/Desktop/tmp/tif_test';
% [filename, pathname] = uigetfile('*.tif','open file path','MultiSelect','on');
% [p q] = size(filename)
% for k = 1:q
%     pic = imread(strcat(pathname,filename{k}))
%     size(pic)
%     for i = 1:16
%         for j = 1:16
%             m_start = 1+(i-1)*128
%             m_end = i*128
%             n_start = 1+(j-1)*128
%             n_end = j*128
%             a(i,j,k) = median(pic(m_start:m_end,n_start:n_end),'all')
%         end
%     end
% end

[x,y,z] = size(a_double)
f = 1/4
t = (0:(z-1))*f

% % % % % matrix transpose
a_trans_double = permute(a_double,[1,3,2])%(x,y,z)->(x,z,y)

% uint16->double for calculate or it will overflow
% for k = 1:y
%     for j = 1:z
%         for i = 1:x
%             a_trans_double(i,j,k) = im2double(a_trans(i,j,k))
%         end
%     end
% end

% % % % deltaF/F0
F0 = squeeze(mean(a_trans_double,2))
for j = 1:z
%     for i = 1:x
%         F0(i,k) = median(a_trans_double(i,:,k),'all')
%         N = a_trans_double(i,:,k)
%         F0(i,k) = median(N(:))
        %         MAD(i,k) = mad(a_trans(i,:,k))
        %         thresh(i,k) = median(a_trans(i,:,k))+3*MAD(i,k) %not double,may
        %         have problem
%         for j = 1:z
% end
%     end 
      s_a = squeeze(a_trans_double(:,j,:))
      F(:,j,:) = (s_a-F0)./F0
end

% % % baseline, SD, noise = baseline + 2*SD.

for k = 1:y
    for i = 1:x
        M = F(i,1:baseline,k)
        baseline_F(i,k) = mean(M(:)) %R2017b matlab or below version
%         baseline_F(i,k) = mean(F(i,1:baseline,k),'all')

%         O = F(i,1:baseline,k)
%         baseline_F(i,k) = mean(O(:))
%         baseline_S(i,k) = std(F(i,1:baseline,k))
%         Event(i,k) = baseline_F(i,k) + 2*baseline_S(i,k)
    end
end

% % plot and findpeaks
% pks(1:x,1:y)=0

for k = 1:y
    for i = 1:x
%         MAD_F(i,k) = mad(F(i,:,k))
%         thesh_F(i,k) = median(F(i,:,k))+3*MAD_F(i,k)
        MAX_F(i,k) = max(F(i,drug_start:drug_end,k))
        [r,c] = find(F(i,drug_start:drug_end,k) == MAX_F(i,k),1) 
        %the index of the first minimum value in F(r_start:r_end)
        col = baseline + c
        T = (col-1)*f
        
%         if isempty(findpeaks(F(i,:,k),'MinPeakDistance',z-2,'MinPeakHeight', thesh_F(i,k))) == 0
%             pks(i,k) = findpeaks(F(i,:,k),'MinPeakDistance',z-2,'MinPeakHeight', thesh_F(i,k))
%         end     % MinPeakDistance = frame-2, so that it will only return one value
%         
        % % % find responding roi
        c(i,k) = MAX_F(i,k)-baseline_F(i,k)

        if c(i,k) >= 0.1
            response(i,k) = 1;
        else
            response(i,k) = 0
        end
%         
        % % % plot
        figure
        plot(t,F(i,:,k))
        if c(i,k) >= 0.1
            text(T,MAX_F(i,k),['(',num2str(T),',',num2str(MAX_F(i,k)),')'],'color','b')
        end
%         plot(t,a_trans(i,:,k))
        ylim([-0.5 0.7])
        hold on
        xlabel('Time (s)'); ylabel('Delta F/F');
        title('20231221_slice1_NE_TTX')
        savefig(fullfile(figures_directory,[num2str(i),'_',num2str(k)]))
        saveas(gcf,fullfile(figures_directory,[num2str(i),'_',num2str(k)]),'png')
    end
end

response_cell = sum(sum(response ~= 0))
csvwrite('/share/home/huhailan/wangjunying/result/20231221/slice1_again/NE_5uM_TTX_response_cell.cvs',response_cell)
filename = '/share/home/huhailan/wangjunying/result/20231221/slice1_again/NE_5uM_TTX_var.mat';
save(filename)
