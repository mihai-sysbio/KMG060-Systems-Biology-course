clear
clc

% Opening from my precise directory. Change if on other machine
T1 = readtable('E:\Users\fabbe\Documents\MATLAB\KMG060\Exercise\Exercise1_RNAseq_Analysis\data\Saccharomyces_RNAseq_RAW_counts.csv');
T2 = standardizeMissing(T1, 0); % The 0 is the value to be replaced with NaN
T3 = rmmissing(T2, 1); % The 2 is for deleting rows; use 1 to delete columns

% 1a) 7126 represented genes 1b) 15 experimental conditions 1c) 5950 non-zero reads

num = T3{:,2:end}; %Gets matrix without gene name column (from non-zero matrix)
mean = geomean(num,2); %Gets one column matrix containing means of each row
norm = (num)./mean; %Matrix with all normalized non-zero reads

Size = median(norm,1);
normCounts = norm./Size;

%% Boxplots
subplot(1,2,1);
boxplot(num)
title('Raw counts')
ylabel('Read counts')
xlabel('Samples')

subplot(1,2,2);
boxplot(normCounts)
title('Normalized counts')
ylabel('Read counts')
xlabel('Samples')

%% Boxplots log scale

lognum = log2(num);
lognormCounts = log2(normCounts);

subplot(1,2,1);
boxplot(lognum)
title('Raw counts')
ylabel('Read counts')
xlabel('Samples')

subplot(1,2,2);
boxplot(lognormCounts)
title('Normalized counts')
ylabel('Read counts')
xlabel('Samples')

%% 2. PCA Analysis

[pc, zscores, pcvars] = pca(normCounts');
colorScheme = colormap(jet(5));
hold on
condStr = {};
for i=1:5
    nameTMP = strsplit(T1.Properties.VariableNames{(i-1)*3+2},'_');
    PC1     = zscores((i-1)*3+1:i*3,1);
    PC2     = zscores((i-1)*3+1:i*3,2);
    condStr = [condStr;nameTMP{1}];
    scatter(PC1,PC2,50,colorScheme(i,:),'fill');
end
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot');
legend(condStr)
hold off
