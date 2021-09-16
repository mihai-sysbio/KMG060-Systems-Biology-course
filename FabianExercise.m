clear
clc
Traw    = readtable('Exercise\Exercise1_RNAseq_Analysis\data\Saccharomyces_RNAseq_RAW_counts.csv');
Factors = Traw(:,2:end).Properties.VariableNames;
ref     = startsWith(Factors,'Control');
Htemp   = startsWith(Factors,'HighT');
LowpH   = startsWith(Factors,'LowPH');
Osmo    = startsWith(Factors,'Osmo');
Anae    = startsWith(Factors,'Anae');

numOnly    = Traw{:,2:end};
Referance  = geomean(numOnly,2);
NaN=Referance > 0;
ratios     = numOnly(NaN,:)./Referance(NaN);
Size       = median(ratios,1);
normCounts = numOnly./Size;

%% Boxplots log scale

lognum = log2(numOnly);
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

[pc,zscores,pcvars] = pca(normCounts');
colorScheme = colormap(turbo(5));
hold on
condStr = {};
for i=1:5
    nameTMP = strsplit(Traw.Properties.VariableNames{(i-1)*3+2},'_');
    PC1     = zscores((i-1)*3+1:i*3,1);
    PC2     = zscores((i-1)*3+1:i*3,2);
    condStr = [condStr;nameTMP{1}];
    scatter(PC1,PC2,50,colorScheme(i,:),'fill');
end
percPC1 = round((abs(sum(PC1))./(abs(sum(PC1))+abs(sum(PC2)))).*100,2);
percPC2 = round((abs(sum(PC2))./(abs(sum(PC1))+abs(sum(PC2)))).*100,2);
xlabel(['First Principal Component ' ,'(',num2str(percPC1) , '%)']);
ylabel(['Second Principal Component ' ,'(',num2str(percPC2) , '%)']);
title('Principal Component Scatter Plot');
legend(condStr)
hold off

%% 3. Differentia gene expression
normCountsStress  = normCounts(:,LowpH);
normCountsref     = normCounts(:,ref);
tLocal = nbintest(normCountsref,normCountsStress,'VarianceLink','LocalRegression');
padj = mafdr(tLocal.pValue,'BHFDR',true);
meanRef    = mean(normCountsref,2);
meanStress = mean(normCountsStress,2);
log2FC     = log2(meanStress./meanRef);
geneTable  = table(meanRef,meanStress,log2FC,tLocal.pValue,padj);
geneTable.Properties.RowNames      = Traw.GeneName;
geneTable.Properties.VariableNames = {'Mean_Ref','Mean_Stress','Log2_FC','pVal','adjPVal'};

x=geneTable.Log2_FC; y=geneTable.adjPVal;
FClim=abs(x)>=2; adjlim=y<=0.01;y=-log10(geneTable.adjPVal);
hold on
scatter(x,y,30,'fill','black')

% For loop to only include values on positions in the matrix where both
% conditions abs(Log2FC)=>2 and adj_pVal<=0.01 are true
for i=1:length(x)
    if adjlim(i)==1 && FClim(i)==1 %BOTH conditions means keep value
        x(i)=x(i);
        y(i)=y(i);
    else                           %Otherwise, put 0 
        x(i)=0;
        y(i)=0;
    end
end

scatter(x,y,30,'fill','red')
line([2 2],[0 60]); line([-2 -2],[0 60]); line([-6 6],-log10([0.01 0.01]));
title('Volcano plot')
xlabel('Log2 FC')
ylabel('-log10(adjPValue)')
hold off

%% 4. 
Desc=readtable('GeneDescriptions.csv'); %Get data as table
Desc=table(Desc{:,1}, Desc{:,3}); %Remove abundant second column, we only want gene names (outerjoin) and description columns
q=geneTable.Properties.RowNames; geneTable{:,6}=q; %Adds row names to 6th column for outerjoin
geneTable = outerjoin(Desc,geneTable,'LeftKeys','Var1','RightKeys','Var6'); %Adds Gene Descriptions to Genetable
geneTable = sortrows(geneTable,'adjPVal','ascend');
geneTable.Var6=[]; %Removes abundand column (which was only used for the purpose of outerjoin)
geneTable.Properties.VariableNames([1 2])={'Gene','Description'}; %Adds headers to Gene names and Description columns
disp(geneTable(1:10,:))
