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
normCountsStress  = normCounts(:,Htemp);
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
%scatter(x,y,30,'fill','black')

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

%scatter(x,y,30,'fill','red')
%line([2 2],[0 60]); line([-2 -2],[0 60]); line([-6 6],-log10([0.01 0.01]));
title('Volcano plot')
xlabel('Log2 FoldChange')
ylabel('-log10(adjPValue)')
hold off

%% 4. First analysis of DE genes
Table4 = geneTable; %In order to leave geneTable untouched
Desc=readtable('GeneDescriptions.csv'); %Get data as table
Desc=table(Desc{:,1}, Desc{:,3}); %Remove abundant second column, we only want gene names (for outerjoin) and description columns
q=Table4.Properties.RowNames; Table4{:,6}=q; %Adds row names to 6th column for outerjoin
Table4 = outerjoin(Desc,Table4,'LeftKeys','Var1','RightKeys','Var6'); %Adds Gene Descriptions to Genetable
Table4 = sortrows(Table4,'adjPVal','ascend');
Table4.Var6 = []; %Removes abundand column (which was only used for the purpose of outerjoin)
Table4.Properties.VariableNames([1 2]) = {'Gene','Description'}; %Adds headers to Gene names and Description columns
disp(Table4(1:10,:))

%% 5. Find Associated GO Terms
GoTerms_table = readtable('../data/GoTermsMapping.txt','delimiter','\t');
GoTermsIDs    = unique(GoTerms_table.GoTerm);

GoTermsGeneMap = containers.Map();
for i = 1:height(GoTerms_table)
    key = GoTerms_table{i,'GeneName'}{1};
    if isKey(GoTermsGeneMap,key)
        GoTermsGeneMap(key) = [GoTermsGeneMap(key),GoTerms_table{i,'GoTerm'}{1}];
    else
        GoTermsGeneMap(key) = {GoTerms_table{i,'GoTerm'}{1}};
    end
end
fprintf('Number of annotated genes related to functional process is %d.\n',GoTermsGeneMap.Count)
fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique(GoTerms_table.GoTerm)))
fprintf('Number of gene-GO term associations is %d.\n',numel(GoTerms_table))

selectedGene      = Table4{1,1}{1,1};
associatedGoTerms = GoTermsGeneMap(selectedGene);
disp(['Associated GO Terms for gene ',selectedGene])
disp(associatedGoTerms)

GO = geneont('File','../data/GoTerms.obo');

for i=1:length(associatedGoTerms)
GoTermID         = associatedGoTerms{i}(4:end);
GoTermName       = GO(str2double(GoTermID)).Terms.Name;
GoTermDefinition = GO(str2double(GoTermID)).Terms.Definition;

disp(['GO Term ',(GoTermID),' - ',GoTermName,' : ',GoTermDefinition])
end

%% Step 6: Find enriched GO Terms in differentially expressed genes

%split gene list into DE and non DE genes
indexDE       = geneTable.adjPVal<=0.01 & abs(geneTable.Log2_FC)>=2;
geneListDE    = geneTable.Properties.RowNames(indexDE);
geneListNonDE = geneTable.Properties.RowNames(~indexDE);
%Select GO Term, in this case the top1 is chosen
rankedPosition = 1;
GoTermID       = str2double(associatedGoTerms{rankedPosition}(4:end));
%GO terms full IDs are of the form GO:XXXXXXX
full_GoTerm_Id = sprintf('GO:%07d',GoTermID);

%################## 'Alternative #1 
%Check enrichment for selected GO Terms
%Count number of occurences for the GO Term in all genes and in all DE
%genes
GoTermCountAll = 0;
GoTermCountDE  = 0;
for key = GoTermsGeneMap.keys
    associatedGoTerms = (GoTermsGeneMap(key{1}));
    %Search the GO term ID in the associated GO terms cell array
    if any(strcmp(associatedGoTerms,full_GoTerm_Id))
        %
        if any(strcmp(geneListDE,key{1}))
            GoTermCountDE = GoTermCountDE+1;
        end
        GoTermCountAll = GoTermCountAll+1;
    end
end
%################## Alternative 2
GoTermsIDs_table  = GoTerms_table.GoTerm;
GoTermCountAll = 0;
GoTermCountDE  = 0;
presence       = find(strcmp(GoTermsIDs_table,full_GoTerm_Id));
if ~isempty(presence)
    relatedGenes  = GoTerms_table.GeneName(presence);
    [~,indexesDE] = intersect(geneListDE,relatedGenes);
    GoTermCountDE  = length(indexesDE);
    GoTermCountAll = length(presence);
end
%##################

disp(['GO Term ',num2str(GoTermID),' is in ',num2str(GoTermCountDE),...
    ' out of ',num2str(numel(geneListDE)),' DE Genes and in total there are ',...
    num2str(GoTermCountAll),' occurences in all ',num2str(numel(GoTermsGeneMap.keys)),' Genes'])

%run hypergeometric test to assess significance
pHyperGeo=hygepdf(GoTermCountDE,numel(GoTermsGeneMap.keys),GoTermCountAll,numel(geneListDE));
disp(['The probablity for this is ',num2str(pHyperGeo)])
%##################
% TASK: 
% Create a loop that goes over all existing GoTerms (GoTermsIDs) and
% calculates the enrichment. Don't forget to use a correction for multiple
% testing, e.g. the Benjamini Hochberg Correction as for the DE detection
% (mafdr).
% Sort GO Terms by adjusted pValue and print out every one that has a
% pValue <=0.01
%##################
pHyperGeos= zeros(length(GoTermsIDs),1);
for i= 1:length(GoTermsIDs)
    GoTermCountAll = 0;
    GoTermCountDE  = 0;
    presence       = find(strcmp(GoTermsIDs_table,GoTermsIDs{i}));
if ~isempty(presence)
    relatedGenes  = GoTerms_table.GeneName(presence);
    [~,indexesDE] = intersect(geneListDE,relatedGenes);
    GoTermCountDE  = length(indexesDE);
    GoTermCountAll = length(presence);
   
end

pHyperGeos(i)=hygepdf(GoTermCountDE,numel(GoTermsGeneMap.keys),GoTermCountAll,numel(geneListDE));
end

pHyperGeosAdj = mafdr(pHyperGeos,'BHFDR',true);
GoTermsTable = table(GoTermsIDs, pHyperGeosAdj);

GoTermsTable.Properties.VariableNames = {'GoTerm', 'adjPVal'};
GoTermsTable = sortrows(GoTermsTable,'adjPVal','ascend');
LimPV = GoTermsTable.adjPVal <= 0.01;
signGO = GoTermsTable(:,LimPV);

for i = 1:size(signGO)
GoTerm = GoTermsTable.GoTerm{i};
GoTermID = str2double(GoTerm(4:end));
GoTermName = GO(GoTermID).Terms.Name;
disp([GoTermsTable.GoTerm{i}, ' - ', GoTermName, ': ', 'is significantly enriched with p_{adj} = ', num2str(GoTermsTable.adjPVal(i))])
end