function output=TumorOnlyWrapper(paramFile,varargin)
%TumorOnlyWrapper - Entry function for tumor only caller
%uses a bayesian framework to call somatic variants and germline variants
%from tumor only exome sequencing
%
% Syntax:  output = TumorOnlyWrapper(paramFile)
%
% Inputs:
%    paramFile - parameter file in yaml format, see configTemplate.yaml 
%
% Outputs:
%    output - returns 0 upon completion
%
% Example: 
%   TumorOnlyWrapper('sampleConfig.yaml')
%
% Other m-files required: callCNA.m, callSNV.m, fitCNA.m, nllCNA.m, 
%   plotTumorOnly.m, preprocessTumorOnly.m, segmentData.m, 
%   writeCloneSummary.m  writeSegVCF.m  writeVCF.m
% Other requirements: parsePileupData.packed.pl, samtools, htslib
% Subfunctions: none
% MAT-files required: cghcbshybridnu.mat

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

inputParam=readInputs(paramFile)
cd(inputParam.workingDirectory);
addpath(genpath(inputParam.workingDirectory));

%%%% If outmat exits, uses outmat, otherwise generates data tables
if(exist([inputParam.outMat],'file'))
    vars={'T','E'};
    load([inputParam.outMat],'-mat',vars{:});
    inputParam=readInputs(paramFile)
else
   [T, E]=preprocessTumorOnly(inputParam,paramFile);
   save([inputParam.outMat]);
end

%%% Filters Exon Data and Segments
exonRD=E{E.MapQC<inputParam.minExonQual & E.perReadPass>inputParam.minPerReadPASS & E.abFrac>inputParam.minABFrac & ~isnan(E.TumorRD) & ~isnan(E.NormalRD),:};
segs=segmentData(exonRD,inputParam.cnaAlpha);

%%%Quality Filtering
[F,postTrust,postArtifact]=qualDiscrim(T,E,inputParam)
filtPos=postTrust(:,2)>inputParam.pGoodThresh & postArtifact(:,1)<inputParam.pGoodThresh;
['PASS positions: ' num2str(sum(filtPos))]

%%%Preliminary Variant Classification
hetPos=min([T.ApopAF T.BpopAF],[],2)>inputParam.minHetPopFreq & filtPos & T.BCountF+T.BCountR>=inputParam.minBCount;
somPos=T.CosmicCount>1 & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq & filtPos;
dataHet=[T.Chr(hetPos) T.Pos(hetPos) T.ControlRD(hetPos) T.ReadDepthPass(hetPos) T.BCountF(hetPos)+T.BCountR(hetPos)];
dataSom=[T.Chr(somPos) T.Pos(somPos) T.ControlRD(somPos) T.ReadDepthPass(somPos) T.BCountF(somPos)+T.BCountR(somPos)];
['Somatic positions: ' num2str(sum(somPos))]
['Het positions: ' num2str(sum(hetPos))]

%%%Make sure segments extend to ends of chromosome
for i=1:22
    idx1=find(segs(:,1)==i,1,'first');
    idx2=find(segs(:,1)==i,1,'last');
    segs(idx1,2)=min([T.Pos(T.Chr==i); exonRD(exonRD(:,1)==i,2)]);
    segs(idx2,3)=max([T.Pos(T.Chr==i); exonRD(exonRD(:,1)==i,3)]);
end

%%%Fit Copy Number Model
[segsTable, W, f, c, nll, pCNAexon]=fitCNA(dataHet,dataSom,exonRD,segs,inputParam);
['clonal fractions: ' num2str(f)]

%%%Add copy number info to data table
idxExon=getPosInRegionSplit([T.Chr T.Pos],exonRD(:,1:3),inputParam.blockSize);
pCNA(~isnan(idxExon),:)=pCNAexon(idxExon(~isnan(idxExon)));
pCNA(isnan(idxExon),:)=NaN;
idx=getPosInRegions([T.Chr T.Pos],segsTable(:,1:3));
T.NumCopies=segsTable(idx,5,1);
T.MinAlCopies=segsTable(idx,6,1);
T.cnaF=segsTable(idx,7,1);
T.W=segsTable(idx,8,1);
T.BmeanBQ(T.BCountF+T.BCountR==0)=inputParam.defaultBQ;

%%%Initial Bayesian Variant Calling
pDataSum(1)=0;
[pSomatic, pGermline, cloneId, pHom, pDataSum(2), pDataSomatic, pDataHet, pDataHom]=callSNV(T, W(1,:), f(1,:), inputParam);

%%%repeat fitting and variant calling until converges
i=1;
while(round(pDataSum(i+1))~=round(pDataSum(i)) && i<inputParam.maxIter)
    hetPos=pGermline(:,i)>inputParam.pGermlineThresh & filtPos;
    somPos=pSomatic(:,i)>inputParam.pSomaticThresh & filtPos & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq;
    dataHet=[T.Chr(hetPos) T.Pos(hetPos) T.ControlRD(hetPos) T.ReadDepthPass(hetPos) T.BCountF(hetPos)+T.BCountR(hetPos)];
    dataSom=[T.Chr(somPos) T.Pos(somPos) T.ControlRD(somPos) T.ReadDepthPass(somPos) T.BCountF(somPos)+T.BCountR(somPos)];
    save([inputParam.outMat]);
    ['Somatic positions: ' num2str(sum(somPos))]
    ['Het positions: ' num2str(sum(hetPos))]
    [segsTable(:,:,i+1), W(i+1,:), f(i+1,:), c(i+1,:), nll(i+1,:), pCNAexon]=fitCNA(dataHet,dataSom,exonRD,segs,inputParam);
    ['clonal fractions: ' num2str(f(i+1,:))]
    pCNA(~isnan(idxExon),i+1)=pCNAexon(idxExon(~isnan(idxExon)));
    pCNA(isnan(idxExon),i+1)=NaN;
    idx=getPosInRegions([T.Chr T.Pos],segsTable(:,1:3,i+1));
    T.NumCopies=segsTable(idx,5,i+1);
    T.MinAlCopies=segsTable(idx,6,i+1);
    T.cnaF=segsTable(idx,7,i+1);
    T.W=segsTable(idx,8,i+1);
    i=i+1
    [pSomatic(:,i), pGermline(:,i), cloneId(:,i), pHom(:,i), pDataSum(i+1), pDataSomatic(:,i), pDataHet(:,i), pDataHom(:,i)]=callSNV(T, W(i,:), f(i,:), inputParam);
    pDataSum
end
T.NumCopies=segsTable(idx,5,end);
T.MinAlCopies=segsTable(idx,6,end);
T.cnaF=segsTable(idx,7,end);
T.W=segsTable(idx,8,end);
['converged at iteration' num2str(i)]

%%%write output files
save([inputParam.outMat]);
writeVCF(T,pSomatic(:,end),postTrust(:,2),postArtifact(:,2),pGermline(:,end),pHom(:,end),cloneId(:,end),f(end,:),W(end,:),inputParam,pDataSomatic(:,end),pDataHet(:,end),pDataHom(:,end),pCNA(:,end),F);
writeSegVCF(segsTable(:,:,end),inputParam);
message='finished writing VCFs'
writeCloneSummary(segsTable(:,:,end),E,T,pSomatic(:,end),postTrust(:,2),postArtifact(:,2),f(end,:),W(end,:),cloneId(:,end),inputParam);
message='finished writing summary table'
plotTumorOnly(exonRD,segsTable(:,:,end),c(end),f(end,:),T,somPos,hetPos,cloneId(:,end),inputParam);

close all;
output=0;
exit;
