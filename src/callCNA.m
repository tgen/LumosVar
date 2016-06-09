function [N, M, F, Wout, log2FC, pCNA] = callCNA(dataHet,exonRD,segs,inputParam,param)
%callCNA - determine most likley copy number state for each segment
%
% Syntax: [N, M, F, Wout, log2FC, pCNA] = callCNA(dataHet,exonRD,segs,inputParam,param)
%
% Inputs:
%   dataHet: data for germline heterozygous positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: cnaPrior, minLik
%   param: vector of paramaters of length 2*numClones+1
%       param(1) - copy number scaling constant
%       param(2:numClones+1)=W (controls width of allele frequency dist)
%       param(numClones+1:2*numClones+1)=f (sample fractions)
%   
% Outputs:
%   N: total copies by segment
%   M: minor allele copies by segment
%   Wout: W parameter by segment
%   log2FC: log2 fold change of tumor/control
%   pCNA: probability of copy number event by exon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, fitCNA, nllCNA

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%% read inputs
CNAscale=param(1);
W=param(2:(length(param)-1)./2+1);
f=param((length(param)-1)./2+2:end);
colHeaders={'Chr' 'Pos' 'ExpReadCount' 'TotalReadCount' 'MinorReadCount'};
D=array2table(dataHet,'VariableNames',colHeaders);
E=array2table(exonRD,'VariableNames',{'Chr','StartPos','EndPos','TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'});

%%% find means accross segments
meanTumorRDexon=getMeanInRegions([E.Chr E.StartPos],E.TumorRD,segs);
meanNormalRDexon=getMeanInRegions([E.Chr E.StartPos],E.NormalRD,segs);
meanTumorRD=getMeanInRegions([D.Chr D.Pos],D.TotalReadCount,segs);
meanMinorRD=getMeanInRegions([D.Chr D.Pos],D.MinorReadCount,segs);

%%% calculate log2FC
log2FC=log2((CNAscale./2).*meanTumorRDexon./(meanNormalRDexon));

%%% find copy number per segment and clone
for i=1:length(f)
    Nseg(:,i)=max(round((CNAscale*(meanTumorRDexon./(meanNormalRDexon))-2*(1-f(i)))/f(i)),0);
    Mseg(:,i)=round((Nseg(:,i)/f(i)).*((meanMinorRD./(meanTumorRD))-0.5*(1-f(i))));
end

%%% lookup copy number per position or exon
idx=getPosInRegions([D.Chr D.Pos], segs);
Mseg(isnan(Mseg))=min(Nseg(isnan(Mseg)),1);
Mseg(Mseg<0)=0;
Nmat=Nseg(idx,:);
Mmat=Mseg(idx,:);
idxExon=getPosInRegions([E.Chr E.StartPos],segs);
NmatExon=Nseg(idxExon,:);

%%% find liklihoods of read counts and depth
for i=1:length(f)
    corr(:,i)=f(i).*Mmat(:,i)./Nmat(:,i)+(1-f(i))*0.5;
    corr(Nmat(:,i)==0,i)=0.5;
    hetlik(:,i)=bbinopdf_ln(D.MinorReadCount,D.TotalReadCount,W(i)*corr(:,i),W(i)*(1-corr(:,i)))+inputParam.minLik;
    expReadCount(:,i)=f(i)*E.NormalRD.*NmatExon(:,i)./CNAscale+(1-f(i))*E.NormalRD*2./CNAscale;
    depthlik(:,i)=poisspdf(round(E.TumorRD),round(expReadCount(:,i)))+inputParam.minLik;
    segLik(:,i)=getMeanInRegions([D.Chr D.Pos],log(hetlik(:,i)),segs)+getMeanInRegions([E.Chr E.StartPos],log(depthlik(:,i)),segs);
end

%%% find which clone contains most likely CNV per segment
[m,cnaIdx]=max(segLik,[],2);
for i=1:length(f)
    N(cnaIdx==i,:)=Nseg(cnaIdx==i,i);
    M(cnaIdx==i,:)=Mseg(cnaIdx==i,i);
end
F(:,1)=f(cnaIdx);
Wout(:,1)=W(cnaIdx);

%%% find probability of copy number per exon
for n=1:length(inputParam.cnaPrior)-1
    expReadCountN(:,n)=F(idxExon).*E.NormalRD.*(n-1)./CNAscale+(1-F(idxExon)).*E.NormalRD*2./CNAscale;
    pDataCNA(:,n)=inputParam.cnaPrior(n)*(poisspdf(round(E.TumorRD),round(expReadCountN(:,n))));
end
nMax=(n)*ones(size(N));
nMax(N>n)=N(N>n);
expReadCountN(:,n+1)=F(idxExon).*E.NormalRD.*nMax(idxExon)./CNAscale+(1-F(idxExon)).*E.NormalRD*2./CNAscale;
pDataCNA(:,n+1)=inputParam.cnaPrior(n+1)*(poisspdf(round(E.TumorRD),round(expReadCountN(:,n+1))));
for n=1:length(inputParam.cnaPrior)-1
    pCNA((n-1)==N(idxExon),:)=pDataCNA((n-1)==N(idxExon),n)./sum(pDataCNA((n-1)==N(idxExon),:),2);
end
pCNA(n<=N(idxExon),:)=pDataCNA(n<=N(idxExon),n)./sum(pDataCNA(n<=N(idxExon),:),2);

return