function [F,pTrust,pArtifact]=qualDiscrim(T,E,inputParam)
%qualDiscrim - finds probability that positions are trusted or artifacts
%using quadratic discriminant analysis on a number of quality metrics
%
% Syntax:  [F,pTrust,pArtifact]=qualDiscrim(T,inputParam)
%
% Inputs:
%    T - table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%       'ControlRD','PosMapQC','perReadPass','abFrac'
%    E - table of exon datawith the following columns: 'Chr','StartPos','EndPos',
%       'TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'
%   inputParam - data structure with the following fields: ReadLength, 
%       minPerReadPASS, minABFrac, minPercentStrand, minMeanBQ, minMeanMQ, 
%       maxPMM, minSeqEndDist, maxStrandDiff, maxBQdiff, maxMQdiff, maxPMMdiff,
%       maxReadPosDiff, minPosQual, minExonQual, perPassReadReject,ABfracReject, 
%       perStrandReject, meanBQReject, meanMQReject, PMMReject,seqEndDistReject, 
%       strandDiffReject, BQdiffReject, MQdiffReject, PMMdiffReject, 
%       ReadPosDiffReject, inputParam.posQualReject, exonQualReject
%   
% Outputs:
%    F - table of quality statistics on positions same order as T
%    pTrust - posterior probability the call should be trusted
%    pArtifact - posterior probabilty the position has an artifact
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, preproccessTumorOnly

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%%set quality metrics of positions with no B allele to NaN
T.BmeanBQ(T.BCountF+T.BCountR==0)=NaN;
T.BmeanMQ(T.BCountF+T.BCountR==0)=NaN;
T.BmeanPMM(T.BCountF+T.BCountR==0)=NaN;
T.BmeanReadPos(T.BCountF+T.BCountR==0)=NaN;

%%%Create Quality Metrics Table
F=table();
F.TumorPerPassReads=T.ReadDepthPass./T.ReadDepth;
F.normalPerReadPass=T.perReadPass;
F.ABfrac=(T.ACountF+T.ACountR+T.BCountF+T.BCountR)./T.ReadDepthPass;
F.normalABfrac=max(T.abFrac,0);
F.minPerStrand=min([T.ACountF./(T.ACountF+T.ACountR) T.ACountR./(T.ACountF+T.ACountR) T.BCountF./(T.BCountF+T.BCountR) T.BCountR./(T.BCountF+T.BCountR)],[],2);
F.minBQ=min([T.AmeanBQ T.BmeanBQ],[],2);
F.minMQ=min([T.AmeanMQ T.BmeanMQ],[],2);
F.maxPMM=max([T.AmeanPMM T.BmeanPMM],[],2);
F.seqEndDist=min([T.AmeanReadPos inputParam.ReadLength-T.AmeanReadPos T.BmeanReadPos inputParam.ReadLength-T.BmeanReadPos],[],2);   
F.strandDiff=max(abs(T.ACountF./(T.ACountF+T.ACountR)-T.BCountF./(T.BCountF+T.BCountR)),0);
F.BQdiff=max(abs(T.AmeanBQ-T.BmeanBQ),0);
F.MQdiff=max(abs(T.AmeanMQ-T.BmeanMQ),0);
F.PMMdiff=max(abs(T.AmeanPMM-T.BmeanPMM),0);
F.ReadPosDiff=max(abs((T.AmeanReadPos-T.BmeanReadPos)),0);
F.posMapQC=-10*log10(T.PosMapQC+1E-6);
idx=getPosInRegions([T.Chr T.Pos],[E.Chr E.StartPos E.EndPos]);
F.exonMapQC=-10*log10(E.MapQC(idx)+1E-6);
F.Properties.VariableDescriptions={'Percent of Reads in Tumor Passing Quality Thresh', ...
    'Percent of Reads in Normals Passing Quality Thresh', ...
    'Fraction of QC Reads in Tumor Supporting A or B Allele', ...
    'Fraction of QC Reads in Normals Supporting A or B Allele', ...
    'Minimum Percentage of Reads from forward or reverse strand supporting A or B allele', ...
    'Minimum of average BQ for reads supporing A or B allele', ...
    'Minimum of average MQ for reads supporting A or B allele', ...
    'Maximum percentage of mismatches in reads supporting A or B allele', ...
    'Minimum average distance to end of sequence of reads supporting A or B allele', ...
    'Phred scaled Fisher test for strand bias', ...
    'Difference in average BQ between bases supporting A and B allele', ...
    'Difference in average MQ between reads supporting A and B allele', ...
    'Difference in average percentage of mismatches in reads supporting A and B allele', ...
    'Difference in average position in sequence between A and B allele', ...
    'Phred scaled postion quality score from Normals', ...
    'mean position quality score in exon'};


indelPos=T.A>4 | T.B>4;

%%%find positions that pass strict quality thresholds
goodPos(:,1)=F.TumorPerPassReads>inputParam.minPerReadPASS;
goodPos(:,2)=F.normalPerReadPass>inputParam.minPerReadPASS;
goodPos(:,3)=F.ABfrac>inputParam.minABFrac;
goodPos(:,4)=F.normalABfrac>inputParam.minABFrac;
goodPos(:,5)=F.minPerStrand>inputParam.minPercentStrand;
goodPos(:,6)=F.minBQ>inputParam.minMeanBQ;
goodPos(:,7)=F.minMQ>inputParam.minMeanMQ;
goodPos(:,8)=F.maxPMM<inputParam.maxPMM;
goodPos(:,9)=F.seqEndDist>inputParam.minSeqEndDist;
goodPos(:,10)=F.strandDiff<inputParam.maxStrandDiff;
goodPos(:,11)=F.BQdiff<inputParam.maxBQdiff;
goodPos(:,12)=F.MQdiff<inputParam.maxMQdiff;
goodPos(:,13)=F.PMMdiff<inputParam.maxPMMdiff;
goodPos(:,14)=F.ReadPosDiff<inputParam.maxReadPosDiff;
goodPos(:,15)=F.posMapQC>inputParam.minPosQual;
goodPos(:,16)=F.exonMapQC>inputParam.minExonQual;

%%%find positions that fail minimal quality thresholds
badPos(:,1)=F.TumorPerPassReads<inputParam.perPassReadReject;
badPos(:,2)=F.normalPerReadPass<inputParam.perPassReadReject;
badPos(:,3)=F.ABfrac<inputParam.ABfracReject;
badPos(:,4)=F.normalABfrac<inputParam.ABfracReject;
badPos(:,5)=F.minPerStrand<inputParam.perStrandReject;
badPos(:,6)=F.minBQ<inputParam.meanBQReject;
badPos(:,7)=F.minMQ<inputParam.meanMQReject;
badPos(:,8)=F.maxPMM>inputParam.PMMReject;
badPos(:,9)=F.seqEndDist<inputParam.seqEndDistReject;
badPos(:,10)=F.strandDiff>inputParam.strandDiffReject;
badPos(:,11)=F.BQdiff>inputParam.BQdiffReject;
badPos(:,12)=F.MQdiff>inputParam.MQdiffReject;
badPos(:,13)=F.PMMdiff>inputParam.PMMdiffReject;
badPos(:,14)=F.ReadPosDiff>inputParam.ReadPosDiffReject;
badPos(:,15)=F.posMapQC<inputParam.posQualReject;
badPos(:,16)=F.exonMapQC<inputParam.exonQualReject;

%%%classify SNVs as variants vs artifacts
sampleSNV=F{~indelPos,:};
trainingSNV=[F{sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos,:}; F{sum(badPos,2)>2 & ~indelPos,:}];
groupSNV=[ones(sum(sum(goodPos,2)>12 & sum(badPos,2)==0 & ~indelPos),1); zeros(sum(sum(badPos,2)>2 & ~indelPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[class(~indelPos),pArtifact(~indelPos,:)] = predict(discrSNV,sampleSNV);

%%%classify SNVs as trusted vs low quality positions
sampleSNV=F{~indelPos,:};
trainingSNV=[F{sum(goodPos,2)==16 & ~indelPos,:}; F{sum(badPos,2)>0 & ~indelPos,:}];
groupSNV=[ones(sum(sum(goodPos,2)==16 & ~indelPos),1); zeros(sum(sum(badPos,2)>0 & ~indelPos),1)];
discrSNV=fitcdiscr(trainingSNV,groupSNV,'DiscrimType', 'pseudoQuadratic');
[class(~indelPos),pTrust(~indelPos,:)] = predict(discrSNV,sampleSNV);

%%%classify indels as variants vs artifacts
finitePos=sum(isfinite(F{:,:}),2)==16;
sampleIndel=F{indelPos & finitePos,:};
trainingIndel=[F{sum(goodPos(:,[1:10 12:16]),2)==15 & indelPos & finitePos,:}; F{sum(badPos(:,[1:10 12:16]),2)>0 & indelPos & finitePos,:}];
groupIndel=[ones(sum(sum(goodPos(:,[1:10 12:16]),2)==15 & indelPos & finitePos),1); zeros(sum(sum(badPos(:,[1:10 12:16]),2)>0 & indelPos & finitePos),1)];
discrIndel=fitcdiscr(trainingIndel(:,[1:10 12:16]),groupIndel,'DiscrimType', 'pseudoQuadratic');
[class(indelPos & finitePos),pTrust(indelPos & finitePos,:)] = predict(discrIndel,sampleIndel(:,[1:10 12:16]));

%%%classify indels as trusted vs low quality positions
finitePos=sum(isfinite(F{:,:}),2)==16;
sampleIndel=F{indelPos & finitePos,:};
trainingIndel=[F{sum(goodPos(:,[1:10 12:16]),2)>11 & sum(badPos,2)==0 &indelPos & finitePos,:}; F{sum(badPos(:,[1:10 12:16]),2)>2 & indelPos & finitePos,:}];
groupIndel=[ones(sum(sum(goodPos(:,[1:10 12:16]),2)>11 & sum(badPos,2)==0 &indelPos & finitePos),1); zeros(sum(sum(badPos(:,[1:10 12:16]),2)>2 & indelPos & finitePos),1)];
discrIndel=fitcdiscr(trainingIndel(:,[1:10 12:16]),groupIndel,'DiscrimType', 'pseudoQuadratic');
[class(indelPos & finitePos),pArtifact(indelPos & finitePos,:)] = predict(discrIndel,sampleIndel(:,[1:10 12:16]));

return;