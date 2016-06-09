function [segsTable, W, f, c, nll, pCNA]=fitCNA(dataHet,dataSom,exonRD,segs,inputParam)
%fitCNA - uses EM to fit copy number parameters and estimate copy number
%
% Syntax: [segsTable, W, f, c, nll, pCNA]=fitCNA(dataHet,dataSom,exonRD,segs,inputParam)
%
% Inputs:
%   dataHet: data for germline heterozygous positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   dataSom: data for somatic positions with columns:
%       1-'Chr',2-'Pos',3-'ControlRD',4-'TumorRD',5-'Bcount'
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%   inputParam: structure with fields: minHetAF, numClones
%   
% Outputs:
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%   W: vector of lenght inputParm.numClones, controls width of allele
%       frequency distributions
%   f: vector of sample fraction of each clone
%   c: centering constant
%   nll: negative log likelihood
%   pCNA: vector of probability of copy number state by exon
%
% Other m-files required: nllCNA.m, callCNA.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%%initial centering on heats
diploidPos=dataHet(:,5)./dataHet(:,4)>inputParam.minHetAF;
cr=median(2*dataHet(diploidPos,3)./dataHet(diploidPos,4));

%%%optimize parameters
[param, nll] = fminsearchbnd(@(param)nllCNA(dataHet,dataSom,exonRD,segs,inputParam,param),[cr nanmedian(exonRD(:,4)*ones(1,inputParam.numClones)) 1/(inputParam.numClones+1):1/(inputParam.numClones+1):1-1/(inputParam.numClones+1)],[0.5*cr inputParam.minW*ones(1,inputParam.numClones) zeros(1,inputParam.numClones)],[2*cr inputParam.maxW*ones(1,inputParam.numClones) ones(1,inputParam.numClones)]);
c=param(1);
W=param(2:(length(param)-1)./2+1);
f=param((length(param)-1)./2+2:end);

%%%use optimized parameters to call copy number
[N, M, F, Wtable, log2FC, pCNA]=callCNA(dataHet,exonRD,segs,inputParam,param);
segsTable=[segs N M F Wtable log2FC];

return;
