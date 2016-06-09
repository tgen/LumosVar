function [pSomatic, pGermline, cloneId, pHom, pDataSum, pDataSomatic, pDataHet, pDataHom, expectedAF]=callSNV(T, W, f, inputParam)
%callSNV - finds probability variants are somatic or germline
%
% Syntax: idx=getPosInRegions(pos,regions)
%
% Inputs:
%   T: table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%       'ControlRD','PosMapQC','perReadPass','abFrac'
%   W: vector of length inputParm.numClones, controls width of allele
%       frequency distributions
%   f: vector of sample fraction of each clone
%   inputParam: structure with fields: priorSomaticSNV, priorSomaticIndel
%
% Outputs:
%   pSomatic: posterior probability somatic variant
%   pGermline: posterior probability germline heterozygous
%   cloneId: most likley clone to have variant if somatic
%   pHom: posterior probability germline homozygous AA
%   pDataSum: sum of marginal
%   pDataSomatic: likliehood of somatic variant
%   pDataHet: likliehood germline heterozygous
%   pDataHom: likliehood of germline homozygous AA
%   expectedAF: expected allele frequency for somatic variant
%
% Other m-files required: none
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

%%% find priors
indelPos=T.A>4 | T.B>4 | T.Ref >4;
priorSomatic(indelPos,:)=(T.CosmicCount(indelPos)+1).*inputParam.priorSomaticIndel;
priorSomatic(~indelPos,:)=(T.CosmicCount(~indelPos)+1).*inputParam.priorSomaticSNV;
priorHet=2*T.ApopAF.*T.BpopAF.*(1-priorSomatic);
priorHom=(T.ApopAF.^2).*(1-priorSomatic);
priorOther=1-priorHet-priorHom-priorSomatic;

%%% calculate expected heterozygous allele frequencies
cnCorr=T.cnaF.*T.MinAlCopies./T.NumCopies+(1-T.cnaF)*0.5;
cnCorr(T.NumCopies==0)=0.5;

%%% calculate likliehoods of germline genotypes
pDataHom=bbinopdf_ln(T.ACountF+T.ACountR,T.ReadDepthPass,T.W.*(1-10.^(-T.BmeanBQ./10)),T.W.*10.^(-T.BmeanBQ./10));
pDataHet=bbinopdf_ln(T.BCountF+T.BCountR,T.ReadDepthPass,cnCorr.*T.W,(1-cnCorr).*T.W);
pDataOther=bbinopdf_ln(max(T.ReadDepthPass-T.ACountF-T.ACountR-T.BCountF-T.BCountR,0),T.ReadDepthPass,T.W.*(1-10.^(-T.AmeanBQ./10)),T.W.*10.^(-T.AmeanBQ./10));

%%% calculate expected somatic allele frequencies
for i=1:length(f)
    matchIdx=round(100.*T.cnaF)==round(100.*f(i));
    expAF(matchIdx,i)=f(i)*(T.NumCopies(matchIdx)-T.MinAlCopies(matchIdx))./(f(i)*T.NumCopies(matchIdx)+(1-f(i))*2);
    if sum(~matchIdx)>0
        expAF(~matchIdx,i)=f(i)./(T.cnaF(~matchIdx).*T.NumCopies(~matchIdx)+(1-T.cnaF(~matchIdx))*2);
        expAF(T.NumCopies==0 & ~matchIdx,i)=min([(1-T.cnaF); f(i)])./2;
    end
end

%%% find likelihood of somatic mutations 
for i=1:length(f)
    alpha(:,i)=min(expAF(:,i),1-expAF(:,i))*W(i);
    beta(:,i)=max(expAF(:,i),1-expAF(:,i))*W(i);
    cloneLik(:,i)=bbinopdf_ln(T.BCountF+T.BCountR,T.ReadDepthPass,alpha(:,i),beta(:,i));
end
[pDataSomatic,cloneId]=max(cloneLik,[],2);

%%% find expected somatic AF for most likely clone
for i=1:length(f)
    expectedAF(cloneId==i,:)=expAF(cloneId==i,i);
end

%%% calculate posteriors
pData=pDataHom.*priorHom+pDataHet.*priorHet+pDataSomatic.*priorSomatic+priorOther.*pDataOther;
pSomatic=pDataSomatic.*priorSomatic./pData;
pGermline=pDataHet.*priorHet./pData;
pHom=pDataHom.*priorHom./pData;
pOther=pDataOther.*priorOther./pData;

pDataSum=sum(pData);
return;
