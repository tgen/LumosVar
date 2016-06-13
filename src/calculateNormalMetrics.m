
function [normalMetrics]=calculateNormalMetrics(NormalData, priorMapError)
%calculateNormalMetrics - calculates mean read depths and 
%position quality scores for a control sample
%calls parsePileupData.packed.pl to parse samtools output
%
% [normalMetrics]=calculateNormalMetrics(NormalData, priorMapError)
%
% Inputs:
%   NormalData - matrix of normal data with the following columns:
%       1-'Chr', 2-'Pos', 3-'ReadDepth', 4-'ReadDepthPass', 5-'Ref', 6-'A',
%       7-'ACountF', 8-'ACountR', 9-'AmeanBQ', 10-'AmeanMQ', 11-'B', 
%       12-'BCountF', 13-'BCountR', 14-'BmeanBQ', 15-'BmeanMQ',
%       16-'ApopAF', 17-'BpopAF'
%   priorMapError - prior probability that the position has poor quality
%   
% Outputs:
%   normalMetrics - matrix with the following columns: 1-'Chr', 2-'Pos', 
%       3-'ReadDepthPass', 4-'pMapError', 5-'perReadPass', 6-'abFrac'
%
% Other m-files required: none
% Other requirements: none
% Subfunctions: none
% MAT-files required: none
%
% See also: printNormalMetrics

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%%% read data into tables
NormalColHeaders={'Chr','Pos','ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ','AmeanMQ','B','BCountF','BCountR','BmeanBQ','BmeanMQ','ApopAF', 'BpopAF'};
N=array2table(NormalData,'VariableNames',NormalColHeaders);
clear NormalColHeaders NormalData;

%%% calculate priors
priorHet=2.*N.ApopAF.*N.BpopAF;
priorHom=N.ApopAF.^2;

%%% calculate liklihoods
pDataMapError=10.^(-min([N.AmeanMQ N.BmeanMQ],[],2)./10);
pDataHom=binopdf(N.BCountF+N.BCountR,N.ReadDepthPass,10.^(-N.BmeanBQ./10));
pDataHet=binopdf(N.BCountF+N.BCountR,N.ReadDepthPass,0.5);

%%% calculate marginal
pData=priorMapError.*pDataMapError+priorHet.*pDataHet+priorHom.*pDataHom;

%%% calculate posteriors
pMapError=priorMapError.*pDataMapError./pData;
pHet=priorHet.*pDataHet./pData;
pHom=priorHom.*pDataHom./pData;

%%% calculate other quality metrics
perReadPass=N.ReadDepthPass./N.ReadDepth;
abFrac=(N.ACountF+N.ACountR+N.BCountF+N.BCountR)./N.ReadDepthPass;

normalMetrics=[N.Chr N.Pos N.ReadDepthPass pMapError perReadPass abFrac];

return;


