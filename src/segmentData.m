function segs=segmentData(exonRD,cnaAlpha)
%segmentData - parellizes segmentation by chromosome
%uses circular binary segmentation from bioinformatics toolbox
%
% Syntax: segs=segmentData(exonRD,cnaAlpha)
%
% Inputs:
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   cnaAlpha: significance cutoff for segmentation
%   
% Outputs:
%   segs: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio'
%
% Other m-files required: none
% Other requirements: bioinformatics toolbox
% Subfunctions: none
% MAT-files required: cghcbshybridnu.mat
%
% See also: TumorOnlyWrapper, cghcbs

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

parfor i=1:22
    clustsegs{i}=cghcbs([exonRD(exonRD(:,1)==i,1) exonRD(exonRD(:,1)==i,2) log((exonRD(exonRD(:,1)==i,4)+1)./(exonRD(exonRD(:,1)==i,5)+1))],'Alpha',cnaAlpha);
end
    
segs=[];
for i=1:22
    segs=[segs; [clustsegs{i}.SegmentData(1).Chromosome*ones(size(clustsegs{i}.SegmentData(1).Start)) clustsegs{i}.SegmentData(1).Start clustsegs{i}.SegmentData(1).End clustsegs{i}.SegmentData(1).Mean]];
end

delete(gcp);
return;
