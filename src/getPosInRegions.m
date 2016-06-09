function idx=getPosInRegions(pos,regions)
%getPosInRegions - finds which region each position falls into
%
% Syntax: idx=getPosInRegions(pos,regions)
%
% Inputs:
%   pos: two column matrix where 1st col is chr and 2nd is pos
%   regions: three column matrix where each row specifies a region
%       1st col is chr, 2nd is start pos, 3rd col is end pos
%
% Outputs:
%   idx: vector same height as pos, number corresponds to row in regions
%       where pos is found, NaN if not in region
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getPosInRegionSplit, getMeanInRegions

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

idx=NaN(size(pos,1),1);
for i=min(pos(:,1)):max(pos(:,1))
    currPos=pos(pos(:,1)==i,:);
    currRegions=regions(regions(:,1)==i,:);
    if isempty(currPos) || isempty(currRegions)
        continue;
    end
    afterStart=ones(size(currRegions,1),1)*currPos(:,2)'>=currRegions(:,2)*ones(1,size(currPos,1));
    beforeEnd=ones(size(currRegions,1),1)*currPos(:,2)'<=currRegions(:,3)*ones(1,size(currPos,1));
    [m,currIdx]=max(afterStart+beforeEnd,[],1);
    currIdx(m<2)=NaN;
    idx(pos(:,1)==i,:)=currIdx+find(regions(:,1)==i,1)-1;
end