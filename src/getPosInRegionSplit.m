function idx=getPosInRegionSplit(pos,regions,split)
%getPosInRegions - finds which region each position falls into
%splits search space to conserve memory usage
%
% Syntax: idx=getPosInRegions(pos,regions,split)
%
% Inputs:
%   pos: two column matrix where 1st col is chr and 2nd is pos
%   regions: three column matrix where each row specifies a region
%       1st col is chr, 2nd is start pos, 3rd col is end pos
%   split: number of regions to examine at once
%
% Outputs:
%   idx: vector same height as pos, number corresponds to row in regions
%       where pos is found, NaN if not in region
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getPosInRegions, getMeanInRegions

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

idx=NaN(length(pos),1);
for i=min(pos(:,1)):max(pos(:,1))
    chrRegions=regions(regions(:,1)==i,:);
    for j=1:split:sum(regions(:,1)==i)
        currRegions=chrRegions(j:min(j+split-1,sum(regions(:,1)==i)),:);
        posIdx=pos(:,1)==i & pos(:,2)>=currRegions(1,2) & pos(:,2)<=currRegions(end,3);
        currPos=pos(posIdx,:);
        if isempty(currPos) || isempty(currRegions)
            continue;
        end
        afterStart=ones(size(currRegions,1),1)*currPos(:,2)'>=currRegions(:,2)*ones(1,size(currPos,1));
        beforeEnd=ones(size(currRegions,1),1)*currPos(:,2)'<=currRegions(:,3)*ones(1,size(currPos,1));
        [m,currIdx]=max(afterStart+beforeEnd,[],1);
        currIdx(m<2)=NaN;
        idx(posIdx,:)=currIdx+find(regions(:,1)==i,1)+j-2;
    end
end