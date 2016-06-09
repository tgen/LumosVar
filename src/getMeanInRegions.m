function regionMean=getMeanInRegions(pos,values,regions)
%getMeanInRegions - finds mean of values within each region
%
% Syntax: regionMean=getMeanInRegions(pos,values,regions)
%
% Inputs:
%   pos: two column matrix where 1st col is chr and 2nd is pos
%   values: values to find mean of, same length as pos
%   regions: three column matrix where each row specifies a region
%       1st col is chr, 2nd is start pos, 3rd col is end pos
%
% Outputs:
%   regionMean: vector same height as regions with mean of value in regions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getPosInRegions

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

regionMean=NaN(size(regions,1),1);

for i=min(pos(:,1)):max(pos(:,1))
    currPos=pos(pos(:,1)==i,:);
    currRegions=regions(regions(:,1)==i,:);
    currValues=values(pos(:,1)==i,:);
    if isempty(currPos) || isempty(currRegions)
        continue;
    end
    afterStart=ones(size(currRegions,1),1)*currPos(:,2)'>=currRegions(:,2)*ones(1,size(currPos,1));
    beforeEnd=ones(size(currRegions,1),1)*currPos(:,2)'<=currRegions(:,3)*ones(1,size(currPos,1));
    regionBool=logical(afterStart+beforeEnd==2);
    idx=sum(regionBool,2)>0;
    currRegionMean=NaN(size(currRegions,1),1);
    currRegionMean(idx,:)=(regionBool(idx,:)*currValues)./sum(regionBool(idx,:),2);
    regionMean(regions(:,1)==i,:)=currRegionMean;
end

return;
