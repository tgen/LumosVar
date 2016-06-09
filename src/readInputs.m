function inputParam=readInputs(paramFile)
%readInputs - reads config file in yaml format
% Syntax:  inputParam=readInputs(paramFile)
%
% Inputs:
%    paramFile - parameter file in yaml format, see configTemplate.yaml 
%
% Outputs:
%    inputParam - data structure with fields from paramFile
%
%
% Other m-files required: none
% Other requirements: none
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

fid=fopen(paramFile);
C=textscan(fid,'%s %s','commentstyle','#');
fclose(fid);

for i=1:length(C{1})
    fieldname=regexprep(C{1}{i},':$','');
    if isempty(str2num(C{2}{i}))
        inputParam.(fieldname)=C{2}{i};
    else
        inputParam.(fieldname)=str2num(C{2}{i});
    end
end