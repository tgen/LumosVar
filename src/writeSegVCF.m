function writeSegVCF(segsTable,inputParam)
%writeSegVCF - writes VCF for copy number alterations
%
% Syntax: writeSegVCF(segsTable,inputParam)
%
% Inputs:
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%  inputParam: structure with all parameters   
%
% Outputs:
%    writes a VCF file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, callCNA, fitCNA

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016
%------------- BEGIN CODE --------------

fout=fopen([inputParam.outName '.cna.seg.vcf'],'w');

%%% print VCF header
fprintf(fout,'%s','##fileformat=VCFv4.2\n');
fprintf(fout,['##fileData=' datestr(clock) '\n']);
inputFields=fieldnames(inputParam);
for i=1:length(inputFields)
    if(isnumeric(inputParam.(inputFields{i})))
        fprintf(fout,['##INPUT=<' inputFields{i} '=' mat2str(inputParam.(inputFields{i})') '>\n']);
    else
        fprintf(fout,['##INPUT=<' inputFields{i} '=' inputParam.(inputFields{i}) '>\n']);
    end
end
fprintf(fout,['##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy Number">\n']);
fprintf(fout,['##INFO=<ID=MACN,Number=1,Type=Integer,Description="Min Allele Copy Number">\n']);
fprintf(fout,['##INFO=<ID=LOG2FC,Number=1,Type=Float,Description="log2 fold change">\n']);
fprintf(fout,['##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of Segment">\n']);
fprintf(fout,['##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural Variant">\n']);
fprintf(fout,['##INFO=<ID=END,Number=1,Type=String,Description="End of Variant">\n']);
fprintf(fout,['##INFO=<ID=CNF,Number=1,Type=Float,Description="Fraction containg Copy Number Alteration">\n']);

%%% determine alteration type
type(segsTable(:,5)>2 & segsTable(:,6)>0,:)={'DUP'};
type(segsTable(:,5)>2 & segsTable(:,6)==0,:)={'DUPLOH'};
type(segsTable(:,5)==2 & segsTable(:,6)==1,:)={'NONE'};
type(segsTable(:,5)==2 & segsTable(:,6)==0,:)={'LOH'};
type(segsTable(:,5)<2,:)={'DEL'};

%%% construct info field
Info=cellstr(strcat('CN=',num2str(segsTable(:,5)),';MACN=',num2str(segsTable(:,6))));
Info=strcat(Info,';LOG2FC=',num2str(segsTable(:,9)),';SVLEN=',num2str(segsTable(:,3)-segsTable(:,2)));
Info=strcat(Info,';SVTYPE=',type,';END=',num2str(segsTable(:,3)),';CNF=',num2str(segsTable(:,7)));

%%% write output
outData=[num2cell(segsTable(:,1)) num2cell(segsTable(:,2)) num2cell(segsTable(:,3)) cellstr(char(ones(size(segsTable,1),1)*78)) strcat('<',type,'>') num2cell(segsTable(:,9)) cellstr(char(ones(size(segsTable,1),1)*46)) Info];
headers={'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILT', 'INFO'};
for i=1:length(headers)
    fprintf(fout,'%s\t',headers{i});
end
for i=1:size(outData,1)
    fprintf(fout,'\n%d\t%d\t%d\t%s\t%s\t%f\t%s\t%s',outData{i,:});
end

