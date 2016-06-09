function writeVCF(T,pSomatic,pTrust,pArtifact,pGermline,pHom,cloneId,f,W,inputParam,pDataSomatic,pDataHet,pDataHom,pCNA,F)
%writeVCF - writes VCF for SNV and indel calls
%
% Syntax: writeVCF(T,pSomatic,posterior,pArtifact,pGermline,pHom,cloneId,f,W,inputParam,pDataSomatic,pDataHet,pDataHom,pCNA,F)
%
% Inputs:
%   T: table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%       'ControlRD','PosMapQC','perReadPass','abFrac'
%   pSomatic: posterior probability somatic variant
%   pTrust: posterior probability call should be trusted
%   pArtifact: posterior probability artificat
%   pGermline: posterior probability germline heterozygous
%   pHom: posterior probability homozygous AA
%   cloneId: most likley clone assuming somatic
%   f: sample fraction of each clone
%   W: w parameter for each clone
%   inputParam: structure with all parameters
%   pDataSomatic: likliehood somatic
%   pDataHet: likliehood germline het
%   pCNA: posterior probability of copy number alteration
%   F: sample fraction containing copy number variant
%   
% Outputs:
%    writes a VCF file
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, callSNV, qualDiscrim

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016
%------------- BEGIN CODE --------------

fout=fopen([inputParam.outName '.tumorOnly.all.vcf'],'w');

%%%print VCF header
fprintf(fout,'##fileformat=VCFv4.2\n');
fprintf(fout,['##fileData=' datestr(clock) '\n']);
inputFields=fieldnames(inputParam);
for i=1:length(inputFields)
    if(isnumeric(inputParam.(inputFields{i})))
        fprintf(fout,['##INPUT=<' inputFields{i} '=' mat2str(inputParam.(inputFields{i})') '>\n']);
    else
        fprintf(fout,['##INPUT=<' inputFields{i} '=' inputParam.(inputFields{i}) '>\n']);
    end
end
for i=1:length(f)
    fprintf(fout,['##CloneID=' num2str(i) ',f=', num2str(f(i)) ',W=' num2str(W(i)) ',PassCount=' num2str(sum(pSomatic(:,end)>inputParam.pSomaticThresh & pTrust>inputParam.pGoodThresh & min([T.ApopAF T.BpopAF],[],2)<inputParam.maxSomPopFreq & cloneId==i))]);
end
fprintf(fout,['\n##INFO=<ID=DP,Number=1,Type=Float,Description="Total Depth">\n']);
fprintf(fout,['##INFO=<ID=DPQC,Number=1,Type=Float,Description="Depth of Reads with MQ and BQ greater than min">\n']);
fprintf(fout,['##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n']);
fprintf(fout,['##INFO=<ID=A,Number=1,Type=String,Description="Major Allele Base">\n']);
fprintf(fout,['##INFO=<ID=B,Number=1,Type=String,Description="Minor Allele Base">\n']);
fprintf(fout,['##INFO=<ID=PT,Number=1,Type=Float,Description="Posterior Probability the Call can be Trusted">\n']);
fprintf(fout,['##INFO=<ID=PV,Number=1,Type=Float,Description="Posterior Probability the Position is not an Artifact">\n']);
fprintf(fout,['##INFO=<ID=PS,Number=1,Type=Float,Description="Posterior Probability of Somatic Mutation">\n']);
fprintf(fout,['##INFO=<ID=PGAB,Number=1,Type=Float,Description="Posterior Probability of No Somatic Mutation and Position is Germline AB">\n']);
fprintf(fout,['##INFO=<ID=PGAA,Number=1,Type=Float,Description="Posterior Probability of No Somatic Mutation and Position is Germline AA">\n']);
fprintf(fout,['##INFO=<ID=LS,Number=1,Type=Float,Description="Likelihood of Somatic Mutation">\n']);
fprintf(fout,['##INFO=<ID=LGAB,Number=1,Type=Float,Description="Likelihood of No Somatic Mutation and Position is Germline AB">\n']);
fprintf(fout,['##INFO=<ID=LGAA,Number=1,Type=Float,Description="Likelihood of No Somatic Mutation and Position is Germline AA">\n']);
fprintf(fout,['##INFO=<ID=ACountF,Number=1,Type=Float,Description="Number of Reads Supporting the A Allele on the Forward Strand">\n']);
fprintf(fout,['##INFO=<ID=ACountR,Number=1,Type=Float,Description="Number of Reads Supporting the A Allele on the Reverse Strand">\n']);
fprintf(fout,['##INFO=<ID=BCountF,Number=1,Type=Float,Description="Number of Reads Supporting the B Allele on the Forward Strand">\n']);
fprintf(fout,['##INFO=<ID=BCountR,Number=1,Type=Float,Description="Number of Reads Supporting the B Allele on the Forward Strand">\n']);
fprintf(fout,['##INFO=<ID=ApopAF,Number=1,Type=Float,Description="population frequency of A allele">\n']);
fprintf(fout,['##INFO=<ID=BpopAF,Number=1,Type=Float,Description="population frequency of B allele">\n']);
fprintf(fout,['##INFO=<ID=CosmicCount,Number=1,Type=Float,Description="Count of mutations at position in COSMIC>\n']);
fprintf(fout,['##INFO=<ID=CloneId,Number=1,Type=Integer,Description="CloneId">\n']);
fprintf(fout,['##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy Number">\n']);
fprintf(fout,['##INFO=<ID=MACN,Number=1,Type=Integer,Description="Min Allele Copy Number">\n']);
fprintf(fout,['##INFO=<ID=CNF,Number=1,Type=Float,Description="Fraction containg Copy Number Alteration">\n']);
fprintf(fout,['##INFO=<ID=PCNA,Number=1,Type=Float,Description="Posterior Probability Copy Number is Correct">\n']);
for i=1:size(F,2)
    fprintf(fout,'%s',char(strcat('##INFO=<ID=',F.Properties.VariableNames(i),',Number=1,Type=Float,Description="',F.Properties.VariableDescriptions(i),'">\n')));
end
fprintf(fout,['##FILTER=<ID=SomaticPASS,Description="PS>pSomaticThresh and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=SomaticLowQC,Description="PS>0.5 and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlinePASS,Description="PG>pGermlineThresh and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineLowQC,Description="PG>0.5 and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=NoCall,Description="PG<0.5 and PS<0.5 and pass filters>\n']);
fprintf(fout,['##FILTER=<ID=REJECT,Description="Below at least one QC filter">\n']);

%%% get reference and alternate bases
RefNT=int2ntIndels(T.Ref);
AltNT(T.Ref~=T.A & T.Ref~=T.B,:)=cellstr(strcat(int2ntIndels(T.A(T.Ref~=T.A & T.Ref~=T.B)),',', int2ntIndels(T.B(T.Ref~=T.A & T.Ref~=T.B))));
AltNT(T.Ref==T.A,:)=cellstr(int2ntIndels(T.B(T.Ref==T.A)));
AltNT(T.Ref==T.B,:)=cellstr(int2ntIndels(T.A(T.Ref==T.B)));

%%% calculate quality
Qual=zeros(size(T.Pos));
Qual(pSomatic(:,end)>0.5,:)=min(min(-10*log10(1-pSomatic(pSomatic(:,end)>0.5,end)),-10*log10(1-pTrust(pSomatic(:,end)>0.5))),-10*log10(1-pArtifact(pSomatic(:,end)>0.5,:)));
Qual(pGermline(:,end)>0.5,:)=min(min(-10*log10(1-pGermline(pGermline(:,end)>0.5,end)),-10*log10(1-pTrust(pGermline(:,end)>0.5))),-10*log10(1-pArtifact(pGermline(:,end)>0.5,:)));
Qual(pHom(:,end)>0.5,:)=min(min(-10*log10(1-pHom(pHom(:,end)>0.5,end)),-10*log10(1-pTrust(pHom(:,end)>0.5))),-10*log10(1-pArtifact(pHom(:,end)>0.5,:)));

%%% assign filters
Filter(pSomatic(:,end)>0.5,:)={'SomaticLowQC'};
Filter(pSomatic(:,end)>inputParam.pSomaticThresh & pTrust>inputParam.pGoodThresh,:)={'SomaticPASS'};
Filter(pSomatic(:,end)>0.5 & min([T.ApopAF T.BpopAF],[],2)>inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
Filter(pGermline(:,end)>0.5,:)={'GermlineHetLowQC'};
Filter(pGermline(:,end)>inputParam.pGermlineThresh & pTrust>inputParam.pGoodThresh,:)={'GermlineHetPASS'};
Filter(pHom(:,end)>0.5,:)={'GermlineHomLowQC'};
Filter(pHom(:,end)>inputParam.pGermlineThresh & pTrust>inputParam.pGoodThresh,:)={'GermlineHomPASS'};
Filter(pSomatic(:,end)<0.5 & pGermline(:,end)<0.5 & pHom(:,end)<0.5,:)={'NoCall'};
Filter(pArtifact<inputParam.pGoodThresh,:)={'REJECT'};

%%% construct info fields
ABfrac=(T.ACountF+T.ACountR+T.BCountF+T.BCountR)./T.ReadDepthPass;
Info=cellstr(strcat('DP=',num2str(T.ReadDepth,'%-d'),';DPQC=',num2str(T.ReadDepthPass,'%-d')));
Info(T.Ref==T.A,:)=strcat(Info(T.Ref==T.A,:),';AF=',num2str((T.BCountF(T.Ref==T.A,:)+T.BCountR(T.Ref==T.A,:))./T.ReadDepthPass(T.Ref==T.A,:),'%-.3f'));
Info(T.Ref==T.B,:)=strcat(Info(T.Ref==T.B,:),';AF=',num2str((T.ACountF(T.Ref==T.B,:)+T.ACountR(T.Ref==T.B,:))./T.ReadDepthPass(T.Ref==T.B,:),'%-.3f'));
Info=strcat(Info,';PT=',num2str(pTrust,'%-.5f'),';PV=',num2str(pArtifact,'%-.5f'),';PS=',num2str(pSomatic(:,end),'%-.5f'),';PGAB=',num2str(pGermline(:,end),'%-.5f'),';PGAA=',num2str(pHom(:,end),'%-.5f'));
Info=strcat(Info,';LS=',num2str(pDataSomatic(:,end),'%-.5f'),';LGAB=',num2str(pDataHet(:,end),'%-.5f'),';LGAA=',num2str(pDataHom(:,end),'%-.5f'));
Info=strcat(Info,';A=',int2ntIndels(T.A),';B=',int2ntIndels(T.B));
Info=strcat(Info,';ACountF=',num2str(T.ACountF,'%-d'),';ACountR=',num2str(T.ACountR,'%-d'),';BCountF=',num2str(T.BCountF,'%-d'),';BCountR=',num2str(T.BCountR,'%-d'));
Info=strcat(Info,';ApopAF=',num2str(T.ApopAF,'%-.5f'),';BpopAF=',num2str(T.BpopAF,'%-.5f'),';CosmicCount=',num2str(T.CosmicCount,'%-d'));
Info=strcat(Info,';CloneId=',num2str(cloneId(:,end)));
Info(T.NumCopies~=2 | T.MinAlCopies~=1)=strcat(Info(T.NumCopies~=2 | T.MinAlCopies~=1),';CN=',num2str(T.NumCopies(T.NumCopies~=2 | T.MinAlCopies~=1)),';MACN=',num2str(T.MinAlCopies(T.NumCopies~=2 | T.MinAlCopies~=1)),';CNF=',num2str(T.cnaF(T.NumCopies~=2 | T.MinAlCopies~=1)));
Info=strcat(Info,';PCNA=',num2str(pCNA,'%-.5f'));
for i=1:size(F,2)
    Info=strcat(Info,';',F.Properties.VariableNames(i),'=',num2str(F{:,i},'%-.3f'));
end

%%% print output
outData=[num2cell(T.Chr) num2cell(T.Pos) cellstr(char(ones(size(T.Pos))*46)) cellstr(RefNT) squeeze(AltNT) num2cell(Qual) Filter Info];
headers={'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILT', 'INFO'};
for i=1:length(headers)
    fprintf(fout,'%s\t',headers{i});
end
for i=1:size(outData,1)
    fprintf(fout,'\n%d\t%d\t%s\t%s\t%s\t%f\t%s\t%s',outData{i,:});
end

