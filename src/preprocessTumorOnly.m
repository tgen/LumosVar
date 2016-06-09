function [T, E]=preprocessTumorOnly(inputParam,paramFile)
%preprocessTumorOnly - creates data structures for tumor only calling
%calls parsePileupData.packed.pl to parse samtools output
%
% Syntax:  [T, E]=preprocessTumorOnly(inputParam,paramFile)
%
% Inputs:
%   inputParam - data structure with the following fields: regionsFile,
%       numCPU, outname, blockSize, snpVCFpath, snpVCFname,
%       workingDirectory, tabixPath, NormalBase
%   
% Outputs:
%   T - table of data by position
%   E - table of data by exon
%
% Other m-files required: none
% Other requirements: parsePileupData.packed.pl, samtools, htslib
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

%%%import bed file
fid=fopen(inputParam.regionsFile);
Regions=cell2mat(textscan(fid,'%d%d%d %*[^\n]'));
fclose(fid);

%%% start parallel pool
delete(gcp('nocreate'));
parpool(inputParam.numCPU);

%%% process by chromosome
chrList=[1:22];
parfor chrIdx=1:length(chrList)
    fout=fopen([inputParam.outName '_' num2str(chrList(chrIdx)) '_log.txt'],'w');
    currRegion=double(Regions(Regions(:,1)==chrList(chrIdx),:));
    largeIdx=find(currRegion(:,3)-currRegion(:,2)>inputParam.blockSize);
    %%%%breakup regions larger than blockSize
    if ~isempty(largeIdx)
        subCount=ceil((currRegion(largeIdx,3)-currRegion(largeIdx,2))./inputParam.blockSize);
        newSize=(currRegion(largeIdx,3)-currRegion(largeIdx,2))./subCount;
        newRegions=nan(sum(subCount),3);
        newRegions(1:subCount(1),:)=[double(chrList(chrIdx))*ones(subCount(1),1) round([currRegion(largeIdx(1),2):newSize(1):currRegion(largeIdx(1),3)-newSize(1)]') round([currRegion(largeIdx(1),2)+newSize(1)-1:newSize(1):currRegion(largeIdx(1),3)]')];
        for i=2:length(largeIdx)
            newRegions(sum(subCount(1:i-1))+1:sum(subCount(1:i)),:)=[double(chrList(chrIdx))*ones(subCount(i),1) round([currRegion(largeIdx(i),2):newSize(i):currRegion(largeIdx(i),3)-newSize(i)]') round([currRegion(largeIdx(i),2)+newSize(i)-1:newSize(i):currRegion(largeIdx(i),3)]')];
        end
        currRegion=[currRegion(currRegion(:,3)-currRegion(:,2)<=inputParam.blockSize,:); newRegions];
        currRegion=sortrows(currRegion,2);
    end
    %%%% make sure snpVCF is valid
    snpVCF=[inputParam.snpVCFpath num2str(chrList(chrIdx)) inputParam.snpVCFname];
    if(~exist(snpVCF,'file'))
        fprintf(fout,'%s\n',['snpVCF file not found in: ' snpVCF]);
        continue;
    else
        fprintf(fout,'%s\n',['snpVCF file found: ' snpVCF]);
    end
    %%% make sure normal data is valid
    NormalPath=[inputParam.NormalBase num2str(chrList(chrIdx)) '.txt.gz'];
    if(~exist(NormalPath,'file'))
        fprintf(fout,'%s\n',['NormalData file not found in: ' NormalPath]);
        continue;
    else
        fprintf(fout,'%s\n',['NormalData file found: ' NormalPath]);
    end
    
    %%%get data by region
    startIdx=1;
    endIdx=1;
    dataIdx=1;
    tempData=zeros(round(sum(currRegion(:,3)-currRegion(:,2))./100),26);
    tempExonRD=zeros(size(currRegion,1),8);
    while(startIdx<=size(currRegion,1))
        %%% define block
        endIdx=find(cumsum(currRegion(startIdx:end,3)-currRegion(startIdx:end,2))>inputParam.blockSize,1)+startIdx-1;
        if isempty(endIdx)
            endIdx=size(currRegion,1);
        end
        block=[num2str(chrList(chrIdx)) ':' num2str(currRegion(startIdx,2)) '-' num2str(currRegion(endIdx,3))]
        fprintf(fout,'\n%s',['Analyzing ' block]);
        %%% get tumor data
        cd(inputParam.workingDirectory);
        output=strsplit(perl('parsePileupData.packed.pl',paramFile,block,'1'),'\n');
        fprintf(fout,'\n%s',output{:});
        idx=~cellfun('isempty',regexp(output,'^\d'));
        TumorData=str2num(char(output(idx)'));
        fprintf(fout,'\n%s',['TumorData has length:' num2str(size(TumorData,1))]);
        idx=~cellfun('isempty',regexp(output,'^\@'));
        temp=char(output(idx)');
        readDepth=str2num(temp(:,2:end));
        fprintf(fout,'\n%s',['readDepth has length:' num2str(size(readDepth,1))]);
        %%% get normal data
        cd(inputParam.tabixPath);
        [status,output]=system(['./tabix ' NormalPath ' ' block]);
        NormalData=str2num(char(output));
        fprintf(fout,'\n%s',['NormalData has length:' num2str(size(NormalData,1))]);
        cd(inputParam.workingDirectory);
        currRegion(startIdx:endIdx,:);
        if isempty(readDepth)
            startIdx=endIdx+1;
            continue
        end
        %%% find exon means
        tempExonRD(startIdx:endIdx,:)=[currRegion(startIdx:endIdx,1) currRegion(startIdx:endIdx,2:3) getMeanInRegions(readDepth(:,1:2),readDepth(:,3),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,3),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,4),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,5),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,6),currRegion(startIdx:endIdx,:))];
        startIdx=endIdx+1;
        %%%% combine tumor and normal data
        if isempty(TumorData)
            continue
        end
        [Lia,Locb]=ismember(TumorData(:,2),NormalData(:,2));
        tempData(dataIdx:dataIdx+sum(Lia)-1,:)=[TumorData(Lia,:) NormalData(Locb(Lia),3:end)];
        dataIdx=dataIdx+sum(Lia);
        fprintf(fout,'\n%s',[' found ' num2str(sum(Lia)) ' canidate positions']);
        waitbar(endIdx./length(currRegion));
    end   
    AllData{chrIdx}=tempData(tempData(:,1)>0,:);
    AllExonData{chrIdx}=tempExonRD(tempExonRD(:,1)>0,:); 
end

message='finished processing chromosomes'

%%%% create position data table
ColHeaders={'Chr','Pos','ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ','AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ','BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount','ControlRD','PosMapQC','perReadPass','abFrac'};
matLen=cellfun(@length,AllData);
dataMat=zeros(sum(matLen),length(ColHeaders));
currIdx=1;
for i=1:length(AllData)
    dataMat(currIdx:currIdx+matLen(i)-1,:)=AllData{i};
    currIdx=currIdx+matLen(i);
end
if isempty(dataMat)
    T=[];
else
    T=array2table(dataMat,'VariableNames',ColHeaders);
end
clear dataMat AllData;
message='finished combining tumor data'

%%% create exon data table
exonColHeaders={'Chr','StartPos','EndPos','TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'};
exonMatLen=cellfun(@length,AllExonData);
exonRD=zeros(sum(exonMatLen),sum(length(exonColHeaders)));
currIdx=1;
for i=1:length(AllExonData)
    exonRD(currIdx:currIdx+exonMatLen(i)-1,:)=AllExonData{i};
    currIdx=currIdx+exonMatLen(i);
end
E=array2table(exonRD,'VariableNames',exonColHeaders);
message='finished combining exon data'

return;
