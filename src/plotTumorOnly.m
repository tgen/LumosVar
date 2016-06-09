function plotTumorOnly(exonRD,segsTable,CNAscale,f,T,somPos,hetPos,cloneId,inputParam)
%plotTumorOnly - plots summary figure of tumor only caller results
%
% Syntax: writeCloneSummary(segsTable,E,T,pSomatic,posterior,f,W,cloneId,inputParam)
%
% Inputs:
%   exonRD: matrix of exon data with columns: 1-'Chr',2-'StartPos',3-'EndPos',
%       4-'TumorRD',5-'NormalRD',6-'MapQC',7-'perReadPass',8-'abFrac'
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%   f: sample fraction of each clone
%   T - table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%   somPos: logical index of postions in T that are somatic
%   hetPos: logical index of positions in T that are germline heterozygous
%   cloneId: most likley clone assuming somatic
%   inputParam: structure with all parameters   
%
% Outputs:
%    writes a png file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, callCNA, callSNV

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016
%------------- BEGIN CODE --------------

%%% transform chromosome coord to linear coord
for i=1:22
    chrLen(i)=max(segsTable(segsTable(:,1)==i,3))+1E7;
end
chrOffset=[0; cumsum(chrLen(1:21))'];
for i=1:22
    segCoord(segsTable(:,1)==i,1)=(segsTable(segsTable(:,1)==i,2)+chrOffset(i))/1E6;
    segCoord(segsTable(:,1)==i,2)=(segsTable(segsTable(:,1)==i,3)+chrOffset(i))/1E6;
    exonCoord(exonRD(:,1)==i,1)=(exonRD(exonRD(:,1)==i,2)+chrOffset(i))/1E6;
    exonCoord(exonRD(:,1)==i,2)=(exonRD(exonRD(:,1)==i,3)+chrOffset(i))/1E6;
    Tcoord(T.Chr==i,1)=(T.Pos(T.Chr==i)+chrOffset(i))/1E6;
end

%%% calculate log2FC
log2FC=log2((CNAscale./2).*exonRD(:,4)./exonRD(:,5));
Nlog2R=log2(median(segsTable(:,7)).*segsTable(:,5)/2+(1-median(segsTable(:,7))));
Mlog2R=log2(median(segsTable(:,7)).*segsTable(:,6)/2+(1-median(segsTable(:,7))));

%%% plot exon log2FC
subplot(2,1,1);
for i=1:2:22
    idx=exonRD(:,1)==i;
    scatter(mean(exonCoord(idx,:),2),log2FC(idx),1,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
    hold on;
end
for i=2:2:22
    idx=exonRD(:,1)==i;
    scatter(mean(exonCoord(idx,:),2),log2FC(idx),1,'.','MarkerFaceColor',[0.8235    0.7059    0.5490],'MarkerEdgeColor',[0.8235    0.7059    0.5490])
    hold on;
end

%%% plot segments
colors='rgbcmy';
for i=1:size(f,2)
    pos=segsTable(:,7)==f(1,i) & (segsTable(:,5)~=2 | segsTable(:,6)~=1);
    plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos)',colors(i),'linewidth',8);
    plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)',colors(i),'linewidth',4);
end
pos=(segsTable(:,5)==2 & segsTable(:,6)==1);
plot(segCoord(pos,:)',ones(2,1)*Nlog2R(pos)','k','linewidth',8);
plot(segCoord(pos,:)',ones(2,1)*Mlog2R(pos)','k','linewidth',4);
ticks=[0 2.^[1:7]];
tickpos=log2(median(segsTable(:,7)).*ticks/2+(1-median(segsTable(:,7))));    
set(gca,'YTick',tickpos,'YTickLabel',ticks,'tickDir','out');
set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',[1:22],'FontSize',8);
axis([0 max(segCoord(:,2)) min(Mlog2R)-1 max(Nlog2R)+1])
title('Copy Number','FontSize',10);

%%% plot het AF
AF=(T.BCountF+T.BCountR)./T.ReadDepthPass;
subplot(2,1,2)
for i=1:2:22
    scatter(Tcoord(hetPos & T.Chr==i),AF(hetPos & T.Chr==i),1,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
    hold on;
end
for i=2:2:22
    scatter(Tcoord(hetPos & T.Chr==i),AF(hetPos & T.Chr==i),1,'.','MarkerFaceColor',[0.8235    0.7059    0.5490],'MarkerEdgeColor',[0.8235    0.7059    0.5490])
end

%%% plot somatic AF
for i=1:size(f,2)
    pos=cloneId(:,1)==i & somPos;
    scatter(Tcoord(pos),AF(pos),2,[colors(i) '*']);
end
axis([0 max(segCoord(:,2)) 0 0.5]);
set(gca,'XTick',(chrOffset'+chrLen./2)/1E6,'XTickLabel',[1:22],'tickDir','out','FontSize',8);
title('Minor Allele Frequencies','FontSize',10);
for i=1:length(f)
    annotation('textbox',[0.8 0.9-0.03*i 0.1 0.1],'String',['Clone ' num2str(i) ', f=',num2str(f(1,i))],'Color',colors(i),'EdgeColor','none','FontSize',10);
end

%%% print plot
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 7 4]);
print(gcf,'-dpng',[inputParam.outName '.png'],'-r300');
close(gcf);
return;