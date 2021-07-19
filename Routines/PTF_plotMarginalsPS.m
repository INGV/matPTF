function fout = PTF_plotMarginalsPS(LongTermInfo,ScenarioProb,ShortTermProbDistr,PreSelection,EarlyEst,POIInfo,UsedFunctions)
disp('----- running: PTF_plotMarginalsPS ------')


% MARGINAL DISTRIBUTIONS
fout=figure();
set(gcf,'position',[216    185   1200    579]);
nc = 2; 
nl = 3;
iSub = 0;

% MAGNITUDE
iSub = iSub+1;
subplot(nc,nl,1)
MargMag = zeros(size(LongTermInfo.Discretizations.PS1_Mag.Val));
for imag = 1:length(LongTermInfo.Discretizations.PS1_Mag.Val)
  isel = find(ScenarioProb.ParScenPS(:,2) == LongTermInfo.Discretizations.PS1_Mag.Val(imag));
  if not(isempty(isel))
    MargMag(imag) = sum(ScenarioProb.ProbScenPS(isel));
  end
end
plot(LongTermInfo.Discretizations.PS1_Mag.Val,ShortTermProbDistr.MagProb,'k:')
hold on
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMag,'r--')
legend('Theoretical','PS')
xlabel('Mw')

% DEPTH
iSub = iSub+1;
subplot(nc,nl,nl+1)
tmpPSDepths = zeros(size(ScenarioProb.ParScenPS(:,7)));
tmpPSDepths=ScenarioProb.ParScenPS(:,7);
maxdepth=ceil(max(tmpPSDepths*0.1))*10;
if isempty(maxdepth)
    maxdepth=1;
end
binEdges = linspace(0, maxdepth, 1+2*maxdepth/5);
[h,whichPSBin] = histc(tmpPSDepths, binEdges);
DepthPSMarg = zeros(length(binEdges)-1,1);
for ibin=1:length(binEdges)-1
    DepthPSMarg(ibin) = sum(ScenarioProb.ProbScenPS(whichPSBin == ibin));
end
binCenters = 0.5*(binEdges(1:end-1)+binEdges(2:end));
bar(binCenters,DepthPSMarg);
xlabel('Depth (center of fault), km')

% POSITION
iSub = iSub+1;
subplot(nc,nl,[2:nl nl+2:2*nl])
if strcmp(LongTermInfo.ProjectName,'ChEESE_PTF')
    deltaDeg=1;
else
    deltaDeg=0.25;
end
if ScenarioProb.BScomputedYN
    binXEdges = floor(min(vertcat(ScenarioProb.ParScenBS(:,3),ScenarioProb.ParScenPS(:,3)))):deltaDeg:ceil(max(vertcat(ScenarioProb.ParScenBS(:,3),ScenarioProb.ParScenPS(:,3))));
    binYEdges = floor(min(vertcat(ScenarioProb.ParScenBS(:,4),ScenarioProb.ParScenPS(:,4)))):deltaDeg:ceil(max(vertcat(ScenarioProb.ParScenBS(:,4),ScenarioProb.ParScenPS(:,4))));
else
    binXEdges = floor(min(ScenarioProb.ParScenPS(:,3))):deltaDeg:ceil(max(ScenarioProb.ParScenPS(:,3)));
    binYEdges = floor(min(ScenarioProb.ParScenPS(:,4))):deltaDeg:ceil(max(ScenarioProb.ParScenPS(:,4)));
end
binXCenters = 0.5*(binXEdges(1:end-1)+binXEdges(2:end));
binYCenters = 0.5*(binYEdges(1:end-1)+binYEdges(2:end));
for x_loc=1:length(binXEdges)-1
  for y_loc=1:length(binYEdges)-1
     tmpIndBar{x_loc, y_loc} = [];
  end
end
tmpX=zeros(length (binXCenters),length (binYCenters));
tmpY=zeros(length (binXCenters),length (binYCenters));
tmpZ=zeros(length (binXCenters),length (binYCenters));
magseltmp=unique(ScenarioProb.ParScenPS(:,2));
tmpXperM= zeros(length (binXCenters),length (binYCenters), length(magseltmp));
tmpYperM= zeros(length (binXCenters),length (binYCenters), length(magseltmp));
tmpZperM = zeros(length (binXCenters),length (binYCenters), length(magseltmp));
for i = 1:length (ScenarioProb.ProbScenPS) %For each data point
  x_loc = find(binXEdges > ScenarioProb.ParScenPS(i,3),1) - 1;
  y_loc = find(binYEdges > ScenarioProb.ParScenPS(i,4),1) - 1;
  if or(x_loc == 0,y_loc==0)
    x_loc
    y_loc
    ParScenPS(i,3:4)
    i
  end        
  
  tmpX(x_loc,y_loc)=binXCenters(x_loc);
  tmpY(x_loc,y_loc)=binYCenters(y_loc);
  tmpZ(x_loc,y_loc)=tmpZ(x_loc,y_loc)+ScenarioProb.ProbScenPS(i);
  
  
  [a,imag]=ismember(ScenarioProb.ParScenPS(i,2),magseltmp);  
  tmpXperM(x_loc,y_loc,imag)=binXCenters(x_loc);
  tmpYperM(x_loc,y_loc,imag)=binYCenters(y_loc);
  tmpZperM(x_loc,y_loc,imag)=tmpZperM(x_loc,y_loc,imag)+ScenarioProb.ProbScenPS(i);
end
hold off
plot(EarlyEst.lon,EarlyEst.lat,'s','markerfacecolor','c','markeredgecolor','b','markersize',10)
hold on
scatter(tmpX(abs(tmpX)>0),tmpY(abs(tmpX)>0),100,tmpZ(abs(tmpX)>0),'filled')

for i=1:size(PreSelection.tmpR_BS_nsigma,1)
  [y(i) x(i)] = utm2ll(PreSelection.tmpR_PS_nsigma(i,2),PreSelection.tmpR_PS_nsigma(i,1),EarlyEst.RefUtmZone);
end
plot(x,y,'k--')
if strcmp(LongTermInfo.ProjectName(1:7),'TSUMAPS')
  plot(POIInfo.lon(POIInfo.Mediterranean),POIInfo.lat(POIInfo.Mediterranean),'.k')
  TSUMAPS_plotRegions(LongTermInfo.Regionalization)
end
colormap(flipud(hot))
colorbar
delta=2.5; % 1.5
set(gca,'xlim',[min(binXEdges)-delta,max(binXEdges)+delta],'ylim',[min(binYEdges)-delta,max(binYEdges)+delta])
set(gca,'DataAspectRatio',ratioMap())
xlabel('Position x,y (center of fault)')

mtit('MARGINAL PS','fontsize',14,'color',[1 0 0],'xoff',0,'yoff',.025)
end

