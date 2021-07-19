function fout = PTF_plotMarginals(LongTermInfo,ScenarioProb,ShortTermProbDistr,EarlyEst,POIInfo,UsedFunctions)
disp('----- running: PTF_plotMarginals ------')

% MARGINAL DISTRIBUTIONS
fout=figure();
set(gcf,'position',[216    185   1112    579]);
nc = 3; 
nl = 3;
iSub = 0;

% MAGNITUDE
iSub = iSub+1;
subplot(nc,nl,1)
MargMagBS = zeros(size(LongTermInfo.Discretizations.PS1_Mag.Val));
MargMagPS = zeros(size(LongTermInfo.Discretizations.PS1_Mag.Val));
for imag = 1:length(LongTermInfo.Discretizations.PS1_Mag.Val)
  isel = find(ScenarioProb.ParScenBS(:,2) == LongTermInfo.Discretizations.PS1_Mag.Val(imag));
  if not(isempty(isel))
    MargMagBS(imag) = sum(ScenarioProb.ProbScenBS(isel));
  end
  isel = find(ScenarioProb.ParScenPS(:,2) == LongTermInfo.Discretizations.PS1_Mag.Val(imag));
  if not(isempty(isel))
    MargMagPS(imag) = sum(ScenarioProb.ProbScenPS(isel));
  end
end
plot(LongTermInfo.Discretizations.PS1_Mag.Val,ShortTermProbDistr.MagProb,'k:')
hold on
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMagBS'+MargMagPS','go--')
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMagBS,'bo--')
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMagPS,'ro--')
legend('Theoretical','BS+PS','BS','PS')
xlabel('Mw')

% MAGNITUDE SEPARATION
iSub = iSub+1;
subplot(nc,nl,4)
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMagBS./(MargMagBS+MargMagPS),'bo-')
hold on
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMagPS./(MargMagBS+MargMagPS),'ro-')
legend('Pr(BS|M)','Pr(PS|M)')
xlabel('Mw')


% DEPTH
iSub = iSub+1;
subplot(nc,nl,7)
tmpBSDepths = zeros(size(ScenarioProb.ParScenBS(:,5)));
for imag = 1:length(LongTermInfo.Discretizations.BS1_Mag.Val)
    isel = find(ScenarioProb.ParScenBS(:,2) == LongTermInfo.Discretizations.BS1_Mag.Val(imag));
    if not(isempty(isel))
        tmpBSDepths(isel) = ScenarioProb.ParScenBS(isel,5) + UsedFunctions.CorrectBSverticalPos(LongTermInfo.Discretizations.BS1_Mag.Val(imag))*1.e-3;
    end
end

tmpPSDepths=ScenarioProb.ParScenPS(:,7);
binEdges = linspace(0, 50, 11);
[h,whichBSBin] = histc(tmpBSDepths, binEdges);
[h,whichPSBin] = histc(tmpPSDepths, binEdges);
DepthMarg = zeros(length(binEdges)-1,1);
for ibin=1:length(binEdges)-1
    DepthMarg(ibin) = sum(ScenarioProb.ProbScenBS(whichBSBin == ibin)) + ...
                      sum(ScenarioProb.ProbScenPS(whichPSBin == ibin));
end
binCenters = 0.5*(binEdges(1:end-1)+binEdges(2:end));
bar(binCenters,DepthMarg);
xlabel('Depth (center of fault), km')


% POSITION
iSub = iSub+1;
subplot(nc,nl,[2 3 5 6 8 9])

if strcmp(LongTermInfo.ProjectName,'ChEESE_PTF')
    deltaDeg=1;
else
    deltaDeg=0.25;
end
if not( ShortTermProbDistr.PScomputedYN)
    binXEdges = floor(min(ScenarioProb.ParScenBS(:,3))):deltaDeg:ceil(max(ScenarioProb.ParScenBS(:,3)));
    binYEdges = floor(min(ScenarioProb.ParScenBS(:,4))):deltaDeg:ceil(max(ScenarioProb.ParScenBS(:,4)));
elseif not( ShortTermProbDistr.BScomputedYN)
    binXEdges = floor(min(ScenarioProb.ParScenPS(:,3))):deltaDeg:ceil(max(ScenarioProb.ParScenPS(:,3)));
    binYEdges = floor(min(ScenarioProb.ParScenPS(:,4))):deltaDeg:ceil(max(ScenarioProb.ParScenPS(:,4)));
else
    binXEdges = floor(min(vertcat(ScenarioProb.ParScenBS(:,3),ScenarioProb.ParScenPS(:,3)))):deltaDeg:ceil(max(vertcat(ScenarioProb.ParScenBS(:,3),ScenarioProb.ParScenPS(:,3))));
    binYEdges = floor(min(vertcat(ScenarioProb.ParScenBS(:,4),ScenarioProb.ParScenPS(:,4)))):deltaDeg:ceil(max(vertcat(ScenarioProb.ParScenBS(:,4),ScenarioProb.ParScenPS(:,4))));
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
for i = 1:length (ScenarioProb.ProbScenPS) %For each data point
  x_loc = find(binXEdges > ScenarioProb.ParScenPS(i,3),1) - 1;
  y_loc = find(binYEdges > ScenarioProb.ParScenPS(i,4),1) - 1;
  if or(x_loc == 0,y_loc==0)
    x_loc
    y_loc
    ScenarioProb.ParScenPS(i,3:4)
    i
  end    

  tmpX(x_loc,y_loc)=binXCenters(x_loc);
  tmpY(x_loc,y_loc)=binYCenters(y_loc);
  tmpZ(x_loc,y_loc)=tmpZ(x_loc,y_loc)+ScenarioProb.ProbScenPS(i);
end
for i = 1:length (ScenarioProb.ProbScenBS) %For each data point
  x_loc = find(binXEdges > ScenarioProb.ParScenBS(i,3),1) - 1;
  y_loc = find(binYEdges > ScenarioProb.ParScenBS(i,4),1) - 1;
  if or(x_loc == 0,y_loc==0)
    x_loc
    y_loc
    ScenarioProb.ParScenBS(i,3:4)
    i
  end    
  
  tmpX(x_loc,y_loc)=binXCenters(x_loc);
  tmpY(x_loc,y_loc)=binYCenters(y_loc);
  tmpZ(x_loc,y_loc)=tmpZ(x_loc,y_loc)+ScenarioProb.ProbScenBS(i);
end
scatter(tmpX(abs(tmpX)>0),tmpY(abs(tmpX)>0),100,tmpZ(abs(tmpX)>0),'filled')
hold on
if strcmp(LongTermInfo.ProjectName(1:7),'TSUMAPS') 
  plot(POIInfo.lon(POIInfo.Mediterranean),POIInfo.lat(POIInfo.Mediterranean),'.b')
  plot(POIInfo.lon([POIInfo.SelectedPOI.Index]),POIInfo.lat([POIInfo.SelectedPOI.Index]),'*b')
  TSUMAPS_plotRegions(LongTermInfo.Regionalization)
end

plot(EarlyEst.lon,EarlyEst.lat,'s','markerfacecolor','c','markeredgecolor','b','markersize',10)
colormap(flipud(hot))
colorbar
delta=2.5;%1.5
set(gca,'xlim',[min(binXEdges)-delta,max(binXEdges)+delta],'ylim',[min(binYEdges)-delta,max(binYEdges)+delta])
set(gca,'DataAspectRatio',ratioMap())

mtit('MARGINALS','fontsize',14,'color',[1 0 0],'xoff',0,'yoff',.025)

end
