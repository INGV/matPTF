function fout = PTF_plotMarginalsBS(LongTermInfo,ScenarioProb,ShortTermProbDistr,PreSelection,EarlyEst,POIInfo,UsedFunctions)
disp('----- running: PTF_plotMarginalsBS ------')


ScenarioProb.ProbScenBS(isnan(ScenarioProb.ProbScenBS))=0;

% MARGINAL DISTRIBUTIONS
fout=figure();
set(gcf,'position',[216    185   1112    579]);
nc = 3;
nl = 3;
iSub = 0;

% MAGNITUDE
iSub = iSub+1;
subplot(nc,nl,1)
MargMag = zeros(size(LongTermInfo.Discretizations.PS1_Mag.Val));
for imag = 1:length(LongTermInfo.Discretizations.PS1_Mag.Val)
    isel = find(ScenarioProb.ParScenBS(:,2) == LongTermInfo.Discretizations.PS1_Mag.Val(imag));
    if not(isempty(isel))
        MargMag(imag) = sum(ScenarioProb.ProbScenBS(isel));
    end
end
plot(LongTermInfo.Discretizations.PS1_Mag.Val,ShortTermProbDistr.MagProb,'k:')
hold on
plot(LongTermInfo.Discretizations.PS1_Mag.Val,MargMag,'r--')
legend('Theoretical','BS')
xlabel('Mw')

% POSITION
iSub = iSub+1;
subplot(nc,nl,[2 3 5 6])
MargPos = zeros(size(PreSelection.PosBSSel));
for i2 = 1:length(PreSelection.PosBSSel)
    isel = find(and(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel(i2),1) == ScenarioProb.ParScenBS(:,3), ...
                    LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel(i2),2) == ScenarioProb.ParScenBS(:,4)));
    if not(isempty(isel))
        MargPos(i2) = sum(ScenarioProb.ProbScenBS(isel));
    end
end
scatter(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,1),LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,2),50,MargPos,'filled')
hold on
if strcmp(LongTermInfo.ProjectName(1:7),'TSUMAPS')
  plot(POIInfo.lon(POIInfo.Mediterranean),POIInfo.lat(POIInfo.Mediterranean),'.b')
  plot(POIInfo.lon([POIInfo.SelectedPOI.Index]),POIInfo.lat([POIInfo.SelectedPOI.Index]),'*b')
  TSUMAPS_plotRegions(LongTermInfo.Regionalization)
end
hold on
for i=1:size(PreSelection.tmpR_BS_nsigma,1)
    [y(i) x(i)] = utm2ll(PreSelection.tmpR_BS_nsigma(i,2),PreSelection.tmpR_BS_nsigma(i,1),EarlyEst.RefUtmZone);
    [yout(i) xout(i)] = utm2ll(PreSelection.tmpR_BS_nsigma_out(i,2),PreSelection.tmpR_BS_nsigma_out(i,1),EarlyEst.RefUtmZone);
end
plot(x,y,'k--')
plot(xout,yout,'k:')
plot(EarlyEst.lon,EarlyEst.lat,'s','markerfacecolor','c','markeredgecolor','b','markersize',10)
colormap(flipud(hot))
xlabel('Position x,y (center of fault)')
colorbar

delta=2.5; % 1.5
xlim =[min(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,1))-delta,max(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,1))+delta];
ylim =[min(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,2))-delta,max(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,2))+delta];
set(gca,'xlim',xlim,'ylim',ylim)
set(gca,'DataAspectRatio',ratioMap())


% DEPTH
iSub = iSub+1;
subplot(nc,nl,4)
tmpBSDepths = zeros(size(ScenarioProb.ParScenBS(:,5)));
for imag = 1:length(LongTermInfo.Discretizations.BS1_Mag.Val)
    isel = find(ScenarioProb.ParScenBS(:,2) == LongTermInfo.Discretizations.BS1_Mag.Val(imag));
    if not(isempty(isel))
        tmpBSDepths(isel) = ScenarioProb.ParScenBS(isel,5) + UsedFunctions.CorrectBSverticalPos(LongTermInfo.Discretizations.BS1_Mag.Val(imag))*1.e-3;
    end
end
binEdges = linspace(0, 50, 11);
[h,whichBSBin] = histc(tmpBSDepths, binEdges);
DepthBSMarg = zeros(length(binEdges)-1,1);
for ibin=1:length(binEdges)-1
    DepthBSMarg(ibin) = sum(ScenarioProb.ProbScenBS(whichBSBin == ibin));
end
binCenters = 0.5*(binEdges(1:end-1)+binEdges(2:end));
bar(binCenters,DepthBSMarg);
xlabel('Depth (center of fault), km')

% ANGLES: STRIKE
iSub = iSub+1;
subplot(nc,nl,7)
tmpStrikes = unique(LongTermInfo.Discretizations.BS4_FocMech.Val(:,1));
AngStrBSMarg = zeros(size(tmpStrikes));
for i = 1:length(tmpStrikes)
    isel = find(tmpStrikes(i) == ScenarioProb.ParScenBS(:,6));
    if not(isempty(isel))
        AngStrBSMarg(i) = sum(ScenarioProb.ProbScenBS(isel));
    end
end
bar(tmpStrikes,AngStrBSMarg)
xlabel('Strike')

% ANGLES: DIP
iSub = iSub+1;
subplot(nc,nl,8)
tmpDip= unique(LongTermInfo.Discretizations.BS4_FocMech.Val(:,2));
AngDipBSMarg = zeros(size(tmpDip));
for i = 1:length(tmpDip)
    isel = find(tmpDip(i) == ScenarioProb.ParScenBS(:,7));
    if not(isempty(isel))
        AngDipBSMarg(i) = sum(ScenarioProb.ProbScenBS(isel));
    end
end
bar(tmpDip,AngDipBSMarg)
xlabel('Dip')

% ANGLES: RAKE
iSub = iSub+1;
subplot(nc,nl,9)
tmpRake = unique(LongTermInfo.Discretizations.BS4_FocMech.Val(:,3));
AngRakBSMarg = zeros(size(tmpRake));
for i = 1:length(tmpRake)
    isel = find(tmpRake(i) == ScenarioProb.ParScenBS(:,8));
    if not(isempty(isel))
        AngRakBSMarg(i) = sum(ScenarioProb.ProbScenBS(isel));
    end
end
bar(tmpRake,AngRakBSMarg)
xlabel('Rake')


mtit('MARGINAL BS','fontsize',14,'color',[1 0 0],'xoff',0,'yoff',.025)
end
