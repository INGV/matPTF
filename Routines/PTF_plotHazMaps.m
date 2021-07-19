function fout=PTF_plotHazMaps(MapThreshold,HazardCurves,EarlyEst,POIInfo,typeMap,PreSettings)
% MapThreshold:  PROB THRESHOLD FOR HAZARD MAP
disp('----- running: PTF_plotHazMaps ------')

% plotting hc_poi
hc_poi = HazardCurves.hc_poi;
[tmp iThSel] = min(abs(MapThreshold - HazardCurves.HazardCurveThreshols));

if (strcmp(PreSettings.ProjectName,'ChEESE_PTF')) 
    POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end)))>180))=POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end)))>180))-360;
    xlim=[min(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))))-10 max(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))))+10];
    ylim=[min(POIInfo.lat(POIInfo.(EarlyEst.ID(11:end))))-10 max(POIInfo.lat(POIInfo.(EarlyEst.ID(11:end))))+10];
else
    xlim=[-10 40];
    ylim=[30 46];
end


fout=figure;
geolimits(ylim,xlim)

tmpZ = zeros(length(POIInfo.SelectedPOI),1);
switch typeMap
    case 'prob'
    for ipt=1:length(POIInfo.SelectedPOI)
       tmpZ(ipt) = interp1(HazardCurves.HazardCurveThreshols,hc_poi(ipt,:),MapThreshold);
    end
    tmpZ = log10(tmpZ);
    texttit = ['log10[Probability Map: Pr(>' num2str(MapThreshold) ' m)] - ' HazardCurves.tsunamiIntensityName];    
    
    
    case 'haz'

    for ipt=1:length(POIInfo.SelectedPOI)
        hc=hc_poi(ipt,:)';
        nonNull = find(hc>0);
        if isempty(nonNull)        
           tmpZ(ipt)=0;
        else
          if length(hc)>length(nonNull)
             nonNull(end+1)=nonNull(end)+1;
          end
          xx = [1;hc(nonNull)];
           yy = [0;HazardCurves.HazardCurveThreshols(nonNull)];
        uniq = [true; diff(xx) ~= 0];
       tmpZ(ipt) = interp1(xx(uniq),yy(uniq),MapThreshold);
        end
    end
    
    texttit = ['Pr(> z m) = ' num2str(MapThreshold) ' - ' HazardCurves.tsunamiIntensityName];
    
%     
    case 'hazmean'

    for ipt=1:length(POIInfo.SelectedPOI)
       tmpZ(ipt) = HazardCurves.hc_poi_mean(ipt);
    end

    texttit = ['E[z] m - ' HazardCurves.tsunamiIntensityName];
    
end

hold off

if (strcmp(PreSettings.ProjectName,'ChEESE_PTF'))  %%% ATTENTION!!! Here we have Maule but WE MUST GENERALIZE for any source in ChEESE!!!
    %plot(POIInfo.lon(POIInfo.maule),POIInfo.lat(POIInfo.maule),'.')
    geoplot(POIInfo.lat(POIInfo.(EarlyEst.ID(11:end))),POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))),'.')
else
    %plot(POIInfo.lon(POIInfo.Mediterranean),POIInfo.lat(POIInfo.Mediterranean),'.')
    geoplot(POIInfo.lat(POIInfo.Mediterranean),POIInfo.lon(POIInfo.Mediterranean),'.')
end

hold on
geoplot(EarlyEst.lat,EarlyEst.lon,'*r')
geoscatter(POIInfo.lat([POIInfo.SelectedPOI(:).Index]),POIInfo.lon([POIInfo.SelectedPOI(:).Index]),30,tmpZ,'filled')
%plot(EarlyEst.lon,EarlyEst.lat,'*r')
%scatter(POIInfo.lon([POIInfo.SelectedPOI(:).Index]),POIInfo.lat([POIInfo.SelectedPOI(:).Index]),30,tmpZ,'filled')
clim=get(gca,'clim');
title(texttit);
grid
a=colorbar;

switch typeMap
    case 'prob'
        set(gca,'clim',[-5 0]);
        txt = ['Log_{10}[Pr(> ' num2str(MapThreshold) ' m ) '];
        a.Label.String = txt;
    case 'haz'
        a.Label.String = HazardCurves.tsunamiIntensityName;
    case 'hazmean'
        a.Label.String = HazardCurves.tsunamiIntensityName;
end

end




