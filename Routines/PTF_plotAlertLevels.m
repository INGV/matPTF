function [f1]=PTF_plotAlertLevels(AlertLevelsInfo,EarlyEst,POIInfo,PreSettings)
disp('----- running: PTF_plotHazMaps ------')

[EpiX EpiY] = ll2utm(EarlyEst.lat,EarlyEst.lon,EarlyEst.RefUtmZone);


if (strcmp(PreSettings.ProjectName,'ChEESE_PTF'))
    POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end)))>180))=POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end)))>180))-360;
    xlim=[min(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))))-10 max(POIInfo.lon(POIInfo.(EarlyEst.ID(11:end))))+10];
    ylim=[min(POIInfo.lat(POIInfo.(EarlyEst.ID(11:end))))-10 max(POIInfo.lat(POIInfo.(EarlyEst.ID(11:end))))+10];
end

npt=1000;
lin100UTM=ellipsedata([1.E5^2 0; 0 1.E5^2],[EpiY,EpiX],npt,1);
lin400UTM=ellipsedata([4.E5^2 0; 0 4.E5^2],[EpiY,EpiX],npt,1);
lin1000UTM=ellipsedata([1.E6^2 0; 0 1.E6^2],[EpiY,EpiX],npt,1);
for i=1:npt
    [lat100(i) lon100(i)] = utm2ll(lin100UTM(i,2),lin100UTM(i,1),EarlyEst.RefUtmZone);
    [lat400(i) lon400(i)] = utm2ll(lin400UTM(i,2),lin400UTM(i,1),EarlyEst.RefUtmZone);
    [lat1000(i) lon1000(i)] = utm2ll(lin1000UTM(i,2),lin1000UTM(i,1),EarlyEst.RefUtmZone);
end


if (strcmp(PreSettings.ProjectName,'pacific-chile'))
    
    selPoints = AlertLevelsInfo.selPoints;
    
    for i=1:length(AlertLevelsInfo.ALtype)
        
        figure();
        
        colormap(AlertLevelsInfo.ALcc);
        scatter(POIInfo.lon(selPoints),POIInfo.lat(selPoints),20,AlertLevelsInfo.AlertLevels(:,i),'filled');
        hold on
        plot(POIInfo.lon(selPoints),POIInfo.lat(selPoints),'ko');
        plot(EarlyEst.lon,EarlyEst.lat,'p','MarkerFaceColor','b','markersize',15)
        plot(lon100,lat100,'k--')
        plot(lon400,lat400,'k--')
        plot(lon1000,lat1000,'k--')
        title(AlertLevelsInfo.ALtypeName{i})
        
        cb=colorbar;
        set(gca,'clim',[-0.5 3.5]);
        set(cb,'ytick',[0:3],'yticklabel',{'Not Defined',AlertLevelsInfo.ALnames{:}})
        set(gca,'xlim',xlim,'ylim',ylim)
    end
    
else
    
    selPoints = AlertLevelsInfo.selPoints;
    
    for i=1:length(AlertLevelsInfo.ALtype)
        
        figure('position',[779   371   900   482]);
        
        Axes1H=subplot(1,1,1);

        colormap(AlertLevelsInfo.ALcc);
        scatter(POIInfo.lon(selPoints),POIInfo.lat(selPoints),20,AlertLevelsInfo.AlertLevels(:,i),'filled');
        hold on
        plot(POIInfo.lon(selPoints),POIInfo.lat(selPoints),'ko');
        plot(EarlyEst.lon,EarlyEst.lat,'p','MarkerFaceColor','b','markersize',15)
        plot(lon100,lat100,'k--')
        plot(lon400,lat400,'k--')
        plot(lon1000,lat1000,'k--')
        title(AlertLevelsInfo.ALtypeName{i})
        
        cb=colorbar;
        set(gca,'clim',[-0.5 3.5]);
        set(cb,'ytick',[0:3],'yticklabel',{'Not Defined',AlertLevelsInfo.ALnames{:}})
        if not(strcmp(PreSettings.ProjectName,'ChEESE_PTF'))
            set(gca,'xlim',[-6 36],'ylim',[30 46])
        end
        grid
        
        pos1H =  get(Axes1H,'position');
        axes('Position',[pos1H(1:2)-0.01 0.3*pos1H(3:4)])
        colormap(AlertLevelsInfo.ALcc);
        scatter(POIInfo.lon(selPoints),POIInfo.lat(selPoints),20,AlertLevelsInfo.AlertLevels(:,i),'filled');
        hold on
        plot(POIInfo.lon(selPoints),POIInfo.lat(selPoints),'ko');   

        plot(EarlyEst.lon,EarlyEst.lat,'p','MarkerFaceColor','b','markersize',15)
        plot(lon100,lat100,'k--')
        set(gca,'clim',[0 3])
        set(gca,'xlim',[EarlyEst.lon-2 EarlyEst.lon+2],'ylim',[EarlyEst.lat-1 EarlyEst.lat+1])
        set(gca,'XTick',[],'YTick',[])
        box on
        
    end
end



end
