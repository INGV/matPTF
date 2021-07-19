function [AlertLevelsInfo] = PTF_AlertLevels(Settings,HazardCurves,POIInfo,EarlyEst)
% evaluate AL corresponding to hazard curves
% ALtype:           method for Alert Level computation:
%                   ='PrNNN': interval including percentile NN.N (for ex. Pr050 is percentile 5)
%                   ='Modal': interval including mode (intevals defined by ALint, with intensity defined by TsuInt)
%                   ='Average': interval including mean (intevals defined by ALint, with intensity defined by TsuInt)     
%                   ='Matrix': Decision Matrix of CAT-INGV (simplified)
%                   ='ENV': Envelop, as defined in Selva et al. 2021, modifying Catalan et al. 2020
%                   ='bestGuessScen': Best Matching Scenario, as defined in Selva et al. 2021
% ALout:            0 nothing (less than 0.25 --> to be updated depending on distance with EE)
%                   1 information 
%                   2 advisory
%                   3 watch   

    disp('----- running: PTF_AlertLevels ------')
    
    AlertLevelsInfo = struct;
    AlertLevelsInfo.DMroundYN = false; 
    % Notation: Pr100 == Percentile 10   
    AlertLevelsInfo.ALtype = {'Matrix','ENV','Pr010','Pr050','Pr100','Pr150','Pr200','Pr300','Pr400','Pr500','Average','bestGuessScen'};   
    AlertLevelsInfo.ALtypeName = {'Decision Matrix (DM)','Envelope (ENV)','Percetile 1','Percetile 5',...
        'Percetile 10','Percetile 15','Percetile 20','Percetile 30','Percetile 40','Percetile 50',...
        'Average','Best Match Scenario (BMS)'}; 
    
    
    AlertLevelsInfo.ALint = [0,0.10,0.5];
    AlertLevelsInfo.ALnames = {'Information','Advisory','Watch'};
    AlertLevelsInfo.ALcc = [1 1 1; 0 1 0; 1 0.8 0; 1 0 0];
    AlertLevelsInfo.MaxDist = 100;
    AlertLevelsInfo.Nnearest = 3;
    AlertLevelsInfo.RefUtmZone = HazardCurves.RefUtmZone;
    AlertLevelsInfo.selPoints = [POIInfo.SelectedPOI.Index];

    %% COMPUTE ALARM LEVEL FOR ALL SELECTED TYPES
    [x y] = ll2utm(POIInfo.lat(AlertLevelsInfo.selPoints),POIInfo.lon(AlertLevelsInfo.selPoints),HazardCurves.RefUtmZone); POIInfo.lonUTM = x; POIInfo.latUTM = y;
    [EpiX EpiY] = ll2utm(EarlyEst.lat,EarlyEst.lon,HazardCurves.RefUtmZone);
    dd = 1.E-3*sqrt((EpiX-POIInfo.lonUTM).^2+(EpiY-POIInfo.latUTM).^2);
    zone1 = dd<=100;
    zone2 = and(dd>100.,dd<=400);
    zone3 = dd>400;
    toll=1.e-3; % for optimizing integration


    % DEFINE INTENSITY MEASURE
    tsuInt = HazardCurves.HazardCurveThreshols;

    % RUN OVER POIs
    ALout = zeros(length(AlertLevelsInfo.selPoints),length(AlertLevelsInfo.ALtype));
    Valout = zeros(length(AlertLevelsInfo.selPoints),length(AlertLevelsInfo.ALtype));
    for itype = 1:length(AlertLevelsInfo.ALtype)
        disp(['...AL from ' AlertLevelsInfo.ALtype{itype}]);
        if strcmpi(AlertLevelsInfo.ALtype{itype},'matrix')
            % implemented assuming dist <= 40 km for M<=6.5 and dist <= 100 km for M>6.5
            
            if AlertLevelsInfo.DMroundYN
               mageff = round(EarlyEst.Mag*10)*0.1;
            else
               mageff = EarlyEst.Mag;
            end
            if and(mageff >= 5.5,mageff <= 6.0)
                ALout(:,itype) = 1;
            elseif and(mageff > 6.0,mageff <= 6.5)  
                ALout(zone1,itype) = 2;
                ALout(zone2,itype) = 1;
                ALout(zone3,itype) = 1;
            elseif and(mageff > 6.5,mageff <= 7.0)
                ALout(zone1,itype) = 3;
                ALout(zone2,itype) = 2;
                ALout(zone3,itype) = 1;
            elseif and(mageff > 7.0,mageff <= 7.5)
                ALout(zone1,itype) = 3;
                ALout(zone2,itype) = 3;
                ALout(zone3,itype) = 2;
            elseif mageff > 7.5
                ALout(zone1,itype) = 3;
                ALout(zone2,itype) = 3;
                ALout(zone3,itype) = 3;
            end
        elseif strcmpi(AlertLevelsInfo.ALtype{itype},'bestGuessScen')
           for ipt=1:length([POIInfo.SelectedPOI.Index])
              valTmp=HazardCurves.bestGuessScen_tsunamiintisity(ipt);
              
              ALout(ipt,itype) = max(find(valTmp<=[AlertLevelsInfo.ALint,inf],1)-1,1);                
              Valout(ipt,itype) = valTmp;                      
           end
        elseif strcmpi(AlertLevelsInfo.ALtype{itype},'ENV')
           for ipt=1:length([POIInfo.SelectedPOI.Index])
              valTmp=HazardCurves.env_tsunamiintisity(ipt);
              
              ALout(ipt,itype) = max(find(valTmp<=[AlertLevelsInfo.ALint,inf],1)-1,1);                
              Valout(ipt,itype) = valTmp;                      
           end
           
        else
            for ipt=1:length([POIInfo.SelectedPOI.Index])
                hc = HazardCurves.hc_poi(ipt,:)';
                iPOI = AlertLevelsInfo.selPoints(ipt);

                if sum(hc) > 0
                        dx = 0.01;
                        tsuMax = 100.;

                        ALint = [0 AlertLevelsInfo.ALint tsuMax];
                        if strcmpi(AlertLevelsInfo.ALtype{itype}(1:2),'Pr')
                           prThr = 0.001*str2double(AlertLevelsInfo.ALtype{itype}(3:5));
                           nonNull = find(hc>0);
			   if length(hc)>length(nonNull)
				   nonNull(end+1)=nonNull(end)+1;
			   end
                           if isempty(find(nonNull, 1))
                               valTmp=0;
                           else
                               xx = [1;hc(nonNull)];
                               yy = [0;tsuInt(nonNull)];
                               uniq        = [true; diff(xx) ~= 0];
                               valTmp = interp1(xx(uniq),yy(uniq),prThr);
                               
                               if (prThr<xx(end))
                                   disp('!!!!Warning: Val>100m')
                                   valTmp=xx(end);
                               end
                           end
                        elseif strcmpi(AlertLevelsInfo.ALtype{itype},'Modal')
                           [xval yval] = hc_decum([AlertLevelsInfo.ALint,2*AlertLevelsInfo.ALint(end)],tsuInt,hc);
                           yy = yval(1:end-1);
                           yy(end)=yval(end-1)+yval(end);
                           
                           [tmp isel]= max(yy);
                           valTmp = xval(isel);

                        elseif strcmpi(AlertLevelsInfo.ALtype{itype},'Average')
                          valTmp = HazardCurves.hc_poi_mean(ipt);
                    
                        else
                           disp('!! Method not recognized')
                           return
                        end
                       
                        ALout(ipt,itype) = find(valTmp<=[AlertLevelsInfo.ALint,inf],1)-1;                
                        Valout(ipt,itype) = valTmp;

                else
                   ALout(ipt,itype) = 1;                                
                end
            end
        end
    end    


    %% STORE IN STRUCTURE
    AlertLevelsInfo.AlertLevels = ALout;
    AlertLevelsInfo.ReferenceIntensity = Valout;
    AlertLevelsInfo.tsunamiIntensity = tsuInt;
 
end
