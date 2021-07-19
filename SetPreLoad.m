function [PreSettings,LongTermInfo,POIInfo] =  SetPreLoad(eqID)
% preload data for specific regions in the domain of interest
% Input: 
% eqID: selected event. Options: 'all_regions','all_events','2003_0521_boumardes', '2015_0416_crete', '2015_1117_lefkada', 
%                                '2016_0125_gibraltar', '2016_1030_norcia','2017_0612_lesbo','2017_0720_kos-bodrum',
%                                '2018_1025_zante','2019_0320_turkey','2019_0921_albania','2019_1126_albania', 
%                                '2020_0502_crete','2020_1030_samos','2018_0000_neamwave17','2010_0227_maule';
% domain: source for long-term info (source and pre-computed scenarios). Options: 'med-tsumaps','pacific-chile'
disp('----- running: SetPreLoad ------')

addpath(genpath('Routines'))

if strcmpi(eqID,'2010_0227_maule')
    domain = 'cheese-chile'; 
else
    domain = 'med-tsumaps';
end
disp(['Domain: ' domain])

PreSettings = struct;
PreSettings.preProcessingYN = false;
PreSettings.domain = domain;
PreSettings.mainFolder =  pwd;
PreSettings.genInputFolder =  [PreSettings.mainFolder '/LocalInput/' PreSettings.domain '/'];

if strcmpi(domain,'med-tsumaps')
   PreSettings.ProjectName = 'TSUMAPSNEAM';
   PreSettings.selectedIntensityMeasure = 'gl';  

   %% REGIONS 
   % Preset only Region Alfeo, one of the smallest. 
   % If other regions are requried, they will be loaded in runtime   
   PreSelectedRegions = [1];    
   % select the regions to lead for each case
   if strcmp(eqID,'all_med_regions')
      % all the regions defined in NEAMTHM18 (Basili et al. 2021)
      PreSelectedRegions = [1:110]; 
   elseif strcmp(eqID,'all_med_events')
      % all regions for 13 earthquakes in paper (about 30 minutes)
      PreSelectedRegions = [2,3,11,15,23,24,25,30,31,32,44,46,48,49,58,59]; 
   elseif strcmp(eqID,'2003_0521_boumardes')
      PreSelectedRegions = [46,58]; 
   elseif strcmp(eqID,'2015_0416_crete')
      PreSelectedRegions = [24,44,49]; 
   elseif strcmp(eqID,'2015_1117_lefkada')
      PreSelectedRegions = [3,11]; % lefkada
   elseif strcmp(eqID,'2016_0125_gibraltar')
      PreSelectedRegions = [59]; %gibraltar
   elseif strcmp(eqID,'2016_1030_norcia')
      PreSelectedRegions = [2,15,23]; %norcia
   elseif strcmp(eqID,'2017_0612_lesbo')
      PreSelectedRegions = [25,32]; %lesbo
   elseif strcmp(eqID,'2017_0720_kos-bodrum')
      PreSelectedRegions = [31,32,49]; 
   elseif strcmp(eqID,'2018_1025_zante')
      PreSelectedRegions = [3,44,48,49]; % zante
   elseif strcmp(eqID,'2019_0320_turkey')
      PreSelectedRegions = [30]; %turkey_2019  
   elseif strcmp(eqID,'2019_0921_albania')
      PreSelectedRegions = [11]; %albania
   elseif strcmp(eqID,'2019_1126_albania')
      PreSelectedRegions = [11]; %albania
   elseif strcmp(eqID,'2020_0502_crete')
      PreSelectedRegions = [44,48];  
   elseif strcmp(eqID,'2020_1030_samos')
      PreSelectedRegions = [30,31,32]; % Samos 2020
   elseif strcmp(eqID,'2018_0000_neamwave17')
      PreSelectedRegions = [3,44,48,49]; %NEAMWAVE
   end

   PreSettings.SelectedRegions=PreSelectedRegions;

   disp(['Regions to load: ' num2str(PreSelectedRegions)])


% call PTF_preLoad
[LongTermInfo,POIInfo] = PTF_preLoad_Med(PreSettings);


   % visualize selected regions
   if true
      TSUMAPS_plotRegions(LongTermInfo.Regionalization)
      for i=1:LongTermInfo.Regionalization.Npoly
         plot(LongTermInfo.Regionalization.Tlon(i,1:LongTermInfo.Regionalization.Tleng(i)+1),LongTermInfo.Regionalization.Tlat(i,1:LongTermInfo.Regionalization.Tleng(i)+1),':r')
         hold on
         text(mean(LongTermInfo.Regionalization.Tlon(i,1:LongTermInfo.Regionalization.Tleng(i)+1)),mean(LongTermInfo.Regionalization.Tlat(i,1:LongTermInfo.Regionalization.Tleng(i)+1)),['#' num2str(i)])
      end
   end


elseif strcmpi(domain,'cheese-chile')
   PreSettings.ProjectName = 'ChEESE_PTF';
   PreSettings.selectedIntensityMeasure = 'gl-sims';   %'os';'gl';'af';'gl-sims'

  [LongTermInfo,POIInfo] = PTF_preLoad_Chile(PreSettings);


else
   stop
end



end
