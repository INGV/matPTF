%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SET INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select event.   Options: '2003_0521_boumardes', '2015_0416_crete', '2015_1117_lefkada', 
%                          '2016_0125_gibraltar', '2016_1030_norcia','2017_0612_lesbo','2017_0720_kos-bodrum',
%                          '2018_1025_zante','2019_0320_turkey','2019_0921_albania','2019_1126_albania', 
%                          '2020_0502_crete','2020_1030_samos','2018_0000_neamwave17','2010_0227_maule';
eqID = '2018_0000_neamwave17';
% select sigma for cutoff
sigmaCutoff = 2.0; 

clc
disp('----- running: runPTF ------')
disp('****************************')
disp(['* eqID: ' eqID])
disp('****************************')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PRELOAD INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preload is required to speed-up PTF computations
% The preload simply load in memory long-term information (list of scenarios, prob of focal mechanisms, precomputed tsunamis)
% If preloads is not perfomed, these data are loaded during the execution of PTF_ncomms.m, slowing its first execution.
% To preload the data relative to one scenario, use:
%  [PreSettings,LongTermInfo,POIInfo] = SetPreLoad(eqID); 
% To preload the data relative to all the mediterranean scenarios of Selva et. al 2021, use:
%    [PreSettings,LongTermInfo,POIInfo] = SetPreLoad('all_med_events'); 
% To preload the data relative to all NEAM regions, use:
%    [PreSettings,LongTermInfo,POIInfo] = SetPreLoad('all_med_regions'); 
if not(exist('LongTermInfo','var'))
    [PreSettings,LongTermInfo,POIInfo] = SetPreLoad(eqID); 
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RUN PTF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PreSettings.nSigmasigmaCO = sigmaCutoff;

% initialize time counter
tic

% load real-time source uncertainty
[EarlyEst] = SetEarlyEst(eqID);

% run PTF
[Settings,HazardCurves,OutTesting,ScenarioProb,ShortTermProbDistr,PreSelection,LongTermInfo,EarlyEst,UsedFunctions] = ...
PTF_1_0(PreSettings,EarlyEst,LongTermInfo,POIInfo);

% compute alert levels out of probabilistic forecast
[AlertLevelsInfo] = PTF_AlertLevels(Settings,HazardCurves,POIInfo,EarlyEst);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  STORING RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results in mat files
save(['Output/' 'HazardCurves_' EarlyEst.ID '_sig' num2str(10*Settings.nSigma)],'HazardCurves')
save(['Output/' 'AlertLevelsInfo_' EarlyEst.ID '_sig' num2str(10*Settings.nSigma)],'AlertLevelsInfo')
save(['Output/' 'OutTesting_' EarlyEst.ID '_sig' num2str(10*Settings.nSigma)],'OutTesting')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PTF_plotMarginals(LongTermInfo,ScenarioProb,ShortTermProbDistr,EarlyEst,POIInfo,UsedFunctions);
if ShortTermProbDistr.BScomputedYN; figures.marginalsBS = PTF_plotMarginalsBS(LongTermInfo,ScenarioProb,ShortTermProbDistr,PreSelection,EarlyEst,POIInfo,UsedFunctions); end
if ShortTermProbDistr.PScomputedYN; figures.marginalsPS = PTF_plotMarginalsPS(LongTermInfo,ScenarioProb,ShortTermProbDistr,PreSelection,EarlyEst,POIInfo,UsedFunctions); end
PTF_plotHazMaps(0.5,HazardCurves,EarlyEst,POIInfo,'hazmean',PreSettings); % hazard map for a probability threshold of 0.5 (median)
PTF_plotHazMaps(0.05,HazardCurves,EarlyEst,POIInfo,'haz',PreSettings);  % hazard map for a probability threshold of 0.05 (5th percentile)
PTF_plotAlertLevels(AlertLevelsInfo,EarlyEst,POIInfo,PreSettings);


