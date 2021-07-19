README.TXT file for code matPTF, version 1.0
date: 16/07/2021
citation: ###

##################
#### RUN PTF #####
##################
This repository includes all the source codes used to run Probabilistic Tsunami Forecast (PTF) for the case studies presented in Selva et al. (2021). To install the required code, follow the instructions in section INSTALLATION.

The code has been developed and tested in MATLAB ver. R2019b (The MathWorks Inc., Natick, Massachusetts, 2019), running in an HP ProLiant DL580 GEN9 four 14-core Intel(R) Xeon(R) with E7-4830 CPUs clocked at 2.0 GHz (56 total compute cores and 3072 GB RAM).

The entire workflow can be run with the script runPTF.m. A complete run of this source code includes: 
1) a pre-processing phase in which the required input is loaded
2) a running phase in which PTF calculations are performed and stored
3) a post-processing phase in which figures are produced.

These three phases are described in the corresponding sections. 


﻿###############################
#### PRE-PROCESSING PHASE #####
﻿﻿###############################
The pre-processing phase consists of defining the main input, and in preloading the required input information.

The main input is defined by two variables to be saved in Matlab workspace. Two variables are required:
1) eqID: unique ID assigned to the earthquake. This should be a text string, and should correspond to a specific input file in folder "LocalInput/EarlyEst". Prefined options are: '2003_0521_boumardes', '2015_0416_crete', '2015_1117_lefkada', '2016_0125_gibraltar', '2016_1030_norcia', '2017_0612_lesbo', '2017_0720_kos-bodrum', '2018_1025_zante','2019_0320_turkey','2019_0921_albania','2019_1126_albania', '2020_0502_crete', '2020_1030_samos','2018_0000_neamwave17', '2010_0227_maule';
2) sigmaCutoff: select the sigma for cutoff for computations. More details about this parameters in Selva et al. (2021).

Preload is required to speed-up PTF computations. Preload is initiated by the script SetPreLoad.m and the preloaded information is stored in the Matlab structures PreSettings, LongTermInfo, and POIInfo. The preload simply loads in memory long-term information (list of scenarios, prob of focal mechanisms, precomputed tsunamis). If preloads is not performed, these data are loaded during the execution of PTF_ncomms.m, slowing its first execution. The following executions (without clearing Matlab memory), will make use of the data loaded during previsous executions.

To preload the data relative to one scenario, use:
>> [PreSettings,LongTermInfo,POIInfo] = SetPreLoad(eqID); 

To preload the data relative to all the mediterranean scenarios of Selva et. al 2021, use:
>> [PreSettings,LongTermInfo,POIInfo] = SetPreLoad('all_med_events'); 

To preload the data relative to the entire Mediterranean region, use:
>> [PreSettings,LongTermInfo,POIInfo] = SetPreLoad('all_med_regions'); 



﻿###############################
#### RUNNING PHASE #####
﻿###############################
The running phase consists of
1) load information from the seismic monitoring;
2) perform the probabilistic forecast
3) evaluate alert levels
4) store the results

The entire workflow is initiated by the script runPTF.m.

The loading of the seismic monitoring data is initiated by the script SetEarlyEst.m, and the loaded information is saved in the Matlab structure EarlyEst. 
Example:
>> [EarlyEst] = SetEarlyEst(eqID);


The quantification of the probabilistic forecast is initiated by the script PTF_1_0.m, and the output information is saved in the Matlab structures Settings, HazardCurves, OutTesting, ScenarioProb, ShortTermProbDistr, PreSelection, LongTermInfo, EarlyEst, and UsedFunctions. The main probabilistic results are reported in HazardCurves, while the main information about the ensemble in ScenarioProb.
Example:
>> [Settings,HazardCurves,OutTesting,ScenarioProb,ShortTermProbDistr,PreSelection,LongTermInfo,EarlyEst,UsedFunctions] = ...
PTF_1_0(PreSettings,EarlyEst,LongTermInfo,POIInfo);

The alert levels estimation is initiated by the script PTF_AlertLevels.m, and the output information is saved in the Matlab structure AlertLevelsInfo.
Example:
>> [AlertLevelsInfo] = PTF_AlertLevels(Settings,HazardCurves,POIInfo,EarlyEst);

The saving phase is performed with Matlab instructions. In script runPTF.m, the only scructures HazardCurves, AlertLevelsInfo, and OutTesting are stored.
Example:
>> save(['Output/' 'HazardCurves_' EarlyEst.ID '_sig' num2str(10*Settings.nSigma)],'HazardCurves')


﻿################################
#### POST-PROCESSING PHASE #####
﻿################################
The post-processing phase consists of producing figures out of PTF results. Figures are produced to represent the ensemble and/or the tsunami forecast results for custom probability thresholds). In script runPTF.m, we reported several examples for producing such figures. 

Examples:
>> PTF_plotMarginals(LongTermInfo,ScenarioProb,ShortTermProbDistr,EarlyEst,POIInfo,UsedFunctions);
>> PTF_plotHazMaps(0.5,HazardCurves,EarlyEst,POIInfo,'hazmean',PreSettings); % hazard map for a probability threshold of 0.5 (median)
>> PTF_plotAlertLevels(AlertLevelsInfo,EarlyEst,POIInfo,PreSettings);

#######################
#### INSTALLATION #####
#######################
All the input files for the scripts in matPTF are included in folder “LocalInput/”. The necessary configuration files are included in Github. 

To run the case studies in the Mediterranean, the input files from TSUMAPS-NEAM project should be additionally downloaded. These data should in included in folder “LocalInput/med-tsumaps”, and can be downloaded from TSUMAPS-NEAM documentation website (http://www.tsumaps-neam.eu/documentation/). These data allow running all the case studies relative to Mediterranean earthquakes, as discussed in Selva et al. (2021). The same data can be used to run whatever else real or hypothetical earthquake originating within the Medirranean area. To do so, follow the instruction in section "MEDITERRANEAN EARTHQUAKES". input source files, 

To run the case studies in the M8.8 Maule case studies, the relative input files should be additionally downloaded. These data should in included in folder “LocalInput/cheese-chile”, and can be downloaded from the public repository figshare (https://figshare.com/). These data allow running the case study relative to the 2010 M8.8 Maule earthquake, as discussed in Selva et al. (2021).

WARNING:
Both Mediterranean and Chilean case studies may be installed in the same folder. However, if the domain has to be changed (from Med to Chile, or viceversa), Matlab memory should be cleared to force a new preLoad phase.


####################################
#### MEDITERRANEAN EARTHQUAKES #####
####################################
Tsumpas-NEAM data allow running the Probabilistic Tsunami Forecast workflow for whatever seismic events occurring within the Mediterranean area. To do so, the Tsumaps-NEAM data should be downloaded and correctly saved (see Section INSTALLATION). The workflow has been tested only the for events listed in Selva et al. (2021). 

To run a different case study, it is required:
1) to define a Eq_ID to be used to run the codes;
2) to create file named "Eq_ID_stat.txt" to define the uncertainty information from the seismic monitoring, which should be placed in folder "LocalInput/EarlyEst". Examples of such files are available through this Github. New files should be formatted in agreement with such files.
3) to run the code, (see RUN PTF)

