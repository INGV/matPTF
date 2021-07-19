function [LongTermInfo,POIInfo] = PTF_preLoad_Med(varargin)
disp('----- running: PTF_PreLoad_Med ------')

PreSettings=varargin{1};
if length(varargin) > 1
   LongTermInfoIn = varargin{2};
   initilizationYN = false;
else
   LongTermInfoIn = struct;
   initilizationYN = true;
end

%% INITIALIZATION
tic
ProjectName=PreSettings.ProjectName;

% CREATE LongTermInfo TO STORE INFO FROM LONG-TERM HAZARD
LongTermInfo = LongTermInfoIn;  % Containing all the information coming from longterm SPTHA
LongTermInfo.ProjectName=ProjectName;
LongTermInfo.vecID = 10.^[8     5     2     0    -2    -4    -6];     % vector for quick search on parameters

% THRESHOLDS FOR HAZARD CURVE
disp('......... reading HCthresholds.txt')
fileThrs = [PreSettings.mainFolder '/LocalInput/HCthresholds.txt'];
LongTermInfo.GenericHazardCurveThresholds = load(fileThrs);   

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT FROM NEAMTHM18 - TSUMAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD FILES AVAILABLE FROM HAZARD
disp('...... loading general input data')

disp('......... reading Regionalization.mat')
load([PreSettings.genInputFolder 'Regionalization.mat']); % INCLUDES REGIONS

disp('......... reading Discretization.mat')
load([PreSettings.genInputFolder 'Discretizations.mat']); % INCLUDES DISCRETIZATIONS

disp('......... reading ModelsWeight.mat')
load([PreSettings.genInputFolder 'ModelsWeight.mat']); % INCLUDES WEIGHTS

% POIs: target points
disp('......... reading POIs.mat')
load([PreSettings.genInputFolder 'POIs.mat']);             
POIInfo = POIs; 
POIInfo.MediterraneanInd = find(POIs.Mediterranean);
POIInfo.MediterraneanIndMed(POIInfo.MediterraneanInd)=1:length(find(POIInfo.MediterraneanInd));
SelectedPOIs = POIs.name(POIs.Mediterranean);
for j = 1:length(SelectedPOIs)
[tmp,iPoi] = ismember(SelectedPOIs{j},POIs.name); if not(tmp); disp(SelectedPOIs{j}); stop; end
POIInfo.SelectedPOI(j).Name = SelectedPOIs{j};
POIInfo.SelectedPOI(j).Index = iPoi;
end

% REGION SELECTION
SelectedRegions=PreSettings.SelectedRegions;
if isempty(SelectedRegions)
SelectedRegions = 1:Regionalization.Npoly;
elseif and(length(SelectedRegions)==1,strcmpi(SelectedRegions,'mediterranean'))
SelectedRegions = find(strcmp(Regionalization.ReferenceCatalog,'EMEC'));
end

%% LOAD PRE-PROCESSED FILES AND CREATE OUTPUT STRUCTURE: LongTermInfo    
LongTermInfo.Regionalization = Regionalization;
LongTermInfo.Discretizations = Discretizations; % DEFINING SCENARIO AVAILABILITY

disp('......... reading PSBarInfo.mat')
load([PreSettings.genInputFolder 'PSBarInfo'],'PSBarInfo'); % DEFINING SCENARIO AVAILABILITY FOR PS (LOCATIONS)
LongTermInfo.PSBarInfo = PSBarInfo;

% Files for MESHES
LongTermInfo.SettingsLambdaBSPS.name_mesh={'mesh/HA_mesh_' ...
'mesh/CA_mesh_' ...
'mesh/Cyprus_mesh_'}';
% Assing to each region to which mesh it owns
% Reg:3 Kefalonia Lefkada
% Reg:10 Calabria
% Reg:16 Calabrian outer wedge
% Reg:24 Hellenic Arc East
% Reg:27 Cyprian Arc South
% Reg:33 Cyprian Arc East
% Reg:35 Cyprian Arc North
% Reg:36 Paphos
% Reg:44 Hellenic Arc South
% Reg:48 Hellenic Arc North
% Reg:49 Aegean South
% Reg:54 Calabrian inner wedge
LongTermInfo.SettingsLambdaBSPS.regionsPerPS = nan(LongTermInfo.Regionalization.Npoly,1);
LongTermInfo.SettingsLambdaBSPS.regionsPerPS([3,24,44,48,49])=1;
LongTermInfo.SettingsLambdaBSPS.regionsPerPS([10,16,54])=2;
LongTermInfo.SettingsLambdaBSPS.regionsPerPS([27,33,35,36])=3;
LongTermInfo.SettingsLambdaBSPS.regionsPerPS([27,33,35,36])=3;    
LongTermInfo.selectedIntensityMeasure = 'gl';        

% STORE SELECTIONS
LongTermInfo.SelectedRegions=SelectedRegions;
LongTermInfo.ModelsWeight = ModelsWeight;
     
    if not(isfield(LongTermInfo,'region'));
        LongTermInfo.region = cell(Regionalization.Npoly,1);
    end
    for i=1:length(SelectedRegions)
        iReg = SelectedRegions(i);
        disp(['......... loading Region #' num2str(iReg) ' (' num2str(i) '/' num2str(length(SelectedRegions)) ')']);
        
        RegionalInfoTmp = struct;
        disp(['............ reading MeanProb_BS4_FocMech_Reg' num2str(iReg,'%03d') '.mat'])
        tmp = load([PreSettings.genInputFolder 'MeanProb_BS4_FocMech_Reg' num2str(iReg,'%03d')],'MeanProb_BS4_FocMech_Reg');
        RegionalInfoTmp = catstruct(RegionalInfoTmp,tmp);
        
        %% READ SCENARIO LISTS
        if sum(Regionalization.Ttypes{iReg}==1)>0
            disp(['............ reading ScenarioListBS_Reg' num2str(iReg,'%03.f') '.mat'])
            fileNamePar = [PreSettings.genInputFolder 'ScenarioListBS_Reg' num2str(iReg,'%03.f') '_' Regionalization.ID{iReg} '_Parameters'];
            tmp=load(fileNamePar);
            tmp2.ScenarioListBSReg.Parameters = tmp.Parameters;
            RegionalInfoTmp = catstruct(RegionalInfoTmp,tmp2);
            
        end
        if sum(Regionalization.Ttypes{iReg}==2)>0
            disp(['............ reading ScenarioListPS_Reg' num2str(iReg,'%03.f') '.mat'])
            fileName = [PreSettings.genInputFolder 'ScenarioListPS_Reg' num2str(iReg,'%03.f') '_' Regionalization.ID{iReg}];
            tmp = load(fileName);
            RegionalInfoTmp = catstruct(RegionalInfoTmp,tmp);
        end                       
                        
        %% READ PRECOMPUTED SCEARIOS
        disp(['............ reading glVal_BS_Reg', num2str(iReg,'%03.f'),'-',LongTermInfo.Regionalization.ID{iReg},'.mat'])
        fileName = strcat(PreSettings.genInputFolder(1:end-1),'/glVal_BS_Reg', num2str(iReg,'%03.f'),'-',LongTermInfo.Regionalization.ID{iReg},'.mat');
        load(fileName);
        RegionalInfoTmp.glVal_BS =  glVal_BS;

        if not(isempty(find(LongTermInfo.Regionalization.Ttypes{iReg}==2)))
           disp(['............ reading glVal_PS_Reg', num2str(iReg,'%03.f'),'-',LongTermInfo.Regionalization.ID{iReg},'.mat'])
           fileName = strcat(PreSettings.genInputFolder(1:end-1),'/glVal_PS_Reg', num2str(iReg,'%03.f'),'-',LongTermInfo.Regionalization.ID{iReg},'.mat');
           load(fileName);
           RegionalInfoTmp.glVal_PS =  glVal_PS;
        end                    
        
        LongTermInfo.region{iReg} = RegionalInfoTmp;
    end

toc
end
