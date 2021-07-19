function [LongTermInfo,POIInfo] = PTF_preLoad_Chile(varargin)
disp('----- running: PTF_PreLoad_Chile ------')


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
% LOAD SIMULATIONS FROM CHEESE PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
disp('...... loading general input data')

disp('......... reading Regionalization.mat')
load([PreSettings.genInputFolder 'Regionalization.mat']); % INCLUDES REGIONS

disp('......... reading Discretization.mat')
load([PreSettings.genInputFolder 'Discretizations.mat']); % INCLUDES DISCRETIZATIONS

disp('......... reading BS3_Depth_ValVec.mat')
load([PreSettings.genInputFolder 'BS3_Depth_ValVec.mat']);
Discretizations.BS3_Depth.ValVec = ValVec;

disp('......... reading PSregions.mat')
load([PreSettings.genInputFolder 'PSregions.mat']); % LOAD subduction names

% POIs: target points
disp('......... reading POIs_chile.mat')
load([PreSettings.genInputFolder 'POIs_chile.mat']); % LOAD POIs
POIInfo = POIs; 
SelectedPOIs = POIs.name;
for j = 1:length(SelectedPOIs)
   [tmp,iPoi] = ismember(SelectedPOIs{j},POIs.name); if not(tmp); disp(SelectedPOIs{j}); stop; end
   POIInfo.SelectedPOI(j).Name = SelectedPOIs{j};
   POIInfo.SelectedPOI(j).Index = iPoi;
end    

% LongTermInfo
LongTermInfo.Regionalization = Regionalization;
LongTermInfo.Discretizations = Discretizations;

disp('......... reading PSBarInfo.mat')
load([PreSettings.genInputFolder 'PSBarInfo'],'PSBarInfo');
LongTermInfo.PSBarInfo = PSBarInfo;

% Files for MESHES
LongTermInfo.SettingsLambdaBSPS.name_mesh={[PreSettings.genInputFolder 'mesh/south-america_mesh_']};

slip_magthr = 7.32;
slip_nslip = 5;

LongTermInfo.SettingsLambdaBSPS.regionsPerPS = ones(1,length(PSregions.Name));

disp('......... reading Probabilities_angles.mat')
load([PreSettings.genInputFolder 'Probabilities_angles.mat'])
LongTermInfo.region =  cell(length(PSregions.Name),1);
    for i=1:length(PSregions.Name)
        k=0;
        kold=1;
        for imag=1:length(LongTermInfo.Discretizations.PS1_Mag.ID)
            if LongTermInfo.Discretizations.PS1_Mag.Val(imag) < slip_magthr
                k=k+1;
                LongTermInfo.region{i}.ScenarioListPSReg.SlipDistribution(kold:k)=0;
            else
                k=k+5;
                LongTermInfo.region{i}.ScenarioListPSReg.SlipDistribution(kold:k)=1:slip_nslip;
            end
            LongTermInfo.region{i}.ScenarioListPSReg.magPSInd(kold:k)=imag;
            kold=k+1;
        end
    end
    LongTermInfo.region{1}.MeanProb_BS4_FocMech_Reg.indPos = 1:length(LongTermInfo.Discretizations.BS2_Pos.Region);
    LongTermInfo.region{1}.MeanProb_BS4_FocMech_Reg.Val = MokN;
    
    
    % computed in prepPSBarInfo.m, MatPTF_PreProcessedFiles
    LongTermInfo.ModelsWeight.PS2_Bar.Wei = PSBarInfo.PSModelWeights;
    
    LongTermInfo.selectedIntensityMeasure = 'gl-sims';   
    % these scenarios are compute ad hoc for the selected events
    % the tsunami results are loaded directly in PTF_ncomms.m
    
    
toc
end
