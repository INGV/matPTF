function [Settings,HazardCurves,OutTesting,ScenarioProb,ShortTermProbDistr,PreSelection,LongTermInfo,EarlyEst,UsedFunctions] = PTF_1_0(PreSettings,EarlyEst,LongTermInfo,POIInfo)
% Quantification of the probabilistic forecast
disp('----- running: PTF_1_0 ------')

%% %%%%%%%%%%%%%
%%% SETTINGS
%%%%%%%%%%%%%%%%
%Settings = struct;
if isfield(PreSettings,'nSigma')
   Settings.nSigma = PreSettings.nSigma;
else
   Settings.nSigma = 2;         % FOR PRE-SELECTION OF POSITION OF INTERESTS
end
disp(['N sigma: ' num2str(Settings.nSigma)]);
Settings.normalizeAngleProb = false; % force normalization of probability of focal mechanisms
Settings.NegligibleProb = 2*normcdf(-Settings.nSigma);  % DEFINE A LOWER LIMIT TO AVOID COMPUTING NOT USEFULL SCENARIOS
Settings.gauss3Dyn = true;   % EVALUATE PROB X,Y AND Z JOINTLY (3D) OR SEPARATED (2D)
Settings.SpaceBin = 2.5E3;  % GRID SIZE IN M FOR 2D OR 3D INTEGRAL (ALONG Z)
Settings.Z2XYfact = 2.5;     % RATION BETWEarlyEst.N GRID SIZE ALONG X,Y AND Z (IF > 1, FINER ON Z)
Settings.magBSmax = 8.1;     % MAXIMUM MODELLED MAGNITUDE FOR BS
Settings.selectedIntensityMeasure = PreSettings.selectedIntensityMeasure; % TO SELECT AMONG AF, GL AND OS
Settings.lambdaFabrizioYN = true;
Settings.lambdaFabruzioPriorYN = false;
% OTHER GENERAL SETTINGS
Settings.commons = 'Routines'; addpath(Settings.commons);
Settings.octaveYN = false;
if Settings.octaveYN
    pkg load statistics
    graphics_toolkit('gnuplot');
end
Settings.writeOutTesting = true;
Settings.verboseYN = false;
Settings.figcheckYN = false;
Settings

%% %%%%%%%%%%%%%
%%% CONVERSION 2 STANDARDIZED UTM
%%%%%%%%%%%%%%%%

% SET AS REFERENCE UTM ZONE THE ONE OF THE EPICENTER
[EarlyEst.lonUTM,EarlyEst.latUTM,EarlyEst.RefUtmZone]=ll2utm(EarlyEst.lat,EarlyEst.lon);
if EarlyEst.RefUtmZone < 0
    EarlyEst.RefUtmZone=-EarlyEst.RefUtmZone;
    [EarlyEst.lonUTM,EarlyEst.latUTM]=ll2utm(EarlyEst.lat,EarlyEst.lon,EarlyEst.RefUtmZone);
end

% PTHA GRIDS
DistAll = distance( LongTermInfo.Discretizations.BS2_Pos.Val(:,2), LongTermInfo.Discretizations.BS2_Pos.Val(:,1),EarlyEst.lat,EarlyEst.lon )*pi*6371/180 ;
PreSelectBS = find(DistAll < 1000);
LongTermInfo.Discretizations.BS2_Pos.ValUTM = nan(size(LongTermInfo.Discretizations.BS2_Pos.Val));
[x y] = ll2utm(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelectBS,2),LongTermInfo.Discretizations.BS2_Pos.Val(PreSelectBS,1),EarlyEst.RefUtmZone);
LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelectBS,:) = [x y];

DistAll = distance( LongTermInfo.Discretizations.PS2_Bar.Val(:,2), LongTermInfo.Discretizations.PS2_Bar.Val(:,1),EarlyEst.lat,EarlyEst.lon )*pi*6371/180 ;
PreSelectPS = find(DistAll < 1000);
LongTermInfo.Discretizations.PS2_Bar.ValUTM = nan(size(LongTermInfo.Discretizations.BS2_Pos.Val));
[x y] = ll2utm(LongTermInfo.Discretizations.PS2_Bar.Val(PreSelectPS,2),LongTermInfo.Discretizations.PS2_Bar.Val(PreSelectPS,1),EarlyEst.RefUtmZone);
LongTermInfo.Discretizations.PS2_Bar.ValUTM = [x y];
[x y] = ll2utm(POIInfo.lat,POIInfo.lon,EarlyEst.RefUtmZone);
POIInfo.lonUTM = x; POIInfo.latUTM = y;

LongTermInfo.PSBarInfo.BarPSperModelUTM = cell(size(LongTermInfo.PSBarInfo.BarPSperModel));
for i1=1:size(LongTermInfo.PSBarInfo.BarPSperModel,1)
    for i2=1:size(LongTermInfo.PSBarInfo.BarPSperModel,2)
        if isempty(LongTermInfo.PSBarInfo.BarPSperModel{i1,i2})
            LongTermInfo.PSBarInfo.BarPSperModelUTM{i1,i2}=[];
        else
            DistAll = distance( LongTermInfo.PSBarInfo.BarPSperModel{i1,i2}(:,2), LongTermInfo.PSBarInfo.BarPSperModel{i1,i2}(:,1),EarlyEst.lat,EarlyEst.lon )*pi*6371/180 ;
            PreSelectTmp = find(DistAll < 1000);
            LongTermInfo.PSBarInfo.BarPSperModelUTM{i1,i2} = nan(size(LongTermInfo.PSBarInfo.BarPSperModel{i1,i2}));
            [x y] = ll2utm(LongTermInfo.PSBarInfo.BarPSperModel{i1,i2}(PreSelectTmp,2),LongTermInfo.PSBarInfo.BarPSperModel{i1,i2}(PreSelectTmp,1),EarlyEst.RefUtmZone);
            LongTermInfo.PSBarInfo.BarPSperModelUTM{i1,i2}(PreSelectTmp,:)=[x y];
        end
    end
end

if Settings.figcheckYN
    figure()
    subplot(2,1,1)
    plot(LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,1),LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,2),'.k')
    hold on
    plot(LongTermInfo.Discretizations.PS2_Bar.ValUTM(:,1),LongTermInfo.Discretizations.PS2_Bar.ValUTM(:,2),'oc')
    plot(EarlyEst.lonUTM,EarlyEst.latUTM,'r*')
    subplot(2,1,2)
    plot(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelectBS,1),LongTermInfo.Discretizations.BS2_Pos.Val(PreSelectBS,2),'.k')
    hold on
    plot(LongTermInfo.Discretizations.PS2_Bar.Val(PreSelectPS,1),LongTermInfo.Discretizations.PS2_Bar.Val(PreSelectPS,2),'oc')
    plot(EarlyEst.lon,EarlyEst.lat,'r*')
    set(gca,'xlim',[round(EarlyEst.lon-10),round(EarlyEst.lon+10)],'ylim',[round(EarlyEst.lat-10),round(EarlyEst.lat+10)])
end

% SET BASIC MEAN & COVARIANCE MATRIX WITHOUT FAULT DIMENSION
EarlyEst.PosCovMat_2D = [EarlyEst.PosSigmaXX EarlyEst.PosSigmaXY; EarlyEst.PosSigmaXY EarlyEst.PosSigmaYY];
EarlyEst.PosMean_2D = [EarlyEst.lonUTM,EarlyEst.latUTM];
EarlyEst.PosCovMat_3D = [EarlyEst.PosSigmaXX EarlyEst.PosSigmaXY EarlyEst.PosSigmaXZ; ...
    EarlyEst.PosSigmaXY EarlyEst.PosSigmaYY EarlyEst.PosSigmaYZ; ...
    EarlyEst.PosSigmaXZ EarlyEst.PosSigmaYZ EarlyEst.PosSigmaZZ];
EarlyEst.PosMean_3D = [EarlyEst.lonUTM,EarlyEst.latUTM,EarlyEst.Dep*1.E3];

clear x y i1 i2

%% %%%%%%%%%%%%%
%%% PREP VARIABLES FOR LAMBDA_BSPS
%%%%%%%%%%%%%%%%
SettingsLambdaBSPS = struct;
SettingsLambdaBSPS.hypo_utm = [EarlyEst.lonUTM,EarlyEst.latUTM,EarlyEst.Dep];
SettingsLambdaBSPS.utmzone_hypo = EarlyEst.RefUtmZone;
SettingsLambdaBSPS.NormCov = EarlyEst.PosCovMat_3D; % covariance matrix (sxx sxy sxz; syx syy syz; szx szy szz) in m^2
SettingsLambdaBSPS.moho_pts=struct('moho_par',{[LongTermInfo.Discretizations.BS2_Pos.Val(:,1),LongTermInfo.Discretizations.BS2_Pos.Val(:,2),LongTermInfo.Discretizations.BS2_Pos.DepthMoho'],LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,1),LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,2),LongTermInfo.Discretizations.BS2_Pos.DepthMoho'});

% confid_lev should be the same as the rest of PTF
SettingsLambdaBSPS.confid_lev= (normcdf(Settings.nSigma)-normcdf(-Settings.nSigma));
SettingsLambdaBSPS.dchi2=chi2inv(SettingsLambdaBSPS.confid_lev,3); % delta-chi squared for conf. lev. = 68.3; can be obtained by chi2inv(0.683,3)
SettingsLambdaBSPS.SD=sqrt(SettingsLambdaBSPS.dchi2);

% the number of volumes depends only on Npts
tmp = load('LocalInput/Lambda_BSPS/tetra_scaling.mat'); SettingsLambdaBSPS=catstruct(SettingsLambdaBSPS,tmp);
SettingsLambdaBSPS.Nvref_vref_ratio=SettingsLambdaBSPS.Ntetra(1,4)./SettingsLambdaBSPS.Vol_tetra(1,4);

if (license('test', 'curve_fitting_toolbox'))
    tmp = load('LocalInput/Lambda_BSPS/Npts_vs_Ntetra.mat'); SettingsLambdaBSPS=catstruct(SettingsLambdaBSPS,tmp);
else
    tmp = load('LocalInput/Lambda_BSPS/Npts_vs_Ntetra_ab.mat'); SettingsLambdaBSPS=catstruct(SettingsLambdaBSPS,tmp);
end
clear tmp

% Settings for Subduction Zones Meshes
SettingsLambdaBSPS.name_mesh= LongTermInfo.SettingsLambdaBSPS.name_mesh;
SettingsLambdaBSPS.regionsPerPS = LongTermInfo.SettingsLambdaBSPS.regionsPerPS;
SettingsLambdaBSPS.Nz=size(SettingsLambdaBSPS.name_mesh,1);
for i=1:SettingsLambdaBSPS.Nz
    if length(SettingsLambdaBSPS.name_mesh{i,1}) > 0
       mesh_faces_file{i,1}=[SettingsLambdaBSPS.name_mesh{i,1} 'faces_x16.dat'];
       mesh_nodes_file{i,1}=[SettingsLambdaBSPS.name_mesh{i,1} 'nodes_x16.dat'];
    
       mnodes{i,1}=load(mesh_nodes_file{i,1});
       mfaces{i,1}=load(mesh_faces_file{i,1});
       [meshx{i,1},meshy{i,1}]=ll2utm(mnodes{i,1}(:,3),mnodes{i,1}(:,2),SettingsLambdaBSPS.utmzone_hypo);
    else
       mnodes{i,1}=[];
       mfaces{i,1}=[];
       meshx{i,1}=[];
       meshy{i,1}=[];
    end
end
field='mesh_par';
value={mesh_faces_file,mesh_nodes_file,mnodes,mfaces,meshx,meshy};
SettingsLambdaBSPS.mesh_subd=struct(field,value);
% Tetrahedra Volumes and Baricenters as output [0/1 --> NO/YES]
SettingsLambdaBSPS.out_tetra_vol=1; 

clear i field mesh_faces_file mesh_nodes_file mnodes mfaces meshx meshy value

%% %%%%%%%%%%%%%
%%% PREDEFINED FUNCTIONS
%%%%%%%%%%%%%%%%
% MULTIVARIATE GAUSSIAN
NormMultiDvec =@(x,mu,Sigma)  (2*pi)^(-length(mu)/2) * (1/sqrt(det(Sigma))) * exp(-.5*diag(((x-repmat(mu,1,size(x,2)))'/Sigma)*(x-repmat(mu,1,size(x,2)))));

% SCALING LAWS, in KM
UsedFunctions=struct;
UsedFunctions.FuncMag2W_BS =@(mag)  1.e3*scalinglaw_WC('M2W',mag);
UsedFunctions.FuncMag2W_PS_Mo =@(mag) 1.e3*scalinglaw_Murotani('M2W',mag);
UsedFunctions.FuncMag2W_PS_St =@(mag) 1.e3*scalinglaw_Murotani('M2W',mag);
UsedFunctions.FuncMag2L_BS =@(mag) 1.e3*scalinglaw_WC('M2L',mag);
UsedFunctions.FuncMag2L_PS_Mo =@(mag) 1.e3*scalinglaw_Murotani('M2L',mag);
UsedFunctions.FuncMag2L_PS_St =@(mag) 1.e3*scalinglaw_Murotani('M2L',mag);
UsedFunctions.FuncMag2L_PS =@(mag) UsedFunctions.FuncMag2L_PS_Mo(mag);  % not differenciated
UsedFunctions.FuncMag2W_PS =@(mag) UsedFunctions.FuncMag2W_PS_Mo(mag);  % not differenciated
UsedFunctions.CorrectBShorizontalPos = @(mag) 0.5*UsedFunctions.FuncMag2L_BS(mag);
UsedFunctions.CorrectBSverticalPos = @(mag) sin(pi/4)*0.5*UsedFunctions.FuncMag2W_BS(mag);
UsedFunctions.CorrectPShorizontalPos = @(mag) 0.5*UsedFunctions.FuncMag2L_PS(mag);
UsedFunctions.CorrectPSverticalPos = @(mag) sin(pi/4)*0.5*UsedFunctions.FuncMag2W_PS(mag);

%% %%%%%%%%%%%%%%
%%% PRE-SELECT SCENARIOS
%%% Output: PosBSSel, PosBSSelOut, tmpR_PS_nsigma
%%%%%%%%%%%%%%%%
PreSelection = struct;
% PRE-SELECTION OF MAGNITUDES FOR PS & BS
PreSelection.MagSel=find(and(LongTermInfo.Discretizations.PS1_Mag.Val > EarlyEst.Mag - Settings.nSigma * EarlyEst.MagSigma , LongTermInfo.Discretizations.PS1_Mag.Val < EarlyEst.Mag+Settings.nSigma*EarlyEst.MagSigma));
% PRE-SELECT MAGNITUDE RANGE FOR BS
PreSelection.MagSelBS = PreSelection.MagSel(LongTermInfo.Discretizations.PS1_Mag.Val(PreSelection.MagSel)<Settings.magBSmax);

% PRE-SELECTION OF POS FOR BS
%    COMPUTE LARGEST SIGMA IN LON AND LAT
PosSigmaLon = sqrt(EarlyEst.PosSigmaXX) + UsedFunctions.CorrectBShorizontalPos(min(Settings.magBSmax,EarlyEst.Mag+Settings.nSigma*EarlyEst.MagSigma));
PosSigmaLat = sqrt(EarlyEst.PosSigmaYY) + UsedFunctions.CorrectBShorizontalPos(min(Settings.magBSmax,EarlyEst.Mag+Settings.nSigma*EarlyEst.MagSigma));
%    SELECT ALL POS WITHIN SELECTED SIGMA
PreSelection.tmpR_BS_nsigma=ellipsedata([PosSigmaLat^2 0; 0 PosSigmaLon^2],[EarlyEst.latUTM,EarlyEst.lonUTM],1000,Settings.nSigma);
PreSelection.PosBSSel=find(inpolygon(LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,1),LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,2),PreSelection.tmpR_BS_nsigma(:,2),PreSelection.tmpR_BS_nsigma(:,1)));
%    SELECT ALL POS WITHIN SELECTED SIGMA + 1
%    --> LARGER AREA TO COMPUTE INTEGRAL WITHOUT BIASES (BORDER EFFECTS)
PreSelection.tmpR_BS_nsigma_out=ellipsedata([PosSigmaLat^2 0; 0 PosSigmaLon^2],[EarlyEst.latUTM,EarlyEst.lonUTM],1000,Settings.nSigma+0.5);
PreSelection.PosBSSelOut=find(inpolygon(LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,1),LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,2),PreSelection.tmpR_BS_nsigma_out(:,2),PreSelection.tmpR_BS_nsigma_out(:,1))); %% LARGER AREA TO COMPUTE INTEGRAL WITHOUT BIASES (BORDER EFFECTS)
% MAP SMALLER IN LARGER POSITIONS
[PreSelection.PosBSSelOutYN PreSelection.PosBSSelOutInd] = ismember(PreSelection.PosBSSelOut,PreSelection.PosBSSel);

% PRE-SELECTION OF BAR FOR PS
PreSelection.nModels = size(LongTermInfo.PSBarInfo.BarPSperModel,2);
PosSigmaLon = sqrt(EarlyEst.PosSigmaXX) + UsedFunctions.CorrectPShorizontalPos(EarlyEst.Mag+Settings.nSigma*EarlyEst.MagSigma);
PosSigmaLat = sqrt(EarlyEst.PosSigmaYY) + UsedFunctions.CorrectPShorizontalPos(EarlyEst.Mag+Settings.nSigma*EarlyEst.MagSigma);
PreSelection.tmpR_PS_nsigma=ellipsedata([PosSigmaLat^2 0; 0 PosSigmaLon^2],[EarlyEst.latUTM,EarlyEst.lonUTM],1000,Settings.nSigma);
PreSelection.BarPSSelPerModMag = cell(length(PreSelection.MagSel),PreSelection.nModels);
PreSelection.MagExistYN = true(length(PreSelection.MagSel),1);
for i1=1:length(PreSelection.MagSel)
    imag=PreSelection.MagSel(i1);
    MagExistYN = true(PreSelection.nModels,1);
    for imod=1:PreSelection.nModels
        if isempty(LongTermInfo.PSBarInfo.BarPSperModelUTM{imag,imod})
            MagExistYN(imod) = false;
            PreSelection.BarPSSelPerModMag{imag,imod} = [];
        else
            PreSelection.BarPSSelPerModMag{imag,imod} = ...
                find(inpolygon(LongTermInfo.PSBarInfo.BarPSperModelUTM{imag,imod}(:,1),LongTermInfo.PSBarInfo.BarPSperModelUTM{imag,imod}(:,2),PreSelection.tmpR_PS_nsigma(:,2),PreSelection.tmpR_PS_nsigma(:,1)));
        end
    end
    PreSelection.MagExistYN(imag)=max(MagExistYN);
end
% UPDATE MAGNITUDE RANGE FOR PRESELECTED PS
PreSelection.MagSelPS = PreSelection.MagSel(PreSelection.MagExistYN(PreSelection.MagSel));

% PLOTS FOR CHECK
if Settings.figcheckYN
    % FOR BS:
    figure()
    isel = not(isnan(LongTermInfo.Discretizations.BS2_Pos.ValUTM(:,1)));
    plot(LongTermInfo.Discretizations.BS2_Pos.Val(isel,1),LongTermInfo.Discretizations.BS2_Pos.Val(isel,2),'.k','markersize',1)
    hold on
    plot(LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,1),LongTermInfo.Discretizations.BS2_Pos.Val(PreSelection.PosBSSel,2),'b.')
    plot(EarlyEst.lon,EarlyEst.lat,'r*')
    % FOR PS:
    for i=1:length(PreSelection.MagSel)
        imag=PreSelection.MagSel(i);
        for imod=1:PreSelection.nModels
            isel=PreSelection.BarPSSelPerModMag{imag,imod};
            
            plot(LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(isel,1),LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(isel,2),'c.')
        end
    end
    set(gca,'xlim',[round(EarlyEst.lon-10),round(EarlyEst.lon+10)],'ylim',[round(EarlyEst.lat-10),round(EarlyEst.lat+10)])
end

clear PosSigmaLat PosSigmaLon tmpR_BS_nsigma tmpR_BS_nsigma



%% %%%%%%%%%%%%%%
%%% COMPUTE PROB DISTR
%%%    Equivalent of shortterm.py with output: node_st_probabilities
%%%    Output: EarlyEst.MagProb, EarlyEst.PosProb, EarlyEst.DepProb, EarlyEst.DepProb, EarlyEst.BarProb, EarlyEst.RatioBSonTot
%%%%%%%%%%%%%%%%
disp('... selecting scenarios and computing probabilities');

% INIZIALIZATIONS
ShortTermProbDistr = struct;
ShortTermProbDistr.PosProb = zeros(length(PreSelection.MagSelBS),length(PreSelection.PosBSSel));
ShortTermProbDistr.DepProb = cell(length(PreSelection.MagSelBS),length(PreSelection.PosBSSel));
ShortTermProbDistr.DepProbPoints = cell(length(PreSelection.MagSelBS),size(PreSelection.PosBSSel,1));
ShortTermProbDistr.BarProb = cell(size(LongTermInfo.PSBarInfo.BarPSperModel,2),length(LongTermInfo.Discretizations.PS1_Mag.Val));
ShortTermProbDistr.RatioBSonTot = zeros(length(LongTermInfo.Discretizations.PS1_Mag.Val),1);

%%%%%%%%%%%%%%%%
% COMPUTE INTEGRAL FOR MAGNITUDES
%%%%%%%%%%%%%%%%
inter_sup = vertcat(0.5*(LongTermInfo.Discretizations.PS1_Mag.Val(1:end-1)+LongTermInfo.Discretizations.PS1_Mag.Val(2:end)),inf); % EVALUATING INTEGRAL
inter_inf = vertcat(-inf,0.5*(LongTermInfo.Discretizations.PS1_Mag.Val(1:end-1)+LongTermInfo.Discretizations.PS1_Mag.Val(2:end)));
prob_sup = normcdf(inter_sup,EarlyEst.Mag,EarlyEst.MagSigma);
prob_inf = normcdf(inter_inf,EarlyEst.Mag,EarlyEst.MagSigma);
ShortTermProbDistr.MagProb = prob_sup - prob_inf;
clear inter_sup inter_inf prob_sup prob_inf

%%%%%%%%%%%%%%%%
% COMPUTE PROBABILITY FOR SEPARATING PS/BS
%%%%%%%%%%%%%%%%
[lambda_BS,lambda_PS,lambda_PS_zone,tetra_vol] = lambdaBSPS_EE(SettingsLambdaBSPS.hypo_utm,SettingsLambdaBSPS.utmzone_hypo,SettingsLambdaBSPS.NormCov,SettingsLambdaBSPS.SD,SettingsLambdaBSPS.Npts_vs_Ntetra,SettingsLambdaBSPS.Nvref_vref_ratio,SettingsLambdaBSPS.moho_pts,SettingsLambdaBSPS.Nz,SettingsLambdaBSPS.mesh_subd,SettingsLambdaBSPS.out_tetra_vol);

ShortTermProbDistr.RatioBSonTot(:) = lambda_BS / (lambda_BS + lambda_PS);
ShortTermProbDistr.RatioPSonTot = 1 - ShortTermProbDistr.RatioBSonTot;
ShortTermProbDistr.RatioPSonPSTot = lambda_PS_zone / lambda_PS;
disp(['Lambdas from lambdaBSPS_EE (BS/PS): ' num2str(lambda_BS) ', ' num2str(lambda_PS)]);
for i=1:length(lambda_PS_zone)
    disp(['  ' SettingsLambdaBSPS.name_mesh{i} ': ' num2str(lambda_PS_zone(i))]);
end
clear lambda_BS lambda_PS tetra_vol
% CORRECT FOR NON EXISTING BS MAGNITUDES
ShortTermProbDistr.RatioBSonTot(LongTermInfo.Discretizations.PS1_Mag.Val>Settings.magBSmax)=0;
ShortTermProbDistr.RatioPSonTot(LongTermInfo.Discretizations.PS1_Mag.Val>Settings.magBSmax)=1;

%%%%%%%%%%%%%%%%
% EVALUATE IF TO COMPUTE OR NOT PS AND BS, BASED ON THE TOTAL PROBABILITY OF PS/BS
%%%%%%%%%%%%%%%%
ShortTermProbDistr.BScomputedYN = sum(ShortTermProbDistr.MagProb(PreSelection.MagSel).*ShortTermProbDistr.RatioBSonTot(PreSelection.MagSel)) > Settings.NegligibleProb;
ShortTermProbDistr.PScomputedYN = sum(ShortTermProbDistr.MagProb(PreSelection.MagSel).*(1-ShortTermProbDistr.RatioBSonTot(PreSelection.MagSel))) > Settings.NegligibleProb;

%%%%%%%%%%%%%%%%
% COMPUTE INTEGRAL FOR BS POS AND DEPTHS
%%%%%%%%%%%%%%%%
%    SET GRID FOR INTEGRATION
ShortTermProbDistr.SpaceGrid  = Settings.Z2XYfact*Settings.SpaceBin; % in m
ShortTermProbDistr.SpaceDepth = Settings.SpaceBin; % in m
xGrid = min(LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,1)):ShortTermProbDistr.SpaceGrid:max(LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,1));
yGrid = min(LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,2)):ShortTermProbDistr.SpaceGrid:max(LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,2));
allDepth = 1.E3*[LongTermInfo.Discretizations.BS3_Depth.ValVec{PreSelection.MagSelBS,PreSelection.PosBSSelOut}];
allDepthMoho = -1.E3*LongTermInfo.Discretizations.BS2_Pos.DepthMoho(PreSelection.PosBSSelOut);
zGrid = min(allDepth):ShortTermProbDistr.SpaceDepth:max(max(allDepth),max(allDepthMoho));
[XX2D,YY2D] = meshgrid(xGrid,yGrid);
GRID2D = [XX2D(:),YY2D(:)];
[XX3D,YY3D,ZZ3D] = meshgrid(xGrid,yGrid,zGrid);
GRID3D = [XX3D(:),YY3D(:),ZZ3D(:)];
clear xGrid yGrid zGrid allDepth allDepthMoho

%    MAP BS POSITIONS IN GRIDS
[tmp ShortTermProbDistr.refPosSelOut2D]=pdist2([LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,2) LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,1)],[YY2D(:) XX2D(:)],'euclidean','smallest',1);
[tmp ShortTermProbDistr.refPosSelOut3D]=pdist2([LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,2) LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSelOut,1)],[YY3D(:) XX3D(:)],'euclidean','smallest',1);
clear XX2D YY2D tmp

% COMPUTE HYPOCENTRAL PROBABILITY DISTRIBUTION FOR BS, IF REQUIRED
tmpBSScen = zeros(length(LongTermInfo.Discretizations.BS1_Mag.Val),length(length(PreSelection.PosBSSel)));
disp('Updated covariance matrices')
if ShortTermProbDistr.BScomputedYN
    for i1 = 1:length(PreSelection.MagSelBS)
        imag = PreSelection.MagSelBS(i1);
        
        % COMPUTE BS HALF WIDTH, TO CORRECT TOP TO CENTER OF THE FAULT
        tmpHalfWidth = UsedFunctions.CorrectBSverticalPos(LongTermInfo.Discretizations.BS1_Mag.Val(imag));
        
        % UPDATE MEAN & COVARIANCE MATRIX WITH FAULT DIMENSION
        tmpMU = EarlyEst.PosMean_3D'; % NO CORRECTIONS FOR CENTER OF FAULTS
        tmpCOV = EarlyEst.PosCovMat_3D;
        tmpCOV(1,1) = tmpCOV(1,1) + (UsedFunctions.CorrectBShorizontalPos(LongTermInfo.Discretizations.BS1_Mag.Val(imag)))^2;
        tmpCOV(2,2) = tmpCOV(2,2) + (UsedFunctions.CorrectBShorizontalPos(LongTermInfo.Discretizations.BS1_Mag.Val(imag)))^2;
        tmpCOV(3,3) = tmpCOV(3,3) + (tmpHalfWidth)^2;
        
        disp(['For M = ' num2str(LongTermInfo.Discretizations.BS1_Mag.Val(imag)) ':'])
        disp(['STD XX = ' num2str(sqrt(EarlyEst.PosCovMat_3D(1,1))) ' -> ' num2str(sqrt(tmpCOV(1,1))) ' m'])
        disp(['STD YY = ' num2str(sqrt(EarlyEst.PosCovMat_3D(2,2))) ' -> ' num2str(sqrt(tmpCOV(2,2))) ' m'])
        disp(['STD ZZ = ' num2str(sqrt(EarlyEst.PosCovMat_3D(3,3))) ' -> ' num2str(sqrt(tmpCOV(3,3))) ' m'])
        
        for i2 = 1:length(PreSelection.PosBSSel)
            ipos = PreSelection.PosBSSel(i2);
            
            % COUNT NUMBER OF SCENARIOS TO TREAT
            tmpBSScen(imag,ipos)=length(LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos})*length(LongTermInfo.Discretizations.BS4_FocMech.ID);
            
            % FIND POSITION OF OUTER GRID (USED TO SOLVE BORDER ISSUE)
            i2out = find(PreSelection.PosBSSelOut == ipos);
            
            % CHECK CONSISTENCY OF INPUT FOR DEPTH
            if isempty(LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos})
                disp(['!!! ERROR IN DATA STORED: NO DEPTH DEFINED ' ...
                    'MAG:' num2str(LongTermInfo.Discretizations.BS1_Mag.Val(imag)) ', ' ...
                    'POS:' num2str(LongTermInfo.Discretizations.BS2_Pos.Val(ipos,:),'%f %f') ',set to 1 km!']);
                LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos} = 1;
            end
            
            % COMPUTE PROBABILITIES OF POSITION AND DEPTH IN 3D OR IN 2D
            if Settings.gauss3Dyn
                
                % SELECT POINTS ABOVE MOHO IN THIS CELL
                cond1=ShortTermProbDistr.refPosSelOut3D == i2out;
                %      cond1=refPosSel3D == i2;
                %                cond2= ZZ3D(:) <= -1.E3*LongTermInfo.Discretizations.BS2_Pos.DepthMoho(ipos);
                cond2= and(ZZ3D(:) > tmpHalfWidth,ZZ3D(:) <= -1.E3*LongTermInfo.Discretizations.BS2_Pos.DepthMoho(ipos) + tmpHalfWidth);
                isel = find(and(cond1,cond2'));
                
                % ASSING EACH INTEGRATION POINT TO SCENARIO DEPTH
                tmpLocalDepths = 1.E3*LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos}+tmpHalfWidth;
                [tmp refDepthSel3D]=pdist2(tmpLocalDepths',ZZ3D(isel)','euclidean','smallest',1);
                
                % COMPUTE PROBABILITY WITHOUT INTEGRAL
                tmpPoints = [repmat([LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSel(i2),1) LongTermInfo.Discretizations.BS2_Pos.ValUTM(PreSelection.PosBSSel(i2),2)],length(tmpLocalDepths),1) tmpLocalDepths'];
                tmpProbPoints = NormMultiDvec(tmpPoints',tmpMU,tmpCOV);
                ShortTermProbDistr.DepProbPoints{i1,i2} = tmpProbPoints';
                
                % COMPUTE PROBABILITY FOR EACH SCENARIO
                tmpProb = NormMultiDvec(GRID3D(isel,:)',tmpMU,tmpCOV);
                for i3=1:length(LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos})
                    ShortTermProbDistr.DepProb{i1,i2}(i3) = sum(tmpProb(refDepthSel3D==i3));
                end
                
                clear cond1 cond2 tmpLocalDepths tmpPoints tmpProbPoints
            else
                % COMPUTING INTEGRAL ON 2D GRID FOR POSITIONS
                selTmp = find(ShortTermProbDistr.refPosSelOut2D == i2out);
                if not(isempty(selTmp))
                    ShortTermProbDistr.PosProb(i1,i2) = sum(NormMultiDvec(GRID2D(selTmp,:)',tmpMU(1:2),tmpCOV(1:2,1:2)));
                end
                
                % COMPUTING INTEGRAL FOR DEPTHS SEPARATELY
                tmpLocalDepths = 1.E3*LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos}+UsedFunctions.CorrectBSverticalPos(LongTermInfo.Discretizations.BS1_Mag.Val(imag));
                % EVALUATING INTEGRAL                
                inter_sup = vertcat(0.5*(tmpLocalDepths(1:end-1)'+tmpLocalDepths(2:end)'),-1.E3*LongTermInfo.Discretizations.BS2_Pos.DepthMoho(ipos)); 
                inter_inf = vertcat(0,0.5*(tmpLocalDepths(1:end-1)'+tmpLocalDepths(2:end)'));
                prob_sup = normcdf(inter_sup,tmpMU(3),sqrt(tmpCOV(3,3)));
                prob_inf = normcdf(inter_inf,tmpMU(3),sqrt(tmpCOV(3,3)));
                ShortTermProbDistr.DepProb{i1,i2} = prob_sup - prob_inf;
                ShortTermProbDistr.DepProb{i1,i2} = ShortTermProbDistr.DepProb{i1,i2} / sum(ShortTermProbDistr.DepProb{i1,i2});
                
                clear inter_sup inter_inf prob_sup prob_inf prob_inf tmpLocalDepths selTmp
            end
        end
        % NORMALIZATION FOR EACH MAGNITUDE SEPARATELY
        if Settings.gauss3Dyn
            NormFact = sum([ShortTermProbDistr.DepProb{i1,:}]);
            NormFactPoints = sum([ShortTermProbDistr.DepProbPoints{i1,:}]);
            for i2 = 1:length(PreSelection.PosBSSel)
                ShortTermProbDistr.DepProb{i1,i2}=ShortTermProbDistr.DepProb{i1,i2}/NormFact;
                ShortTermProbDistr.DepProbPoints{i1,i2}=ShortTermProbDistr.DepProbPoints{i1,i2}/NormFactPoints;
            end
            % SET SPATIAL PROBABILITY TO 1, SINCE ALREADY INCLUDED IN 3D INTEGRATION
            ShortTermProbDistr.PosProb(i1,:) = ones(1,length(PreSelection.PosBSSel));
        else
            ShortTermProbDistr.PosProb(i1,:) = ShortTermProbDistr.PosProb(i1,:) / sum(ShortTermProbDistr.PosProb(i1,:));
        end
        
    end
end
ShortTermProbDistr.TotBSScen = sum(tmpBSScen(:));
clear XX3D YY3D ZZ3D GRID3D GRID2D NormFact NormFactPoints refDepthSel3D
clear tmp tmpMU tmpCOV tmpBSScen tmpProb tmpW2
clear i1 i2 i2out i3 imag ipos isel

%%%%%%%%%%%%%%%%
% COMPUTE PS BAR PROBABILITY DISTRIBUTION FOR EACH PS MODEL, IF REQUIRED
%%%%%%%%%%%%%%%%
tmpPSScen = zeros(length(LongTermInfo.Discretizations.PS1_Mag.Val),size(LongTermInfo.PSBarInfo.BarPSperModel,2));
if ShortTermProbDistr.PScomputedYN
    ShortTermProbDistr.PSmodelYN=ones(length(LongTermInfo.Discretizations.PS1_Mag.Val),size(LongTermInfo.PSBarInfo.BarPSperModel,2));
    
    if strcmp(PreSettings.ProjectName(1:7),'TSUMAPS')
        ShortTermProbDistr.nModels = 8; % IN NEAMTHM18, THE FIRST 8 MODELS ARE EQUAL TO THE SECOND 8
        if not(length(LongTermInfo.ModelsWeight.PS2_Bar.Wei)==16)
            stop
        end
    else
        ShortTermProbDistr.nModels = size(LongTermInfo.PSBarInfo.BarPSperModel,2);
    end

    for imod=1:ShortTermProbDistr.nModels
        for i1 = 1:length(PreSelection.MagSel)
            imag = PreSelection.MagSel(i1);
            
            if isempty(LongTermInfo.PSBarInfo.BarPSperModel{imag,imod})
                ShortTermProbDistr.PSmodelYN(imag,imod) = 0;
            else
                % COUNT THE NUMBER OF SCENARIOS
                tmpPSScen(imag,imod) = size(PreSelection.BarPSSelPerModMag{imag,imod},1);
                
                if Settings.figcheckYN
                    x=LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(PreSelection.BarPSSelPerModMag{imag,imod},1);
                    y=LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(PreSelection.BarPSSelPerModMag{imag,imod},2);
                    plot(x,y,'.')
                    hold on
                end
                
                % COMPUTING PROBABILITY OF EACH BARICENTER
                EarlyEst.BarProb{imod,imag} = zeros(size(PreSelection.BarPSSelPerModMag{imag,imod},1),1);
                % UPDATE COVARIANCE MATRIX WITH FAULT DIMENSION
                tmpMU = EarlyEst.PosMean_2D';
                tmpCOV = EarlyEst.PosCovMat_2D;
                tmpCOV(1,1) = tmpCOV(1,1) + (UsedFunctions.CorrectPShorizontalPos(LongTermInfo.Discretizations.PS1_Mag.Val(imag)))^2;
                tmpCOV(2,2) = tmpCOV(2,2) + (UsedFunctions.CorrectPShorizontalPos(LongTermInfo.Discretizations.PS1_Mag.Val(imag)))^2;
                % COMPUTE PROBABILITY OF EACH BARICENTER
                for i2=1:size(PreSelection.BarPSSelPerModMag{imag,imod},1)
                    ibar=PreSelection.BarPSSelPerModMag{imag,imod}(i2);
                    ireg = LongTermInfo.PSBarInfo.BarPSperModelReg{imag,imod}(ibar);
                    tmpVAR = LongTermInfo.PSBarInfo.BarPSperModelUTM{imag,imod}(ibar,:)';
                    
                    % if the subduction is not selected, also the bar is discharged
                    if ShortTermProbDistr.RatioPSonPSTot(SettingsLambdaBSPS.regionsPerPS(ireg)) > 0
                        ShortTermProbDistr.BarProb{imod,imag}(i2) = NormMultiDvec(tmpVAR,tmpMU,tmpCOV);
                    else
                        ShortTermProbDistr.BarProb{imod,imag}(i2) = 0;
                    end
                end
                
                % NORMALIZE OVER BARICENTRES WITH EQUAL MAGNITUDE AND MODEL
                ShortTermProbDistr.BarProb{imod,imag} = ShortTermProbDistr.BarProb{imod,imag} / sum(ShortTermProbDistr.BarProb{imod,imag});
                
                if Settings.figcheckYN
                    figure()
                    tmpCoord=LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(PreSelection.BarPSSelPerModMag{imag,imod},:);
                    tmpDepth=LongTermInfo.PSBarInfo.BarPSperModelDepth{imag,imod}(PreSelection.BarPSSelPerModMag{imag,imod},:);
                    tmpProb = ShortTermProbDistr.BarProb{imod,imag};
                    plot(EarlyEst.lon,EarlyEst.lat,'s','markerfacecolor','c','markeredgecolor','b','markersize',10)
                    hold on
                    scatter(tmpCoord(:,1),tmpCoord(:,2),50,tmpProb,'filled');
                    hold on
                    colorbar;
                    delta = 5;
                    set(gca,'xlim',[min(tmpCoord(:,1))-delta,max(tmpCoord(:,1))+delta],'ylim',[min(tmpCoord(:,2))-delta,max(tmpCoord(:,2))+delta])
                    
                    clear tmpCoord tmpDepth tmpProb
                end
            end
        end
    end
end
ShortTermProbDistr.TotPSScen = sum(tmpPSScen(:));
clear tmpPSScen tmpMU tmpCOV tmpVAR
clear i1 i2 ibar imag imod

%% %%%%%%%%%%%%%%
%%% COMPUTE PROBABILITIES SCENARIOS
%%% Output: ParScenBS, ProbScenBSFact, ProbScenBS, ParScenPS, ProbScenBSFact, ProbScenPS
%%%         ScenarioProb contains the list of scenarios of the ensemble, including parameters and probabilities 
%%%%%%%%%%%%%%%%
% INITIALIZATIONS
ScenarioProb = struct;

ScenarioProb.ParScenBS = zeros(ShortTermProbDistr.TotBSScen,11);
ScenarioProb.ParScenBS_ID = zeros(ShortTermProbDistr.TotBSScen,1);
ScenarioProb.ProbScenBSFact = zeros(ShortTermProbDistr.TotBSScen,5);
ScenarioProb.ProbScenBS = [];
iScenBS=0;

ScenarioProb.ParScenPS = zeros(ShortTermProbDistr.TotPSScen,7);
ScenarioProb.ParScenPS_ID = zeros(ShortTermProbDistr.TotPSScen,1);
ScenarioProb.ProbScenPSFact = zeros(ShortTermProbDistr.TotPSScen,5);
ScenarioProb.ProbScenPS = [];
iScenPS=0;

ScenarioProb.BScomputedYN = sum(ShortTermProbDistr.MagProb.*ShortTermProbDistr.RatioBSonTot) > Settings.NegligibleProb;
ScenarioProb.PScomputedYN = sum(ShortTermProbDistr.MagProb.*(1-ShortTermProbDistr.RatioBSonTot)) > Settings.NegligibleProb;

% RUN OVER MAGNITUDES
if ShortTermProbDistr.BScomputedYN
    for i1 = 1:length(PreSelection.MagSel)
        imag = PreSelection.MagSel(i1);
        
        if ismember(imag,PreSelection.MagSelBS)
            for i2 = 1:length(PreSelection.PosBSSel)
                ipos = PreSelection.PosBSSel(i2);
                ireg = LongTermInfo.Discretizations.BS2_Pos.Region(ipos);
                
                if strcmp(PreSettings.ProjectName,'ChEESE_PTF')
                    iregFocMec = 1;
                    iregParam = 1; 
                    iposFocMec = ireg;
                else
                    iregFocMec = ireg;
                    iregParam = ireg;
                    iposFocMec = ipos;
                end
                
                if isempty(LongTermInfo.region{iregFocMec})
                    disp(['!!!! region not yet loaded: ' num2str(iregFocMec)]);
                    PreSettings.SelectedRegions=[iregFocMec];
                    [LongTermInfo,tmp] = PTF_preLoad_Med(PreSettings,LongTermInfo);
                end
                
                if Settings.normalizeAngleProb
                    RegMeanProb_BS4=LongTermInfo.region{iregFocMec}.MeanProb_BS4_FocMech_Reg.ValNorm;
                else
                    RegMeanProb_BS4=LongTermInfo.region{iregFocMec}.MeanProb_BS4_FocMech_Reg.Val;
                end
                if isempty(RegMeanProb_BS4)
                    disp(['!!!! angles not found in region ' num2str(ireg) ' and position ' num2str(ipos)]);
                else
                    tmpProbAngles = ...
                        RegMeanProb_BS4(iposFocMec == LongTermInfo.region{iregFocMec}.MeanProb_BS4_FocMech_Reg.indPos,:);
                    
                    % ENUMERATE ALL RELEVANT SCENARIOS FOR EACH MAG AND POS
                    for i3 = 1:length(LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos})
                        for i4 = 1:length(LongTermInfo.Discretizations.BS4_FocMech.ID)
                            iang = i4;
                            iScenBS = iScenBS+1;
                            
                            AreaLength = LongTermInfo.Discretizations.BS5_AreaLength.Val{iregParam,imag,iang}; % AREA, L in KM
                            Slip = LongTermInfo.Discretizations.BS6_Slip.Val{iregParam,imag,iang}; % SLIP jn M
                            
                            % SAVE PARAMETERS
                            ScenarioProb.ParScenBS_ID(iScenBS) = iScenBS;
                            ScenarioProb.ParScenBS(iScenBS,1) = ireg;
                            ScenarioProb.ParScenBS(iScenBS,2) = ...
                                LongTermInfo.Discretizations.BS1_Mag.Val(imag);
                            ScenarioProb.ParScenBS(iScenBS,3:4) = ...
                                LongTermInfo.Discretizations.BS2_Pos.Val(ipos,:); % long, lat
                            ScenarioProb.ParScenBS(iScenBS,5) = ...
                                LongTermInfo.Discretizations.BS3_Depth.ValVec{imag,ipos}(i3);
                            ScenarioProb.ParScenBS(iScenBS,6:8) =  ...
                                LongTermInfo.Discretizations.BS4_FocMech.Val(i4,:);
                            ScenarioProb.ParScenBS(iScenBS,9:10) = AreaLength; % AREA, L in KM
                            ScenarioProb.ParScenBS(iScenBS,11) = Slip; % SLIP jn M
                            
                            % SAVE PROBABILITIES (ALL FACTORS)
                            ScenarioProb.ProbScenBSFact(iScenBS,1)=ShortTermProbDistr.MagProb(imag);
                            ScenarioProb.ProbScenBSFact(iScenBS,2)=ShortTermProbDistr.PosProb(i1,i2);
                            ScenarioProb.ProbScenBSFact(iScenBS,3)=ShortTermProbDistr.RatioBSonTot(imag);
                            ScenarioProb.ProbScenBSFact(iScenBS,4)=ShortTermProbDistr.DepProb{i1,i2}(i3);
                            ScenarioProb.ProbScenBSFact(iScenBS,5)=tmpProbAngles(iang);
                        end
                    end
                end
            end
        end
    end
    ScenarioProb.ProbScenBS = prod(ScenarioProb.ProbScenBSFact,2);

end

% PS SCENARIOS, IF REQUIRED
if ShortTermProbDistr.PScomputedYN
    for i1 = 1:length(PreSelection.MagSel)
        imag = PreSelection.MagSel(i1);
        if strcmp(PreSettings.ProjectName(1:7),'TSUMAPS')
            nmod=8;
        else
            nmod=size(LongTermInfo.PSBarInfo.BarPSperModel,2);
        end
        for imod = 1:nmod
            if strcmp(PreSettings.ProjectName(1:7),'TSUMAPS')
                % IN NEAMTHM18, THE FIRST 8 MODELS ARE EQUAL TO THE SECOND 8
                probwei = LongTermInfo.ModelsWeight.PS2_Bar.Wei(imod)+LongTermInfo.ModelsWeight.PS2_Bar.Wei(imod+8);
            else
                probwei = LongTermInfo.ModelsWeight.PS2_Bar.Wei(imod);
            end
            
            
            if and(imod<=ShortTermProbDistr.nModels,not(isempty(PreSelection.BarPSSelPerModMag{imag,imod})))
                for i2=1:size(PreSelection.BarPSSelPerModMag{imag,imod},1)
                    ibar = PreSelection.BarPSSelPerModMag{imag,imod}(i2);
                    ireg = LongTermInfo.PSBarInfo.BarPSperModelReg{imag,imod}(ibar);
                    
                    if ShortTermProbDistr.RatioPSonPSTot(SettingsLambdaBSPS.regionsPerPS(ireg)) > 0 % only if the subduction is activated
                        % vectmp (slip) ENUMERATE ALL RELEVANT SCENARIOS FOR EACH MAG, MODEL AND BARICENTER
                        if isempty(LongTermInfo.region{ireg}) && not(strcmp(PreSettings.ProjectName,'ChEESE_PTF'))
                            disp(['!!!! region not yet loaded: ' num2str(ireg)]);
                           PreSettings.SelectedRegions=[ireg];                            
                            [LongTermInfo,tmp] = PTF_preLoad_Med(PreSettings,LongTermInfo);
                        end
                        
                        slipVal =  unique(LongTermInfo.region{ireg}.ScenarioListPSReg.SlipDistribution(LongTermInfo.region{ireg}.ScenarioListPSReg.magPSInd==imag));
                        nScen = length(slipVal);
                        locScen = iScenPS+1:iScenPS+nScen;
                        vectmp = ones(nScen,1);
                        
                        % COUNT SCENARIOS
                        iScenPS = iScenPS+nScen;
                        
                        % SAVE PARAMETERS
                        ScenarioProb.ParScenPS(locScen,1) = vectmp*ireg;
                        ScenarioProb.ParScenPS(locScen,2) = vectmp*LongTermInfo.Discretizations.PS1_Mag.Val(imag);
                        ScenarioProb.ParScenPS(locScen,3) = vectmp*LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(ibar,1);
                        ScenarioProb.ParScenPS(locScen,4) = vectmp*LongTermInfo.PSBarInfo.BarPSperModel{imag,imod}(ibar,2);
                        ScenarioProb.ParScenPS(locScen,5) = vectmp*imod; % Models 1/2 Strasser, 3/4 Mourotani
                        ScenarioProb.ParScenPS(locScen,6) = slipVal;
                        ScenarioProb.ParScenPS(locScen,7) = vectmp*abs(LongTermInfo.PSBarInfo.BarPSperModelDepth{imag,imod}(ibar));
                        
                        ScenarioProb.ProbScenPSFact(locScen,1)=vectmp*ShortTermProbDistr.MagProb(imag);
                        ScenarioProb.ProbScenPSFact(locScen,2)=vectmp*ShortTermProbDistr.BarProb{imod,imag}(i2);
                        ScenarioProb.ProbScenPSFact(locScen,3)= vectmp*ShortTermProbDistr.RatioPSonTot(imag)* ...
                            ShortTermProbDistr.RatioPSonPSTot(SettingsLambdaBSPS.regionsPerPS(ireg)); % SEPARATING DIFFERENT MESHES
                        
                        % NORMALIZING CONSIDERING ONLY THE MODELS THAT HAVE THE MAGNITUDE
                        totwei = sum(LongTermInfo.ModelsWeight.PS2_Bar.Wei(ShortTermProbDistr.PSmodelYN(imag,:)==1));
                        ScenarioProb.ProbScenPSFact(locScen,4)=vectmp*probwei / totwei ;
                        
                        ScenarioProb.ProbScenPSFact(locScen,5)=vectmp./nScen;
                        if nScen == 0
                            ireg
                            imag
                        end
                    end
                end
            end
        end
    end
end
ScenarioProb.ProbScenPS = prod(ScenarioProb.ProbScenPSFact,2);

disp(['N. BS scenarios: ' num2str(length(ScenarioProb.ProbScenBS))]);
disp(['N. PS scenarios:	' num2str(length(ScenarioProb.ProbScenPS))]);

% RENORMALIZING (TO MANAGE EVENTS OUTSIDE nSigma, BS OF LARGE MAGNITUDES, ETC.
ScenarioProb.TotProbBS_preNorm = sum(ScenarioProb.ProbScenBS);
ScenarioProb.TotProbPS_preNorm = sum(ScenarioProb.ProbScenPS);
ScenarioProb.TotProb_preNorm = ScenarioProb.TotProbBS_preNorm + ScenarioProb.TotProbPS_preNorm;
if ScenarioProb.TotProb_preNorm < 1.0
    disp([' ---> ScenarioProb.TotProbBS_preNorm = ' num2str(ScenarioProb.TotProbBS_preNorm)]);
    disp([' ---> ScenarioProb.TotProbPS_preNorm = ' num2str(ScenarioProb.TotProbPS_preNorm)]);
    disp(' ---> Total probabilty renormalized to 1');
    ScenarioProb.ProbScenBS = ScenarioProb.ProbScenBS / ScenarioProb.TotProb_preNorm;
    ScenarioProb.ProbScenPS = ScenarioProb.ProbScenPS / ScenarioProb.TotProb_preNorm;
end

clear iScenBS iScenPS
clear i1 i2 i3 i4 iang ibar imag imod ipos iposAng ireg
clear vectmp slipVal nScen locScen tmpProbAngles
toc

% the scenario list produced at this point can be used to set external simulations, to be read below to
%    quantify hazard curves


%% %%%%%%%%%%%%%%
%%% COMPUTING HAZARD CURVES
%%%%%%%%%%%%%%%%%
disp('... computing hazard curves at POIs');
HazardCurves = struct;
BStsunamiIntensity = zeros(size(ScenarioProb.ParScenBS,1),length([POIInfo.SelectedPOI]));
PStsunamiIntensity = zeros(size(ScenarioProb.ParScenPS,1),length([POIInfo.SelectedPOI]));
if Settings.writeOutTesting
    OutTesting = struct;
    OutTesting.BSprob = nan(size(ScenarioProb.ParScenBS,1),1);
    OutTesting.PSprob = nan(size(ScenarioProb.ParScenPS,1),1);
end
HazardCurves.RefUtmZone = EarlyEst.RefUtmZone;

% INITIALIZATION OF X-AXIS OF HAZARD CURVES
HazardCurves.OriginalHazardCurveThreshols = LongTermInfo.GenericHazardCurveThresholds;
HazardCurves.tsunamiIntensityName = 'Near-coast wave amplitude (m)';
HazardCurves.tsunamiIntensityRunupAmpFact = 1;
HazardCurves.HazardCurveThreshols = HazardCurves.OriginalHazardCurveThreshols * HazardCurves.tsunamiIntensityRunupAmpFact;
HazardCurves.nThresholds = length(HazardCurves.HazardCurveThreshols);

% HAZARD CURVES FOR BS, IF REQUIRED
HazardCurves.regBS=unique(ScenarioProb.ParScenBS(:,1));
HazardCurves.regBS=HazardCurves.regBS(HazardCurves.regBS>0);
HazardCurves.hc_poiBS = zeros(length(POIInfo.SelectedPOI),HazardCurves.nThresholds);
HazardCurves.hc_poiBS_mean = zeros(length(POIInfo.SelectedPOI),1);
if ShortTermProbDistr.BScomputedYN
    if strcmpi(Settings.selectedIntensityMeasure,'gl-sims')
        %LOAD SIMULATIONS
        simPS_hmax=importdata([ PreSettings.genInputFolder EarlyEst.ID '_BS_hmax_ALL.mat']);
        simPS_hmax=simBS_hmax.maxts_gl_p2t;
    end
    for k=1:length(HazardCurves.regBS)
        % index within the selected list of BS scenarios
        isel = find(ScenarioProb.ParScenBS(:,1) == HazardCurves.regBS(k)); 
        probsel = ScenarioProb.ProbScenBS(isel);
        
        % CREATE IDs FOR THE SELECTED SCENARIOS IN REGIONS
        tmpsel = ScenarioProb.ParScenBS(isel,2:8);
        IDsel_BS = tmpsel * LongTermInfo.vecID';
        
        % CREATE IDs FOR ALL THE SCENARIOS IN REGIONS
        tmpreg = LongTermInfo.region{HazardCurves.regBS(k)}.ScenarioListBSReg.Parameters(:,1:7);
        IDreg_BS = tmpreg * LongTermInfo.vecID';
        
        % SELECT RELEVANT SCENARIOS
        [a,b]=ismember(IDsel_BS,IDreg_BS);
        if not(isempty(find(b==0,1)))
            itmp=find(b==0); %itmp=itmp(4);
            disp(['!!! BS simulations not found (' num2str(length(itmp)) ...
                ') in Region #' num2str(regBS(k)) ' - Slip ' num2str(islip) ...
                '   ---> max bias: ' num2str(sum(probsel(itmp))) ])
        else
            disp([' ---> all BS simulations found in Region #' num2str(HazardCurves.regBS(k))])
        end
        
        % SAVE PROBABILITY IN OUTTESTING
        if Settings.writeOutTesting
            OutTesting.BSprob(isel) = probsel;
        end
        
        % CUMULATE OVER HC AT EACH POI
        bfound = find(b>0);
        tmpmax=0;
        for ip = 1:length(POIInfo.SelectedPOI)
            % manage the different intensity measures
            % it sets the simulated intensity for all scenarios at the ip forecast point (InputIntensityMeasure)
            % it sets indeces (itmpInd) for all scenarios at the ip forecast point. Index points to precoputed hazard thresholds (ConditionalHazardCurve)
            switch Settings.selectedIntensityMeasure
                case 'gl'
                    InputIntensityMeasure =  LongTermInfo.region{HazardCurves.regBS(k)}.glVal_BS(ip,b(bfound))';
                    itmpInd = InputIntensityMeasure > 0;
                case 'gl-sims'
                    itmpInd = LongTermInfo.region{HazardCurves.regBS(k)}.POI{POIInfo.SelectedPOI(ip).Index}.BS.gl.GenericIndexHazardCurves(b(bfound));
                    
                    InputIntensityMeasure = simBS_hmax(isel,ip);
                    ConditionalHazardCurve = LongTermInfo.GenericConditionalHazardCurves_GL;
            end
            tmpmax=max(max(InputIntensityMeasure(:)),tmpmax);
            if not(isempty(find(itmpInd>0,1)))
                itmp2 = itmpInd(itmpInd>0);
                probTmpScen =    probsel(bfound(itmpInd>0));
                probTmp = repmat(probTmpScen,1,size(HazardCurves.OriginalHazardCurveThreshols,1))';
                
                % it recomputes the conditional hazard curve based on modeled intensity
                mu=log(InputIntensityMeasure(itmpInd>0));sigma=1;                
                xxmat=repmat(HazardCurves.OriginalHazardCurveThreshols',length(mu),1);
                mumat=repmat(mu,1,length(HazardCurves.OriginalHazardCurveThreshols));
                condhcTmp = 1-logncdf(xxmat,mumat,sigma)';               
                
                condhcTmpMean = exp(mu+0.5*sigma^2);
                
                % SAVE TSUNAMI INTENSITY IN OUTTESTING
                BStsunamiIntensity(isel,ip) = InputIntensityMeasure;
                
                if max(condhcTmp(:))>1; stop; end
                HazardCurves.hc_poiBS(ip,:) = HazardCurves.hc_poiBS(ip,:) + sum(probTmp.*condhcTmp,2)';
                HazardCurves.hc_poiBS_mean(ip) = HazardCurves.hc_poiBS_mean(ip) + sum(probTmpScen.*condhcTmpMean);
            end
        end
    end
    clear a b IDreg_BS IDsel_BS
    clear k ip isel
    clear tmpreg probsel bfound indSel tmpsel
    clear itmp itmp2 probTmp condhcTmp condhcTmp
    clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7
end

% HAZARD CURVES FOR PS, IF REQUIRED
HazardCurves.regPS=unique(ScenarioProb.ParScenPS(ScenarioProb.ParScenPS(:,1)>0,1));
HazardCurves.hc_poiPS = zeros(length(POIInfo.SelectedPOI),HazardCurves.nThresholds);
HazardCurves.hc_poiPS_mean = zeros(length(POIInfo.SelectedPOI),1);

if ShortTermProbDistr.PScomputedYN
    if strcmpi(Settings.selectedIntensityMeasure,'gl-sims')
        %LOAD SIMULATIONS
        simPS_hmax=importdata([ PreSettings.genInputFolder EarlyEst.ID '_PS_hmax_ALL.mat']);
        simPS_hmax=simPS_hmax.maxts_gl_p2t;
    end
    for k=1:length(HazardCurves.regPS)
        % CUMULATE SCENARIOS WITH OR WITHOUT SLIP DISTRIBUTION
        cond1 = ScenarioProb.ParScenPS(:,1) == HazardCurves.regPS(k);
        isel = find(cond1);
        probsel = ScenarioProb.ProbScenPS(isel);
        
        % CREATE IDs FOR THE SELECTED SCENARIOS IN REGIONS
        tmpsel = ScenarioProb.ParScenPS(isel,2:6);
        tmpselCorrMod=mod(tmpsel(:,4),8); tmpselCorrMod(tmpselCorrMod==0)=8;
        tmpsel(:,4)=tmpselCorrMod;
        IDsel_PS = tmpsel * LongTermInfo.vecID(1:5)';
        
        % CREATE IDs FOR ALL THE SCENARIOS IN REGIONS
        if strcmp(PreSettings.ProjectName,'ChEESE_PTF') 
            % !!!! ATTENTION !!! it works only for 2 sigma because it assumes that the scenarios simualted coincide with the ones generated here
            tmpreg = tmpsel;
            IDreg_PS = IDsel_PS;
        else
            tmp=vertcat(LongTermInfo.region{HazardCurves.regPS(k)}.ScenarioListPSReg.modelVal,LongTermInfo.region{HazardCurves.regPS(k)}.ScenarioListPSReg.SlipDistribution)';
            tmpreg = horzcat(LongTermInfo.region{HazardCurves.regPS(k)}.ScenarioListPSReg.Parameters(:,1:3),tmp);
            IDreg_PS = tmpreg * LongTermInfo.vecID(1:5)';
        end
        
        % SELECT RELEVANT SCENARIOS
        [a,b]=ismember(IDsel_PS,IDreg_PS);
        if not(isempty(find(b==0,1)))
            itmp=find(b==0);
            disp('!!! PS simulations not found')
            if Settings.verboseYN
                tmpsel(itmp,:)
                ScenarioProb.ParScenPS(isel(itmp),:)
                itmp2 = itmp(1);tmpsel(itmp2,:)
                ScenarioProb.ParScenPS(isel(itmp2),:)
                isell = and(tmpreg(:,1)==tmpsel(itmp2,1),and(abs(tmpreg(:,2)-tmpsel(itmp2,2))<10,abs(tmpreg(:,3)-tmpsel(itmp2,3))<10));
                tmpreg(isell,:)
                tmpreg(and(and(tmpreg(:,2)==tmpsel(itmp2,2),tmpreg(:,3)==tmpsel(itmp2,3)),tmpreg(:,1)==tmpsel(itmp2,1)),:)
            end
            stop
        else
            disp([' ---> all PS simulations found in Region #' num2str(HazardCurves.regPS(k))])
        end
        
        % SAVE PROBABILITY IN OUTTESTING
        if Settings.writeOutTesting
            OutTesting.PSprob(isel) = probsel;
        end
        
        
        % CUMULATE OVER HC AT EACH POI
        bfound = find(b>0);
        
        for ip = 1:length(POIInfo.SelectedPOI)
            % manage the different intensity measures
            switch Settings.selectedIntensityMeasure
                case 'gl'
                    InputIntensityMeasure =  LongTermInfo.region{HazardCurves.regPS(k)}.glVal_PS(ip,b(bfound));
                    itmpInd = InputIntensityMeasure > 0;
                case 'gl-sims'
                    ipp = POIInfo.SelectedPOI(ip).Index;
                    InputIntensityMeasure = simPS_hmax(isel,ipp)';
                    itmpInd = InputIntensityMeasure > 0;

            end
            if not(isempty(find(itmpInd>0,1)))
                itmp2 = itmpInd(itmpInd>0);
                
                probTmpScen =    probsel(bfound(itmpInd>0));
                probTmp = repmat(probTmpScen,1,length(HazardCurves.OriginalHazardCurveThreshols))';
                
                mu=log(InputIntensityMeasure(itmpInd>0));sigma=1;
                xxmat=repmat(HazardCurves.OriginalHazardCurveThreshols',length(mu),1);
                mumat=repmat(mu',1,length(HazardCurves.OriginalHazardCurveThreshols));
                condhcTmp = 1-logncdf(xxmat,mumat,sigma)';                                               
                condhcTmpMean = exp(mu+0.5*sigma^2);
                
                % SAVE TSUNAMI INTENSITY IN OUTTESTING
                PStsunamiIntensity(isel,ip) = InputIntensityMeasure;
                HazardCurves.hc_poiPS(ip,:) = HazardCurves.hc_poiPS(ip,:) + sum(probTmp.*condhcTmp,2)';
                HazardCurves.hc_poiPS_mean(ip) = HazardCurves.hc_poiPS_mean(ip) + sum(probTmpScen'.*condhcTmpMean);
            end
        end
        
    end
end
HazardCurves.hc_poi =  HazardCurves.hc_poiBS + HazardCurves.hc_poiPS;
HazardCurves.hc_poi_mean =  HazardCurves.hc_poiBS_mean + HazardCurves.hc_poiPS_mean;


% save intensity of scenario best-guess
pmax1=0;
pmax2=0;
if ShortTermProbDistr.BScomputedYN; [pmax1,isel1]=max(ScenarioProb.ProbScenBS);end
if ShortTermProbDistr.PScomputedYN; [pmax2,isel2]=max(ScenarioProb.ProbScenPS); end
if pmax1 > pmax2
   HazardCurves.bestGuessScen_tsunamiintisity = squeeze(BStsunamiIntensity(isel1,:));
else
   HazardCurves.bestGuessScen_tsunamiintisity = squeeze(PStsunamiIntensity(isel2,:));
end

% save intensity of envelope, similar to Catalan et al. 2020
HazardCurves.env_tsunamiintisity = zeros(size(HazardCurves.bestGuessScen_tsunamiintisity));
if ScenarioProb.BScomputedYN
   % spatial selection
   DistAll = distance( ScenarioProb.ParScenBS(:,4),ScenarioProb.ParScenBS(:,3),EarlyEst.lat,EarlyEst.lon )*pi*6371/180;
   DistMax = UsedFunctions.FuncMag2L_BS(EarlyEst.Mag)/2;
   % magnitude selection selection   
   mags = unique(ScenarioProb.ParScenBS(:,2));
   [tmp,iselmag]=min(abs(EarlyEst.Mag + 0.5 - mags));
   magref = mags(iselmag);
   % envelope
   isel=and(DistAll<DistMax,ScenarioProb.ParScenBS(:,2)==magref);
   HazardCurves.env_tsunamiintisity = max(BStsunamiIntensity(isel,:));
end
if ScenarioProb.PScomputedYN
   % spatial selection
   DistAll = distance( ScenarioProb.ParScenPS(:,4),ScenarioProb.ParScenPS(:,3),EarlyEst.lat,EarlyEst.lon )*pi*6371/180;
   DistMax = UsedFunctions.FuncMag2L_PS(EarlyEst.Mag)/2;
   % magnitude selection selection   
   mags = unique(ScenarioProb.ParScenPS(:,2));
   [tmp,iselmag]=min(abs(EarlyEst.Mag + 0.5 - mags));
   magref = mags(iselmag);
   isel=and(DistAll<DistMax,ScenarioProb.ParScenPS(:,2)==magref);
   % envelope
   if sum(HazardCurves.env_tsunamiintisity) == 0
      HazardCurves.env_tsunamiintisity = max(PStsunamiIntensity(isel,:));
   else
      tmpint = max(PStsunamiIntensity(isel,:));
      for i=1:length(tmpint)
         HazardCurves.env_tsunamiintisity(i)=max(HazardCurves.env_tsunamiintisity(i),tmpint(i));
      end
   end
end


clear a b cond1 IDreg_PS IDsel_PS
clear k isel ip itmp itmp2 probTmp scenSelTmp condhcTmp
clear tmpreg probsel bfound indSel tmpsel
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7

toc


%% %%%%%%%%%%%%%
%%% SAVINGS
%%%%%%%%%%%%%%%%
if Settings.writeOutTesting
    disp('...storing results for testing')
    OutTesting.ScenarioProb = ScenarioProb;
    OutTesting.POIlat = POIInfo.lat([POIInfo.SelectedPOI.Index]);
    OutTesting.POIlatUTM = POIInfo.latUTM([POIInfo.SelectedPOI.Index]);
    OutTesting.POIlon = POIInfo.lon([POIInfo.SelectedPOI.Index]);
    OutTesting.POIlonUTM = POIInfo.lonUTM([POIInfo.SelectedPOI.Index]);
    OutTesting.POIdep = POIInfo.dep([POIInfo.SelectedPOI.Index]);
    OutTesting.tsunamiIntensityName = HazardCurves.tsunamiIntensityName;
    OutTesting.tsunamiIntensityRunupAmpFact = HazardCurves.tsunamiIntensityRunupAmpFact;
    OutTesting.EarlyEst = EarlyEst;
    OutTesting.BStsunamiIntensity=BStsunamiIntensity;
    OutTesting.PStsunamiIntensity=PStsunamiIntensity;
        
end


%% %%%%%%%%%%%%%
%%% PLOTS
%%%%%%%%%%%%%%%%
if Settings.octaveYN
    fflush(stdout);
end

end




