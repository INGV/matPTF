function [EarlyEst] = SetEarlyEst(varargin)
disp('----- running: SetEarlyEst ------')

%% %%%%%%%%%%%%%
%%% INPUT PTF from Early Est
%%%%%%%%%%%%%%%%
EarlyEst = struct;
EarlyEst.ID = '';
EarlyEst.Folder = 'LocalInput/EarlyEst/';
EarlyEst.ID = varargin{1};


%% SET INPUT FROM REAL-TIME MONITORING
if strcmp(EarlyEst.ID,'2010_0227_maule')
  % set central values
    EarlyEst.lon = -72.898; EarlyEst.lat = -36.122; EarlyEst.Dep = 30.5; EarlyEst.Mag = 8.8; % 2020/10/19 MAULE CHILE (from USGS, depth W-PHASE DUPUTEL)

  % set uncertinaty
    EarlyEst.PosSigmaXX = 10000.^2;      % in m^2
    EarlyEst.PosSigmaXY = 0.; 
    EarlyEst.PosSigmaXZ = 0.; 
    EarlyEst.PosSigmaYY = 10000.^2; 
    EarlyEst.PosSigmaYZ = 0.; 
    EarlyEst.PosSigmaZZ = 10000.^2;     
    EarlyEst.MagSigma = 0.2;    

elseif strcmp(EarlyEst.ID,'2018_0000_neamwave17')        
  % set central values
    EarlyEst.lon = 21.00; EarlyEst.lat = 37.50; EarlyEst.Dep = 12.; EarlyEst.Mag = 8.5; 

  % set uncertinaty
    EarlyEst.PosSigmaXX = 10000.^2;      % in m^2
    EarlyEst.PosSigmaXY = 0.; 
    EarlyEst.PosSigmaXZ = 0.; 
    EarlyEst.PosSigmaYY = 10000.^2; 
    EarlyEst.PosSigmaYZ = 0.; 
    EarlyEst.PosSigmaZZ = 10000.^2;     
    EarlyEst.MagSigma = 0.2;    

else  

  % read from EarlyEst file (from LocalInput)
  fid=fopen([EarlyEst.Folder EarlyEst.ID '_stat.txt']);
  for i=1:10
     strLine{i} = fgetl(fid);
  end
  fclose(fid);

  % extract required information from EE in hindcasting mode
  tmp = regexp(strLine{10},' ','Split');
  EarlyEst.Mag = str2num(tmp{6});  %p50
  EarlyEst.MagSigma = 0.5*(str2num(tmp{7})-str2num(tmp{5})); %p84-p16
  tmp = regexp(strLine{4},'\s','Split');
  EarlyEst.lon =  str2num(tmp{5}); 
  EarlyEst.lat =  str2num(tmp{4}); 
  EarlyEst.Dep =  str2num(tmp{6}); 
  tmp = regexp(strLine{5},'\s','Split'); 
  fact = 1.E6;
  EarlyEst.PosSigmaXX =  str2num(tmp{2})*fact; 
  EarlyEst.PosSigmaXY =  str2num(tmp{4})*fact;  
  EarlyEst.PosSigmaXZ =  str2num(tmp{6})*fact;  
  EarlyEst.PosSigmaYY =  str2num(tmp{8})*fact;  
  EarlyEst.PosSigmaYZ =  str2num(tmp{10})*fact;  
  EarlyEst.PosSigmaZZ =  str2num(tmp{12})*fact;    
  clear InData tmp fact strLine fid i
end


end
