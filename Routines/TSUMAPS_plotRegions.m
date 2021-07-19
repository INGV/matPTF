function TSUMAPS_plotRegions(varargin)
% TSUMAPS_plotRegions(LongTermInfo.regions,HazardCurves.regBS(k))
   regions = varargin{1};
   if length(varargin) == 2;
      regsel = varargin{2};
   else
      regsel = [];       
   end
               

    for i=1:regions.Npoly
       plot(regions.Tlon(i,1:regions.Tleng(i)+1),regions.Tlat(i,1:regions.Tleng(i)+1),':r')
       hold on
%       text(mean(regions.Tlon(i,1:regions.Tleng(i)+1)),mean(regions.Tlat(i,1:regions.Tleng(i)+1)),regions.ID{i})
%       text(mean(regions.Tlon(i,1:regions.Tleng(i)+1)),mean(regions.Tlat(i,1:regions.Tleng(i)+1)),['#' num2str(i)])
    end
    if not(isempty(regsel))
       for j=1:length(regsel)
           i=regsel(j);
           plot(regions.Tlon(i,1:regions.Tleng(i)+1),regions.Tlat(i,1:regions.Tleng(i)+1),'-g')
           hold on
    %       text(mean(regions.Tlon(i,1:regions.Tleng(i)+1)),mean(regions.Tlat(i,1:regions.Tleng(i)+1)),regions.ID{i})
    %       text(mean(regions.Tlon(i,1:regions.Tleng(i)+1)),mean(regions.Tlat(i,1:regions.Tleng(i)+1)),['#' num2str(i)])
        end     
    end
    
end