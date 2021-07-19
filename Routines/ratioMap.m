function [ratioVec]=ratioMap()
%  USE:  set(gca,'DataAspectRatio',ratioMap())
    xlim= get(gca,'XLim');ylim= get(gca,'YLim');
    lenX = distance(ylim(1),xlim(1),ylim(1),xlim(2));
    lenY = distance(ylim(1),xlim(1),ylim(2),xlim(1));
    ratioX = lenX/max(lenX,lenY);
    ratioY = lenY/max(lenX,lenY); 
    ratioZ = 1;
    ratioVec = [ratioX,ratioY,ratioZ];
end