function outVal = scalinglaw_Murotani(typeScal,M)
    a =-3.806;
    b =1.000;         
    areaTmp = 10.^(a+b*M); 
    if strcmp(typeScal,'M2W')
      outVal = sqrt(2.5*areaTmp)/2.5;
    elseif strcmp(typeScal,'M2L') 
      outVal = sqrt(2.5*areaTmp);
    else
      disp(['Type ' typeScal ' not recongnized!!']);
    end  
end




