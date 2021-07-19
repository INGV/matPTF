function outVal = scalinglaw_WC(typeScal,M)
    
    if strcmp(typeScal,'M2L')
      a =-2.440;
      b =0.590;
      outVal = 10.^(a+b*M);
    elseif strcmp(typeScal,'M2W') 
      a=-1.010;
      b=0.320;
      outVal = 10.^(a+b*M);
    else
      disp(['Type ' typeScal ' not recongnized!!']);
    end  

end

    
    

