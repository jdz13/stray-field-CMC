function [outputArg1,outputArg2] = MagH(inputArg1,inputArg2)
%MAGH Summary of this function goes here
%   Detailed explanation goes here
             for k = 1:2
                for l = 1:2
                    for m = 1:2
                        
                        Xfac = a + ((-1)^k)*xb;
                        Yfac = c + ((-1)^l)*zb;
                        Zfac = b + ((-1)^m)*yb;
                        SqrtXYZ = sqrt(Xfac^2+Yfac^2+Zfac^2);
                        
                        dHxMagH = (-1)^(k+l+m)*log(Zfac+SqrtXYZ);
                      
                        dHzAMagH = (Yfac*Xfac)/((abs(Yfac))*(abs(Xfac)));
                        dHzBMagH = atan((dot(abs(Xfac),Zfac))/(dot(abs(Yfac),SqrtXYZ)));
                        dHzMagH = ((-1)^(k+l+m))*dHzAMagH*dHzBMagH;
                        
                        dHyMagH = (-1)^(k+l+m)*log(Xfac+SqrtXYZ);
                        
                        unit1 = round((a+xa)*(div/2/xa))+1;
                        unit2 = round((c+za)*(div/2/za))+1;
                        unit3 = round((b+ya)*(div/2/ya))+1;
                        
                        HxMagH (unit1,unit2,unit3) = HxMagH(unit1,unit2,unit3) + 1/4/pi()*dHxMagH*M;
                        HzMagH (unit1,unit2,unit3) = HzMagH(unit1,unit2,unit3) - 1/4/pi()*dHzMagH*M;
                        HyMagH (unit1,unit2,unit3) = HyMagH(unit1,unit2,unit3) + 1/4/pi()*dHyMagH*M;
                     
                    end
                end 
            end
end

