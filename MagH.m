function [HxMagH, HyMagH, HzMagH] = MagH(Rx,Ry,Rz,cell_size)
%MAGH Summary of this function goes here
%   Detailed explanation goes here
             for k = 1:2
                for l = 1:2
                    for m = 1:2
                        
                        Xfac = Rx + ((-1)^k)*cell_size(1);
                        Yfac = Rz + ((-1)^l)*cell_size(3);
                        Zfac = Ry + ((-1)^m)*cell_size(2);
                        SqrtXYZ = sqrt(Xfac^2+Yfac^2+Zfac^2);
                        
                        dHxMagH = (-1)^(k+l+m)*log(Zfac+SqrtXYZ);
                      
                        dHzAMagH = (Yfac*Xfac)/((abs(Yfac))*(abs(Xfac)));
                        dHzBMagH = atan((dot(abs(Xfac),Zfac))/(dot(abs(Yfac),SqrtXYZ)));
                        dHzMagH = ((-1)^(k+l+m))*dHzAMagH*dHzBMagH;
                        
                        dHyMagH = (-1)^(k+l+m)*log(Xfac+SqrtXYZ);

                        
                        HxMagH  = HxMagH + 1/4/pi()*dHxMagH;
                        HzMagH  = HzMagH - 1/4/pi()*dHzMagH;
                        HyMagH  = HyMagH + 1/4/pi()*dHyMagH;
                     
                    end
                end 
            end
end

