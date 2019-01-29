function [HxAkoun, HyAkoun, HzAkoun] = Jannsen(Rx,Ry,Rz,cell_size)

cell_size = cell_size./2;

HxAkoun  = 0;
HyAkoun  = 0;
HzAkoun  = 0;

             for k = 0:1
                for l = 0:1
                    for m = 0:1
                        
                        S = Rx - ((-1)^k)*cell_size(1);
                        T = Ry - ((-1)^l)*cell_size(2);
                        U = Rz - ((-1)^m)*cell_size(3);
                        R = sqrt(S^2+T^2+U^2);
                 
                        dHxAkoun = ((-1)^(k+l+m))* log(R-T);
                        dHyAkoun = ((-1)^(k+l+m))* log(R-S);
                        dHzAkoun = ((-1)^(k+l+m))* atan((S*T)/(U*R));
                                                                   
                        HxAkoun  = HxAkoun + (dHxAkoun);
                        HyAkoun  = HyAkoun + (dHyAkoun);
                        HzAkoun  = HzAkoun + (dHzAkoun);
                     
                    end
                end 
             end

end

