function [Mask] = plane_mask(px,py,radius)
%PLANE_MASK Summary of this function goes here
%   Detailed explanation goes here

    Mask = zeros(length(px), length(py));
    for nX = 1:length(px)
        for nY = 1:length(py)
               
                if (px(nX)^2)+(py(nY)^2) <= radius^2 
                    Mask(nX,nY) = 1;
                end
            
        end
    end
    
end

