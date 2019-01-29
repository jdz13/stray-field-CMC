function [Mask] = shapes (magdim, cell_size, shape,px,py,pz,disp)

magdim = magdim/2;

%error code if the magnet is smaller than any cell dimension
if any(lt((magdim - cell_size), 0)) == 1
    error('A dimension of the particle is smaller than the cell')
end 
    

% error code for if input is the wrong size 
if any(ne(rem(magdim,cell_size),[0,0,0])) == 1
    error('Particle is not an integer number of cells. Check your inputs.')
end

% error code for if the input is in the wrong position (no cell overlap). 



squares
 
if shape == 'square'
    
  Mask = zeros(length(px), length(py),length(pz));
    for nX = 1:length(px)
        for nY = 1:length(py)
            for nZ = 1:length(pz)     
                if (px(nX)) <= disp(1)+magdim(1) && (px(nX)) >= disp(1)-magdim(1)
                    if (py(nY)) <= disp(2)+magdim(2) && (py(nY)) >= disp(2)-magdim(2)
                        if (pz(nZ)) <= disp(3)+magdim(3) && (pz(nZ)) >= disp(3)-magdim(3)
                            
                            Mask(nX,nY,nZ) = 1;
                        
                        end 
                    end
                end
            end
        end 
    end 
end


%-------------------------------------------------------------------------

% circles 
if shape == 'sphere'
   
    
    Mask = zeros(length(px), length(py),length(pz));
    for nX = 1:length(px)
        for nY = 1:length(py)
            for nZ = 1:length(pz)     
                if (px(nX)^2)+(py(nY)^2) +(pz(nZ)^2)<= magdim(1)^2
                    Mask(nX,nY) = 1;
                end
            end
        end
    end
    
end 

%------------------------------------------------------------------------

% discs 
if shape == 'disc-z'
   
    
    Mask = zeros(length(px), length(py),length(pz));
    for nX = 1:length(px)
        for nY = 1:length(py)
            for nZ = 1:length(pz)     
                if (px(nX)^2)+(py(nY)^2) <= magdim(1)^2 && (pz(nZ)) <= magdim(2)
                    Mask(nX,nY) = 1;
                end
            end
        end
    end
    
end 


end 
