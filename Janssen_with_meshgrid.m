function [Akoun] = Janssen_with_meshgrid(xprime,yprime,zprime,mag_size)
%JANSSEN_WITH_MESHGRID Summary of this function goes here
%   Detailed explanation goes here

for nX = 1:size(xprime,1)
    for nY = 1:size(yprime,2)
        for nZ = 1:size(zprime,3)
            % Running the function for each variable
            [Akoun.HxAkoun(nX,nY,nZ), Akoun.HyAkoun(nX,nY,nZ), ...
                Akoun.HzAkoun(nX,nY,nZ)] = Jannsen(xprime(nX,nY,nZ),yprime(nX,nY,nZ),...
                zprime(nX,nY,nZ),mag_size);
        end 
    end 
end

Akoun.modBAkoun = sqrt(Akoun.HxAkoun.^2 + Akoun.HyAkoun.^2 + Akoun.HzAkoun.^2);
end

