function [middledipole] = MiddleVariableLine(Bzdipole)
%MIDDLEVARIABLELINE Summary of this function goes here
%   Detailed explanation goes here
middledipole = zeros(1,size(Bzdipole,3));
sz = size(Bzdipole);

for o = 1:size(Bzdipole,3)  
    go = Bzdipole(:,:,o);
    go = reshape(go,sz(1),sz(2));
    
    middledipole(o) = go(round(sz(1)/2),round(sz(2)/2));
    
end 
end

