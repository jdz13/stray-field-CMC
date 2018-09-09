function [out1,out2,out3] = quiver_filter(filter,argX,argY,argZ)
%QUIVER_FILTER Summary of this function goes here
%   Detailed explanation goes here


    filterx = argX<= filter(1);
    filtery = argY<= filter(2);
    filterz = argZ<= filter(3);
    filter = filterx.*filtery.*filterz;

    out1 = filter.*argX;
    out2 = filter.*argY;
    out3 = filter.*argZ;

end

