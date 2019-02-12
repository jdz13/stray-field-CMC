function [rotxx,rotyy,rotzz] = rotation_matricies(alpha,beta,gamma)
%ROTATION_MATRICIES Summary of this function goes here
%   Detailed explanation goes here
rotxx = [1,0,0; 0, cos(alpha),sin(alpha); 0, -sin(alpha), cos(alpha)];
rotyy = [cos(beta), 0 ,-sin(beta); 0,1,0; sin(beta), 0, cos(beta)];
rotzz = [cos(gamma),sin(gamma),0; -sin(gamma),cos(gamma),0; 0,0,1];
end

