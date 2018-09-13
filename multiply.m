function [varargout] = multiply(multi,varargin)
%TIMES Summary of this function goes here
%   Detailed explanation goes here

numin = max(size(varargin));

    for no = 1:numin
        varargout{no} = varargin{no}*multi;
    end
end

