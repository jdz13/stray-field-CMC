function [numin,numout, varargout] = field_comps(varargin)
%FIELDS Summary of this function goes here
%   Detailed explanation goes here
numin = max(size(varargin));
numout = nargout;

    for no = 1:numin
        varargout{no} = no;
    end
end

