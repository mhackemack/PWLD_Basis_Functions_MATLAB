function varargout = emptyFunction(varargin)
if nargout > 0
    for i=1:nargout
        varargout{i} = 0;
    end
end