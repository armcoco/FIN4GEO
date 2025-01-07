function sxt = xt__(xi)

global c_parameterx m_parx ISINF;

if ~exist('c_parameterx','var') || isempty(c_parameterx)
    c_parameterx = 1;
end

if ISINF
    sxt = xi./(1.0-xi.^2).^m_parx;
    sxt=c_parameterx*sxt;
else
    sxt=c_parameterx*xi;
end

if ISINF
    sxt = xi./(1.0-xi.^2).^m_parx;
    sxt=c_parameterx*sxt;
else
    sxt=c_parameterx*xi;
end

