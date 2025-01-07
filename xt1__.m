function out1 = xt1__(xi)

global c_parameterx m_parx ISINF;

if ISINF
    out1=(xi.^2*(2*m_parx-1)+1)./(1-xi.^2).^(m_parx+1);
    out1=c_parameterx*out1;
else
    out1=c_parameterx*ones(size(xi));
end
