function out1 = yt1__(eta)

global c_parameter m_par ISINF;

if ISINF
    out1=(eta.^2*(2*m_par-1)+1)./(1-eta.^2).^(m_par+1);
    out1=c_parameter*out1;
else
    out1=c_parameter*ones(size(eta));
end
