function syt = yt__(eta)

global c_parameterz m_par ISINF;

if ~exist('c_parameterz','var') || isempty(c_parameterz)
    c_parameterz = 1;
end

if ISINF
    syt = eta./(1.0-eta.^2).^m_par;
    syt=c_parameterz*syt;
%     syt=c_parametery*syt-6500; %for Gottsmann 23 Jan 2018
else
    syt=c_parameterz*eta;
end

