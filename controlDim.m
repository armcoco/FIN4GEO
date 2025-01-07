function [varargout]=controlDim(nn,mm,varargin)

for i=1:nargin-2
    if length(varargin{i})==1
        varargout(i)={varargin{i}*ones(nn,mm)};
    else
        varargout(i)={varargin{i}};
    end
end

% function [Phi1,Phi2,sol_i,sol_e,Gamma_i,Gamma_e,F_i,F_e]=controlDim(Phi1,Phi2,sol_i,sol_e,Gamma_i,Gamma_e,F_i,F_e,nn,mm)
% 
% if length(Phi1)==1
%     Phi1=Phi1*ones(nn,mm);
% end
% 
% if length(Phi2)==1
%     Phi2=Phi2*ones(nn,mm);
% end
% 
% if length(sol_i)==1
%     sol_i=sol_i*ones(nn,mm);
% end
% 
% if length(sol_e)==1
%     sol_e=sol_e*ones(nn,mm);
% end
% 
% if length(Gamma_i)==1
%     Gamma_i=Gamma_i*ones(nn,mm);
% end
% 
% if length(Gamma_e)==1
%     Gamma_e=Gamma_e*ones(nn,mm);
% end
% 
% if length(F_i)==1
%     F_i=F_i*ones(nn,mm);
% end
% 
% if length(F_e)==1
%     F_e=F_e*ones(nn,mm);
% end
