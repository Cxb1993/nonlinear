%===============
% Editor: Filippo Gatti
% Ecole Centrale Paris - Laboratoire MSSMat
% Copyright 2015
% NOTES
% get_1d_stress_state: function to evaluate effective vertical/horizontal/mean
% stress given a horizontal multi-layer soil profile.
% INPUT:  H_strata (vector containing heights of each stratum)
%         rho_strata (vector containing unit weight of each strata)
%         z_stress (depth at which vertical stress is evaluated)
%         z_water (groundwater level (optional))
%         k0 (coefficient of lateral stress (optional))
% OUTPUT: p_eff (effective mean stress)
%         p_tot (total mean stress)
%         sigma_v_eff (effective vertical stress state)
%         sigma_h_eff (effecgive horizontal stress state)
%         sigma_v_tot (total vertical stress state)
%         sigma_h_tot (total horizontal stress state)
%===============

function [varargout]=get_1d_stress_state(varargin)
    
    H_strata    = [0 varargin{1}(:).'];                                         % heights [m] of each strata (row column)
    rho_strata  = [0 varargin{2}(:).'];                                         % unit weights [kg/m3] of each strata (row column)
    z_strata    = cumsum(H_strata);                                             % depth of each interface between strata [m]
    z_stress    = varargin{3};                                                  % depth [m] at each evaluate stress state
    
    g=10;                                                                       % gravity acceleration [m/s2]
    % groundwater depth [m]
    rho_water=1000;                                                             % water unit weight [kg/m3]
    if nargin>=4
        z_water=varargin{4};
    else
        z_water=0;
    end
    
    % k0=sigma_H/sigma_V
    if nargin>=5
        k0=varargin{5};
    else
        k0=.7;
    end
    
    % locating z_stress
    idx=find(z_strata<=z_stress,1,'last');
    
    sigma_v_tot=sum(g*rho_strata(1:idx).*H_strata(1:idx));
    
    if z_strata(idx)<z_stress
        sigma_v_tot=sigma_v_tot+g*rho_strata(idx+1)*(z_stress-z_strata(idx));   % total vertical stress state
    end
    
    sigma_v_eff=sigma_v_tot-g*rho_water*(z_stress-z_water)*(z_stress>z_water);  % effective vertical stress state
    
    sigma_h_eff=sigma_v_eff*k0;
    
    p_eff=(sigma_v_eff+2*sigma_h_eff)/3;
    
    sigma_h_tot=sigma_h_eff+(z_stress>z_water)*g*rho_water*(z_stress-z_water);
    
    p_tot=p_eff+(z_stress>z_water)*g*rho_water*(z_stress-z_water);
    
    varargout{1}=p_eff/1000;
    varargout{2}=p_tot/1000;
    varargout{3}=sigma_v_eff/1000;
    varargout{4}=sigma_h_eff/1000;
    varargout{5}=sigma_v_tot/1000;
    varargout{6}=sigma_h_tot/1000;
    
end

