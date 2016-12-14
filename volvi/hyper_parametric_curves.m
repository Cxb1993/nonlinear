function [varargout]=hyper_parametric_curves(varargin)
    % NOTES
    % hyper_parametric_curve: function returning G/Gmax-gamma-D curve according
    % to HYPERBOLIC model (1/2 parameters: alpha-beta). Stress
    % correction is applied.
    %
    % G/Gmax = 1 / ( 1 + alpha * abs( gamma ) ^ beta )
    % D/Dinf = 1 - G/Gmax
    %
    % INPUT:  name (reference string for selected curve pair)
    %         alpha,beta (coefficients for hyperbolic model)
    %
    % OUTPUT: G_Gmax (G/Gmax curve)
    %         D (D-D0 model)
    %
    % REFERENCE:
    % "Nakagawa, K. and Soga, K. (1995). Nonlinear cyclic strees - strain relations of soils.
    % In 3rd International Conference on recent Advances in Geotechnical
    % Earthquake Engineering and Soil Dynamics, volume I. St. Louis, Missouri.
    % Eds. Prakash."
    %===============
    
    name  = varargin{1};
    alpha = varargin{2};
    beta  = varargin{3};
    gi    = varargin{4};
    Dinf  = varargin{5};
    if nargin>5
        sv    = varargin{6};
    end
    pa = 101.325;
    
    switch(name)
        
        case 'Nakagawa-Soga' % Nakagawa & Soga curves
            
            G_Gmax=1./(1+alpha.*abs(gi).^beta); % G/Gmax
            D=Dinf.*(1-G_Gmax);                 % D-D0
            
        case 'Mod. Nakagawa-Soga'
            
            G_Gmax=1./(1+alpha.*(abs(gi)./sqrt(sv/pa)).^beta); % G/Gmax
            D=Dinf.*(1-G_Gmax);                               % D-D0
            
    end
    
    varargout{1}=G_Gmax;
    varargout{2}=D;
    return
end