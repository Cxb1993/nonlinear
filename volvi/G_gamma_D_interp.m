%% *G-gamma-D curve interpolation*
% Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% Politecnico di Milano - DICA
% Copyright 2015
%% NOTES
% _G_gamma_D_interp_: function to interpolate G/Gmax-gamma-D curves. 
% G/Gmax is interpolated according to Hyperbolic or Nakagawa-Soga models. 
% D is interpolated via spline method.
% INPUT
% * _g (gamma sampled values)_
% * _G (G/G0 sampled values)_
% * _D (D sampled values)_
% * _gi (interpolating gamma values)_
% * _D0 (small strain damping ratio)_
% * _Dinf (large strain damping ratio)_
% * _ctrl_gi (strain range limits)_
%% OUTPUT
% * _Gi (interpolated G/Gmax values)_
% * _Di (interpolated D values)_
%% REFERENCE
% Nakagawa, K. and Soga, K. (1995). Nonlinear cyclic strees - strain relations of soils.
% In 3rd International Conference on recent Advances in Geotechnical
% Earthquake Engineering and Soil Dynamics, volume I. St. Louis, Missouri. Eds. Prakash.

function [varargout] = G_gamma_D_interp(varargin)
    
    %% *SET-UP*
    g       = varargin{1}(:).';   % gamma sampled values [1]
    G       = varargin{2}(:).';   % G/G0 sampled values [1]
    D       = varargin{3}(:).';   % D sampled values [1]
    D       = D(~isnan(D));
    gi      = varargin{4}(:).';   % interpolating gamma values)
    if nargin == 7
        D0      = varargin{5};         % small strain D [1]
        Dinf    = varargin{6};         % large strain D [1]
        ctrl_gi = varargin{7}(:).';    % gamma range [1]
    else
        D0      = D(1);
        Dinf    = D(end);
        ctrl_gi = g(end) + g(end)/100.*[10 20];
    end
    
    method = 'NS95';
    if nargin == 5
        method = upper(varargin{5});
    end
    %% *INTERPOLATION G/G0 CURVE
    % G_Gmax curve interpolated with original curve points (no 0/1 added at
    % the beginning or at the end)
    [alpha_G,beta_G] = alphabeta(g,G,method);
    Gi = 1./(1+alpha_G*(abs(gi).^beta_G));
    
    %% *INTERPOLATION D CURVE
    % D curve interpolated with adding 0/Dinf values at
    % the beginning or at the end)
    temp = ~isnan(D);
    try
        gD=g(temp); D = D(temp);
    catch
        keyboard
    end
    if g(end) < ctrl_gi(1)
        gD = [gD ctrl_gi];
        D  = [D Dinf*ones(size(ctrl_gi))];
    end
    
    if g(1) > gi(1)
        gD = [gi(1) gD];
        D  = [D0     D];
    end
    
    Di = pchip(log10(gD),D,log10(gi));
    
    %% *OUTPUT*
    varargout{1} = Gi;
    varargout{2} = Di;
    
    return
end