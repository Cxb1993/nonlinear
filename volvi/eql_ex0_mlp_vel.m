function [varargout] = eql_ex0_mlp_vel
    %% *SET-UP SOIL GEOLOGY*
    Vs0_layers  = 100; % layers shear wave velocities [m/s]
    H_layers    = 1000;
    f_max       = 2000;  % maximum frequency to propagate
    n           = 10;
    dH_max      = 10;%Vs0_layers./n./f_max;
    N_layers    = floor(H_layers./dH_max);   % number of layers
    dH_layers   = H_layers./N_layers;
    temp_H = [];
    temp_V = [];
    for i = 1:numel(H_layers)
        temp_H = [temp_H;eql_ex0_multi_layers(dH_layers(i),N_layers(i),1)];
        temp_V = [temp_V;eql_ex0_multi_layers(Vs0_layers(i),N_layers(i),1)];
        if sum(temp_H)<sum(H_layers(1:i))
            temp_H(end)=temp_H(end)+sum(H_layers(1:i))-sum(temp_H);
        end
    end
    H_layers    = temp_H; clear temp_H;
    N_layers    = numel(H_layers);
    Vs0_layers  = temp_V; clear temp_V;
     plot_vs_profile({H_layers},{Vs0_layers},{'VS-PROFILE'});
    
    % layers density
    rho_layers  = 1800*ones(N_layers,1); 
    % layers critical damping ratios
    xsi0_layers = 0.0001*ones(N_layers,1);                     
%     plot_vs_profile({H_layers},{Vs0_layers},{'D-PROFILE'});
    
    Dinf        = .20;  % large strain damping              [1]
    z_water     =  10;  % ground water depth                [m]
    
    %% *SET-UP RIGID BEDROCK GEOLOGY (WITHIN CONDITION)*
    CI        = 1;
    Vs0_rock  = CI*Vs0_layers(end);  % bedrock Vs value                  [m/s]
    rho_rock  = rho_layers(end);     % bedrock unit weight               [kg/m3]
    xsi0_rock = .0001;               % bedrock damping coefficient       [1]
    
    %% SET-UP G/G0-GAMMA-D CURVES
    
    gi=(1e-5:1e-4:10)./100;    
    G_Gmax = ones(N_layers,numel(gi));
    D      = xsi0_layers*ones(1,numel(gi));    
    gamma=repmat(gi,N_layers,1);
    
    mat = struct('H_layers',H_layers,'N_layers',N_layers,...
        'rho_layers',rho_layers,'rho_rock',rho_rock,...
        'Vs0_layers',Vs0_layers,'xsi0_layers',xsi0_layers,...
        'Vs0_rock',Vs0_rock,'xsi0_rock',xsi0_rock,...
        'G_Gmax',G_Gmax,'D',D,'gamma',gamma);
    
    %% *OUTPUT*
    varargout{1} = mat;
    
    return
end