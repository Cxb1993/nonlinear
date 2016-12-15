function [varargout] = eql_ex1_mlp_vel
    %% *SET-UP SOIL GEOLOGY*
    Vs0_layers  = [130;200;300;450;600;700;1250;2600]; % layers shear wave velocities [m/s]
    rho_layers  = [2050;2150;2075;2100;2155;2200;2500;2600]; 
    H_layers    = [7;13;34;23.5;50;59;10;103.5];
    f_max       = 1;  % maximum frequency to propagate
    n           = 1;
    dH_max      = H_layers/10;%Vs0_layers./n./f_max;
    N_layers    = floor(H_layers./dH_max);   % number of layers
    dH_layers   = H_layers./N_layers;
    temp_H = [];
    temp_V = [];
    temp_R = [];
    G_Gmax = [];
    D      = [];
    lab    = [];
    gamma = importdata('ggmax_0.txt');
    count = 0;
    for i = 1:numel(H_layers)
        temp_H = [temp_H;eql_ex0_multi_layers(dH_layers(i),N_layers(i),1)];
        temp_V = [temp_V;eql_ex0_multi_layers(Vs0_layers(i),N_layers(i),1)];
        temp_R = [temp_R;eql_ex0_multi_layers(rho_layers(i),N_layers(i),1)];
        if sum(temp_H)<sum(H_layers(1:i))
            temp_H(end)=temp_H(end)+sum(H_layers(1:i))-sum(temp_H);
        end
        if i<numel(H_layers)
            temp_G = importdata(sprintf('ggmax_%u.txt',i-1));
            temp_D = importdata(sprintf('ddmax_%u.txt',i-1));
            G_Gmax = [G_Gmax;repmat(temp_G(:,2)',[N_layers(i),1])];
            D      = [D;repmat(temp_D(:,2)',[N_layers(i),1])];
        end
        if i==numel(H_layers)
            G_Gmax = [G_Gmax;repmat(ones(size(gamma(:,1)))',[N_layers(i),1])];
            D = [D;repmat(0.0001*ones(size(gamma(:,1)))',[N_layers(i),1])];
        end
        for j=1:N_layers(i)
            lab = [lab;{sprintf('%u',count+j)}];
        end
        count = count+N_layers(i);
    end
    H_layers    = temp_H; clear temp_H;
    N_layers    = numel(H_layers);
    Vs0_layers  = temp_V; clear temp_V;
    rho_layers  = temp_R; clear temp_R;
    gamma = repmat(gamma(:,1)',N_layers,1);
%     G_Gmax = ones(size(G_Gmax));
%     D = 0.0001*ones(size(G_Gmax));
    xsi0_layers = 0.0001*ones(N_layers,1); 
    plot_vs_profile({H_layers},{Vs0_layers},{'$V_S$'},'Volvi','[Kg/m^3]');
    plot_vs_profile({H_layers},{rho_layers},{'$\rho_S$'},'Volvi','[Kg/m^3]');
    plot_ggd(gamma(1,:),G_Gmax,D,lab);
    
    
    Dinf        = .20;  % large strain damping              [1]
    z_water     =  300;  % ground water depth                [m]
    
    %% *SET-UP RIGID BEDROCK GEOLOGY (WITHIN CONDITION)*
    CI        = 1;
    Vs0_rock  = CI*Vs0_layers(end);  % bedrock Vs value                  [m/s]
    rho_rock  = rho_layers(end);     % bedrock unit weight               [kg/m3]
    xsi0_rock = .0001;               % bedrock damping coefficient       [1]
    
    
    mat = struct('H_layers',H_layers,'N_layers',N_layers,...
        'rho_layers',rho_layers,'rho_rock',rho_rock,...
        'Vs0_layers',Vs0_layers,'xsi0_layers',xsi0_layers,...
        'Vs0_rock',Vs0_rock,'xsi0_rock',xsi0_rock,...
        'G_Gmax',G_Gmax,'D',D,'gamma',gamma);
    
    %% *OUTPUT*
    varargout{1} = mat;
    
    return
end