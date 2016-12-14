%===============
% Editor: Filippo Gatti
% Ecole Centrale Paris - Laboratoire MSSMat
% Copyright 2015
% NOTES
% eql_acc_strain: function to perform Equivalent Linear Analysis (EQL)
% (calculations in frequency domain)
% INPUT:  H_layers (row vector containing the thicknesses of half-layers set)
%         Vs_layers (row vector containing shear wave velocity values for half-layers set)
%         rho_layers (row vector containing half-layers densities)
%         xsi_layers (row vector containing half-layers damping)
%         Vs_rock (shear wave velocity value)
%         rho_rock (bedrock density)
%         xsi0_rock (bedrock initial damping)
%         inacc_th (row vecotr containing input acceleration time history)

% OUTPUT: outacc_th (matrix of output acceleration response for all half-layers => row=nr. half-layers - column=time step)
%         gamma_th (matrix of strain time-histories for all half-layers => row=nr. half-layers - column=time step)
%         inacc_th (row vector containing processed input signal (CyberQuake damping))
%         freq (row vector containing output frequencies)
%         TF_to_base (matrix containing impendance matrix)
% References: Kramer S.J. - Geotechnical Earthquake Engineering (1996) -
%               pagg. 268-269
%             Kwok et al. - Use of Exact Solutions of Wave Propagation Problems to Guide
%               Implementation of Nonlinear Seismic Ground Response Analysis Procedures (2007)
%===============

function [varargout] = eql_ex0_acc_strain(varargin)
    
    %======================================================================
    % SET-UP
    %======================================================================
    H_layers   = varargin{1}(:);                                            % layer thicknesses
    Vs0_layers = varargin{2}(:);                                            % layer INITIAL shear wave velocities
    Vs_layers  = varargin{3};                                               % layer CURRENT (ITERATIVE) shear wave velocities
    xsi_layers = varargin{4};                                               % layer CURRENT (ITERATIVE) damping
    rho_layers = varargin{5}(:);                                            % layer densities
    %
    Vs0_rock   = varargin{6};                                               % bedrock shear wave velocity
    xsi0_rock  = varargin{7};                                               % bedrock initial damping
    rho_rock   = varargin{8};                                               % bedrock density
    %
    gamma_G    = varargin{9};                                               % shear strain for G/Gmax-gamma curves
    GGmax      = varargin{10};                                              % G/Gmax curves
    gamma_D    = varargin{11};                                              % shear strain for D-gamma curves
    Damp       = varargin{12};                                              % D curves
    %
    inacc_th  = varargin{13};                                               % input acceleration time history
    inacc_fft = varargin{14};                                               % input acceleration fft
    d_freq    = varargin{15};                                               % input acceleration frequency-step
    %
    alpha_eff = varargin{16};                                               % coefficient to compute gamma_eff
    flag_freq = varargin{17};                                               % flag for frequency dependent properties
    %
    N_layers  = numel(H_layers);                                            % number of layers
    N_t       = numel(inacc_th);                                            % number of time points
    N_f       = numel(inacc_fft);                                           % number of frequency points
    freq      = (0:(N_f-1))*d_freq;                                         % frequency
    %======================================================================
    % COMPLEX IMPEDANCE FUNCTIONS
    %======================================================================
    
    if flag_freq~=0  % check for frequency dependency
        max_freq   = 25;
        N_omega    = round(max_freq/d_freq)+1;
        omega      = 2*pi*d_freq*(0:N_omega-1);
        if size(xsi_layers,2)==1 % initialize xsi vector
            xsi_layers = repmat(xsi_layers,1,N_f);
        end
    else
        N_omega    = 1;
    end
    
    %% *IMPEDANCE COEFFICIENTS*
    alpha_z = zeros(N_layers,N_f); % N_layers * N_f
    
    for ii = 1 : size(xsi_layers,2)
        
        NOM     = (rho_layers(1:end-1).*Vs_layers(1:end-1).*(1+1i*xsi_layers(1:end-1,ii)));
        DENOM   = (rho_layers(2:end).*Vs_layers(2:end).*(1+1i*xsi_layers(2:end,ii)));
        alpha_z_temp = NOM./DENOM;
        
        NOM     = (rho_layers(end).*Vs_layers(end).*(1+1i*xsi_layers(end,ii)));
        DENOM   = (rho_rock.*Vs0_rock.*(1+1i*xsi0_rock));
        
        alpha_z(:,ii) = [alpha_z_temp;NOM./DENOM];
        
    end
    
    % N.B.:
    % rigid bedrock "within": Vs0_rock >> and rho_rock >>
    % alpha_z(end)~=0
    
    %======================================================================
    % COMPLEX WAVE NUMBERS
    %======================================================================
    K  = zeros(N_layers+1,N_f);
    KH = zeros(N_layers+1,N_f);
    
    NOM   = repmat(2.*pi.*freq(:).',N_layers,1); % 1 * N_f
    NOMH  = (2.*pi.*H_layers(:)*freq(:).');      % 1 * N_f
    
    if flag_freq~=0
        for ii = 1 : N_f
            DENOM = Vs_layers(:).*(1+xsi_layers(:,ii).*1i);
            K(1:end-1,ii)  = NOM(:,ii)./DENOM;
            KH(1:end-1,ii) = NOMH(:,ii)./DENOM;
        end
    else
        DENOM = repmat(Vs_layers(:).*(1+xsi_layers(:).*1i),1,N_f);
        K(1:end-1,:)  = NOM./DENOM;
        KH(1:end-1,:) = NOMH./DENOM;
    end
    %======================================================================
    
    %======================================================================
    % LAYER'S AMPLITUDES
    %======================================================================
    % defining amplitude coefficients En and Fn for each soil layer
    % free surface condition E1=F1
    
    E = ones(1,N_f);                                                          % upgoing wave amplitude
    F = ones(1,N_f);                                                          % downgoing wave amplitude
    
    for j = 2 : N_layers+1
        if flag_freq~=0
            for ii = 1 : N_f
                E(j,ii_) = .5.*(E(j-1,:).*(1+alpha_z(j-1)).*exp( 1i.*KH(j-1,:)) + ...
                    F(j-1,:).*(1-alpha_z(j-1)).*exp(-1i.*KH(j-1,:)));
                F(j,:) = .5.*(E(j-1,:).*(1-alpha_z(j-1)).*exp( 1i.*KH(j-1,:)) + ...
                    F(j-1,:).*(1+alpha_z(j-1)).*exp(-1i.*KH(j-1,:)));
            end
        else
            E(j,:) = .5.*(E(j-1,:).*(1+alpha_z(j-1)).*exp( 1i.*KH(j-1,:)) + ...
                F(j-1,:).*(1-alpha_z(j-1)).*exp(-1i.*KH(j-1,:)));
            F(j,:) = .5.*(E(j-1,:).*(1-alpha_z(j-1)).*exp( 1i.*KH(j-1,:)) + ...
                F(j-1,:).*(1+alpha_z(j-1)).*exp(-1i.*KH(j-1,:)));
        end
    end
    
    U = E + F;                                                              % total wave field amplitudes in each soil layer
    gamma_tf = 1i.*K.*(E.*exp(1i.*KH)-F.*exp(-1i.*KH));                     % strain fourier spectra in each soil layer
    %======================================================================
    
    %======================================================================
    % OUTCROPPING -> BEDROCK
    %======================================================================
    TF_to_base = U./repmat(U(end,:),N_layers+1,1);                          % layer-to-bedrock
    %======================================================================
    % EERA METHOD
    %======================================================================
    TF_BO = 0.5*repmat(U(end,:)./E(end,:),N_layers+1,1);                    % outcrop-to-bedrock
    TF_to_base=TF_to_base.*TF_BO;
    %======================================================================
    % FREQUENCY RESPONSE
    %======================================================================
    outacc_fft = repmat(inacc_fft,N_layers+1,1).*TF_to_base;
    %======================================================================
    
    %======================================================================
    % TIME RESPONSE
    %======================================================================
    outacc_th = real(ifft(outacc_fft,[],2,'symmetric'));
    outacc_th = outacc_th(:,1:N_t);
    %======================================================================
    
    %======================================================================
    % GAMMA FIELD
    %======================================================================
    %% NOM   = ? [TODO]
    %% DENOM = ? [TODO]
    
    gamma_fft = NOM./DENOM;
    
    gamma_th = real(ifft(gamma_fft(1:N_layers,:),[],2,'symmetric'));
    gamma_th = gamma_th(:,1:N_t).*9.81;
    %======================================================================
    
    %======================================================================
    % COMPUTE GAMMA EFF
    %======================================================================
    gamma_eff = -1*ones(N_layers/2,size(xsi_layers,2));
    
    for ii = 1 : N_layers/2
        if flag_freq~=0
            gamma_spectrum  = abs(gamma_fft(ii,1:N_omega));
            geff = eql_ex0_gamma_eff(flag_freq,omega,...
                gamma_spectrum,alpha_eff);
            
            gamma_eff(ii,:) = [geff zeros(1,N_f-2*N_omega) flip(geff)];
            
        else
            gamma_eff(ii,:) = eql_ex0_gamma_eff(flag_freq,gamma_G(ii,end),...
                gamma_th(2*ii,:),alpha_eff);
        end
        
    end
    %======================================================================
    
    %======================================================================
    % CHECK CONVERGENCE
    %======================================================================
    
    [sumerr,Vs_layers,xsi_layers,err,status] = ...
        eql_ex0_check_conv(gamma_eff,Vs0_layers(:),Vs_layers(:),...
        gamma_G,GGmax,...
        gamma_D,Damp,flag_freq);
    %======================================================================
    
    %======================================================================
    % OUTPUT
    %======================================================================
    varargout{1} = sumerr;
    varargout{2} = outacc_th;
    varargout{3} = outacc_fft;
    varargout{4} = gamma_th;
    varargout{5} = Vs_layers;
    varargout{6} = xsi_layers;
    varargout{7} = gamma_eff;
    varargout{8} = err;
    varargout{9} = status;
    %======================================================================
    return
end