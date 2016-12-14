% ===============
% Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% Politecnico di Milano - DICA
% Copyright 2015
% Notes
% eql_multi_reflection: function to perform Equivalent Linear Analysis (EQL)
% INPUT:  H_layers (row vector containing the thicknesses of half-layers set)
%         Vs0_layers (row vector containing shear wave velocity values for half-layers set)
%         rho_layers (row vector containing half-layers densities)
%         xsi0_layers (row vector containing half-layers initial damping)
%         Vs0_rock (shear wave velocity value)
%         rho_rock (bedrock density)
%         xsi0_rock (bedrock initial damping)
%         gamma_G (matrix containing shear strain values for G/Gmax-gamma degradation curves (one per row))
%         G (matrix containing G/Gmax values for G/Gmax-gamma degradation curves (one per row))
%         gamma_D (matrix containing shear strain values for D-gamma degradation curves (one per row))
%         D (matrix containing D values for D-gamma degradation curves (one per row))
%         inacc_th (row vector containing input acceleration time history)
%         d_time (input accelerogram time-step)
%         alpha_eff (coefficient for gamma_eff)
%         flag_freq (flag for frequency dependent properties)
% OUTPUT: output_acc (matrix of output acceleration response for all layers)
%         gamma_th   (matrix of strain time-histories for all layers)
%         Vs_layers  (final values of shear wave velocity)
%         xsi_layers (final values of damping ratio)
% ===============

function [varargout] = eql_ex0_multi_reflection(varargin)
    
    %======================================================================
    % INPUT PARAMETERS
    %======================================================================
    if nargin==1
        varargin = struct2cell(varargin{1});
    end
    
    H_layers    = varargin{1}(:);                                           % layers thicknesses
    Vs0_layers  = varargin{2}(:);                                           % layers INITIAL shear wave velocities
    rho_layers  = varargin{3}(:);                                           % layers unit weight
    xsi0_layers = varargin{4}(:);                                           % layers INITIAL damping
    Vs0_rock    = varargin{5};                                              % bedrock shear wave velocity
    rho_rock    = varargin{6};                                              % bedrock unit weight
    xsi0_rock   = varargin{7};                                              % bedrock initial damping
    gamma_G     = varargin{8};                                              % shear strain values for G/Gmax-gamma curves
    GGmax       = varargin{9};                                              % G/Gmax curves
    gamma_D     = varargin{10};                                             % shear strain values for D-gamma curves
    Damp        = varargin{11};                                             % D curves
    inacc_th    = varargin{12}(:).';                                        % input acceleration time history [g]
    d_time      = varargin{13};                                             % input acceleration time-step
    alpha_eff   = varargin{14};                                             % coefficient for gamma_eff
    flag_freq   = varargin{15};                                             % flag for frequency dependency
    fig_name    = varargin{16};
    N_layers    = numel(H_layers);
    z_layers    = [0;-eql_ex0_multi_layers(cumsum(H_layers),2,1)];
    
    % doubling vectors: 2*j-1 refers to stratum's surface, 2*j refers to
    % stratum mid-thickness
    
    % initial shear wave velocities and damping
    Vs_layers = -1*ones(N_layers,1);
    
    for m = 1 : N_layers
        Vs_layers(m) = G_gamma_D_extrap(gamma_G(m,:),GGmax(m,:),gamma_G(m,1))*...
            Vs0_layers(m);
    end
    
    H_layers   = eql_ex0_multi_layers(H_layers,2,2);
    Vs0_layers = eql_ex0_multi_layers(Vs0_layers,2,1);
    Vs_layers  = eql_ex0_multi_layers(Vs_layers,2,1);
    rho_layers = eql_ex0_multi_layers(rho_layers,2,1);
    xsi_layers = eql_ex0_multi_layers(xsi0_layers,2,1);
    %======================================================================
    
    %======================================================================
    % PLOT PROFILE
    %======================================================================
    
    %     fprintf('NUMBER OF LAYERS: %.0f\n\n',N_layers);
    %     fprintf('GEOLOGY: \n\n');
    %
    %     for j=1:N_layers
    %         fprintf('STRATUM %u : z = [%.3f;%.3f] m - H = %.3f m - Vs = %.3f m/s - xsi = %.3f - rho = %.3f kg/m3\n',...
    %             j,z_layers(2*j-1),z_layers(2*j),H_layers(j),Vs0_layers(j),xsi0_layers(j),rho_layers(j));
    %     end
    
    %     fprintf('\nBEDROCK: \n\n');
    %     fprintf('z < %.3f - Vs = %.3f m/s - xsi = %.3f - rho = %.3f kg/m3\n\n',...
    %         z_layers(end),Vs0_rock,xsi0_rock,rho_rock);
    %
    %     plot_set_up;
    %     fig_geo=figure('name',sprintf('profile_%s',fig_name),'visible','on');
    %     hax_geo=axes('parent',fig_geo,'xtick',100:100:600);
    %     hold(hax_geo,'all'); box(hax_geo,'on');
    %     plot(hax_geo,[Vs0_layers;Vs0_layers(end)],z_layers,'b','displayname','V_S_0');
    %     xlabel(hax_geo,'V_S [m/s]'); ylabel(hax_geo,'depth [m]');
    %     title(hax_geo,fig_name);
    %======================================================================
    
    %======================================================================
    % PROCESS INPUT MOTION
    %======================================================================
    mod_inacc_th = [0 inacc_th(2:end)];                                     % processed input acceleration
    N_t          = numel(mod_inacc_th);
    time         = (0:N_t-1)*d_time;
    %     coeff        = log10(5)/(max(time));
    %     mod_inacc_th = mod_inacc_th.*exp(-coeff*time);                        % applying CyberQuake filter
    N_f          = 2^nextpow2(numel(mod_inacc_th))+1;                       % number of points for FFT computation
    d_freq       = 1/(N_f-1)/d_time;                                        % frequency increment [Hz]
    inacc_fft = fft(mod_inacc_th,N_f);                                      % modified input FFT
    % freq = (0:(N_f-1))*d_freq;
    
    %======================================================================
    % SET UP ITERATIVE PROCEDURE
    %======================================================================
    MAX_ITER = 8*(~flag_freq)+100*(flag_freq);
    iter     = 0;
    sumerr   = 1;
    
    %     fig_err=figure('name',sprintf('error_map_%s',fig_name),'visible','on');
    %     hax_err=axes('parent',fig_err,'xtick',1:N_layers,'position',[.15 .125 .8 .775]);
    %     hold(hax_err,'all');    box(hax_err,'on');
    %     title(hax_err,fig_name);
    %     xlabel(hax_err,'N^o layers'); ylabel(hax_err,'error [%]');
    %
    %     fig_geff=figure('name',sprintf('geff_map_%s',fig_name),'visible','on');
    %     hax_geff=axes('parent',fig_geff,'xtick',1:N_layers,'position',[.15 .125 .8 .775]);
    %     hold(hax_geff,'all');    box(hax_geff,'on');
    %     title(hax_geff,fig_name);
    %     xlabel(hax_geff,'N^o layers'); ylabel(hax_geff,'\gamma_E_F_F [%]');
    %======================================================================
    
    %======================================================================
    % SOLVE TILL CONVERGENCE
    %======================================================================
    status = 1;
    while (sumerr~=0) && (iter<MAX_ITER) && status
        [sumerr,output_acc,output_fft,gamma_th,Vs_layers,xsi_layers,gamma_eff,err,status] = ...
            eql_ex0_acc_strain(...
            H_layers,Vs0_layers,Vs_layers,xsi_layers,rho_layers,...
            Vs0_rock,xsi0_rock,rho_rock,...
            gamma_G,GGmax,gamma_D,Damp,...
            mod_inacc_th,inacc_fft,d_freq,...
            alpha_eff,flag_freq);
        % %
        %         plot(hax_err,1:N_layers,err,'displayname',sprintf('it-%u',iter));
        %         plot(hax_geff,1:N_layers,gamma_eff(:,1),'displayname',sprintf('it-%u',iter));
        %
        iter = iter + 1;
        
    end
    
    %     plot(hax_geo,[Vs_layers;Vs_layers(end)],z_layers,'r--','displayname','V_S');
    %     leg = legend(hax_geo,'show'); set(leg,'box','off');
    %     leg = legend(hax_err,'show'); set(leg,'box','off');
    %     leg = legend(hax_geff,'show');set(leg,'box','off');
    %     xlim(hax_err,[0;N_layers+1]); xlim(hax_geff,[0;N_layers+1]);
    %     ylim(hax_geo,[z_layers(end);z_layers(1)]); xlim(hax_geo,[100,600]);
    %     format_figures(hax_err); format_figures(hax_geff); format_figures(hax_geo);
    %     rule_fig(fig_err); rule_fig(fig_geff); rule_fig(fig_geo);
    %     saveas(fig_err,get(fig_err,'name')); close(fig_err);
    %     saveas(fig_geff,get(fig_geff,'name')); close(fig_geff);
    %     saveas(fig_geo,get(fig_geo,'name')); close(fig_geo);
    %
    %======================================================================
    % OUTPUT
    %======================================================================
    varargout{1} = real(output_acc);
    varargout{2} = output_fft;
    varargout{3} = gamma_th;
    varargout{4} = Vs_layers(1:2:end,:);
    varargout{5} = xsi_layers(1:2:end,:);
    varargout{6} = status;
    %======================================================================
    
    return
end
