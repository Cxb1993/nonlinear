function [varargout] = eql_ex0_check_conv(varargin)
    
    gamma_eff  = varargin{1};
    Vs0_strata = varargin{2};
    Vs_strata  = varargin{3};
    gamma_G    = varargin{4};
    GGmax      = varargin{5};
    gamma_D    = varargin{6};
    D          = varargin{7};
    flag_freq  = logical(varargin{8}~=0);
    
    N_strata   = size(Vs_strata,1)/2;
    N_f        = size(gamma_eff,2);
    
    Vs_iter  = -1*ones(2*N_strata,1);
    xsi_iter = -1*ones(2*N_strata,N_f);
    err      = -1*ones(2*N_strata,1);
    sumerr   = 0;
    
    % iteration on Vs and convergence check
    status = 1;
    for j = 1 : N_strata
        try
            Vs_iter(2*j-1:2*j) = Vs0_strata(2*j-1:2*j).*...
            G_gamma_D_extrap(gamma_G(j,:),GGmax(j,:),gamma_eff(j,1));
        catch
            err = -1;
            Vs_iter=-1;
            xsi_iter=-1;
            status = 0;
            break;
        end
        for ii = 1 : N_f
            xsi_iter(2*j-1:2*j,ii) = ...
                repmat(G_gamma_D_extrap(gamma_D(j,:),D(j,:),gamma_eff(j,ii)),2,1);
        end
        
        err(j) = abs((Vs_iter(2*j)-Vs_strata(2*j))./Vs_strata(2*j));
        
%         if err(j) > 0.05*flag_freq+0.01*(~flag_freq)
%             sumerr = sumerr + 1;
%         end
        
    end
    sumerr = (sum(err)>0.05*flag_freq+0.01*(~flag_freq));
    
    varargout{1} = sumerr;
    varargout{2} = Vs_iter;
    varargout{3} = xsi_iter;
    varargout{4} = err;
    varargout{5} = status;
    return
end
