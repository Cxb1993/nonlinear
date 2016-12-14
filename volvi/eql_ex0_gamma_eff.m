function [varargout] = eql_gamma_eff(varargin)
    global omega omega0 fft_gamma gamma0 ab0
    %======================================================================
    % INPUT PARAMETERS
    %======================================================================
    flag_freq = varargin{1};
    alpha_eff = varargin{4};
    
    if flag_freq~=0
        omega     = varargin{2}(:);
        fft_gamma = varargin{3}(:);
        
        switch(flag_freq)
            case 1
                %==========================================================
                % GAMMA EFF - FREQUENCY (ASSIMAKI-KAUSEL 2002)
                %==========================================================
                d_omega = mean(diff(omega));
                idx0 = round(2*pi*.1/d_omega)+1;
                
                Tm1 = fft_gamma(idx0:end).^2./omega(idx0:end);
                Tm2 = fft_gamma(idx0:end).^2;
                omega0 = d_omega*round(sum(Tm2)/sum(Tm1)/d_omega);
                gamma0 = fft_gamma(omega==omega0);
                
                X0  = [1,1];
                ab0 = [.23;2.3];
                
                options = optimset('Display','Iter','TolFun',1e-15,'MaxFunEvals',200,...
                    'MaxIter',200,'TolCon',1e-10,'DiffMinChange',1e-4);
                A = []; B = []; Aeq = []; Beq =[]; Ub = [1.5,1.5]; Lb = [.5 .5]; c =[]; ceq = [];
                [X_opt,~,~,~,~,~,~]=...
                    fmincon(@(X) min2_gamma_spectra(X),X0,A,B,Aeq,Beq,Lb,Ub,[c,ceq],options);
                
                gamma_eff = smooth_gamma_spectrum(omega,omega0,gamma0,X_opt(1)*.23,X_opt(2)*2.1);
                
            case 2
                %==========================================================
                % GAMMA EFF - FREQUENCY (YOSHIDA 2002)
                %==========================================================
                N_f = numel(omega);
                [gamma_max,idx_0] = max(fft_gamma); 
                fp = omega(idx_0)/2/pi; fp=fp(1);
                fe = 6; idx_1 = find(omega/2/pi > fe,1,'first');
                
                if idx_1 > idx_0
                    % f < fp
                    gamma_eff = ones(idx_0-1,1);
                    % fp <= f <= fe
                    gamma_eff=[gamma_eff; ...
                        (1-(log10(omega(idx_0:idx_1)./fp/2/pi)./log10(fe/fp)).^2)];
                    % f > fe
                    gamma_eff=[gamma_eff; ...
                        zeros(N_f-idx_1,1)];
                    
                    gamma_eff=gamma_max.*gamma_eff;
                else
                    gamma_eff=gamma_max*ones(N_f,1);
                end
                
        end
    else
        gamma_lim = varargin{2};
        gamma_th  = varargin{3}(:);
        
        gamma_eff = alpha_eff * max(abs(gamma_th));
        if gamma_eff > gamma_lim
            gamma_eff = 1e-4;
        end
        
    end
    
    varargout{1} = gamma_eff(:).';
    
    return
end