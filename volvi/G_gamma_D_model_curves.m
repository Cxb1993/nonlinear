% ================
% Editor: Filippo Gatti
% Ecole Centrale Paris - Laboratoire MSSMat
% Copyright 2015
% NOTES
% G_gamma_D_model_curves: function to interpolate and plot several
% G/Gmax-gamma_int-D curves from literature
% N.B.: need for check_transp
% INPUT: curve_name (label name of the curve to be computed)
%        gamma_int (interpolating gamam values NOT IN PERCENTAGE!)
%        p_eff (effective stress (optional))
%        PI (Plastic Index (optional))
%        OCR_soil (Over Consolidation Ratio (optional))
%        f_exc (excitation frequency (optional))
%        N_cyc (number of cycles (optional))
% OUTPUT: N_CURVES (number of curves)
%         G (vector/matrix containing G/Gmax values interpolated along gamma_int)
%         D (vector/matrix containing D values interpolated along gamma_int)
%         ref_label (string containing the bibliogrpahy reference for the
%         evaluated curves)
%         std_G (standard deviation on G/Gmax curve (optional))
%         std_D (standard deviation on D curve (optional))
% ================

function [varargout]=G_gamma_D_model_curves(varargin)
    
    pa=101.325;                                                                 % atmospheric pressure [kPa]
    
    curve_name=upper(char(varargin{1}));                                        % curve name (uppercase)
    
    %fprintf('CURVE TYPE: %s\n',curve_name);
    
    gamma_int=varargin{2};                                                      % interpolating gamma_int values
    
    %fprintf('STRAIN RANGE %.6f-%.6f\n',gamma_int(1),gamma_int(end));
    
    if nargin>=3
        
        p_eff=varargin{3};                                                      % effective stress (kPa)
        
    end
    
    
    if nargin>=4
        
        if strcmpi(curve_name,'DARE')
            
            PI=varargin{4};                                                      % plasticity index
            
        elseif strcmpi(curve_name,'MENQ')
            
            CU=varargin{4};                                                      % uniformity coefficient
            
        end
        
    end
    
    if nargin>=5
        
        OCR_soil=varargin{5};                                                        % Over-consolidation Ratio
        
    end
    
    if nargin>=6
        
        f_exc=varargin{6};                                                      % excitation frequency [Hz]
        
    end
    
    if nargin>=7
        
        N_cyc=varargin{7};                                                      % cycles of loading
        
    end
    
    if ~exist('p_eff','var')
        
        p_eff=pa;
        
        %fprintf('EFFECTIVE STRESS NOT DEFINED - DEFAULT VALUE: %.2f kPa\n',p_eff)
        
    else
        
        %fprintf('EFFECTIVE STRESS: %.2f kPa\n',p_eff)
        
    end
    
    if strcmpi(curve_name,'DARE')
        
        if ~exist('PI','var')
            
            PI=0;
            
            %fprintf('PLASTICITY INDEX NOT DEFINED - DEFAULT VALUE: %.0f\n',PI)
            
        else
            
            %fprintf('PLASTICITY INDEX: %.0f\n',PI)
            
        end
        
        if ~exist('OCR_soil','var')
            
            OCR_soil=1;
            
            %fprintf('OCR NOT DEFINED - DEFAULT VALUE: %.2f\n',OCR_soil)
            
        else
            
            %fprintf('OCR: %.2f\n',OCR_soil)
            
        end
        
    end
    
    if  strcmpi(curve_name,'MENQ') || strcmpi(curve_name,'YEE')
        
        if ~exist('CU','var')
            
            CU=10;
            
            %fprintf('UNIFORMITY COEFFICIENT NOT DEFINED - DEFAULT VALUE: %.2f\n',CU)
            
        else
            
            %fprintf('UNIFORMITY COEFFICIENT: %.2f\n',CU)
            
        end
        
    end
    
    if ~exist('f_exc','var')
        
        f_exc=1;
        
        %fprintf('FREQUENCY NOT DEFINED - DEFAULT VALUE: %.2f\n',f_exc)
        
    else
        
        %fprintf('FREQUENCY: %.2f\n',f_exc)
        
    end
    
    if ~exist('N_cyc','var')
        
        N_cyc=10;
        
        %fprintf('Nr CYCLES NOT DEFINED - DEFAULT VALUE: %.0f\n',N_cyc)
        
    else
        
        %fprintf('Nr CYCLES: %.0f\n',N_cyc)
        
    end
    
    %% G/Gmax-GAMMA-D CURVES
    
    
    switch curve_name
        
        case 'SEED'                                                             % Seed and Idriss (1970)
            
            ref_label=sprintf('Seed_&Idriss (1970)');
            
            N_CURVES=2;
            
            % Seed and Idriss digitalized curves
            
            % curve 1
            
            gamma_G1=[0.000001,5.9936E-06,8.9508E-06,0.000015211,0.000026925,0.000045474,0.000067485,0.000088737,0.00015712,0.00022322,0.0003301,0.00045685,0.00060196,0.00081985,0.0011549,0.0019129,0.0024926,0.0033454,0.0046625,0.0068603,0.0088777,0.01];
            G1=[1,0.99875,0.99469,0.98544,0.96027,0.92144,0.88636,0.85503,0.77358,0.71186,0.62556,0.54799,0.48563,0.41818,0.35502,0.27509,0.23367,0.19654,0.16009,0.12214,0.10165,0.095782];
            gamma_G2=[0.000001,1.0749E-06,1.8503E-06,3.4744E-06,0.000008165,0.00001031,0.000019562,0.000035565,0.00006218,0.00009204,0.00014488,0.00027921,0.000472,0.00071937,0.001139,0.0019052,0.0032296,0.0048364,0.0071855,0.0083496,0.0096228];
            G2=[1,0.9949,0.98996,0.97266,0.93491,0.92166,0.87332,0.81204,0.72844,0.65585,0.56515,0.43887,0.35242,0.27907,0.2129,0.15532,0.10928,0.08429,0.065808,0.060589,0.057544];
            
            % curve 2
            
            gamma_D1=[1.0194E-06,0.000001676,3.2607E-06,6.9667E-06,0.000011544,0.000031951,0.000047526,0.000081817,0.00013624,0.00021294,0.00035637,0.00049446,0.00074008,0.0011206,0.0013849,0.0019387,0.0027229,0.0036295,0.0050396,0.0061106,0.0080475,0.0092514];
            D1=[0.2397,0.28185,0.36467,0.5556,0.77381,1.4023,1.8338,2.4352,3.2144,4.1087,5.5453,6.7015,8.3674,10.241,11.435,13.312,14.949,16.62,18.145,19.179,20.417,20.973];
            gamma_D2=[1.0356E-06,1.4718E-06,2.8331E-06,4.5922E-06,0.000007997,0.000011552,0.000016818,0.000033339,0.000052773,0.00011252,0.00018936,0.000271,0.00049022,0.00082246,0.0012791,0.0022521,0.0030893,0.0043814,0.0061418,0.0076619,0.0090009,0.009415];
            D2=[0.65595,0.75265,1.1086,1.5041,2.2494,2.9708,3.9003,5.8746,7.506,10.423,12.774,14.458,17.414,20.101,22.247,24.772,25.816,26.715,27.389,27.684,27.853,27.916];
            
            % interpolating curve 1
            
            D00=D1(1);
            
            Dinf=max(D1);
            
            [G1_SI,~]=G_gamma_D_interp(gamma_G1,G1,G1,gamma_int,D00,Dinf,[1 1.1]/10);
            
            [~,D1_SI]=G_gamma_D_interp(gamma_D1,D1,D1,gamma_int,D00,Dinf,[1 1.1]/10);
            
            D00=D2(1);
            
            Dinf=max(D2);
            
            [G2_SI,~]=G_gamma_D_interp(gamma_G2,G2,G2,gamma_int,D00,Dinf,[1 1.1]/10);
            
            [~,D2_SI]=G_gamma_D_interp(gamma_D2,D2,D2,gamma_int,D00,Dinf,[1 1.1]/10);
            
            G=[G1_SI(:) G2_SI(:)];
            
            D=[D1_SI(:) D2_SI(:)];
            
            
        case 'ISZH'                                                             % Ishibashi & Zhang (1993)
            
            ref_label=sprintf('Ishibashi_&Zhang (1993)');
            
            N_CURVES=1;
            
            % coefficient n_IZ(PI)
            n_IZ=0+...                                                          % PI=0
                3.37*10^-6*PI^1.404*(PI>0 && PI<=15)+...                        % 0<PI<=15
                7.0*10^-7*PI^1.976*(PI>15 && PI<=70)+...                        % 15<PI<=70
                2.7*10^-5*PI^1.115*(PI>70);                                     % PI>70
            
            K_IZ=.5.*(ones(size(gamma_int))+...                                 % coefficient K_IZ(gamma_int,PI)
                tanh((log((.000102+n_IZ)./(gamma_int))).^.492));
            
            m_m0_IZ=.272.*(ones(size(gamma_int))-...                            % coefficient m(gamma_int,PI)-m0
                tanh((log((.000556)./(gamma_int))).^.4)).*...
                exp(-.0145*PI^1.3);
            
            G=K_IZ.*((p_eff).^m_m0_IZ);
            
            D=100*(ones(size(G))-G);
            
        case 'DARE'
            
            N_CURVES=1;
            
            gamma_int=gamma_int*100;
            
            ref_label=sprintf('Darendeli (2001)');
            
            a=.919;
            
            gamma_ref=((p_eff/pa)^.3483)*(.0352+.001*PI*OCR_soil^.3246);        % not in percent!
            
            G=1./(1+(gamma_int./gamma_ref).^a);
            
            Dmin=(p_eff^(-.2889))*(.8005+.0129*PI*OCR_soil^(-.1069))*(1+.2919*log(f_exc)); % minimum damping [%]
            
            D_MAS_1=100/pi*(4*((gamma_int-gamma_ref.*...
                log((gamma_int+gamma_ref)./gamma_ref))./...
                (gamma_int.^2./(gamma_int+gamma_ref)))-2);
            
            c1=-1.1143*a.^2+1.8618*a+.2533;
            
            c2=.0805*a.^2-0.071*a-.0095;
            
            c3=-.0005*a.^2+.0002*a+.0003;
            
            D_MAS=c1.*D_MAS_1+c2.*D_MAS_1.^2+c3.*D_MAS_1.^3;
            
            b=.6329-.0057*log(N_cyc);
            
            D=b.*(G.^.1).*D_MAS+Dmin;
            
            std_G = exp(-4.23)+(1/sqrt(exp(3.62)))*sqrt(.25-(G-.5).^2);
            
            std_D = exp(-5) + .78.*sqrt(D);
            
        case 'MENQ'
            
            N_CURVES=1;
            
            gamma_int=gamma_int*100;
            
            ref_label=sprintf('Menq (2003)');
            
            a=.86+.1*log10(p_eff/pa);
            
            gamma_ref=.12*(CU^-.6)*((p_eff/pa)^(.5*CU^(-.15)));        % not in percent!
            
            G=1./(ones(size(gamma_int))+(gamma_int./gamma_ref).^a);
            
            D50=.5;                                                             % D50 [mm]
            
            Dmin=((p_eff/pa)^(-.08))*.55*CU^(.1)*D50^(-.3);                     % minimum damping [%]
            
            D_MAS_1=100/pi*(4*((gamma_int-gamma_ref.*...
                log((gamma_int+ones(size(gamma_int))*gamma_ref)./gamma_ref))./...
                (gamma_int.^2./(gamma_int+ones(size(gamma_int))*gamma_ref)))-2*ones(size(gamma_int)));
            
            c1=-1.1143*a.^2+1.8618*a+.2533;
            
            c2=.0805*a.^2-0.071*a-.0095;
            
            c3=-.0005*a.^2+.0002*a+.0003;
            
            D_MAS=c1.*D_MAS_1+c2.*D_MAS_1.^2+c3.*D_MAS_1.^3;
            
            b=.6329-.0057*log(N_cyc);
            
            D=b.*((G).^.1).*D_MAS+Dmin;
            
        case 'YEE'
            
            N_CURVES=1;
            
            gamma_int=gamma_int*100;
            
            ref_label=sprintf('Yee et al. (2011)');
            
            a=.82+.34*log10(p_eff/pa);
            
            % gamma_ref for upper sandy layers 70m at KSH
            
            gamma_ref=.0904*((p_eff/pa)^.4345); % 0.0904 corresponds to CU=1.6 in Menq's model
            
            G=1./(ones(size(gamma_int))+(gamma_int./gamma_ref).^a);
            
            D50=.25; % D50 [mm] corresponds to CU of 2.55
            
            Dmin=((p_eff/pa)^(-.08))*.55*CU^(.1)*D50^(-.3);                     % minimum damping [%]
            
            D_MAS_1=100/pi*(4*((gamma_int-gamma_ref.*...
                log((gamma_int+ones(size(gamma_int))*gamma_ref)./gamma_ref))./...
                (gamma_int.^2./(gamma_int+ones(size(gamma_int))*gamma_ref)))-2*ones(size(gamma_int)));
            
            c1=-1.1143*a.^2+1.8618*a+.2533;
            
            c2=.0805*a.^2-0.071*a-.0095;
            
            c3=-.0005*a.^2+.0002*a+.0003;
            
            D_MAS=c1.*D_MAS_1+c2.*D_MAS_1.^2+c3.*D_MAS_1.^3;
            
            b=.6329-.0057*log(N_cyc);
            
            D=1.4*(b.*((G).^.1).*D_MAS)+Dmin; % to fit experimental data
            
    end
    
    varargout{1}=G;
    
    varargout{2}=D./100;
    
    varargout{3}=ref_label;
    
    varargout{4}=N_CURVES;
    
    if exist('std_G','var')
        
        varargout{5}=std_G;
        
    else
        
        varargout{5}=[];
        
    end
    
    if exist('std_D','var')
        
        varargout{6}=std_D/100;
        
    else
        
        varargout{6}=[];
        
    end
    
    return
end