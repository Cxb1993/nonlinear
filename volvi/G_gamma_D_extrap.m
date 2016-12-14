%===============
% G-gamma-D curve extrapolation
% Editor: Filippo Gatti
% Ecole Centrale Paris - Laboratoire MSSMat
% Copyright 2015
% NOTES
% G_gamma_D_extrap: function to extrapolate (via linear interpolation)
% values from G/Gmax-gamma-D curves
% INPUT:  X (row vector curves gamma values)
%         Y (row vector G/Gmax or D values)
%         XE (desired value of gamma)
% OUTPUT: YE (extrapolated value at XE)
%===============

function YE=G_gamma_D_extrap(X,Y,XE)
    
    X=X(:).';
    Y=Y(:).';
    
    N_X=length(X);
    
    if X(1)>=XE      % desired value is outside the curve range (too small value)
        YE=Y(1);
    elseif X(end)<XE % desired value is outside the curve range (too high value)
        YE=X(end);
        % looking for X value next to the
    else
        for j=2:N_X
            if X(j)>XE
                ab=polyfit(log10(X(j-1:j)),Y(j-1:j),1); % coefficient of linear interpolating polynomial
                YE=polyval(ab,log10(XE));               % interpolated value
                break;
            end
        end
    end
    
    return
end