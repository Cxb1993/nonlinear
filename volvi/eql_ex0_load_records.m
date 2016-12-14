function [varargout] = eql_ex0_load_records(varargin)
    %% *SET-UP*
    inp.fnm = varargin{1};
    inp.dtm = 0.01;
    inp.vtd = 10.;
    inp.vtm = (0:inp.dtm:inp.vtd)';
    inp.ntm = numel(inp.vtm);
    
    switch inp.fnm
        case 'gabor'
            inp.thv = importdata('gabor.txt');
            inp.vtm = inp.thv(:,1);
            inp.dtm = mean(diff(inp.vtm));
            inp.ntm = numel(inp.vtm);
            inp.thv = 2*inp.thv(:,2)*0.008;
            inp.tha = [0;diff(inp.thv)./inp.dtm];
        case 'gaussian'
            inp.tau   = 1;
            inp.t0    = 2;
            inp.amp   = 5;
            inp.tha=-inp.amp.*2.*(inp.vtm-inp.t0).*exp(-((inp.vtm-inp.t0)./inp.tau).^2);
            inp.tha(inp.vtm-inp.t0>=8) = 0.;
            
    end
    
    inp.psa = max(abs(inp.tha));
    
    %% *OUTPUT*
    varargout{1} = inp;
    
    return
end