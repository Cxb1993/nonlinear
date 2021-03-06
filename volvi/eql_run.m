function [varargout] = eqlrun(varargin)
    %% *SET-UP*
    inp = varargin{1};
    mat = varargin{2};
    
    %% *EQL ANALYSIS*
    [eql.tha,~,eql.gam,eql.tau,~,~,~] = eql_multi_reflection(...
        mat.H_layers,...
        mat.Vs0_layers,...
        mat.rho_layers,...
        mat.xsi0_layers,...
        mat.Vs0_rock,...
        mat.rho_rock,...
        mat.xsi0_rock,...
        mat.gamma,...
        mat.G_Gmax,...
        mat.gamma,...
        mat.D,...
        inp.tha,inp.dtm,0.65,0,...
        sprintf('argostoli-s-motion'));
    eql.vtm = inp.vtm;
    
    for i_=1:2*mat.N_layers+1
        eql.thv(i_,:) = cumtrapz(eql.tha(i_,:))*inp.dtm;
    end
    
    %% *OUTPUT*
    varargout{1} = eql;
    return
end