function Y=multi_layers(X,Nr,Nd)
    
    dimX=size(X,1);
    dimY=size(X,2);
    
    dim=(max(dimX,dimY)==dimX)+2*(max(dimX,dimY)==dimY);
    
    X=X(:)./Nd;
    
    Y=[];
    
    for i=1:numel(X);
        
        Y=[Y;repmat(X(i),Nr,1)];
        
    end
    
    if dimX && dimY
        Y = Y(:);
    else
        Y=check_transp(Y,dim);
    end
    
    return
end