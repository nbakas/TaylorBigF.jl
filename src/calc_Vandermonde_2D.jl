function calc_Vandermonde_2D(x,y,ni,nj)

    N=1
    NN=2
    pows=reverse.(Iterators.product(fill(0:N-1,NN)...))[:]
    while length(pows)<nj
        N+=1
        pows=reverse.(Iterators.product(fill(0:N-1,NN)...))[:]
    end
    
    V=zeros(BigFloat, ni, nj)
    for i=1:ni
        for j=1:nj
            V[i,j] = (x[i]^BigInt(pows[j][1]))*(y[i]^BigInt(pows[j][2]))
        end
    end

    return V

end