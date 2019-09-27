function calc_Vandermonde_xxxx_2D(x,y,ni,nj)

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
            if pows[j][1]-4>=0
                # a1=BigInt(pows[j][1]-1)*BigInt(pows[j][1]-2)*BigInt(pows[j][1]-3)*BigInt(pows[j][1]-4)
                a1=BigInt(pows[j][1]-0)*BigInt(pows[j][1]-1)*BigInt(pows[j][1]-2)*BigInt(pows[j][1]-3)
                V[i,j] = a1*(x[i]^BigInt(pows[j][1]-4))*(y[i]^BigInt(pows[j][2]))
            end
        end
    end

    return V

end