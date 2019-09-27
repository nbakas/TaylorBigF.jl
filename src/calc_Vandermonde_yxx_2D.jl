function calc_Vandermonde_yxx_2D(x,y,ni,nj)

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
            if pows[j][1]-2>=0 && pows[j][2]-1>=0
                a1=BigInt(pows[j][1]-0)*BigInt(pows[j][1]-1)
                a2=BigInt(pows[j][2]-0)
                V[i,j] = a1*a2*(x[i]^BigInt(pows[j][1]-2))*(y[i]^BigInt(pows[j][2]-1))
            end
        end
    end

    return V

end