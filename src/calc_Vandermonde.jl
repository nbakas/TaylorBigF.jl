function calc_Vandermonde(x,ni,nj)

    V=zeros(BigFloat, ni, nj)
    for i=1:ni
        for j=1:nj
            V[i,j] = (x[i])^(BigInt(j)-BigInt(1))
        end
    end

    return V

end