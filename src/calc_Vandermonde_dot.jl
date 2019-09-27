function calc_Vandermonde_dot(x,ni,nj)

    V=zeros(BigFloat, ni, nj)
    for i=1:ni
        for j=2:nj
            V[i,j] = (BigInt(j)-BigInt(1))*x[i]^(BigInt(j)-BigInt(2))
        end
    end

    return V

end