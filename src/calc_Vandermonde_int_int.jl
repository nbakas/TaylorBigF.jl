function calc_Vandermonde_int_int(x,ni,nj)

    V=zeros(BigFloat, ni, nj)
    for i=1:ni
        for j=1:nj
            V[i,j] = (x[i]^BigInt(j+1))/BigInt(j+1)/BigInt(j)
        end
    end

    return V

end