function calc_Vandermonde_dot_dot_dot_dot(x,ni,nj)

    V=zeros(BigFloat, ni, nj)
    for i=1:ni
        for j=5:nj
            V[i,j] = (BigInt(j)-BigInt(1))*(BigInt(j)-BigInt(2)*(BigInt(j)-BigInt(3))*(BigInt(j)-BigInt(4)))*x[i]^(BigInt(j)-BigInt(5))
        end
    end

    return V

end