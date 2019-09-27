function calc_Vandermonde_inv_formula(x,n)

    detX=BigInt(1)
    for j=1:n
        for i=1:j-1
            detX*=(x[j]-x[i])
        end
    end

    L_1=zeros(BigFloat, n, n)
    for i=1:n
        for j=1:n
            if i>=j
                lij=BigInt(1)
                for k=1:i
                    if k!=j
                        lij/=(x[j]-x[k])
                    end
                end
                L_1[i,j]=lij
            end
        end
    end
    L_1[1,1]=BigInt(1)

    U_1=zeros(BigFloat, n, n)
    for i=1:n
        for j=1:n
            if i==j
                U_1[i,j]=BigInt(1)
            elseif j>i
                if i==1
                    U_1[i,j]=BigInt(0)-U_1[i,j-1]*x[j-1]
                else
                    U_1[i,j]=U_1[i-1,j-1]-U_1[i,j-1]*x[j-1]
                end
            end
        end
    end

    U_1L_1=U_1*L_1

    return U_1L_1,detX

end
