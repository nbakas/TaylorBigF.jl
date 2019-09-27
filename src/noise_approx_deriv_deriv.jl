

p=50;setprecision(p); setrounding(BigFloat, RoundUp); 
L=BigFloat("10"); n=BigInt(100); iN=convert(Int64,n); 
dx=BigFloat("2")*L/(n-1); x=convert(Array{Float64},-L:dx:L)
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

f=sin.(x)+sin.(x./BigFloat("2"))
fi=sin.(xi)+sin.(xi./BigFloat("2"))
f+=convert(Array{BigFloat},(rand(n).-1/2)./1)

ff=3.2
typeof(ff)
ff16=convert(Float16,ff)
ff32=convert(Float32,ff)
ff16-ff32

scatter(x,f)
plot!(x,f)

i_dd=Array{Int64}(undef,0)
for i=2:length(d)-1
    if (f[i+1]-f[i])*(f[i-1]-f[i])>0
        push!(i_dd,i) 
    end 
end
scatter!(x[i_dd],f[i_dd])



V=TaylorBigF.calc_Vandermonde(x,n,n+BigInt(length(i_dd)))
xi=(x.+dx/BigFloat("2"))[1:end-1]
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,n+BigInt(length(i_dd)))
ddV=TaylorBigF.calc_Vandermonde_dot_dot(x,n,n+BigInt(length(i_dd)))

Eq1=ddV[i_dd,:]
q1=zeros(BigFloat,length(i_dd))

VV=[V;Eq1]
q_all=[f;q1]
a=VV\q_all
plot(a)

r_root=Array{Float64}(undef,0); for i=1:Int64(n-1) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
plot(r_rootF)
error1=maximum(abs.(f-V*a))
plot(x,f)
pred=V*a
scatter!(x[20:80],pred[20:80])
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,BigInt(150));
fi_comp=Vi*a
error2=maximum(abs.(fi_comp-fi))
scatter!(xi,Vi*a)

