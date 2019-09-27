

p=50;setprecision(p); setrounding(BigFloat, RoundUp); 
L=BigFloat("10"); n=BigInt(100); iN=convert(Int64,n); 
dx=BigFloat("2")*L/(n-1); x=convert(Array{Float64},-L:dx:L)
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

V=TaylorBigF.calc_Vandermonde(x,n,n)
xi=(x.+dx/BigFloat("2"))[1:end-1]
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,n)
dV=TaylorBigF.calc_Vandermonde_dot(x,n,n)
ddV=TaylorBigF.calc_Vandermonde_dot_dot(x,n,n)
dddV=TaylorBigF.calc_Vandermonde_dot_dot_dot(x,n,n)
ddddV=TaylorBigF.calc_Vandermonde_dot_dot_dot_dot(x,n,n)
intV=TaylorBigF.calc_Vandermonde_int(x,n,n)

f=sin.(x)+sin.(x./BigFloat("2"))
fi=sin.(xi)+sin.(xi./BigFloat("2"))
f+=convert(Array{BigFloat},(rand(n).-1/2)./1)
int_f=zeros(BigFloat,n)
for i=2:n 
    int_f[i]=int_f[i-1]+(f[i]+f[i-1])*dx/2 
end


d1=f[2:end]-f[1:end-1]
d2=(d1[2:end]+d1[1:end-1])/2
d=[d1[1];d2;d1[end]]

scatter(x,f)
plot!(x,f)


Eq1=V[1,:]
q1=f[1]
Eq2=V[end,:]
q2=f[end]

VV=[intV;Eq1';Eq2']
q_all=[f;q1;q2]
a=VV\q_all


r_root=Array{Float64}(undef,0); for i=1:Int64(n-1) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
plot(r_rootF)
error1=maximum(abs.(f-V*a))
plot(x,f)
scatter!(x,V*a)
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,BigInt(150));
fi_comp=Vi*a
error2=maximum(abs.(fi_comp-fi))
scatter!(xi,Vi*a)

