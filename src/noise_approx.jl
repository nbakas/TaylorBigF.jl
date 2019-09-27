
using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf, JLD2, FileIO

# noise1=convert(Array{BigFloat},(rand(n).-1/2)./BigFloat("1e1000"))
# f=sin.(25x)
# f=f+noise1


p=50;setprecision(p); setrounding(BigFloat, RoundUp); 
L=BigFloat("10"); n=BigInt(100); iN=convert(Int64,n); 
dx=BigFloat("2")*L/(n-1); x=convert(Array{Float64},-L:dx:L)
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

V=TaylorBigF.calc_Vandermonde(x,n,n)
xi=(x.+dx/BigFloat("2"))[1:end-1]
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,n)


# f=convert(Array{BigFloat},(rand(n).-1/2))
# tansf1=BigFloat("1e100")
# f_transform=exp.(-f.^2)
f=sin.(x)+sin.(x./BigFloat("2"))

# fi=(f[2:end]+f[1:end-1])/2
fi=sin.(xi)+sin.(xi./BigFloat("2"))
# fi_transform=exp.(-fi.^2)



# f[75]+=BigFloat("1e-20")
f+=convert(Array{BigFloat},(rand(n).-1/2)./1)
# plot(x,f)

a=V\f
# a.*=(1 .+(rand(n).-1/2)/1e3)
r_root=Array{Float64}(undef,0); for i=1:Int64(n-1) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
# plot(r_rootF)
# plot(sign.(a).*log10.(abs.(a)))
error1=maximum(abs.(f-V*a))
fi_comp=Vi*a
error2=maximum(abs.(fi_comp-fi))

plot(x,f)
scatter!(x,V*a,markersize=1)
scatter!(xi,fi,markersize=1)
scatter!(xi,fi_comp,markersize=1)
# fi_comp=sqrt.(log.(-fi_comp_transform))
# fi_comp=fi_comp_transform*tansf1
# fi_analyt=sin.(25xi);



# dV=TaylorBigF.calc_Vandermonde_dot(x,n,n)

# Eq1=V[1,:]
# q1=f[1]
# Eq2=dV[1,:]
# q2=(f[2]-f[1])/dx
# Eq3=V[end,:]
# q3=f[end]
# Eq4=dV[end,:]
# q4=(f[end]-f[end-1])/dx

# VV=[V;Eq1';Eq2';Eq3';Eq4']
# q_all=[f;q1;q2;q3;q4]
# a=VV[1:end,:]\q_all[1:end]
# r_root=Array{Float64}(undef,0); for i=1:Int64(n-1) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
# plot(r_rootF)
# error1=maximum(abs.(f-V*a))
# plot(x,f)
# scatter!(x,V*a)
# Vi=TaylorBigF.calc_Vandermonde(xi,n-1,n);
# fi_comp=Vi*a
# error2=maximum(abs.(fi_comp-fi))
# scatter!(xi,Vi*a)




