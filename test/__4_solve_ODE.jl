
using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 3.6 Solution of Ordinary Differential Equations


p=50;setprecision(p); setrounding(BigFloat, RoundUp);
L=BigFloat("1"); n=BigInt(100); iN=convert(Int64,n);
dx=BigFloat("1")*L/(n-1); x=0:dx:(n-1)*dx;
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

V=TaylorBigF.calc_Vandermonde(x,n,n+BigInt(4))
# Vi=calc_Vandermonde(xi,n-1,n+BigInt(4))
dV=TaylorBigF.calc_Vandermonde_dot(x,n,n+BigInt(4))
ddV=TaylorBigF.calc_Vandermonde_dot_dot(x,n,n+BigInt(4))
dddV=TaylorBigF.calc_Vandermonde_dot_dot_dot(x,n,n+BigInt(4))
ddddV=TaylorBigF.calc_Vandermonde_dot_dot_dot_dot(x,n,n+BigInt(4))

xi=(x.+dx/BigFloat("2"))[1:end-1];
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,n+BigInt(4));

q=zeros(BigFloat,n)
# q.-=1e2sin.(exp.(3x))
# q.+=convert(Array{BigFloat},(rand(n).-1/2).*10)
# plot(q)
Eq1=V[1,:]
q1=BigFloat("0")
Eq2=dV[1,:]
q2=BigFloat("0")
Eq3=V[end,:]
q3=BigFloat("1")
Eq4=dV[end,:]
q4=BigFloat("0")

VV=[ddddV;Eq1';Eq2';Eq3';Eq4']
q_all=[q;q1;q2;q3;q4]
a=VV[5:end,:]\q_all[5:end]
# a=convert(Array{Float64},VV[5:end,:])\convert(Array{Float64},q_all[5:end])
# plot(a,xlabel="n",ylabel="a",color=:black,width=2,
#         label=string("p=",p),linestyle=:dash)
# plot!(a,xlabel="n",ylabel="a",color=:black,width=2,
#         label=string("p=",p),linestyle=:dashdotdot)
# savefig(string("ode_beam.pdf"))




# r_root=Array{Float64}(undef,0); for i=1:Int64(n-1) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
# plot(r_rootF)
sol1=V*a
scatter(convert(Array{Float64},x),sol1)


exact1=-2x.^3+3x.^2
plot!(convert(Array{Float64},x),exact1)

error1=exact1-sol1
# plot(error1)
@sprintf("error max=%1.3e",maximum(abs.(error1)))


exacti1=-2xi.^3+3xi.^2
soli1=Vi*a
errori1=exacti1-soli1
@sprintf("error i max=%1.3e",maximum(abs.(errori1)))
plot!(convert(Array{Float64},xi),soli1)



M=ddV*a
plot(convert(Array{Float64},x),M)

precision(BigFloat)
m1=string(M[1])
