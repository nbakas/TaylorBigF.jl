
using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf, FileIO
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 3.7 System Identification


p=2000;setprecision(p); setrounding(BigFloat, RoundUp);
L=BigFloat("1"); n=BigInt(100); iN=convert(Int64,n);
dx=BigFloat("1")*L/(n-1); x=convert(Array{BigFloat},0:dx:(n-1)*dx)
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

V=TaylorBigF.calc_Vandermonde(x,n,n)
dV=TaylorBigF.calc_Vandermonde_dot(x,n,n)
ddV=TaylorBigF.calc_Vandermonde_dot_dot(x,n,n)
dddV=TaylorBigF.calc_Vandermonde_dot_dot_dot(x,n,n)
ddddV=TaylorBigF.calc_Vandermonde_dot_dot_dot_dot(x,n,n)


# x=time (t in manuscript) s=space dds=a=F/m
s=x.^BigInt(2)
# plot(x,s)


a=V\s
# plot(a)

# s_dot=s[2:end]-s[1:end-1]
# s_dot_dot=s_dot[2:end]-s_dot[1:end-1]
# VV=[ones(n-2) s[2:end-1] s_dot[1:end-1] s_dot_dot]
# aa=VV\ones(n-2)
# plot(aa)

# ss=V*a;ssd=dV*a;ssdd=ddV*a;
# VV=[ss ssd ssdd ss.^2 ss.*ssd ss.*ssdd ssd.^2 ssd.*ssdd ssdd.^2]
VV=[V*a dV*a ddV*a]
b=VV\ones(n)
# plot(b)


# xi=x[end]+dx:dx:x[end]+100dx
xi=x[end]+BigFloat("9525000000000")dx:dx:x[end]+BigFloat("9525000000010")*dx
Vi=TaylorBigF.calc_Vandermonde(xi,length(xi),n)
IVi=TaylorBigF.calc_Vandermonde_int(xi,length(xi),n)
IIVi=TaylorBigF.calc_Vandermonde_int_int(xi,length(xi),n)
IV=TaylorBigF.calc_Vandermonde_int(x,length(x),n)
IIV=TaylorBigF.calc_Vandermonde_int_int(x,length(x),n)


# t2=(x.^BigInt(2))/BigInt(2)
# mat1=(-t2+b[1]*IIV*a+b[2]*IV*a+b[3]*V*a)
# cc=[s[1] BigInt(1);s[end] BigInt(1)]\[mat1[1,:];mat1[end,:]]
# s_comp=(t2+cc[1].*x+cc[2].*ones(BigFloat,length(x))-b[1]*IIV*a-b[2]*IV*a)/b[3]
# # plot(s_comp)#confirmation
# t2=(xi.^BigInt(2))/BigInt(2)
# s_comp=(t2+cc[1].*xi+cc[2].*ones(BigFloat,length(xi))-b[1]*IIVi*a-b[2]*IVi*a)/b[3] # Equation 5
# plot(x,s,legend=:topleft,xaxis="t",yaxis="s",label="given s")
# plot!(xi,s_comp,label=string("identified s, p=",p))
# # savefig("syst_id.pdf")



t2=(xi.^BigInt(2))/BigInt(2)
s_comp=(t2+-b[1]*IIVi*a-b[2]*IVi*a)/b[3] # Equation 5
# plot(x,s,legend=:topleft,xaxis="t",yaxis="s",label="given s")
# plot!(xi,s_comp,label=string("identified s, p=",p))
# savefig("syst_id.pdf")



s_analyt=xi.^BigInt(2)
error1=s_comp-s_analyt
iok=abs.(error1).<1
xi[iok][end]
