
using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 3.4 Numerical Integration
p=50;setprecision(p); setrounding(BigFloat, RoundUp);
L=BigFloat("1"); n=BigInt(201); iN=convert(Int64,n); dx=BigFloat("2")*L/(n-1); x=-L:dx:L;
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end
V=TaylorBigF.calc_Vandermonde(x,n,n)

f=sin.(x); a=V\f;

L_n=zeros(BigFloat,n)
for i=1:n
    L_n[i]=(L^i)/i
    if iseven(i) L_n[i]=-L_n[i] end
end
c=a'*L_n

Vintegral=TaylorBigF.calc_Vandermonde_int(x,n,n)
F_all_x=Vintegral*a.+c

F_comp=F_all_x[end]
@sprintf("F_comp=%1.3e",F_comp)
