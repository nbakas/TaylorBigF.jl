using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf, LaTeXStrings
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 3.3 Extrapolation

p=500;setprecision(p); setrounding(BigFloat, RoundUp);
L=BigFloat("1"); n=BigInt(201); iN=convert(Int64,n); dx=BigFloat("2")*L/(n-1); x=-L:dx:L;
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end
V=TaylorBigF.calc_Vandermonde(x,n,n)
xi=x[end]+dx:dx:100x[end]
Vi=TaylorBigF.calc_Vandermonde(xi,length(xi),n)

f=sin.(x); a=V\f;
fi_comp=Vi*a
fi_analyt=sin.(xi)
error1=convert(Array{Float64},fi_comp-fi_analyt)
iok=abs.(error1).<1


# plot(x,f,color=:black,width=5,xlabel= "x" ,ylabel="f(x)",
#     label="given f")
# plot!(xi[iok],fi_analyt[iok],color=:black,width=2,label="exact f")
# plot!(xi[iok],fi_comp[iok],color=:black,width=2,
#         label=string("extrapolated f, p=",p),linestyle=:dash)
# plot!(xi[iok],fi_comp[iok],color=:black,width=2,
#         label=string("extrapolated f, p=",p),linestyle=:dashdot)
# plot!(xi[iok],fi_comp[iok],color=:black,width=2,
#         label=string("extrapolated f, p=",p),linestyle=:dashdotdot)

# savefig(string("extrap.pdf"))


r_root=Array{Float64}(undef,0); for i=1:Int64(n) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
1/r_root[end]
xi[iok][end]
