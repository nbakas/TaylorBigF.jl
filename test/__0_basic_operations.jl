using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 3.1 Basic Operations

p=1000;setprecision(p); setrounding(BigFloat, RoundUp);
L=BigFloat("1"); n=BigInt(201); iN=convert(Int64,n); dx=BigFloat("2")*L/(n-1); x=-L:dx:L;
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

V=TaylorBigF.calc_Vandermonde(x,n,n)
inv_V_comp = V\I; det_V_comp=det(V);
inv_V_analyt,det_V_analyt=TaylorBigF.calc_Vandermonde_inv_formula(x,n);
diff_V=maximum(abs.(inv_V_comp.-inv_V_analyt));@sprintf("diff_V=%1.3e",diff_V)
diff_det=det_V_analyt-det_V_comp;@sprintf("det_V_comp=%1.3e det_V_analyt=%1.3e diff_det=%1.3e",det_V_comp,det_V_analyt,diff_det)


f=sin.(x); a=V\f;
# plot(a)
r_root=Array{Float64}(undef,0); for i=1:Int64(n) push!(r_root,abs(a[i])^(1/(n-1))) end; r_rootF=convert(Array{Float64},r_root);
# ticks = round.(range(minimum(r_rootF),stop=maximum(r_rootF),length = 4),digits=2);
# plot(1:Int64(n),r_rootF,color=:black,legend=false,xlabel="n",ylabel="1/r",yticks=ticks)
# scatter!(1:Int64(n),r_rootF,color=:black)
# savefig(string("r_sin",p,".pdf"))
# r_ratio=Array{Float64}(undef,0); for i=1:Int64(n)-1 push!(r_ratio,abs(a[i+1]/a[i])) end; plot(r_ratio[2:2:end])

t = Taylor1(Float64, Int64(n)-1); ts=sin(t)
a_analyt=Array{Float64}(undef,0); for i=1:Int64(n) push!(a_analyt,getcoeff(ts,i-1)) end
diff_a=maximum(abs.(a_analyt-a));@sprintf("diff_a=%1.3e",diff_a)

f_comp=V*a
diff_f=maximum(abs.(f-f_comp));@sprintf("diff_f=%1.3e",diff_f)
@sprintf("Rn=%1.3e",1/factorial(n))
