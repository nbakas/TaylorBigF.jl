


using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf, JLD2, FileIO
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 4.1 Multidimensional Interpolation

p=2000;setprecision(p); setrounding(BigFloat, RoundUp);
L=BigFloat("1"); n=BigInt(300); dx=BigFloat("2")*L/(n-1);
x=convert(Array{BigFloat},(rand(n).-1/2)*1)
y=convert(Array{BigFloat},(rand(n).-1/2)*1)
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end
V=TaylorBigF.calc_Vandermonde_2D(x,y,n,n)
@time Vinv=V\I

# @save string("ws",p,".jld2") p,L,n,dx,x,y,I,V,Vinv
# p=1000;setprecision(p); setrounding(BigFloat, RoundUp);
# @load string("ws",p,".jld2") p,L,n,dx,x,y,I,V,Vinv

# sos pigaine poli kalitera me tria matrixes gia x y xy me psiles dinameis, alla mono se eidikes sinartisis ph sinx + cosy
# V=calc_Vandermonde_2D_highO(x,y,n,BigInt(34))

function func(x,y) sin.(BigFloat("5")x).+cos.(exp.(BigFloat("2")y)) end

f=func(x,y)
# V=V[:,end:-1:1]
# a=V\f
a=Vinv*f

# plot(a)
# scatter3d(x,y,f); scatter3d!(x,y,V*a)
# maximum(abs.(f-V*a))
# xi=x.+BigFloat("0.01")
# yi=y.+BigFloat("0.01")
# Vi=calc_Vandermonde_2D_highO(xi,yi,n,BigInt(34))
xi=convert(Array{BigFloat},(rand(n).-1/2)*0.7)
yi=convert(Array{BigFloat},(rand(n).-1/2)*0.7)
Vi=TaylorBigF.calc_Vandermonde_2D(xi,yi,n,n)
fi_comp=Vi*a
fi_analyt=func(xi,yi)
diff_fi=maximum(abs.(fi_analyt-fi_comp));@sprintf("diff_fi=%1.3e",diff_fi)
# scatter3d(xi,yi,fi_analyt,label="exact values"
#     ,markershape=:rtriangle,color=:black,)
# scatter3d!(xi,yi,fi_comp,label=string("approximated, p=",p)
#     ,markershape=:cross,legend=:topleft,color=:black)
# scatter3d!(xi,yi,fi_comp,label=string("approximated, p=",p)
#     ,markershape=:xcross,legend=:topleft,color=:black,xaxis="x")
# zlims!(-1.4,2.0)
# maximum(fi_comp)
# minimum(fi_comp)
# savefig("interp3Dp2000.pdf")
# mean(abs.(V))-mean(abs.(Vi))
# mean(abs.(x))-mean(abs.(xi))
