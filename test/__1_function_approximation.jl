using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf

path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

# 3.2 Function Approximation

p=2000;setprecision(p); setrounding(BigFloat, RoundUp);
precision(BigFloat)==p
L=BigFloat("1"); n=BigInt(201); iN=convert(Int64,n); dx=BigFloat("2")*L/(n-BigFloat("1"));
x=convert(Array{BigFloat},-L:dx:L)
length(x)==n
# x=-L:dx:L
# x=convert(Array{BigFloat},(rand(n).-1/2)/10)
I=zeros(BigFloat, n, n);for i=1:n I[i,i]=BigFloat("1") end

V=TaylorBigF.calc_Vandermonde(x,n,n)


function func(x) sin(x) end

f=func.(x)
a=V\f
# plot(a)
# inv_V_comp = inv(V)
# a=inv_V_comp*f

# inv_V_analyt,det_V_analyt=TaylorBigF.calc_Vandermonde_inv_formula(x,n);
# a=inv_V_analyt*f

# inv_V_comp = V\I
# a=inv_V_comp*f

# a=convert(Array{BigFloat},convert(Array{Float64},V)\convert(Array{Float64},f))

xi=(x.+dx/BigFloat("2"))[1:end-1];
Vi=TaylorBigF.calc_Vandermonde(xi,n-1,n);

fi_comp=Vi*a;
fi_analyt=func.(xi)
diff_f=maximum(abs.(f-V*a));@sprintf("diff_f=%1.3e",diff_f)
diff_fi=maximum(abs.(fi_analyt-fi_comp));@sprintf("diff_fi=%1.3e",diff_fi)


# l1=30;fF=convert(Array{Float64},fi_comp);xF=convert(Array{Float64},x);xiF=convert(Array{Float64},xi);
# ticksY = round.(range(minimum(fF[1:l1]),stop=maximum(fF[1:l1]),length = 5),digits=2);
# ticksX = round.(range(minimum(xiF[1:l1]),stop=maximum(xiF[1:l1]),length = 5),digits=2);
# scatter(xF[1:l1],f[1:l1],color=:black,markersize=3,xlabel="x",
# ylabel="f(x)",yticks=ticksY,xticks=ticksX,label="f(x) analytical",
# legend=:bottomright,markershape=:rtriangle)
# scatter!(xiF[1:l1],fi_comp[1:l1],color=:black,markersize=3,
# label=string("f(x) numerical p=",p),markershape=:circle)# circle,rect,cross

# savefig(string("f__.pdf"))
