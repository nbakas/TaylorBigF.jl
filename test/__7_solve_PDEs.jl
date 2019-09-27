using LinearAlgebra, Plots, Statistics, TaylorSeries, Printf, JLD2, FileIO
# if !(pwd() in LOAD_PATH) push!(LOAD_PATH, pwd()) end
# using TaylorBigF


path1=realpath(dirname(@__FILE__)*"/..")
include(string(path1,"/src/TaylorBigF.jl"))

p=2000;setprecision(p); setrounding(BigFloat, RoundUp);
n=BigInt(20);# zigos lliws i prwti grammi tou Vall einai =0 sto x=0, (to idio k sta ODES)
nb=BigInt(2)*n+BigInt(2)*(n-BigInt(2)) # boundary points
dx=BigInt(1)/BigInt(99)
L=(n-BigInt(1))*dx
x0=zeros(BigFloat,n).-L/BigFloat("2")
x=Array{BigFloat}(undef,0)
for i=BigInt(1):n
    global x0,x,dx
    x=[x;x0.+(i-BigInt(1))*dx]
end
y=convert(Array{BigFloat},range(-L/BigFloat("2"),L/BigFloat("2"),length=n))
y=repeat(y,n)
# scatter(x,y)

V=TaylorBigF.calc_Vandermonde_2D(x,y,n*n,n*n+nb)
Vxxxx=TaylorBigF.calc_Vandermonde_xxxx_2D(x,y,n*n,n*n+nb)
Vyyyy=TaylorBigF.calc_Vandermonde_yyyy_2D(x,y,n*n,n*n+nb)
Vxxyy=TaylorBigF.calc_Vandermonde_xxyy_2D(x,y,n*n,n*n+nb)

Vall0=Vxxxx+2Vxxyy+Vyyyy
# boundaries
i_keep1=x.==-L/BigFloat("2")
i_keep2=x.==L/BigFloat("2")
i_keep3=y.==-L/BigFloat("2")
i_keep4=y.==L/BigFloat("2")
i_keep=i_keep1+i_keep2+i_keep3+i_keep4
i_b=i_keep.>0
# scatter!(x[i_b],y[i_b])
BigInt(sum(i_b))==nb
Vall1=TaylorBigF.calc_Vandermonde_2D(x[i_b],y[i_b],nb,n*n+nb)

Vall=[Vall0;Vall1]
q=[-ones(n*n);zeros(BigFloat,nb)]
@time a=Vall\q
maximum(abs.(Vall*a-q))
# scatter(x,y,zcolor=V*a)
# minimum(V*a)
# scatter3d(x,y,V*a,label=string("deflection (w), p=",p),legend=:topleft)



Vxxx=TaylorBigF.calc_Vandermonde_xxx_2D(x[i_b],y[i_b],nb,n*n+nb)
Vyyy=TaylorBigF.calc_Vandermonde_yyy_2D(x[i_b],y[i_b],nb,n*n+nb)
Vxyy=TaylorBigF.calc_Vandermonde_xyy_2D(x[i_b],y[i_b],nb,n*n+nb)
Vyxx=TaylorBigF.calc_Vandermonde_yxx_2D(x[i_b],y[i_b],nb,n*n+nb)
Qx=-(Vxxx+Vxyy)*a
Qy=-(Vyxx+Vyyy)*a
# scatter(x[i_b],y[i_b],zcolor=Qx)



i1=x.==-L/BigFloat("2")
i2=x.==L/BigFloat("2")
i12=i1+i2
i12=i12.>0
Vxxx=TaylorBigF.calc_Vandermonde_xxx_2D(x[i12],y[i12],2n,n*n+nb)
Vxyy=TaylorBigF.calc_Vandermonde_xyy_2D(x[i12],y[i12],2n,n*n+nb)
Qx=-(Vxxx+Vxyy)*a
# scatter(x[i12],y[i12],zcolor=Qx)
i3=y.==-L/BigFloat("2")
i4=y.==L/BigFloat("2")
i34=i3+i4
i34=i34.>0
Vyyy=TaylorBigF.calc_Vandermonde_yyy_2D(x[i34],y[i34],2n,n*n+nb)
Vyxx=TaylorBigF.calc_Vandermonde_yxx_2D(x[i34],y[i34],2n,n*n+nb)
Qy=-(Vyxx+Vyyy)*a
# scatter!(x[i34],y[i34],zcolor=Qy)

mean(abs.(Qx))*2*L+mean(abs.(Qy))*2*L-1*L*L






# point load
q=zeros(BigFloat,n*n+nb)
dd=zeros(BigFloat,n*n)
for i=1:n*n
    dd[i]=sqrt(x[i].^2+y[i].^2)
end
iso=sortperm(dd)[1:4]
# scatter(x[iso],y[iso])
q[iso].=-BigFloat("1")

@time a=Vall\q
maximum(abs.(Vall*a-q))
# scatter(x,y,zcolor=V*a)
minimum(V*a)
# scatter3d(x,y,V*a,label=string("deflection (w), p=",p),legend=:topleft)
# savefig(string("pdes_p_",p,".pdf"))

Qx=-(Vxxx+Vxyy)*a
Qy=-(Vyxx+Vyyy)*a
# scatter(x[i34],y[i34],zcolor=Qy)

mean(abs.(Qx))*2*L+mean(abs.(Qy))*2*L-1*dx*dx
