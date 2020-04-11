
using Profile #@profiler

#test allele_var_mtx
 # x=rand(5)
 # B=allele_var_mtx(x)
 # round.(B,digits=4)

pwd()
cd("D:\\Julia\\func")#pathway to your working directory
include("2_1allele_var_mtx.jl")
#test Gaussian_multinomial()
function simu_p(n)
    x=randn(n)
    x=x/sum(x)
    return x
 end
 xx=simu_p(100)
 N=1e6
include("2_2gaussian_multinomial.jl") #only need to include the last file
 #as former functions are called in this file by include("")
@time z=gaussian_multinomial(N,xx) #sum(xx)=1, input xx is the frequency of all alleles

sum(z)
println(xx)
plot(z)#same output as Matlab code

include("2_1allele_var_mtx.jl")
 include("2_2gaussian_multinomial.jl")
 include("2_3allele_frequency_gaussian.jl")
 using Profile
 # allele_frequency_Gaussian(K,N,twoNmu,s0,T)
@profiler x,X,X_deter,K=allele_frequency_Gaussian(1,1e6,10,0.05,2000)
println(size(x),size(X),size(X_deter),size(K))
@time x,X,X_deter,K=allele_frequency_Gaussian(1,1e6,10,0.05,2000)
,X,X_deter,K


#test2_1 to 2_4
include("2_1allele_var_mtx.jl")
 include("2_2gaussian_multinomial.jl")
 include("2_3allele_frequency_gaussian.jl")
 include("2_4num_origin_gaussian.jl")
 using Profile

@time ndx,simu_mtx,simu_mtx_se,semi_mtx,semi_mtx_se=num_origin_Gaussian([0.05; 0.005],[1],1e8,2)
# ndx,simu_mtx,simu_mtx_se,semi_mtx,semi_mtx_se,ewen_mtx,ewen_mtx_se=
display(ndx);display(simu_mtx);display(semi_mtx)
plot(semi_mtx)
plot!(simu_mtx)#one case is distict
pwd()
savefig("..\\photo\\distinct_semi_simu")#why?

@time ndx,simu_mtx1,simu_mtx_se1,semi_mtx1,semi_mtx_se1=num_origin_Gaussian([0.05; 0.005; 0.0005],[1; 10],1e6,50)#50
@time ndx,simu_mtx2,simu_mtx_se2,semi_mtx2,semi_mtx_se2=num_origin_Gaussian([0.05; 0.005; 0.0005],[1 10],1e7,10)#10
@time ndx,simu_mtx3,simu_mtx_se3,semi_mtx3,semi_mtx_se3=num_origin_Gaussian([0.05; 0.005; 0.0005],[1 10],1e8,100)#100

# 8185.972500 seconds (20.32 G allocations: 13.483 TiB, 11.48% gc time)
# 1876.964082 seconds (4.40 G allocations: 3.089 TiB, 11.45% gc time)
# 18864.345870 seconds (45.42 G allocations: 31.684 TiB, 10.74% gc time)
#8.05hr

using JLD #save(), load()
save("simu_semi_3test.jld","ndx",ndx,
"simu_mtx1",simu_mtx1,"simu_mtx_se1",simu_mtx_se1,
"semi_mtx1",semi_mtx1,"semi_mtx_se1",semi_mtx_se1,

"simu_mtx2",simu_mtx2,"simu_mtx_se2",simu_mtx_se2,
"semi_mtx2",semi_mtx2,"semi_mtx_se2",semi_mtx_se2,

"simu_mtx3",simu_mtx3,"simu_mtx_se3",simu_mtx_se3,
"semi_mtx3",semi_mtx3,"semi_mtx_se3",semi_mtx_se3)


pwd()
cd("D:\\Julia\\func")
@load "simu_semi_3test.jld"
ndx

using Plots
 #plot setting & subplot println("final plotting")
 legend=String[]
 s0=Float64.(ndx[1,:])
 twoNmu=Int64.(ndx[2,:])
 for i in 1:6
  push!(legend,"s0=$(s0[i]), 2N*mu=$(twoNmu[i])")
 end

 symbol=[:rect :rect :utriangle :utriangle :circle :circle]
 colors=[:green :cornflowerblue :red :orange :purple :steelblue]
 a=collect(201:200:2001)
 #plot1
 default(lab="")
 p1=plot(semi_mtx1,color=colors,linewidth=2.5)
 scatter!(a,simu_mtx1[a,:], markersize=5, color=color,markershape=symbol,yerror=simu_mtx_se1[a,:])
 xlabel!("Generation(t)")
 ylabel!("Average number of origins 𝛍(t) ")
 #plot2
 default(lab="")
 p2=plot(semi_mtx2,color=color,linewidth=2.5)
 scatter!(a,simu_mtx2[a,:], markersize=5, color=color,markershape=symbol,yerror=simu_mtx_se2[a,:])
 xlabel!("Generation(t)")
 #plot3
 default(lab="")
 p3=plot(semi_mtx3,color=color,linewidth=2.5)
 scatter!(a,simu_mtx3[a,:], markersize=5, color=color,markershape=symbol,
  yerror=simu_mtx_se3[a,:],label=legend, legend=:outertopright,
  foreground_color_legend = nothing,legendfontsize=9,grid=false)
 xlabel!("Generation(t)")

#final plot
gr()
 l = @layout [a{0.265w} b{0.265w} c{0.47w}]#use layout to creat a fake plot with legend
 using Plots.PlotMeasures
 plot(p1,p2,p3,
 layout=l,
 size=(1200,300),
 left_margin=7mm,bottom_margin=5mm,right_margin=-3mm,
 title=["N=10^6" "N=10^7" "N=10^8"],
 tickfont   = font("Times new roman", 6),
 yticks=0:5:50,
 titlefont  = font("Times new roman", 12),
 dpi=300,
 grid=false,
 ylims=(0,50),
 legendfontsize=6,
 background_color_legend=:transparent)
savefig("../photo/simu_semi_3")











savefig(p3,"../photo/10^8_simu_semi")
 #why simu almost twice of semi-Semideterministic prediction? count(>1/ns),ns as sample size

plot(p1,p2,p3, layout=(1,3))

 p.attr[:title_style][:fontsize] = 10.0
