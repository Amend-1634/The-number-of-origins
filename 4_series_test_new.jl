#call the functions in "4_4_figjoejyn_series.jl" to get the summary
 #matrix and plot

#mainly two series data: data at fixation and time-series data
 #for data at fixation you may need to estimate and provide a
  #slightly redundant T unlike time sampling data

##Fix: series as nDemes
cd("D:\\Julia\\func\\mig")

# twoNDememig_vec=[1e-4, 1e-3, 1e-2, 1e-1, 1]
# N=1e7;nDemes=100;NDeme=1e5
# mig_vec=twoNDememig_vec/2/NDeme
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]


##Spatial structure paramter collection
#3 influential factor: nDemes, mig, nu
##Spatial structure paramter collection
#nu-subplots: nDemes~mig
nu=1 #@fig6b_nu1 #tmux name in @clifford
T=15000; nDemes_vec=[500]
include("4_4_figjoejyn_series.jl")
include("4_3_1d_fixsamp.jl")
include("4_1_1nonlocal_frequency_sparse.jl")
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; M=10; Gauss=1;
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
  onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
using JLD
save("4_fig6b_nu1.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
"eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)

#will not do due to time limit---test whether it is approximating the binomial distribution
nu=30 #@fig6b_nu30
T=10000; nDemes_vec=[50,500]#8000 for the one omits jump, with jump the time could be longer
include("4_4_figjoejyn_series_threads.jl")
include("4_3_1d_fixsamp.jl")
# normal size
twoNDememig_vec=[1e-4, 1e-3, 1e-2, 1e-1, 1]
N=1e7;

mig_vec=twoNDememig_vec/N

twoNmu=1; s0=0.1; mig=1e-2; M=10; Gauss=1;
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
  onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,twoNDememig_vec,M,T,Gauss, nu)
using JLD
save("4_fig6b_nu30.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
"eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)




###############################


#fig6b. x dispersal, y time to fixation, series: nDemes
#test
N=5e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1; M=replicates=2;
T=3000; Gauss=1; nu=2
nDemes_vec=[5,10]; twoNDememig_vec=[5e-3,5e-4]
include("4_4_figjoejyn_series_threads.jl")
include("4_3_1d_fixsamp.jl")
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
  onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)


#tmux @4fig6b_100
T=15000; nDemes_vec=[10,50,100]#8000 for the one omits jump, with jump the time could be longer
include("4_4_figjoejyn_series.jl")
include("4_3_1d_fixsamp.jl")
# normal size
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=10; Gauss=1; nu=2
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
  onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
using JLD
save("4_fig6b_till100.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
"eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)

#tmux @4fig6b_200
T=20000; nDemes_vec=[10,50,100,200]#ralph-Coop model
include("4_4_figjoejyn_series.jl")
include("4_3_1d_fixsamp.jl")
# normal size
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=10; Gauss=1; nu=2
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
  onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
using JLD
save("4_fig6b_till200.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
"eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


#previous: 75901s=20hr 5min


using JLD
cd("D:\\Julia\\func\\mig")
readdir()
@load "4_fig6b_till100.jld"
@load "4_fig6b_till200.jld"

#time
using Plots.PlotMeasures #for mm
 mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2] #x
 colors=[:steelblue :red :orange]
 scatter(twoNDememig_vec, t_avn_mtx, markersize=4, label=legend,
 markershape=:rect, legend=:outertopright,legendfontsize=8,dpi=300,
 grid=false,fg_legend = :transparent,xscale=:log10,ylims=(0,8000),
 margin=8mm, size=(800,600), yerror=t_se_mtx,color=colors)
 plot!(twoNDememig_vec, t_avn_mtx,color=colors,lab="")
 xlabel!("2*NDeme*mig")
 ylabel!("Time to fixation (generation)")

savefig("..//..//photo//4_fig6b_t")


#eta
using Plots.PlotMeasures #for mm
  mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2] #x
  colors=[:steelblue :red :orange]
  scatter(twoNDememig_vec, eta_avn_mtx, markersize=4, label=legend,
  markershape=:rect, legend=:outertopright,legendfontsize=8,
  grid=false,fg_legend = :transparent,scale=:log10,dpi=300,
  margin=8mm, size=(800,600), yerror=eta_se_mtx,color=colors,
    guidefont   = font("Times new roman", 12))
  plot!(twoNDememig_vec, eta_avn_mtx,color=colors,lab="")
  xlabel!("2*NDeme*mig")
  ylabel!("Asymptotic number of origins")

  for i in 1:length(twoNDememig_vec), j in 1:length(colors)
    annotate!([(twoNDememig_vec[i], eta_avn_mtx[i,j], text(string(fix_M_mtx[i,j]),:bottom,12))]) #all fixed
  end #unlike joejyn's, all fixed

fix_M_mtx

#power law (related to Ralph-Coop)
a=float.(-4:0)
b=10 .^ a #as x, which is d in the equation
b
# c = (1/b) .^ (1/4)#1/b not functional here
e = 10^1.7 * b .^ (-1/4)
size(c)
# plot(b, c,scale=:log10)
plot!(b, e,scale=:log10,lab="",color=:black)
annotate!(0.1,10^2.25,text("eta=A*(1/mig)^1/4",:purple))

plot(c)
plot(b)

10^(-4)

#Ralph-Coop model
eta=A*(s/d)^(1/4)
#as s is always kept constant here
eta=A*(1/d)^(1/4)


savefig("..//..//photo//4_fig6b_eta")

###fix  sereise as nu
##nDemes=10    and     nDemes=50

#353s no need for test parameter set
nDemes=10;T=4000 #@nu_fix_10
nu_vec=[1 2 30]#so 95% range are [6.3 2.9 1.69]
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=replicates=50; Gauss=1; nu=2
# include("4_4_figjoejyn_series.jl")
include("4_4_figjoejyn_series_threads.jl")
include("4_3_1d_fixsamp.jl")

@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx =
 onedim_xdispersal_yfix_nuseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss, nu_vec)

using JLD
save("nu_fix_10.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
 "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


nDemes=50;T=10000 #@nu_fix_50
nu_vec=[1 2 30]
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=replicates=10;Gauss=1; nu=2
include("4_4_figjoejyn_series.jl")
include("4_3_1d_fixsamp.jl")

@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx =
 onedim_xdispersal_yfix_nuseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss, nu_vec)

using JLD
save("nu_fix_50.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
 "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


using JLD
cd("D:\\Julia\\func\\mig")
@load "nu_fix_10.jld"
@load "nu_fix_50.jld"

t_avn_mtx


#time
using Plots.PlotMeasures #for mm
 using Plots
 mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2] #x
 colors=[:steelblue :red :orange]
 scatter(twoNDememig_vec, t_avn_mtx, markersize=4, label=legend,
 markershape=:rect, legend=:outertopright,legendfontsize=8,dpi=300,
 grid=false,fg_legend = :transparent,xscale=:log10,
 margin=8mm, size=(800,600), yerror=t_se_mtx,color=colors)
 plot!(twoNDememig_vec, t_avn_mtx,color=colors,lab="")
 xlabel!("2*NDeme*mig")
 ylabel!("Time to fixation (generation)")

savefig("..//..//photo//nu_fix_t_10nDemes")


#eta
using Plots.PlotMeasures #for mm
  mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2] #x
  colors=[:steelblue :red :orange]
  scatter(twoNDememig_vec, eta_avn_mtx, markersize=4, label=legend,
  markershape=:rect, legend=:outertopright,legendfontsize=8,
  grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
  margin=8mm, size=(800,600), yerror=eta_se_mtx,color=colors)
  plot!(twoNDememig_vec, eta_avn_mtx,color=colors,lab="")
  xlabel!("2*NDeme*mig")
  ylabel!("Asymptotic number of origins")

savefig("..//..//photo//nu_fix_eta_10nDemes")


###y  as fixation eta and time, 𝝂 (nu) as series
##nDemes=10        and       nDemes=50

# #test
# nDemes=10
#  N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=replicates=2
#  T=1000; Gauss=1; nu=2
#  nu_vec=[5 20]
#  twoNDememig_vec=[1e-4, 1e-3]

nDemes=10;T=900 #@nu_fix_10
N=1e7; twoNmu=1; s0=0.1;
M=50; Gauss=1 #replicate number raise as it vary a lot in replicate=10
nu_vec=[1 2 30]#so 95% range are [6.3 2.9 1.69]
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
include("4_4_figjoejyn_series.jl")
include("4_3_1d_fixsamp.jl")
include("4_1_1nonlocal_frequency.jl")

@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx =
  onedim_xdispersal_yfix_nuseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss, nu_vec)
using JLD
save("nu_fix_10_2.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
 "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


nDemes=50;T=8000 #@nu_fix_50
N=1e7; twoNmu=1; s0=0.1;  M=replicates=10; Gauss=1;
nu_vec=[1 2 30]#so 95% range are [6.3 2.9 1.69]
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
include("4_4_figjoejyn_series.jl")
include("4_3_1d_fixsamp.jl")
include("4_1_1nonlocal_frequency_sparse.jl")

@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx =
  onedim_xdispersal_yfix_nuseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss, nu_vec)
using JLD
 save("nu_fix_50_2.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
  "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


using JLD
cd("D:\\Julia\\func\\mig")
@load "nu_fix_10.jld"
@load "nu_fix_50.jld"

t_avn_mtx


#time
cd("D:\\Julia\\func")

@load "nu_fix_10_2.jld"
 legend=["nDemes = 10, nu = 1"  "nDemes = 10, nu = 2"  "nDemes = 10, nu = 30"]
 using Plots.PlotMeasures #for mm
 using Plots
 mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
 colors=[:steelblue :red :orange]
 p2=scatter(mig_vec, t_avn_mtx, markersize=4, label=legend,
 markershape=:rect, legend=:topright,legendfontsize=12,dpi=300,
 grid=false,fg_legend = :transparent,xscale=:log10,
   guidefont   = font("Times new roman", 16),
 margin=8mm, size=(700,700), yerror=t_se_mtx,color=colors)
 plot!(mig_vec, t_avn_mtx,color=colors,lab="")
 xlabel!("Migration rate (mig)")
 ylabel!("Time to fixation (generation)")
 title!("(b) Time to fixation")

 @load "nu_fix_50_2.jld"
 legend=["nDemes = 50, nu = 1"  "nDemes = 50, nu = 2"  "nDemes = 50, nu = 30"]
 scatter!(mig_vec, t_avn_mtx, markersize=4, label=legend,
 markershape=:circle, legend=:topright,legendfontsize=14,dpi=300,
 grid=false,fg_legend = :transparent,xscale=:log10,
   guidefont   = font("Times new roman", 16),
 margin=8mm, size=(700,700), yerror=t_se_mtx,color=colors)
 plot!(mig_vec, t_avn_mtx,color=colors,lab="",line=:dash)

savefig("..//..//photo//nu_fix_10nDemes")


#eta
@load "nu_fix_10_2.jld"
  legend=["nDemes = 10, nu = 1"  "nDemes = 10, nu = 2"  "nDemes = 10, nu = 30"]
  using Plots.PlotMeasures #for mm
  mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
  colors=[:steelblue :red :orange]
  p1=scatter(mig_vec, eta_avn_mtx, markersize=4, label="",
  markershape=:rect, legend=:topright,legendfontsize=11,
  grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
    guidefont   = font("Times new roman", 16),
  margin=10mm, size=(600,700), yerror=eta_se_mtx,color=colors)
  plot!(mig_vec, eta_avn_mtx,color=colors,lab="")
  xlabel!("Migration rate (mig)")
  ylabel!("Asymptotic number of origins")
  title!("(a) Asymptotic number of origins")

  xlabel!("Migration rate (mig)")
  ylabel!("Time to fixation (generation)")
  title!("(b) Time to fixation")


  @load "nu_fix_50_2.jld"
  legend=["nDemes = 50, nu = 1"  "nDemes = 50, nu = 2"  "nDemes = 50, nu = 30"]
  scatter!(mig_vec, t_avn_mtx, markersize=4, label=legend,
  markershape=:circle, legend=:topright,legendfontsize=14,dpi=300,
  grid=false,fg_legend = :transparent,xscale=:log10,
    guidefont   = font("Times new roman", 16),
  margin=8mm, size=(700,700), yerror=t_se_mtx,color=colors)
  plot!(mig_vec, t_avn_mtx,color=colors,lab="",line=:dash)

 savefig("..//..//photo//nu_fix_10nDemes")

  @load "nu_fix_50_2.jld"
  legend=["nDemes = 50, nu = 1"  "nDemes = 50, nu = 2"  "nDemes = 50, nu = 30"]
  using Plots.PlotMeasures #for mm
    mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
    colors=[:steelblue :red :orange]
    scatter!(mig_vec, eta_avn_mtx, markersize=4, label="",
    markershape=:circle,legendfontsize=15,
    grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
      guidefont   = font("Times new roman", 16),
    margin=8mm, size=(600,700), yerror=eta_se_mtx,color=colors)
    plot!(mig_vec, eta_avn_mtx,color=colors,lab="",line=:dash)
    xlabel!("Migration rate (mig)")
    ylabel!("Asymptotic number of origins")

plot(p1,p2,layout=(1,2),size=(1300,600))


  #power law (related to Ralph-Coop)
a=float.(-5.5:0.1:-5)
  b=10 .^ a #as x, which is d in the equation
  e =   b .^ (-1/4)
  plot!(b, e,scale=:log10,lab="",color=:black)

savefig("..//..//photo//nu_fix_50nDemes")
###tsamp, y  Average eta, nDemes as series
##nDemes*mig[5e-8, 5e-4]

cd("D:\\Julia\\func\\mig")


mig_vec=[5e-8, 5e-4]#@3_t_mig8
T=2500
nDemes_vec=[10,50,100]
N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1;
M=replicates=10;Gauss=1
include("2_4num_origin_gaussian.jl") #non-spatial
include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
include("3_4_figjoejyn_series.jl")
@time tsamp, legend, avn_mtx, se_mtx =
  onedim_series_mig_nDemes(N,nDemes_vec,twoNmu,s0,mig_vec,M,T)

using JLD
save("3_fig6a_mig_8_4.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
      "legend", legend)


 T=3000;  #@4_fig6a_3000_nu1 #@3000_nu1
 nDemes_vec=[10,50,100]
 N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1;mig=1e-2;
 M=replicates=10;Gauss=1; nu=1
 include("2_4num_origin_gaussian.jl")
 include("4_2_1d_tsamp.jl")
 include("4_4_figjoejyn_series.jl")
 @time tsamp, legend, avn_mtx, se_mtx = onedim_nDemes_series(nDemes_vec,N,twoNmu,s0,mig, M, T, Gauss, nu)
 #6266s
 using JLD
  save("4_fig6a_3000_nu1.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
  "legend", legend)


  T=600;  #3_fig6a
    nDemes_vec=[10,50,100]
    N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1;mig=1e-2;
     M=replicates=50;Gauss=1
    include("2_4num_origin_gaussian.jl")
    include("3_2_1d_tsamp.jl")
    include("3_4_figjoejyn_series.jl")
    @time tsamp, legend, avn_mtx, se_mtx = onedim_nDemes_series(nDemes_vec,N,twoNmu,
                                twoNmu_seq,s0,s0_seq,mig, M, replicates, T, Gauss)

using JLD
   save("3_fig6a_600_2.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
   "legend", legend)

 T=600;  #@4_fig6a_600
   nDemes_vec=[10,50,100]
   N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1;mig=1e-2;
    M=replicates=10;Gauss=1; nu=1
   include("2_4num_origin_gaussian.jl")
   include("4_2_1d_tsamp.jl")
   include("4_4_figjoejyn_series.jl")
 @time tsamp, legend, avn_mtx, se_mtx = onedim_nDemes_series(nDemes_vec,N,twoNmu,s0,mig, M, T, Gauss, nu)
 using JLD
  save("4_fig6a_600.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
  "legend", legend)

T=600;  #@4_fig6a_600
  nDemes_vec=[10,50,100]
  N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1;mig=1e-2;
   M=replicates=50;Gauss=1; nu=1
  include("2_4num_origin_gaussian.jl")
  include("4_2_1d_tsamp.jl")
  include("4_4_figjoejyn_series.jl")
@time tsamp, legend, avn_mtx, se_mtx = onedim_nDemes_series(nDemes_vec,N,twoNmu,s0,mig, M, T, Gauss, nu)
using JLD
 save("4_fig6a_600_nu1_50M.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
 "legend", legend)

T=600;  #@4_fig6a_600
  nDemes_vec=[10,50,100]
  N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1;mig=1e-2;
   M=replicates=50;Gauss=1; nu=4
  include("2_4num_origin_gaussian.jl")
  include("4_2_1d_tsamp.jl")
  include("4_4_figjoejyn_series.jl")
@time tsamp, legend, avn_mtx, se_mtx = onedim_nDemes_series(nDemes_vec,N,twoNmu,s0,mig, M, T, Gauss, nu)
using JLD
 save("4_fig6a_600_nu4.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
 "legend", legend)


using JLD
@load "3_fig6a_3000.jld"
@load "3_fig6a_3000_2.jld"

@load "3_fig6a_600.jld"



using Plots
 using Plots.PlotMeasures #for mm
 colors=[:steelblue :red :orange :black]
 scatter!(tsamp, avn_mtx, markersize=2, label=legend,
 color=colors,markershape=:circle, legend=:outertopright,
 legendfontsize=8, grid=false,fg_legend = :transparent,
 size=(600,600), yerror=se_mtx,margin=5mm,dpi=300)
 plot!(tsamp, avn_mtx,color=colors,lab="",line=:dash)
 xlabel!("Generation(t)")
 ylabel!("Average number of origins 𝛍(t) ")
 # scatter!(tsamp[1], avn_mtx[1,1],color=colors[1],markersize=:circle,label="local migration")

 # title!("Local migration")

@load "4_fig6a_3000.jld"
using JLD
cd("D:\\Julia\\func\\mig")
@load "4_fig6a_3000_nu1.jld"

cd("D:\\Julia\\func\\mig")
using JLD
@load "4_fig6a_600.jld"#nu=4
@load "4_fig6a_600_nu1.jld"#nu=1  #2330s

#M=50
@load "4_fig6a_600_nu1_50M.jld"#nu=1 #differ from last file #employed
@load "4_fig6a_600_nu4.jld"#M=50  #2330s


@load "4_fig6a_3000_nu1.jld"




using Plots.PlotMeasures #for mm
  using Plots
  colors=[:steelblue :red :orange :black]
  scatter(tsamp, avn_mtx, markersize=2,label="",
  color=colors,markershape=:rect, legend=:outertopright,
  legendfontsize=8, grid=false,fg_legend = :transparent,
  size=(500,500), yerror=se_mtx, margin=5mm, dpi=300)
  plot!(tsamp, avn_mtx,color=colors,lab="")
  xlabel!("Generation(t)")
  ylabel!("Average number of origins 𝛍(t) ")

savefig("..//..//photo//both_fig6a")


  # annotate!("▯ t distribution")
  # scatter!([tsamp[1]], [avn_mtx[1,1]],color=colors[1],markersize=:rect,label="t distribution")
  # title!("long-distance dispersal")
  # title!("Non-local migration")

##tsamp, oscilalting population peak-to-trough ration (𝛟)phi as series

# #test
# twoNDememig_vec=[1e-4, 1e-3]
# phi_vec=[10 50] #peak-to-trough ratio
# s0=0.1;nDemes=5;twoNmu=1;M=10;T=3000;Gauss=1;N=1e7;nu=2
# include("4_3_1d_fixsamp.jl")
# include("4_3_1d_fixsamp_osci.jl")
# include("4_4_figjoejyn_series.jl")
# @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx = dispersal_oscillating(N,
#   phi_vec, har_geo,nDemes,twoNmu,s0,M,T,Gauss, twoNDememig_vec, nu)

# phi=10; har_geo="H" #constrained harmonic mean of populaiton size #@osci_H
# phi=10; har_geo="G" #constrained geometric mean #@osci_G

har_geo="H" #osci_h
include("4_3_1d_fixsamp.jl")
include("4_3_1d_fixsamp_osci.jl")
include("4_4_figjoejyn_series.jl")
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
phi_vec=[10 50 100] #peak-to-trough ratio
s0=0.1;twoNmu=1;M=10;T=3000;Gauss=1;N=1e7;nu=2;nDemes=10
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx = dispersal_oscillating(N,
  phi_vec, har_geo,nDemes,twoNmu,s0,M,T,Gauss, twoNDememig_vec, nu)
using JLD
save("4_fig10_h.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
    "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx)

har_geo="G" #osci_g
include("4_3_1d_fixsamp.jl")
include("4_3_1d_fixsamp_osci.jl")
include("4_4_figjoejyn_series.jl")
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
phi_vec=[10 50 100] #peak-to-trough ratio
s0=0.1;twoNmu=1;M=10;T=3000;Gauss=1;N=1e7;nu=2;nDemes=10
@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx = dispersal_oscillating(N,
  phi_vec, har_geo,nDemes,twoNmu,s0,M,T,Gauss, twoNDememig_vec, nu)
using JLD
save("4_fig10_g.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
    "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx)



using JLD
@load "4_fig10.jld"

using Plots.PlotMeasures #for mm
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
colors=[:steelblue :red :orange :black] #
eta_avn_mtx#6*4 (mig,phi)
legend=string.(legend)
typeof(legend)
legend=reshape(legend,1,4)

using Plots
scatter(mig_vec,eta_avn_mtx,label=legend,color=colors,xscale=:log10,dpi=300,
        grid=false,fg_legend = :transparent,legend=:bottomleft)
plot!(mig_vec,eta_avn_mtx,color=colors,lab="")
xlabel!("Dispersal rate")
ylabel!("Asymptotic number of origins")

savefig("..//photo//fig10")


using Plots
scatter(mig_vec, t_avn_mtx,label=legend,color=colors,xscale=:log10,dpi=300,
        grid=false,fg_legend = :transparent,legend=:bottomleft)
plot!(mig_vec,t_avn_mtx,color=colors,lab="")
xlabel!("Dispersal rate")
ylabel!("Mean time to fixation")

savefig("..//photo//fig10_t")

###tsamp y total frequency
##nDemes=10           nDemes=50
cd("D:\\Julia\\func\\mig")


#normal size
nDemes=10; T=600
include("4_2_1d_tsamp.jl")
include("4_4_figjoejyn_series.jl")
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=replicates=10;Gauss=1; nu=1
@time tsamp, legend, X_avn_mtx, X_se_mtx, X_medn_mtx, X_lqrn_mtx, X_uqrn_mtx =
  xtime_yX_twoNDememigseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss,nu)

using JLD
save("tot_freq10.jld","tsamp", tsamp, "legend", legend, "X_avn_mtx", X_avn_mtx,
 "X_se_mtx", X_se_mtx, "X_medn_mtx", X_medn_mtx, "X_lqrn_mtx", X_lqrn_mtx, "X_uqrn_mtx",X_uqrn_mtx)

nDemes=50; T=3500 #10903s=3hr
include("4_2_1d_tsamp.jl")
include("4_4_figjoejyn_series.jl")
mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
N=1e7; twoNmu=1; s0=0.1; mig=1e-2; M=replicates=10;Gauss=1; nu=1
@time tsamp, legend, X_avn_mtx, X_se_mtx, X_medn_mtx, X_lqrn_mtx, X_uqrn_mtx =
  xtime_yX_twoNDememigseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss,nu)

using JLD
 save("tot_freq50.jld","tsamp", tsamp, "legend", legend, "X_avn_mtx", X_avn_mtx,
  "X_se_mtx", X_se_mtx, "X_medn_mtx", X_medn_mtx, "X_lqrn_mtx", X_lqrn_mtx, "X_uqrn_mtx",X_uqrn_mtx)


# #test
# twoNDememig_vec=[1e-4, 1e-3]
# N=1e4; twoNmu=1; s0=0.1; mig=1e-2; M=2
# T=1000; Gauss=1; nu=4
# nDemes=10

@time tsamp, legend, X_avn_mtx, X_se_mtx, X_medn_mtx, X_lqrn_mtx, X_uqrn_mtx =
  xtime_yX_twoNDememigseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss,nu)

legend


pwd()
cd("D:\\Julia\\func\\mig")
using JLD
@load "tot_freq10.jld"
@load "tot_freq50.jld"
tsamp

using Plots.PlotMeasures #for mm
 using Plots
 colors=[:steelblue :red :orange :purple :green]
 scatter(tsamp, X_avn_mtx, markersize=2, label=legend,
 markershape=:rect, legend=:outertopright,legendfontsize=8,
 grid=false,fg_legend = :transparent,dpi=300,
 margin=15mm, size=(800,600), yerror=X_se_mtx,color=colors)
 plot!(tsamp, X_avn_mtx,color=colors,lab="")
 xlabel!("Generations")
 ylabel!("The mean total frequency of mutants")

savefig("..//..//photo//4_total_freq10")



X_avn_mtx
X_medn_mtx
#checked, all right


##
@load "nu_fix_50_2.jld"
  legend= ["v = 1"  "v = 2"  "v = 30"]
  p1=scatter(mig_vec, t_avn_mtx, markersize=3, label="",
  markershape=:rect, legend=:topright,legendfontsize=14,dpi=300,
  grid=false,fg_legend = :transparent,xscale=:log10,
    guidefont   = font("Times new roman", 16),
    titlefont   = font("Times new roman", 18),
  margin=14mm, size=(700,700), yerror=t_se_mtx,color=colors)
  plot!(mig_vec, t_avn_mtx,color=colors,lab="",linewidth=2,alpha=0.5)
  xlabel!("Migration rate (mig)")
  ylabel!("Time to fixation (generation)")
  title!("(a) Time to fixation")
  # savefig("..//..//photo//nu_fix_10nDemes")

  @load "nu_fix_50_2.jld"
  legend= ["v = 1"  "v = 2"  "v = 30"]
  using Plots.PlotMeasures #for mm
  mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
  colors=[:steelblue :red :orange]
  p2=scatter(mig_vec, eta_avn_mtx, markersize=2, label=legend,
  markershape=:rect,legendfontsize=15,
  titlefont   = font("Times new roman", 18),
  grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
    guidefont   = font("Times new roman", 16),
  margin=0mm, size=(600,700), yerror=eta_se_mtx,color=colors)
  plot!(mig_vec, eta_avn_mtx,color=colors,lab="",linewidth=2,alpha=0.5)
  xlabel!("Migration rate (mig)")
  ylabel!("Asymptotic number of origins")
  title!("(b) Asymptotic number of origins")

  plot(p1,p2,layout=(1,2),size=(1300,600))


@load "nu_fix_10_2.jld"
    legend= ["v = 1"  "v = 2"  "v = 30"]
    p1=scatter(mig_vec, t_avn_mtx, markersize=3, label="",
    markershape=:rect, legend=:topright,legendfontsize=14,dpi=300,
    grid=false,fg_legend = :transparent,xscale=:log10,
      guidefont   = font("Times new roman", 16),
      titlefont   = font("Times new roman", 18),
    margin=14mm, size=(700,700), yerror=t_se_mtx,color=colors)
    plot!(mig_vec, t_avn_mtx,color=colors,lab="",linewidth=2,alpha=0.5)
    xlabel!("Migration rate (mig)")
    ylabel!("Time to fixation (generation)")
    title!("(a) Time to fixation")
    # savefig("..//..//photo//nu_fix_10nDemes")

    @load "nu_fix_10_2.jld"
    legend= ["v = 1"  "v = 2"  "v = 30"]
    using Plots.PlotMeasures #for mm
    mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
    colors=[:steelblue :red :orange]
    p2=scatter(mig_vec, eta_avn_mtx, markersize=2, label=legend,
    markershape=:rect,legendfontsize=15,
    titlefont   = font("Times new roman", 18),
    grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
      guidefont   = font("Times new roman", 16),
    margin=0mm, size=(600,700), yerror=eta_se_mtx,color=colors)
    plot!(mig_vec, eta_avn_mtx,color=colors,lab="",linewidth=2,alpha=0.5)
    xlabel!("Migration rate (mig)")
    ylabel!("Asymptotic number of origins")
    title!("(b) Asymptotic number of origins")

    plot(p1,p2,layout=(1,2),size=(1300,600))


    #fixation

    ###nu=1
    ##1. #@nest1
    nu=1;s0=0.05;nDemes_vec=[10,50,100]
    T=7500 #as s lower=>longer time(according to joejyn's graph, 2000 generation more for 0.05 than 0.1)
    mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
    N=1e7; twoNmu=1; M=10; Gauss=1;
    include("4_4_figjoejyn_series_threads.jl")
    include("4_3_1d_fixsamp.jl")

    @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
      onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
    using JLD
    save("fix_nest1.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
    "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


    ##2.#@nest2
    nu=1;s=0.1;nDemes_vec=[10,50,100]
    T=10000 #as s lower=>longer time(according to joejyn's graph, 2000 generation more for 0.05 than 0.1)
    mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
    N=1e7; twoNmu=1; s0=0.1; M=10; Gauss=1;
    include("4_4_figjoejyn_series_threads.jl")
    include("4_3_1d_fixsamp.jl")

    @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
      onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
    using JLD
    save("fix_nest2.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
    "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)

    #2_2 #@nest_2_2
    nu=1;s0=0.1;nDemes_vec=[500]
    T=15000 #as s lower=>longer time(according to joejyn's graph, 2000 generation more for 0.05 than 0.1)
    mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
    N=1e7; twoNmu=1; M=10; Gauss=1;
    include("4_4_figjoejyn_series.jl")
    include("4_3_1d_fixsamp.jl")
    include("4_1_1nonlocal_frequency_sparse.jl")

    @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
      onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
    using JLD
    save("fix_nest2_2.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
    "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)
    #2_2

    #2_3
    s=0.01

    ###local
    ##3          #@nest3
    s0=0.05;nDemes_vec=[10,50,100]
    # T=7500 #no fixation

    mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
    N=1e7; twoNmu=1; M=10; Gauss=1;
    include("3_1_2notrack_cons_sparse.jl")
    include("3_3_1d_fixsamp.jl")
    include("3_4_figjoejyn_series.jl")

    @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
          onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
    using JLD
    save("fix_nest3.jld","legend",legend,"t_avn_mtx",t_avn_mtx,
      "t_se_mtx", t_se_mtx, "eta_avn_mtx",eta_avn_mtx,
      "eta_se_mtx",eta_se_mtx,"fix_M_mtx", fix_M_mtx)

    ##4    #@nest4 #16.6hr
    s0=0.1;nDemes_vec=[10,50,100]
    T=7000
    mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
    N=1e7; twoNmu=1; M=10; Gauss=1;
    include("3_1_2notrack_cons.jl")
    include("3_3_1d_fixsamp.jl")
    include("3_4_figjoejyn_series.jl")

    @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
          onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
    using JLD
    save("fix_nest4.jld","legend",legend,"t_avn_mtx",t_avn_mtx,
      "t_se_mtx", t_se_mtx, "eta_avn_mtx",eta_avn_mtx,
      "eta_se_mtx",eta_se_mtx,"fix_M_mtx", fix_M_mtx)

    #@nest4_2
    s0=0.1;nDemes_vec=[500]
      T=15000
      mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2]
      N=1e7; twoNmu=1; M=10; Gauss=1;
      include("3_1_2notrack_cons_sparse.jl")
      include("3_3_1d_fixsamp.jl")
      include("3_4_figjoejyn_series.jl")

      @time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
            onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
      using JLD
      save("fix_nest4_2.jld","legend",legend,"t_avn_mtx",t_avn_mtx,
        "t_se_mtx", t_se_mtx, "eta_avn_mtx",eta_avn_mtx,
        "eta_se_mtx",eta_se_mtx,"fix_M_mtx", fix_M_mtx)




    using JLD
    cd("D:\\Julia\\func")
    readdir()
    @load "fix_nest4.jld"
    legend
    t_avn_mtx
    mig_vec

    fix_M_mtx

    #time
    using Plots.PlotMeasures #for mm
     using Plots
     mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]#x
     colors=[:steelblue :red :orange]
     scatter(mig_vec, t_avn_mtx, markersize=4, label=legend,
     markershape=:rect, legend=:outertopright,legendfontsize=8,dpi=300,
     grid=false,fg_legend = :transparent,xscale=:log10,ylims=(0,6000),
     margin=8mm, size=(800,600), yerror=t_se_mtx,color=colors)
     plot!(mig_vec, t_avn_mtx,color=colors,lab="")
     xlabel!("migration rate (mig)")
     ylabel!("Time to fixation (generation)")
     title!("s=0.1")


    scatter(mig_vec, t_avn_mtx)

    using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]#x
      colors=[:steelblue :red :orange]
      scatter(mig_vec, t_avn_mtx, markersize=4, label=legend,
      markershape=:circle, legend=:outertopright,legendfontsize=8,dpi=300,
      grid=false,fg_legend = :transparent,xscale=:log10,ylims=(0,8000),
      margin=8mm, size=(800,600), yerror=t_se_mtx,color=colors)
      plot!(mig_vec, t_avn_mtx,color=colors,lab="",line=:dash)
      xlabel!("migration rate (mig)")
      ylabel!("Time to fixation (generation)")




    savefig("..//..//photo//4_fig6b_t")


    #eta
    using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]#x
      colors=[:steelblue :red :orange]
      scatter(mig_vec, eta_avn_mtx, markersize=4, label=legend,
      markershape=:rect, legend=:outertopright,legendfontsize=8,
      grid=false,fg_legend = :transparent,scale=:log10,dpi=300,
      margin=8mm, size=(800,600), yerror=eta_se_mtx,color=colors,
        guidefont   = font("Times new roman", 12))
      plot!(mig_vec, eta_avn_mtx,color=colors,lab="")
      xlabel!("migration rate (mig)")
      ylabel!("Asymptotic number of origins")

    # fix_M_mtx

    #how twoNDememig change with nDemes
    mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]#x
    nDemes_vec=[10,50,100]
    N=1e7

    twoNDememig=N/nDemes_vec .* mig_vec #5*3
    legend

    using Plots.PlotMeasures #for mm
     colors=[:steelblue :red :orange]
     plot(mig_vec,twoNDememig,scale=:log10,label="",grid=false,
      fg_legend = :transparent,margin=8mm,legendfontsize=12,linewidth=2,alpha=0.5,
      color=colors)
      xlabel!("Migration rate (mig)")
      ylabel!("Migration from one deme")
    savefig("..\\photo\\twoNDememig")
    #so higher deme number=>higher migration cases between deme





    #power law (related to Ralph-Coop)
    a=float.(mig_vec)
     e =  a .^ (-1/4)
     plot!(mig_vec, e,scale=:log10,lab="",color=:green,label="power law")
     # annotate!(1e-7,10^2.15,text("eta=A*(1/mig)^1/4",:purple,:right))

    plot(c)
    plot(b)

    10^(-4)

    #Ralph-Coop model
    eta=A*(s/d)^(1/4)
    #as s is always kept constant here
    eta=A*(1/d)^(1/4)


    savefig("..//..//photo//4_fig6b_eta")

    NDeme=N/nDemes
    twoND


    ##
    s0_vec=[1e-5, 1e-4, 1e-3, 1e-2, 0.1]
    mig_vec=[5e-8, 5e-4]
    N_vec=[1e7, 1e9]
    nDemes=100 #in large deme number


    #compare s
    @load "fix_nest1.jld" #s=0.05, local
    @load "fix_nest2.jld" #s=0.05, local
    fix_M_mtx

    @load "fix_nest3.jld" #s=0.05, local
    @load "fix_nest4.jld" #s=0.05, local
    fix_M_mtx2

    #all are fixed

    #
    legend1=["s = 0.05, nDemes = 10"  "s = 0.05, nDemes = 50"  "s = 0.05, nDemes = 100"]
     legend2=["s = 0.1, nDemes = 10"  "s = 0.1, nDemes = 50"  "s = 0.1, nDemes = 100"]
    pwd()
    cd("D:\\Julia\\func")
     #time
    using JLD
     @load "fix_nest3.jld" #s=0.05, local
     using Plots.PlotMeasures #for mm
     using Plots
     mig_vec=[5e-10, 5e-9, 5e-8, 5e-7, 5e-6] #x
     colors=[:steelblue :red :orange]
     scatter(mig_vec, t_avn_mtx2, markersize=4, label=legend1,
     markershape=:rect, legend=:topright,legendfontsize=14,dpi=300,
     grid=false,fg_legend = :transparent,xscale=:log10,ylims=(0,7000),
     margin=8mm, size=(600,600), yerror=t_se_mtx2,color=colors)
     plot!(mig_vec, t_avn_mtx2,color=colors,lab="")
     xlabel!("migration rate (mig)")
     ylabel!("Time to fixation (generation)")

     @load "fix_nest4.jld" #s=0.1, local
      scatter!(mig_vec, t_avn_mtx2, markersize=4, label=legend2,
      markershape=:circle,
      grid=false,fg_legend = :transparent,xscale=:log10, yerror=t_se_mtx2,color=colors)
      plot!(mig_vec, t_avn_mtx2,color=colors,lab="",line=:dash)



    savefig("..//..//photo//4_fig6b_t")


    #eta
    @load "fix_nest3.jld" #s=0.05, local
     using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-9, 5e-8, 5e-7, 5e-6] #x
      colors=[:steelblue :red :orange]
      scatter(mig_vec, eta_avn_mtx2, markersize=4, label=legend1,
      markershape=:rect, legend=:topright,legendfontsize=12,
      grid=false,fg_legend = :transparent,scale=:log10,dpi=300,
      margin=8mm, size=(700,600), yerror=eta_se_mtx2,color=colors,
        guidefont   = font("Times new roman", 12))
      plot!(mig_vec, eta_avn_mtx2,color=colors,lab="")
      xlabel!("migration rate (mig)")
      ylabel!("Asymptotic number of origins")

    @load "fix_nest4.jld" #s=0.1, local
      using Plots.PlotMeasures #for mm
        using Plots
        mig_vec=[5e-10, 5e-9, 5e-8, 5e-7, 5e-6] #x
        colors=[:steelblue :red :orange]
        scatter(mig_vec, eta_avn_mtx2, markersize=4, label=legend2,
        markershape=:circle,grid=false,fg_legend = :transparent,scale=:log10,dpi=300,
        margin=8mm,yerror=eta_se_mtx2,color=colors,size=(600,600),legendfontsize=12,
          guidefont   = font("Times new roman", 12))
        plot!(mig_vec, eta_avn_mtx2,color=colors,lab="",line=:dash)
        xlabel!("migration rate (mig)")
        ylabel!("Asymptotic number of origins")

        #power law (related to Ralph-Coop)
        a=float.(mig_vec)
         e =  a .^ (-1/4)
        plot!(mig_vec, e,scale=:log10,lab="",color=:green,label="power law",linewidth=2,
         alpha=0.5)
         # annotate!(1e-7,10^2.15,text("eta=A*(1/mig)^1/4",:purple,:right))


    #but can this prove that s isn't influential? As 0.05 and 0.1 slight different
    #from each other---even in time


    ##non-local

    #time
    cd("D://Julia/func")

    legend1=["s = 0.05, nDemes = 10"  "s = 0.05, nDemes = 50"  "s = 0.05, nDemes = 100"]
     legend2=["s = 0.1, nDemes = 10"  "s = 0.1, nDemes = 50"  "s = 0.1, nDemes = 100"]
     @load "fix_nest1.jld" #problematic
     #guess:
    #as stopped in the 11th simulation=>could be lower


     using Plots.PlotMeasures #for mm
     using Plots
     mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]
     colors=[:steelblue :red :orange]
     scatter(mig_vec, t_avn_mtx, markersize=4, label=legend1,
       markershape=:rect, legend=:topright,legendfontsize=12,dpi=300,
       grid=false,fg_legend = :transparent,xscale=:log10,
       margin=8mm, size=(800,600), yerror=t_se_mtx,color=colors)
     plot!(mig_vec, t_avn_mtx,color=colors,lab="")
     xlabel!("migration rate (mig)")
     ylabel!("Time to fixation (generation)")
     # title!("s=")



     @load "fix_nest2.jld"
     using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]
      colors=[:steelblue :red :orange]
      scatter!(mig_vec, t_avn_mtx, markersize=4, label=legend2,
      markershape=:circle,dpi=300,
      grid=false,fg_legend = :transparent,xscale=:log10,ylims=(0,7000),
      margin=8mm, size=(600,600), yerror=t_se_mtx,color=colors)
      plot!(mig_vec, t_avn_mtx,color=colors,lab="",line=:dash)


    # savefig("..//..//photo//4_fig6b_t")

    legend
    #eta
    @load "fix_nest1.jld"
      using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]
      colors=[:steelblue :red :orange]
      scatter(mig_vec, eta_avn_mtx, markersize=4, label=legend1,
        markershape=:rect, legend=:topright,legendfontsize=12,
        grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
        margin=8mm, size=(600,600), yerror=eta_se_mtx,color=colors,
        guidefont   = font("Times new roman", 12))
      plot!(mig_vec, eta_avn_mtx,color=colors,lab="")

      @load "fix_nest2.jld"
      scatter!(mig_vec, eta_avn_mtx, markersize=4, label=legend2,
      markershape=:circle, grid=false,fg_legend = :transparent,dpi=300,
      margin=8mm, yerror=eta_se_mtx,color=colors,
        guidefont   = font("Times new roman", 12))
      plot!(mig_vec, eta_avn_mtx,color=colors,lab="",line=:dash)

      xlabel!("migration rate (mig)")
      ylabel!("Asymptotic number of origins")



    #power law (related to Ralph-Coop)
    a=float.(mig_vec)
     e =  a .^ (-1/4)
    plot!(mig_vec, e,scale=:log10,lab="",color=:green,label="power law",linewidth=2,
     alpha=0.5)
     # annotate!(1e-7,10^2.15,text("eta=A*(1/mig)^1/4",:purple,:right))

    plot(c)
    plot(b)

    10^(-4)

    #Ralph-Coop model
    eta=A*(s/d)^(1/4)
    #as s is always kept constant here
    eta=A*(1/d)^(1/4)


    ##compare local and non-local
    #time
    #s=0.05
    using JLD
     legend1=["s = 0.05, nDemes = 10"  "s = 0.05, nDemes = 50"  "s = 0.05, nDemes = 100"]
     @load "fix_nest1.jld" #s=0.05, non-local
     using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-9, 5e-8, 5e-7, 5e-6] #x
      colors=[:steelblue :red :orange]
      scatter(mig_vec, t_avn_mtx, markersize=4, label=legend1,xticks=mig_vec,
      markershape=:rect, legend=:topright,legendfontsize=12,dpi=300,
      grid=false,fg_legend = :transparent,xscale=:log10,ylims=(0,7000),
      margin=8mm, size=(600,600), yerror=t_se_mtx,color=colors)
      plot!(mig_vec, t_avn_mtx,color=colors,lab="")
      xlabel!("migration rate (mig)")
      ylabel!("Time to fixation (generation)")

     @load "fix_nest3.jld" #s=0.05, local
     scatter!(mig_vec, t_avn_mtx2, markersize=4, label="",
       markershape=:circle,yerror=t_se_mtx2,color=colors)
       plot!(mig_vec, t_avn_mtx2,color=colors,lab="",line=:dash)


    #s=0.1
    cd("D://Julia//func")

    #compare local and non-local
    cd("D://Julia//func")
    using JLD
     legend1=["non-local, nDemes = 10"  "non-local, nDemes = 50"  "non-local, nDemes = 100"]
     legend2=["local, nDemes = 10"  "local, nDemes = 50"  "local, nDemes = 100"]
     @load "fix_nest2.jld" #non-local
     using Plots.PlotMeasures #for mm
      using Plots
      mig_vec=[5e-10, 5e-9, 5e-8, 5e-7, 5e-6] #x
      colors=[:steelblue :red :orange]
      scatter(mig_vec, t_avn_mtx, markersize=4, label=legend1,
      markershape=:rect, legend=:topright,legendfontsize=12,dpi=300,
      grid=false,fg_legend = :transparent,scale=:log10,ylims=(10^2.4,1e4),
      margin=8mm, size=(600,600), yerror=t_se_mtx,color=colors)
      plot!(mig_vec, t_avn_mtx,color=colors,lab="")
      xlabel!("migration rate (mig)")
      ylabel!("Time to fixation (generation)")

     @load "fix_nest4.jld" #local
     scatter!(mig_vec, t_avn_mtx2, markersize=4, label=legend2,xticks=mig_vec,
       markershape=:circle,yerror=t_se_mtx2,color=colors)
       plot!(mig_vec, t_avn_mtx2,color=colors,lab="",line=:dash)

    ##eta

    #s=0.1
    legend1=["non-local, nDemes = 10"  "non-local, nDemes = 50"  "non-local, nDemes = 100"]
        legend2=["local, nDemes = 10"  "local, nDemes = 50"  "local, nDemes = 100"]
        using JLD
        mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2]
        @load "fix_nest2.jld" # non-local
        using Plots.PlotMeasures #for mm
         using Plots
         colors=[:steelblue :red :orange]
         scatter(mig_vec, eta_avn_mtx, markersize=4, label=legend1,
         markershape=:rect, legend=:topright,legendfontsize=13,dpi=300,
         grid=false,fg_legend = :transparent,xscale=:log10,xlims=(1e-10,50),
         margin=8mm, size=(600,600), yerror=eta_se_mtx,color=colors)
         plot!(mig_vec, eta_avn_mtx,color=colors,lab="")
         xlabel!("migration rate (mig)")
         ylabel!("Asymptotic number of origins")

        @load "fix_nest4.jld" #s=0.1, local
        scatter!(mig_vec, eta_avn_mtx2, markersize=4, label=legend2,
          markershape=:circle,yerror=eta_se_mtx2,color=colors)
          plot!(mig_vec, eta_avn_mtx2,color=colors,lab="",line=:dash)

      #power law (related to Ralph-Coop)
     a=float.(mig_vec)
      e =  10*a .^ (-1/4)
      plot!(mig_vec, e,scale=:log10,lab="",color=:black,label="power law",linewidth=2,
      alpha=0.5)
      # annotate!(1e-7,10^2.15,text("eta=A*(1/mig)^1/4",:purple,:right))

    mig_vec
    eta_avn_mtx

    ##Compare s in non-local
    nest1 and nest2
    @load
