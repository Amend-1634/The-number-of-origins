#call the functions in "3_4_figjoejyn_series.jl" to get the summary
 #matrix and plot

#mainly two series data: data at fixation and time-series data
 #for data at fixation you may need to estimate and provide a
  #slightly redundant T unlike time sampling data

#common parameter in all comparison unless specified
N=1e7
s0=0.1
twoNmu=1
M=10 #replicate

##figure4: timeline, compare dispersal rate series and non-spatial
include("2_4num_origin_gaussian.jl")
# include("3_1_1track_cons.jl")
include("3_1_2notrack_cons.jl") #to compare with Joejyn's codes
  #to define tracking population or oscillating population size
include("3_2_1d_tsamp.jl")
include("3_4_figjoejyn_series.jl")

mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #shorter time to fixation
N=1e7;nDemes=20;twoNmu=1;s0=0.1;M=20;T=400

@time tsamp, legend, avn_mtx_g, se_mtx_g, avn_mtx_m, se_mtx_m =
        onedim_dispersal_series(N,nDemes,twoNmu,s0,mig_vec,M,T)
        #tmux fig4 #823s

using JLD
save("fig4.jld", "tsamp", tsamp, "avn_mtx_g", avn_mtx_g, "se_mtx_g", se_mtx_g, "avn_mtx_m",
 avn_mtx_m, "se_mtx_m", se_mtx_m, "legend", legend) #not included legend

using JLD
cd("D:\\Julia\\func")
pwd()
readdir()
@load "fig4.jld"
legend
avn_mtx_g #20*7
tsamp
legend

#Gaussian
using Plots
 using Plots.PlotMeasures #for mm
 colors=[:cornflowerblue :red :orange :purple :green :steelblue :black]
 scatter(tsamp, avn_mtx_g, markersize=4, label=legend,
 color=colors,markershape=:rect, legend=:outertopright,
 legendfontsize=8, grid=false,fg_legend = :transparent,
 size=(600,600), yerror=se_mtx_g,margin=5mm,dpi=300)
 plot!(tsamp, avn_mtx_g,color=colors,lab="")
 xlabel!("Generation(t)")
 ylabel!("Average number of origins 𝛍(t) ")

savefig("..//photo//fig4_gaussian")


#Multinomial
using Plots.PlotMeasures #for mm
 scatter(tsamp, avn_mtx_m, markersize=4, label=legend,
 color=colors,markershape=:rect, legend=:outertopright,
 legendfontsize=8, grid=false,fg_legend = :transparent,
 size=(600,600), yerror=se_mtx_m,margin=5mm,dpi=300)
 plot!(tsamp, avn_mtx_m,color=colors,lab="")
 xlabel!("Generation(t)")
 ylabel!("Average number of origins 𝛍(t) ")

savefig("..//photo//fig4_multinomial")
#
# #gaussian
# using Plots.PlotMeasures #for mm
# for i in 1:length(legend)
#       if i==1
#           scatter(tsamp, avn_mtx_g[:,i], markersize=4, label=legend[i],
#           color="$(colors[i])",markershape=:rect, legend=:outertopright,
#            legendfontsize=8, grid=false,fg_legend = :transparent,
#            margin=5mm, size=(800,800), yerror=se_mtx_g[:,i])
#       else
#           scatter!(tsamp, avn_mtx_g[:,i], markersize=4, label=legend[i],
#           color="$(colors[i])",markershape=:rect, legend=:outertopright,
#            legendfontsize=8, grid=false,fg_legend = :transparent,
#             yerror=se_mtx_g[:,i])
#       end
#       plot!(tsamp, avn_mtx_g[:,i],color="$(colors[i])",lab="")
#  end
# xlabel!("Generation(t)")
# ylabel!("Average number of origins 𝛍(t) ")
#
# #multinomial
# using Plots.PlotMeasures #for mm
# for i in 1:length(legend)
#       if i==1
#           scatter(tsamp, avn_mtx_m[:,i], markersize=4, label=legend[i],
#           color="$(colors[i])",markershape=:rect, legend=:outertopright,
#            legendfontsize=8, grid=false,fg_legend = :transparent,
#            margin=5mm, size=(800,800), yerror=se_mtx_m[:,i])
#       else
#           scatter!(tsamp, avn_mtx_m[:,i], label=legend[i],
#           color="$(colors[i])",yerror=se_mtx_m[:,i])
#       end
#       plot!(tsamp, avn_mtx_m[:,i],color="$(colors[i])",lab="")
#  end
# xlabel!("Generation(t)")
# ylabel!("Average number of origins 𝛍(t) ")
#
#

savefig("..\\photo\\1d_local_fig4",dpi=300)
#d=1 and d=5e-4 not following the d higher num of origins lower rule

##figure6.A: timeline, compare nDemes series and non-spatial
# cd("D:\\Julia\\func")
include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
include("2_4num_origin_gaussian.jl") #non-spatial
include("3_4_figjoejyn_series.jl")

N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1; mig=0.005; M=replicates=10
T=3000; Gauss=1
nDemes_vec=[10,20,30,40]

# #test
# N=5e4; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1; mig=0.000005; M=replicates=3
# T=3000;Gauss=1
# nDemes_vec=[5,10]

@time tsamp, legend, avn_mtx, se_mtx = onedim_nDemes_series(nDemes_vec,N,twoNmu,
                            twoNmu_seq,s0,s0_seq, mig, M, replicates, T, Gauss)
                            #tmux fig6a  #353s
                            #T=3000 2hr

#through ssh(@Clifford)
using JLD
save("3_fig6a_3000.jld","tsamp", tsamp, "avn_mtx", avn_mtx, "se_mtx", se_mtx,
      "legend", legend)

using JLD
@load "figure6a.jld"
legend
colors=[:steelblue :red :yellow :purple :black]

#plotting
using Plots.PlotMeasures #for mm
scatter(tsamp, avn_mtx, markersize=4, label=legend,
 color=colors,markershape=:rect, legend=:outertopright,
 legendfontsize=8, grid=false,fg_legend = :transparent,
 size=(600,600), yerror=se_mtx,margin=5mm,dpi=300)
 plot!(tsamp, avn_mtx,color=colors,lab="")
 xlabel!("Generation(t)")
 ylabel!("Average number of origins 𝛍(t) ")
savefig("..//photo//fig6a")


##6.B. x dispersal, y time to fixation, series: nDemes
#test
# include("3_3_1d_fixsamp.jl")
# include("3_4_figjoejyn_series.jl")
# N=5e4; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1; mig=0.000005; M=replicates=3;T=3000;Gauss=1;
# nDemes_vec=[5,10];mig_vec=[5e-3,5e-4]

#normal size
include("3_4_figjoejyn_series.jl")
include("3_3_1d_fixsamp.jl")
N=1e7; twoNmu=twoNmu_seq=1; s0=s0_seq=0.1; mig=0.005; M=replicates=10;T=3000; Gauss=1
nDemes_vec=[10,20,30,40] #series: nDemes
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
#row: mig; col: nDemes

@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx, fix_M_mtx =
onedim_xdispersal_yfix_nDemesseries(mig_vec,nDemes_vec) #2019s

using JLD
save("fig6b.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
 "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)

using JLD
pwd()
cd("D:\\Julia\\func")
@load "fig6b.jld"

@load "fig8_largenDemes.jld"
@load "fig8_largenDemes2.jld"
legend

#plot
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
 using Plots.PlotMeasures #for mm
 using Plots
 colors=[:steelblue :red :yellow :purple]
 scatter(mig_vec, t_avn_mtx, markersize=4, label=legend,
 markershape=:rect, legend=:outertopright,legendfontsize=8,
 grid=false,fg_legend = :transparent,xscale=:log10,
 margin=8mm, size=(800,600), yerror=t_se_mtx,color=colors)
 plot!(mig_vec, t_avn_mtx,color=colors,lab="")
 xlabel!("dispersal rate")
 ylabel!("time to fixation (generation)")

savefig("..//photo//fig6b_t")


##8. higher nDemes higher symptotic η and longer time to fixation
 ##8.A
  # x: dispersal rate
  # y: symptotic η
  # series: nDemes, s

#difer from fig6B 1. y as eta_fix 2. series: nDemes+s
 # 1) s=0.05
  # #test
  # include("3_3_1d_fixsamp.jl")
  # include("3_4_figjoejyn_series.jl")
  # s0=s0_seq=0.05
  # nDemes_vec=[5,10];mig_vec=[5e-3,5e-4]
  # N=5e4; twoNmu=twoNmu_seq=1; mig=0.000005; M=replicates=3;T=3000;Gauss=1;

  #row: mig; col: nDemes
t_avn_mtx
eta_avn_mtx

  #normal size
  include("3_4_figjoejyn_series.jl")
  include("3_3_1d_fixsamp.jl")
  s0=s0_seq=0.05
  N=1e7; twoNmu=twoNmu_seq=1; mig=0.005; M=replicates=10;T=3000; Gauss=1
  nDemes_vec=[10,20,30,40] #series: nDemes
  mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x

  @time legend, t_avn_mtx1, t_se_mtx1, eta_avn_mtx1, eta_se_mtx1, fix_M_mtx1 =
                      onedim_xdispersal_yfix_nDemesseries(mig_vec,nDemes_vec)
                      #4953s
  #2)s=0.1

    # #test
    # include("3_4_figjoejyn_series.jl")
    # s0=s0_seq=0.1
    # nDemes_vec=[5,10];mig_vec=[5e-3,5e-4]
    # N=5e4; twoNmu=twoNmu_seq=1; mig=0.000005; M=replicates=3;T=3000;Gauss=1;

    #row: mig; col: nDemes
    t_avn_mtx
    eta_avn_mtx

    #normal size
    include("3_4_figjoejyn_series.jl")
    include("3_3_1d_fixsamp.jl")
    s0=s0_seq=0.1
    N=1e7; twoNmu=1; M=10;T=3000; Gauss=1
    nDemes_vec=[10,20,30,40] #series: nDemes
    mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x


    @time legend, t_avn_mtx2, t_se_mtx2, eta_avn_mtx2, eta_se_mtx2, fix_M_mtx2 =
                   onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
                  #1741s

using JLD
save("fig8.jld","legend",legend,"t_avn_mtx1",t_avn_mtx1,"t_se_mtx1", t_se_mtx1,
"eta_avn_mtx1",eta_avn_mtx1,"eta_se_mtx1",eta_se_mtx1,"fix_M_mtx1", fix_M_mtx1,
"t_avn_mtx2",t_avn_mtx2,"t_se_mtx2", t_se_mtx2, "eta_avn_mtx2",eta_avn_mtx2,
"eta_se_mtx2",eta_se_mtx2,"fix_M_mtx2", fix_M_mtx2)


save("t-1.jld", "t_avn_mtx2",t_avn_mtx2,"t_se_mtx2", t_se_mtx2, "eta_avn_mtx2",eta_avn_mtx2,
"eta_se_mtx2",eta_se_mtx2,"fix_M_mtx2", fix_M_mtx2)

using JLD
        pwd()
        cd("D:\\Julia\\func")
        @load "fig8.jld"
        legend

#8.A############################################y as eta_fix
@load "fig8_largenDemes2.jld"
  using Plots.PlotMeasures #for mm
  using Plots
  mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
  colors=[:steelblue :red :orange :purple]
     #s=0.05
        scatter(mig_vec, eta_avn_mtx1, markersize=4, label=legend,
        markershape=:rect, legend=:outertopright,xscale=:log10,
        legendfontsize=8, grid=false,fg_legend = :transparent,
        margin=5mm, size=(800,800), yerror=eta_se_mtx1, color=colors,
        series_annotations = text.(fix_M_mtx, :bottom))
        plot!(mig_vec, eta_avn_mtx1,lab="", color=colors,linewidth=3,alpha=0.6)

      #s=0.1
        scatter!(mig_vec, eta_avn_mtx2, markersize=4, label=legend,
        color=colors,markershape=:circle,yerror=eta_se_mtx2)
        plot!(mig_vec, eta_avn_mtx2, lab="",line=:dash, color=colors,linewidth=3,alpha=0.6)
   xlabel!("Dispersal rate")
   ylabel!("Asymptotic number of origins")

savefig("..//photo//fig8b_eta")


#8.B######################################## y as t_fix
using JLD
@load "fig8.jld"
  legend1=["s=0.05, nDemes = 10"  "s=0.05, nDemes = 20"  "s=0.05, nDemes = 30"  "s=0.05, nDemes = 40"]
  legend2=["s=0.1, nDemes = 10"  "s=0.1, nDemes = 20"  "s=0.1, nDemes = 30"  "s=0.1, nDemes = 40"]
  using Plots.PlotMeasures #for mm
  using Plots
  mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
  colors=[:steelblue :red :orange :purple]
     #s=0.05
        scatter(mig_vec, t_avn_mtx1, markersize=2, label=legend1,
        markershape=:rect, legend=:topright,xscale=:log10,
        legendfontsize=12, grid=false,fg_legend = :transparent,
        margin=10mm, size=(600,800), yerror=eta_se_mtx1, color=colors)
        plot!(mig_vec, t_avn_mtx1,lab="", color=colors,linewidth=3,alpha=0.6)

      #s=0.1
        scatter!(mig_vec, t_avn_mtx2, markersize=2, label=legend2,guidefont=font(20),
        color=colors,markershape=:circle,yerror=eta_se_mtx2)
        plot!(mig_vec, t_avn_mtx2, lab="",line=:dash, color=colors,linewidth=3,alpha=0.6)
   xlabel!("Dispersal rate")
   ylabel!("Time to fixation (generations)")

@load "fig8.jld"
     legend1=["s=0.05, nDemes = 10"  "s=0.05, nDemes = 20"  "s=0.05, nDemes = 30"  "s=0.05, nDemes = 40"]
     legend2=["s=0.1, nDemes = 10"  "s=0.1, nDemes = 20"  "s=0.1, nDemes = 30"  "s=0.1, nDemes = 40"]
     using Plots.PlotMeasures #for mm
     using Plots
     mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
     colors=[:steelblue :red :orange :purple]
        #s=0.05
           scatter(mig_vec, eta_avn_mtx1, markersize=2, label="",
           markershape=:rect, legend=:topright,xscale=:log10,
           legendfontsize=12, grid=false,fg_legend = :transparent,
           margin=10mm, size=(600,800), yerror=eta_se_mtx1, color=colors)
           plot!(mig_vec, eta_avn_mtx1,lab="", color=colors,linewidth=3,alpha=0.6)

         #s=0.1
           scatter!(mig_vec, eta_avn_mtx2, markersize=2, label="",guidefont=font(20),
           color=colors,markershape=:circle,yerror=eta_se_mtx2)
           plot!(mig_vec, eta_avn_mtx2, lab="",line=:dash, color=colors,linewidth=3,alpha=0.6)
      xlabel!("Dispersal rate")
      ylabel!("Asymptotic number of origins")



savefig("..//photo//fig8b_t")


#how to notate the series s?

##Fig10.Oscillating population:

#test
# include("3_3_1d_fixsamp_osci.jl")
# include("3_4_figjoejyn_series.jl")
# s0=0.1;nDemes=5;twoNmu=1; mig=0.000005; M=3;T=3000;Gauss=1;
#   mig_vec=[5e-8 5e-7] #x
#   phi_vec=[10 50 100] #peak-to-trough ratio
#   N=1e7; phi=10; har_geo="H" #constrained harmonic mean of populaiton size
#   N=1e7; phi=10; har_geo="G" #constrained geometric mean

#Joejyn's setting
include("3_3_1d_fixsamp.jl")
include("3_3_1d_fixsamp_osci.jl")
include("3_4_figjoejyn_series.jl")
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
phi_vec=[10 50 100] #peak-to-trough ratio
N=1e7; phi=10; har_geo="H" #constrained harmonic mean of populaiton size #@osci_H
N=1e7; phi=10; har_geo="G" #constrained geometric mean #@osci_G
s0=0.1;nDemes=5;twoNmu=1; mig=0.000005; M=10;T=3000;Gauss=1;


@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx=
  dispersal_oscillating(N, phi, har_geo,nDemes,twoNmu,s0,M,T,Gauss,mig_vec,phi_vec) #1381s

using JLD
save("fig10.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
    "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx)
using JLD
@load "fig10.jld"

using Plots.PlotMeasures #for mm
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
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




#x-column vector
#legend,color--all row vector


##Large nDemes
# cd("D:\\Julia\\func")
include("3_1_2notrack_cons.jl")
include("3_3_1d_fixsamp.jl")
include("3_4_figjoejyn_series.jl")
s0=s0_seq=0.1
N=1e7; twoNmu=1; M=10;T=3000;Gauss=0 # Gauss=1
nDemes_vec=[10,50, 100, 500] #series: nDemes
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x

@time legend, t_avn_mtx2, t_se_mtx2, eta_avn_mtx2, eta_se_mtx2, fix_M_mtx2 =
      onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
using JLD
save("fig8_largenDemes2.jld","legend",legend,"t_avn_mtx2",t_avn_mtx2,
  "t_se_mtx2", t_se_mtx2, "eta_avn_mtx2",eta_avn_mtx2,
  "eta_se_mtx2",eta_se_mtx2,"fix_M_mtx2", fix_M_mtx2)



include("3_1_2notrack_cons.jl")
include("3_3_1d_fixsamp.jl")
include("3_4_figjoejyn_series.jl")
s0=s0_seq=0.1
N=1e7; twoNmu=1; M=10;T=3000;Gauss=0 # Gauss=1
nDemes_vec=[1000] #to plug in  nDemes_vec=[10,50, 100, 500]
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x

@time legend, t_avn_mtx2, t_se_mtx2, eta_avn_mtx2, eta_se_mtx2, fix_M_mtx2 =
      onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
using JLD
save("3_fig8_1000nDemes.jld","legend",legend,"t_avn_mtx2",t_avn_mtx2,
  "t_se_mtx2", t_se_mtx2, "eta_avn_mtx2",eta_avn_mtx2,
  "eta_se_mtx2",eta_se_mtx2,"fix_M_mtx2", fix_M_mtx2)


  #Runnig time:
    # mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
        #nDemes_vec=[10,50, 100]  #7337s = 2hr
        #nDemes_vec=[10,50, 100, 500]  #47957s, 13 hr
        #nDemes_vec=[10,50, 100, 1000]#no, overload

   # mig_vec=[5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3] #
       #nDemes_vec=[10,50, 100]  #
       #nDemes_vec=[10,50, 100, 500]  #

using JLD
pwd()
cd("D:\\Julia\\func")
@load "fig8_largenDemes.jld" # mig_vec=[5e-5,..];  nDemes_vec=[10,50, 100]
@load "fig8_largenDemes2.jld" # mig_vec=[5e-5,..];  nDemes_vec=[10,50, 100, 500]


using Plots.PlotMeasures #for mm
  # legend=["nDemes = 10"  "nDemes = 50"  "nDemes = 100"  "nDemes = 500"]
  using Plots
  mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
  colors=[:steelblue :red :orange :green]
    #Asymptotic number of origins
    p1=Plots.scatter(mig_vec, eta_avn_mtx2, markersize=2, label="",
    markershape=:rect, legend=:topright,xscale=:log10,xlims=(1e-5,10),
    legendfontsize=12, grid=false,fg_legend = :transparent,
    guidefont  = font("Times new roman", 18), #for axis
    tickfont  = font("Times new roman", 10), #for axis
    titlefont  = font("Times new roman", 18), #for axis
    margin=10mm, size=(500,600), yerror=eta_se_mtx2, color=colors,
    dpi=300)
    # for i in 1:length(mig_vec), j in 1:length(colors)
    #   annotate!([(mig_vec[i], eta_avn_mtx2[i,j], text(string(fix_M_mtx2[i,j]),:bottom,12))]) #all fixed
    # end #unlike joejyn's, all fixed
    Plots.plot!(mig_vec, eta_avn_mtx2,lab="", color=colors,)
    xlabel!("Migration rate (mig)")
    ylabel!("Asymptotic number of origins")
    title!("(a) x in log10 scale")

p2=Plots.scatter(mig_vec, eta_avn_mtx2, markersize=2, label=legend,
    markershape=:rect, legend=:topright,xscale=:log10,xlims=(1e-5,10),
    legendfontsize=11, grid=false,fg_legend = :transparent,
    guidefont  = font("Times new roman", 18), #for axis
    tickfont  = font("Times new roman", 10), #for axis
    titlefont  = font("Times new roman", 18), #for axis
    margin=10mm, size=(500,600), yerror=eta_se_mtx2, color=colors,
    dpi=300)
    # for i in 1:length(mig_vec), j in 1:length(colors)
    #   annotate!([(mig_vec[i], eta_avn_mtx2[i,j], text(string(fix_M_mtx2[i,j]),:bottom,12))]) #all fixed
    # end #unlike joejyn's, all fixed
    Plots.plot!(mig_vec, eta_avn_mtx2,lab="", color=colors)
    xlabel!("Migration rate (mig)")
    ylabel!("Asymptotic number of origins")
    title!("(b) x and y in log10 scale")


#power law (related to Ralph-Coop)
 a=float.(-4:0)
  b=10 .^ a #as x, which is d in the equation
  b
  # c = (1/b) .^ (1/4)#1/b not functional here
  e = 10^1.3 * b .^ (-1/4)
  # size(c)
  # plot(b, c,scale=:log10)
  Plots.plot!(b, e,scale=:log10,lab="",color=:black)
  Plots.plot(p1,p2,layout=(1,2),size=(1400,600))

  # savefig("..\\photo\\local_largenDemes500_eta")
  # annotate!(0.1,10^2,text("eta=A*(1/mig)^1/4",:black))


        #time to fixation
scatter(mig_vec, t_avn_mtx2, markersize=4, label=legend,
    markershape=:rect, legend=:outertopright,xscale=:log10,
    legendfontsize=8, grid=false,fg_legend = :transparent,
    margin=5mm, size=(600,600), yerror=t_se_mtx2, color=colors,
    dpi=300)
    for i in 1:length(mig_vec), j in 1:length(colors)
      annotate!([(mig_vec[i], t_avn_mtx2[i,j], text(string(fix_M_mtx2[i,j]),:bottom,12))]) #all fixed
    end
    plot!(mig_vec, t_avn_mtx2,lab="", color=colors,linewidth=3,alpha=0.6)
    xlabel!("Dispersal rate")
    ylabel!("Time to fixation")
savefig("..\\photo\\largenDemes500_t")


##Is t-1 -> t in Joejyn's code and t->t in my code cause my mean η half of
#Joejyn's?
#not much difference
pwd()
cd("D:\\Julia\\func")
include("3_1_2notrack_cons_t-1.jl") #from t-1
# include("3_1_2notrack_cons.jl")
include("3_3_1d_fixsamp.jl")
include("3_4_figjoejyn_series.jl")
s0=s0_seq=0.1
N=1e7; twoNmu=1; M=10;T=3000; Gauss=1
nDemes_vec=[10,20,30,40] #series: nDemes
mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x

@time legend, t_avn_mtx2, t_se_mtx2, eta_avn_mtx2, eta_se_mtx2, fix_M_mtx2 =
      onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
              #nDemes_vec=[10,20,30,40]
              #3640s

using JLD
save("t-1.jld","legend",legend, "t_avn_mtx2", t_avn_mtx2, "t_se_mtx2", t_se_mtx2,
 "eta_avn_mtx2", eta_avn_mtx2, "eta_se_mtx2", eta_se_mtx2, "fix_M_mtx2", fix_M_mtx2)


using JLD
pwd()
cd("D:\\Julia\\func")
eta_avn_mtx2

using Plots.PlotMeasures #for mm
  using Plots
  mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
  colors=[:steelblue :red :orange :purple]
    #Asymptotic number of origins
    #t-1->t
    @load "t-1.jld"
    scatter(mig_vec, eta_avn_mtx2, markersize=4, label=legend,
    markershape=:rect, legend=:outertopright,xscale=:log10,
    legendfontsize=8, grid=false,fg_legend = :transparent,
    margin=5mm, size=(600,600), yerror=eta_se_mtx2, color=colors,
    dpi=300)
    plot!(mig_vec, eta_avn_mtx2,lab="", color=colors,linewidth=3,alpha=0.6)

    #t->t  #dash line and circle
    @load "fig8.jld"
    scatter!(mig_vec, eta_avn_mtx2, markersize=4, markershape=:circle,
     yerror=eta_se_mtx2, color=colors, lab="")
    plot!(mig_vec, eta_avn_mtx2,lab="", color=colors,linewidth=3,alpha=0.6,line=:dash)

    xlabel!("Dispersal rate")
    ylabel!("Asymptotic number of origins")
    #for notating fixation probability
    # for i in 1:length(mig_vec), j in 1:length(colors)
    #   annotate!([(mig_vec[i], eta_avn_mtx2[i,j], text(string(fix_M_mtx2[i,j]),:bottom,12))]) #all fixed
    # end
savefig("..\\photo\\t-1_tcircle_eta")



#time to fixation
using Plots.PlotMeasures #for mm
  using Plots
  mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1] #x
  colors=[:steelblue :red :orange :purple]
    #Asymptotic number of origins
    #t-1->t
    @load "t-1.jld"
    scatter(mig_vec, t_avn_mtx2, markersize=4, label=legend,
    markershape=:rect, legend=:outertopright,xscale=:log10,
    legendfontsize=8, grid=false,fg_legend = :transparent,
    margin=5mm, size=(600,600), yerror=t_se_mtx2, color=colors,
    dpi=300)
    plot!(mig_vec, t_avn_mtx2,lab="", color=colors,linewidth=3,alpha=0.6)

    #t->t  #dash line and circle
    @load "fig8.jld"
    scatter!(mig_vec, t_avn_mtx2, markersize=4, markershape=:circle,
     yerror=t_se_mtx2, color=colors, lab="")
    plot!(mig_vec, t_avn_mtx2,lab="", color=colors,linewidth=3,alpha=0.6,line=:dash)

    xlabel!("Dispersal rate")
    ylabel!("Time to fixation (generation)")
    #for notating fixation probability
    # for i in 1:length(mig_vec), j in 1:length(colors)
    #   annotate!([(mig_vec[i], t_avn_mtx2[i,j], text(string(fix_M_mtx2[i,j]),:bottom,12))]) #all fixed
    # end
savefig("..\\photo\\t-1_tcircle_eta")

##track migrants and not tracking just scaling
#fig4: time-series

include("2_4num_origin_gaussian.jl")
include("3_2_1d_tsamp.jl")
include("3_4_figjoejyn_series.jl")

mig_vec=[5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 1]
N=1e7;nDemes=20;twoNmu=1;s0=0.1;M=20;T=3000

include("3_1_1track_cons.jl") #to compare with not tracking #769s
# include("3_1_2notrack_cons.jl") #to compare with Joejyn's codes

@time tsamp, legend, avn_mtx_g, se_mtx_g, avn_mtx_m, se_mtx_m =
        onedim_dispersal_series(N,nDemes,twoNmu,s0,mig_vec,M,T)
        #tmux fig4 #823s

using JLD
save("fig4_tracking.jld", "tsamp", tsamp, "avn_mtx_g", avn_mtx_g, "se_mtx_g", se_mtx_g, "avn_mtx_m",
 avn_mtx_m, "se_mtx_m", se_mtx_m, "legend", legend) #not included legend

using JLD
cd("D:\\Julia\\func")
pwd()
readdir()
@load "fig4_tracking.jld"
legend=["d = 5e-5" "d = 5e-4" "d = 5e-3" "d = 5e-2" "d = 5e-1" "d = 1" "Non-spatial"]
avn_mtx_g #20*7
tsamp
legend

#Gaussian
using Plots
colors=[:cornflowerblue :red :orange :purple :green :steelblue :black]

#Gaussian
using Plots.PlotMeasures #for mm
scatter(tsamp, avn_mtx_g, markersize=4, label=legend,
 color=colors,markershape=:rect, legend=:outertopright,
 legendfontsize=8, grid=false,fg_legend = :transparent,
 size=(600,600), yerror=se_mtx_g,margin=5mm,dpi=300)
 plot!(tsamp, avn_mtx_g,color=colors,lab="")
 xlabel!("Generation(t)")
 ylabel!("Average number of origins 𝛍(t) ")
savefig("..//photo//fig4_gaussian")


##tsamp, ave_origins #@3_s0_twoNmu_N
#local--as local not significantly differ from non-local (to save time)
N_vec=[1e7 1e9]

N=1e7
s0_vec=[0.01 0.1]
twoNmu_vec=[1 10]
nDemes=100
mig=0.005
Gauss=1
M=10
T=1000 #should already reach fixation
include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
include("3_4_figjoejyn_series_threads.jl")
@time tsamp, legend, avn_mtx, se_mtx = s0_twoNmu_N_series(nDemes,N,twoNmu_vec,s0_vec, mig, M, T, Gauss)
using JLD
save("3_s0_twoNmu_N1.jld","tsamp", tsamp, "legend", legend, "avn_mtx", avn_mtx, "se_mtx", se_mtx)


N=1e9
s0_vec=[0.01 0.1]
twoNmu_vec=[1 10]
nDemes=100
mig=0.005
Gauss=1
M=10
T=1000 #should already reach fixation
include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
include("3_4_figjoejyn_series_threads.jl")
@time tsamp, legend, avn_mtx, se_mtx = s0_twoNmu_N_series(nDemes,N,twoNmu_vec,s0_vec, mig, M, T, Gauss)
using JLD
save("3_s0_twoNmu_N2.jld","tsamp", tsamp, "legend", legend, "avn_mtx", avn_mtx, "se_mtx", se_mtx)

cd("D:\\Julia\\func")

@load "3_s0_twoNmu_N1.jld"
legend

Plots.GRBackend()
using JLD
  using Plots
  colors=[:yellow :orange :lightgreen :green]
  @load "3_s0_twoNmu_N1.jld" #N=1e7
 using Plots.PlotMeasures #for mm
 scatter(tsamp, avn_mtx, markersize=2, label=legend,
 color=colors,markershape=:rect, legend=:legend,
 guidefont  = font("Times new roman", 24), #for axis
 legendfontsize=16, grid=false,fg_legend = :transparent,
 size=(700,600), yerror=se_mtx,margin=5mm,dpi=300)
 plot!(tsamp, avn_mtx,color=colors,lab="")
 ylabel!("Average number of origins")

 @load "3_s0_twoNmu_N2.jld"#N=1e9
 using Plots.PlotMeasures #for mm
 legend=["N=1e9, s=0.01, twoNmu=0.1" "N=1e9, s=0.01, twoNmu=1.0" "N=1e9, s=0.1,   twoNmu=0.1" "N=1e9, s=0.1,   twoNmu=1.0"]
 scatter!(tsamp, avn_mtx, markersize=2, label=legend,
 color=colors,markershape=:circle, legend=:outerright,
 legendfontsize=24, grid=false,fg_legend = :transparent,
 guidefont  = font("Times new roman", 24),
 tickfont  = font("Times new roman", 18),
 size=(1400,900), yerror=se_mtx,margin=16mm,dpi=300)
 plot!(tsamp, avn_mtx,color=colors,lab="",line=:dash)
 xlabel!("Generation (t)")

savefig("..//photo//s_twoNmu_N_1000generation")

#whether the two series are equal at fixation:
using GLM
ols = lm(@formula(Y ~ X), data)

cd("D:\\Julia\\func")
@load "3_s0_twoNmu_N1.jld" #N=1e7
x=avn_mtx[10:end]
@load "3_s0_twoNmu_N2.jld"#N=1e9
y=avn_mtx[10:end]




##
N=1e7
T=3000
s0_vec=[1e-2, 5e-2, 0.1, 0.5]
mig_vec=[5e-8, 5e-4]
nDemes_vec=[100,200]
twoNmu=10
M=10;Gauss=1
include("3_4_figjoejyn_series.jl")
include("3_3_1d_fixsamp.jl")
include("3_s0_mig_nDemes.jl")
include("3_1_2notrack_cons.jl")

@time legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx=
    onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0_vec,mig_vec,M,T,Gauss)
using JLD
save("s0_mig_nDemes.jld","legend",legend, "t_avn_mtx", t_avn_mtx, "t_se_mtx", t_se_mtx,
 "eta_avn_mtx", eta_avn_mtx, "eta_se_mtx", eta_se_mtx, "fix_M_mtx", fix_M_mtx)


##
#Plot total frequency, average number of origins in the same plot

s0_vec=[0.1 0.05]
mig=5e-4
nDemes=100
N=1e7; M=10; Gauss=1; twoNmu=10
T=2000 #should be 600

# cd("D:\\Julia\\func")

include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
include("3_X_eta.jl")

@time tsamp, legend, eta_avn_mtx, eta_se_mtx, X_avn_mtx, X_se_mtx =
                        eta_X_tsamp(N,nDemes,twoNmu,s0_vec,mig,M,T)

using JLD
save("3_eta_X.jld", "tsamp", tsamp, "legend",legend,  "eta_avn_mtx", eta_avn_mtx,
        "eta_se_mtx", eta_se_mtx,"X_avn_mtx", X_avn_mtx, "X_se_mtx", X_se_mtx)

#@compare_s

using JLD
cd("D:\\Julia\\func")
readdir()
@load "3_eta_X.jld"
 using Plots.PlotMeasures #for mm
 using Plots
 colors=[:green :lightgreen]
 scatter(tsamp,eta_avn_mtx, markersize=3, label=legend,
 markershape=:rect, legend=:right,legendfontsize=14,dpi=300,
 grid=false,fg_legend = :transparent,
 guidefont  = font("Times new roman", 18), #for axis
 legendfont = font("Times new roman", 18),
 tickfont = font("Times new roman", 14),
 margin=8mm, size=(800,600), yerror=eta_se_mtx,color=colors)
 plot!(tsamp,eta_avn_mtx, color=colors,lab="")
 xlabel!("Generation (t)")
 ylabel!("Average number of origins")


#slight difference between s series---
tsamp
