pwd()
cd("D:\\Julia\\func\\twoNDememig")

nDemes=10;N=1e7;
  NDeme=N/nDemes
  twoNDememig_vec=[1e-4, 1e-3, 1e-2, 1e-1, 1]
  mig_vec=twoNDememig_vec/NDeme

using JLD
# @load "nu_fix_10_2.jld" #M=50
  @load "nu_fix_10.jld"
 #time
 using Plots.PlotMeasures #for mm
 using Plots
 colors=[:steelblue :red :orange]
 p1=scatter(mig_vec, t_avn_mtx, markersize=4, label="",
   markershape=:rect, legend=:topright,legendfontsize=6,dpi=300,
   grid=false,fg_legend = :transparent,xscale=:log10,
   margin=4mm, size=(1200,600), yerror=t_se_mtx,color=colors)
 plot!(mig_vec, t_avn_mtx,color=colors,lab="")
 xlabel!("Migration rate (mig)")
 ylabel!("Time to fixation (generation)")
 title!("(a) nDemes = 10")

 # @load "nu_fix_50.jld"
 # scatter!(mig_vec, t_avn_mtx, markersize=4, label=legend,
 # markershape=:circle, legend=:topright, yerror=t_se_mtx,color=colors)
 # plot!(mig_vec, t_avn_mtx,color=colors,lab="",line=:dash)

 using JLD
 #time
 using Plots.PlotMeasures #for mm
 using Plots
 colors=[:steelblue :red :orange]
 @load "nu_fix_50.jld"
 p2=scatter(mig_vec, t_avn_mtx, markersize=4, label=legend,
   markershape=:rect, legend=:topright,dpi=300,
   grid=false,fg_legend = :transparent,xscale=:log10,
   yerror=t_se_mtx,color=colors)
 plot!(mig_vec, t_avn_mtx,color=colors,lab="")
 xlabel!("Migration rate (mig)")
 ylabel!("Time to fixation (generation)")
 title!("(b) nDemes = 50")

 plot(p1,p2,layout=(1,2))
savefig("..//..//photo//nu_fix_t_10_50")




cd("D:\\Julia\\func")

#eta
@load "nu_fix_10.jld"
  using Plots.PlotMeasures #for mm
  colors=[:steelblue :red :orange]
  p1=scatter(mig_vec, eta_avn_mtx, markersize=3, label="",
  guidefont  = font("Times new roman", 16), #for axis
  markershape=:rect,
  grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
  titlefont  = font("Times new roman", 20),
  margin=10mm, yerror=eta_se_mtx,color=colors)
  plot!(mig_vec, eta_avn_mtx,color=colors,lab="")
  xlabel!("Migration rate (mig)")
  ylabel!("Asymptotic number of origins")
  title!("(a) nDemes = 10")

  @load "nu_fix_50.jld"
  using Plots.PlotMeasures #for mm
    # mig_vec=[5e-10 5e-8 5e-6 5e-4 5e-2] #x
    colors=[:steelblue :red :orange]
    p2=scatter(mig_vec, eta_avn_mtx, markersize=3, label=legend,
    guidefont  = font("Times new roman", 16),
    titlefont  = font("Times new roman", 20),
    markershape=:rect, size=(1000,600),legend=:bottomleft,legendfontsize=16,
    grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
    margin=6mm,yerror=eta_se_mtx,color=colors)
    plot!(mig_vec, eta_avn_mtx,color=colors,lab="")
    xlabel!("Migration rate (mig)")
    ylabel!("Asymptotic number of origins")
    title!("(b) nDemes = 50")


    plot(p1,p2,layout=(1,2))

savefig("..//..//photo//nu_fix_t_10_50")

##
#eta
cd("D:\\Julia\\func")


  #
  # b=mig_vec[2:4]
  # a=float.(b)
  # e =  a .^ (-1/4)
  # plot!(b, e,xscale=:log10,lab="",color=:black,label="power law",linewidth=2,
  # alpha=0.5)

cd("D:\\Julia\\func")

  @load "nu_fix_50_2.jld"
      mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
      using Plots.PlotMeasures #for mm
      colors=[:steelblue :red :orange]
      p2=scatter(mig_vec, t_avn_mtx, markersize=3, label=legend,
      markershape=:circle,size=(500,600),
      legendfontsize=12,
      guidefont  = font("Times new roman", 16), #for axis
      titlefont  = font("Times new roman", 18),
      tickfont  = font("Times new roman", 13),
      grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
      margin=10mm, yerror=t_se_mtx,color=colors)
      plot!(mig_vec, t_avn_mtx,color=colors,lab="",line=:dash,linewidth=2)
      xlabel!("Migration rate (mig)")
      ylabel!("Time to fixation")
      title!("(b) Asymptotic number of origins")


  @load "nu_fix_50_2.jld"
    mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
    using Plots.PlotMeasures #for mm
    colors=[:steelblue :red :orange]
    p1=scatter(mig_vec, eta_avn_mtx, markersize=3, label=legend,
    markershape=:circle,size=(500,600),
    grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
    legendfontsize=12,
    guidefont  = font("Times new roman", 16), #for axis
    titlefont  = font("Times new roman", 18),
    tickfont  = font("Times new roman", 13),
    margin=10mm, yerror=eta_se_mtx,color=colors)
    plot!(mig_vec, eta_avn_mtx,color=colors,lab="",line=:dash)
    xlabel!("Migration rate (mig)")
    ylabel!("Asymptotic number of origins")
    title!("(a) Time to fixation")

    # plot!(mig_vec, e,xscale=:log10,lab="",color=:black,label="power law",linewidth=2,alpha=0.5)


    plot(p1,p2,layout=(1,2),size=(1000,600))
savefig("D://Julia//photo//compare_nu_50")



##t; M=50
legend

cd("D:\\Julia\\func")
  using JLD
  using Plots
  @load "nu_fix_10_2.jld"
  mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
  using Plots.PlotMeasures #for mm
  colors=[:steelblue :red :orange]
  p2=scatter(mig_vec, eta_avn_mtx, markersize=3, label=legend,
  guidefont  = font("Times new roman", 16), #for axis
  markershape=:rect,size=(500,600),legendfontsize=12,
  grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
  titlefont  = font("Times new roman", 18),
  tickfont  = font("Times new roman", 13),

  margin=10mm, yerror=eta_se_mtx,color=colors)
  plot!(mig_vec, eta_avn_mtx,color=colors,lab="")
  xlabel!("Migration rate (mig)")
  ylabel!("Asymptotic number of origins")
  title!("(b) Asymptotic number of origins")

 @load "nu_fix_10_2.jld"
    mig_vec=[5e-10, 5e-8, 5e-6, 5e-4, 5e-2] #x
    using Plots.PlotMeasures #for mm
    colors=[:steelblue :red :orange]
    p1=scatter(mig_vec, t_avn_mtx, markersize=3, label=legend,
    guidefont  = font("Times new roman", 16), #for axis
    grid=false,fg_legend = :transparent,xscale=:log10,dpi=300,
    markershape=:rect,size=(500,600),
    legendfontsize=12,
    titlefont  = font("Times new roman", 18),
    tickfont  = font("Times new roman", 13),
    margin=10mm, yerror=t_se_mtx,color=colors)
    plot!(mig_vec, t_avn_mtx,color=colors,lab="",linewidth=2)
    xlabel!("Migration rate (mig)")
    ylabel!("Time to fixation")
    title!("(a) Time to fixation")

    plot(p1,p2,layout=(1,2),size=(1000,600))
    savefig("D://Julia//photo//compare_nu_10")
