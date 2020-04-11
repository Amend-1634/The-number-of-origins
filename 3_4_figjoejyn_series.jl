#the functions in this file is for collecting data in matrix for
 #comparison in parameters and will be plotted in *_test.jl later
#the interested parameter (y or series in legend) is inputed in the form of
  #*_vec
#all figure number are corresponding to the figure in Joejyn Wan's figure

using Base.Threads
function onedim_series_mig_nDemes(N,nDemes_vec,twoNmu,s0,mig_vec,M,T) #,Gauss
 len_mig=length(mig_vec)+1 #last one for non-spatial
 len_nDemes=length(nDemes_vec)
 legend=String[]
 len_tsamp=20 #20 generations
 avn_mtx_g=zeros(len_tsamp,0)
 se_mtx_g=zeros(len_tsamp,0)

 avn_mtx_m=zeros(len_tsamp,0)
 se_mtx_m=zeros(len_tsamp,0)
 for j in 1:len_nDemes
     for i in 1:(len_mig-1)
        avn_g, se_g = num_origin_1d_t(N,nDemes_vec[j],twoNmu,s0,mig_vec[i],M,T,1)[1:2]

        avn_m, se_m= num_origin_1d_t(N,nDemes_vec[j],twoNmu,s0,mig_vec[i],M,T,0)[1:2]
        #                  # num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
        #*_g in gaussian sampling while *_m in multinomial sampling
        avn_mtx_g=hcat(avn_mtx_g,avn_g)
        se_mtx_g=hcat(se_mtx_g,se_g)
        avn_mtx_m=hcat(avn_mtx_m,avn_m)
        se_mtx_m=hcat(se_mtx_m,se_m)
        push!(legend,"mig=$(mig_vec[i]), nDemes=$(nDemes_vec[j])")
    end
 end

 #non-spatial
 simu_mtx,simu_mtx_se=num_origin_gaussian(s0,   twoNmu,N,M,T)[2:3]
                            #num_origin_gaussian(s0_seq,twoNmu_seq,N,replicates,T)
 push!(legend,"Non-spatial")

 len_tsamp=20 #20 generations
 tsamp=range(1,T+1,length=len_tsamp)
 tsamp=round.(Int,tsamp)

 avn_mtx_g=hcat(avn_mtx_g, simu_mtx[tsamp,:])
 se_mtx_g=hcat(se_mtx_g, simu_mtx_se[tsamp,:])
 avn_mtx_m=hcat(avn_mtx_m, simu_mtx[tsamp,:])
 se_mtx_m=hcat(se_mtx_m, simu_mtx_se[tsamp,:])

 return tsamp, legend, avn_mtx_g, se_mtx_g, avn_mtx_m, se_mtx_m
 end






#figure4 in Joejyn Wan's papaer
function onedim_dispersal_series(N,nDemes,twoNmu,s0,mig_vec,M,T) #,Gauss
 len_mig=length(mig_vec)+1 #last one for non-spatial
 legend=Array{String,2}(undef,1,len_mig)
 len_tsamp=20 #20 generations
 avn_mtx_g=zeros(len_tsamp,0)
 se_mtx_g=zeros(len_tsamp,0)

 avn_mtx_m=zeros(len_tsamp,0)
 se_mtx_m=zeros(len_tsamp,0)
 for i in 1:(len_mig-1)
    # avn, se, tsamp = num_origin_1d_t(1e5,5,1,0.1,mig_vec[i],10,2000,1)[1:3]
    # @time avn_g, se_g, tsamp = num_origin_1d_t(5e5,5,1,0.1,mig_vec[i],3,200,1)[1:3]
    avn_g, se_g = num_origin_1d_t(N,nDemes,twoNmu,s0,mig_vec[i],M,T,1)[1:2]

    # # @time avn_m, se_m= num_origin_1d_t(5e5,5,1,0.1,mig_vec[i],3,200,0)[1:2]
    avn_m, se_m= num_origin_1d_t(N,nDemes,twoNmu,s0,mig_vec[i],M,T,0)[1:2]
    #                  # num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
    #_g for gaussian while _m for multinomial sampling
    avn_mtx_g=hcat(avn_mtx_g,avn_g)
    se_mtx_g=hcat(se_mtx_g,se_g)
    avn_mtx_m=hcat(avn_mtx_m,avn_m)
    se_mtx_m=hcat(se_mtx_m,se_m)
    legend[1,i]="d = $(mig_vec[i])" #dispersal rate
 end

 #non-spatial
 simu_mtx,simu_mtx_se=num_origin_gaussian(0.1,    1,         1e7,20,3000)[2:3]
                            #num_origin_gaussian(s0_seq,twoNmu_seq,N,replicates,T)

 len_tsamp=20 #20 generations
 tsamp=range(1,T+1,length=len_tsamp)
 tsamp=round.(Int,tsamp)

 avn_mtx_g=hcat(avn_mtx_g, simu_mtx[tsamp,:])
 se_mtx_g=hcat(se_mtx_g, simu_mtx_se[tsamp,:])
 avn_mtx_m=hcat(avn_mtx_m, simu_mtx[tsamp,:])
 se_mtx_m=hcat(se_mtx_m, simu_mtx_se[tsamp,:])

 legend[1,end]="Non-spatial"

 return tsamp, legend, avn_mtx_g, se_mtx_g, avn_mtx_m, se_mtx_m
 end



##figure6.A:timeline, compare nDemes series and non-spatial

#figure6a
function onedim_nDemes_series(nDemes_vec,N,twoNmu,twoNmu_seq,s0,s0_seq, mig, M, replicates, T, Gauss)
    len_nDemes=length(nDemes_vec)
    len_tsamp=20 #20 generations
     avn_mtx=zeros(len_tsamp,0)
     se_mtx=zeros(len_tsamp,0)
     tsamp_mtx=zeros(len_tsamp,0)
     for i in 1:len_nDemes
        # avn, se, tsamp = num_origin_1d_t(1e7,nDemes_vec[i],1, 0.1,0.005,10,3000,1)[1:3]
        avn, se, tsamp = num_origin_1d_t(N,nDemes_vec[i],twoNmu,s0,mig,M, T, Gauss)
        avn_mtx=hcat(avn_mtx,avn)
        se_mtx=hcat(se_mtx,se)
     end

    #non-spatial
    simu_mtx,simu_mtx_se=num_origin_gaussian(0.1,1,1e7,10,3000)[2:3]
                        #num_origin_Gaussian(s0_seq,twoNmu_seq,N,replicates,T)
    len_tsamp=20 #20 generations
    T=3000
    tsamp=range(1,T+1,length=len_tsamp)
    tsamp=round.(Int,tsamp)

    avn_mtx=hcat(avn_mtx,simu_mtx[tsamp,:])
    se_mtx=hcat(se_mtx,simu_mtx_se[tsamp,:])

    legend=Array{String,2}(undef,1,len_nDemes+1)
    for i in 1:len_nDemes
        legend[1,i] = "nDemes = $(nDemes_vec[i])"
    end
    legend[1,end]="Non-spatial"
 return tsamp, legend, avn_mtx, se_mtx
 end


##Fig6.B. x dispersal, y time to fixation, series: nDemes
#Fig8. difer from fig6B 1. y as eta_fix 2. series: nDemes+s

include("3_3_1d_fixsamp.jl")
function onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss)
    len_mig=length(mig_vec)
    len_nDemes=length(nDemes_vec)
    #y: time to fixation
    t_avn_mtx=zeros(len_mig,len_nDemes)
    t_se_mtx=zeros(len_mig,len_nDemes)
    eta_avn_mtx=zeros(len_mig,len_nDemes)
    eta_se_mtx=zeros(len_mig,len_nDemes)
    fix_M_mtx=zeros(len_mig,len_nDemes)

    for i in 1:len_nDemes
        for j in 1:len_mig
            t_avn, t_se, eta_avn, eta_se, fix_M=
            num_origin_1d_fix(N,nDemes_vec[i],twoNmu,s0,mig_vec[j],M,T,Gauss)
            t_avn_mtx[j,i]=t_avn
            t_se_mtx[j,i]=t_se
            eta_avn_mtx[j,i]=eta_avn
            eta_se_mtx[j,i]=eta_se
            fix_M_mtx[j,i]=fix_M
        end
    end

    legend=Array{String,2}(undef,1,len_nDemes)
    for i in 1:len_nDemes
        legend[1,i] = "nDemes = $(nDemes_vec[i])"
    end
    fix_M_mtx=fix_M_mtx/M #fixation probability
    return legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx
end

##figure7.
# x: dispersal rate
# y: symptotic number of origin
# series: nDemes & nonspatial & Ralph-Coop model


#  t_avn_mtx1=zeros(len_mig,len_series)
#  t_se_mtx1=zeros(len_mig,len_series)
#  eta_avn_mtx1=zeros(len_mig,len_series)
#  eta_se_mtx1=zeros(len_mig,len_series)
#      for nDemes in nDemes_vec
#          t_avn, t_se, eta_avn, eta_se=
#          num_origin_1d_fix(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
#          t_avn_mtx1[mig]=t_avn
#          t_se_mtx1[mig]=t_se
#          eta_avn_mtx1[mig]=eta_avn
#          eta_se_mtx1[mig]=eta_se
#      end
#
# s0=0.1
# t_avn_mtx2=zeros(len_mig,len_series)
# t_se_mtx2=zeros(len_mig,len_series)
# eta_avn_mtx2=zeros(len_mig,len_series)
# eta_se_mtx2=zeros(len_mig,len_series)
#
# for nDemes in nDemes_vec
#     t_avn, t_se, eta_avn, eta_se=
#     num_origin_1d_fix(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
#     t_avn_mtx2[mig]=t_avn
#     t_se_mtx2[mig]=t_se
#     eta_avn_mtx2[mig]=eta_avn
#     eta_se_mtx2[mig]=eta_se
# end
#
# legend=Array{String,2}(undef,1,len_nDemes)
# for i in 1:len_nDemes
#     legend[1,i] = "nDemes = $(nDemes_vec[i])"
# end
#
# scatter(mig_vec,t_avn_mtx1,yerror=t_se_mtx1,label=legend_nDemes,markershape=:rect)
# xlabel!("dispersal rate")
# ylabel!("time to fixation (generation)")
# title!(s=0.05)
#
# scatter(mig_vec,t_avn_mtx2,yerror=t_se_mtx2,label=legend_nDemes,markershape=:circle)
# xlabel!("dispersal rate")
# ylabel!("time to fixation (generation)")
# title!(s=0.1)


##Fig10: x dispersal rate, y asymptotic number of origins, series peak-to-trough ratio
function dispersal_oscillating(Nmax,Nmin,har_geo,nDemes,twoNmu,s0,M,T,Gauss,
                                                            mig_vec,phi_vec)
    #phi as peak-to-trough ratio
    len_mig=length(mig_vec)
    len_phi=length(phi_vec)+1 #last +1 for constant population size

    t_avn_mtx=zeros(len_mig,len_phi)
    t_se_mtx=zeros(len_mig,len_phi)
    eta_avn_mtx=zeros(len_mig,len_phi)
    eta_se_mtx=zeros(len_mig,len_phi)
    legend=Array{String,2}(undef,1,len_phi)

    #series
    for i in 1:len_mig
        for j in 1:len_phi
            if j!=len_phi #Oscillating population size
                t_avn, t_se, eta_avn, eta_se=
                num_origin_1d_fix_osci(N, phi_vec[j], har_geo,nDemes,twoNmu,s0,mig_vec[i],M,T,Gauss)
            else #Constant  population size
                t_avn, t_se, eta_avn, eta_se, fix_M =
                num_origin_1d_fix(N,nDemes,twoNmu,s0,mig_vec[i],M,T,Gauss)
            end
            t_avn_mtx[i,j]=t_avn
            t_se_mtx[i,j]=t_se
            eta_avn_mtx[i,j]=eta_avn
            eta_se_mtx[i,j]=eta_se
        end
     end

    for k in 1:(len_phi-1)
        if har_geo=="H" #harmonic
            legend[1,k]= "HOP = $(phi_vec[k])"
        elseif har_geo=="G" #geometric
            legend[1,k]= "GOP = $(phi_vec[k])"
        end
    end
    legend[1,end]= "CP" #constant population size

    return legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx
 end


 #sampling of time-series total mutant frequency (X) and the number of
  #origin(eta, which is η)

 using Base.Threads

 function eta_X_tsamp(N,nDemes,twoNmu,s0_vec,mig,M,T) #Gauss=1
  len_s0=length(s0_vec) #last one for non-spatial
  legend=Array{String,2}(undef,1,len_s0)
  len_tsamp=30 #30 generations

  eta_avn_mtx=zeros(len_tsamp,0)
  eta_se_mtx=zeros(len_tsamp,0)

  X_avn_mtx=zeros(len_tsamp,0)
  X_se_mtx=zeros(len_tsamp,0)

  @threads for i in 1:len_s0
     tsamp, eta_avn, eta_se, X_avn, X_se = num_origin_1d_t(N,nDemes,twoNmu,s0_vec[i],mig,M,T,1)

     eta_avn_mtx=hcat(eta_avn_mtx,eta_avn)
     eta_se_mtx=hcat(eta_se_mtx,eta_se)

     X_avn_mtx=hcat(X_avn_mtx,X_avn)
     X_se_mtx=hcat(X_se_mtx,X_se)

     legend[1,i]="s = $(s0_vec[i])" #dispersal rate
  end

  len_tsamp=30 #30 generations
  tsamp=range(1,T+1,length=len_tsamp)
  tsamp=round.(Int,tsamp)

  return tsamp, legend, eta_avn_mtx, eta_se_mtx, X_avn_mtx, X_se_mtx
  end
