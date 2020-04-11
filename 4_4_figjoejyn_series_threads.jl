#only difference from "4_4_figjoejyn_series.jl" is
 #@threads applied in this file
#should include("4_4_figjoejyn_series_threads.jl") instead of
 #include("4_4_figjoejyn_series.jl") when:
   #a. the parameter vector '*_vec' is surely a vector instead of one-element object
      #as @threads not applied to one-element obj)
   #b. when the data set if not large (in other words your parameter range
    #will not cause too much haplotypes & longer time to fixation=> larger
     #data set; as @threads tend to increase the memory allocation when reducing
      #the computational time)

using Base.Threads

#figure4
#x: time series; y: eta; series: mig
function onedim_dispersal_series(N,nDemes,twoNmu,s0,mig_vec,M,T,nu) #,Gauss
 len_mig=length(mig_vec)+1 #last one for non-spatial
 legend=Array{String,2}(undef,1,len_mig)
 len_tsamp=20 #20 generations
 avn_mtx_g=zeros(len_tsamp,0)
 se_mtx_g=zeros(len_tsamp,0)

 avn_mtx_m=zeros(len_tsamp,0)
 se_mtx_m=zeros(len_tsamp,0)

 @threads for i in 1:(len_mig-1)
    # avn, se, tsamp = num_origin_1d_t(1e5,5,1,0.1,mig_vec[i],10,2000,1)[1:3]
    # @time avn_g, se_g, tsamp = num_origin_1d_t(5e5,5,1,0.1,mig_vec[i],3,200,1)[1:3]
    avn_g, se_g = num_origin_1d_t(N,nDemes,twoNmu,s0,mig_vec[i],M,T,1,nu)[1:2]

    # # @time avn_m, se_m= num_origin_1d_t(5e5,5,1,0.1,mig_vec[i],3,200,0)[1:2]
    avn_m, se_m= num_origin_1d_t(N,nDemes,twoNmu,s0,mig_vec[i],M,T,0,nu)[1:2]
    #                  # num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
    #_g for gaussian while _m for multinomial sampling
    avn_mtx_g=hcat(avn_mtx_g,avn_g)
    se_mtx_g=hcat(se_mtx_g,se_g)
    avn_mtx_m=hcat(avn_mtx_m,avn_m)
    se_mtx_m=hcat(se_mtx_m,se_m)
    legend[1,i]="2*NDeme*mig = $(mig_vec[i])" #dispersal rate
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
function onedim_nDemes_series(nDemes_vec,N,twoNmu,s0,mig, M, T, Gauss, nu)
    len_nDemes=length(nDemes_vec)
    len_tsamp=20 #20 generations
     avn_mtx=zeros(len_tsamp,0)
     se_mtx=zeros(len_tsamp,0)
     tsamp_mtx=zeros(len_tsamp,0)
     # @threads
     for i in 1:len_nDemes
        # avn, se, tsamp = num_origin_1d_t(1e7,nDemes_vec[i],1, 0.1,0.005,10,3000,1)[1:3]
        avn, se, tsamp = num_origin_1d_t(N,nDemes_vec[i],twoNmu,s0,mig,M, T, Gauss, nu)
        avn_mtx=hcat(avn_mtx,avn)
        se_mtx=hcat(se_mtx,se)
     end

    #non-spatial
    simu_mtx,simu_mtx_se = num_origin_gaussian(s0,twoNmu,N,M,T)[2:3]
                         # num_origin_gaussian(s0_seq,twoNmu_seq,N,replicates,T)[2:3]

                          #num_origin_gaussian(0.1,1,1e7,10,3000)[2:3]
    len_tsamp=20 #20 generations
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

using JuliaDB
function onedim_xdispersal_yfix_nDemesseries(N,nDemes_vec,twoNmu,s0,mig_vec,M,T,Gauss, nu)
    len_mig=length(mig_vec)
    len_nDemes=length(nDemes_vec)
    #y: time to fixation
    t_avn_mtx=zeros(len_mig,len_nDemes)
    t_se_mtx=zeros(len_mig,len_nDemes)
    eta_avn_mtx=zeros(len_mig,len_nDemes)
    eta_se_mtx=zeros(len_mig,len_nDemes)
    fix_M_mtx=zeros(len_mig,len_nDemes)

     for i in 1:len_nDemes
          @threads for j in 1:len_mig
            t_avn, t_se, eta_avn, eta_se, fix_M=
            num_origin_1d_fix(N,nDemes_vec[i],twoNmu,s0,mig_vec[j],M,T,Gauss, nu)
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


##Fig10: x dispersal rate, y asymptotic number of origins, series peak-to-trough ratio

function dispersal_oscillating(N, phi_vec, har_geo,nDemes,twoNmu,s0,mig_vec,M,T,Gauss,nu)
    #phi as peak-to-trough ratio
    len_mig=length(mig_vec)
    len_phi=length(phi_vec)+1 #last +1 for constant population size

    t_avn_mtx=zeros(len_mig,len_phi)
    t_se_mtx=zeros(len_mig,len_phi)
    eta_avn_mtx=zeros(len_mig,len_phi)
    eta_se_mtx=zeros(len_mig,len_phi)
    legend=Array{String,2}(undef,1,len_phi)

    #series
    @threads for i in 1:len_mig
        for j in 1:len_phi
            if j!=len_phi #Oscillating population size
                t_avn, t_se, eta_avn, eta_se=
                num_origin_1d_fix_osci(N, phi_vec[j], har_geo,nDemes,twoNmu,s0,mig_vec[i],M,T,Gauss, nu)
            else #Constant  population size
                t_avn, t_se, eta_avn, eta_se, fix_M =
                num_origin_1d_fix(N,nDemes,twoNmu,s0,mig_vec[i],M,T,Gauss, nu)
            end
            t_avn_mtx[i,j]=t_avn
            t_se_mtx[i,j]=t_se
            eta_avn_mtx[i,j]=eta_avn
            eta_se_mtx[i,j]=eta_se
        end#for j
    end#for i

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

##Average total frequency of all mutants
#time-series, compare mig as series
function xtime_yX_migseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss,nu)
    len_mig=length(mig_vec)
    len_tsamp=20 #20 generations

    tsamp=range(1,T+1,length=len_tsamp)
    tsamp=round.(Int,tsamp)

    X_avn_mtx=zeros(len_tsamp,0)
    X_se_mtx=zeros(len_tsamp,0)
    X_medn_mtx=zeros(len_tsamp,0)
    X_lqrn_mtx=zeros(len_tsamp,0)
    X_uqrn_mtx=zeros(len_tsamp,0)
    legend=Array{String,2}(undef,1,len_mig)

    @threads for i in 1:length(mig_vec)
        X_avn, X_se, X_medn, X_lqrn, X_uqrn =
         num_origin_1d_t(N,nDemes,twoNmu,s0,mig_vec[i],M,T,Gauss,nu)[7:end]

        X_avn_mtx=hcat(X_avn_mtx,X_avn)
        X_se_mtx=hcat(X_se_mtx,X_se)
        X_medn_mtx=hcat(X_medn_mtx,X_medn)
        X_lqrn_mtx=hcat(X_lqrn_mtx,X_lqrn)
        X_uqrn_mtx=hcat(X_uqrn_mtx,X_uqrn)
        legend[1,i]="2*NDeme*mig = $(mig_vec[i])"

    end

   return tsamp, legend, X_avn_mtx, X_se_mtx, X_medn_mtx, X_lqrn_mtx, X_uqrn_mtx
end

##Compare 𝝂 (called nu) the degree of freedom between student t distributions
#time-series
function series_nu(N,nDemes,twoNmu,s0,mig,M,T,nu_vec) #,Gauss
 len_nu=length(nu_vec) #last one for non-spatial
 legend=Array{String,2}(undef,1,len_nu)
 len_tsamp=20 #20 generations
 avn_mtx_g=zeros(len_tsamp,0)
 se_mtx_g=zeros(len_tsamp,0)

 avn_mtx_m=zeros(len_tsamp,0)
 se_mtx_m=zeros(len_tsamp,0)

 @threads for i in 1:len_nu
    # avn, se, tsamp = num_origin_1d_t(1e5,5,1,0.1,mig_vec[i],10,2000,1)[1:3]
    # @time avn_g, se_g, tsamp = num_origin_1d_t(5e5,5,1,0.1,mig_vec[i],3,200,1)[1:3]
    avn_g, se_g = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,1,nu_vec[i])[1:2]

    # # @time avn_m, se_m= num_origin_1d_t(5e5,5,1,0.1,mig_vec[i],3,200,0)[1:2]
    avn_m, se_m= num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,0,nu_vec[i])[1:2]
    #                  # num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
    #_g for gaussian while _m for multinomial sampling
    avn_mtx_g=hcat(avn_mtx_g,avn_g)
    se_mtx_g=hcat(se_mtx_g,se_g)
    avn_mtx_m=hcat(avn_mtx_m,avn_m)
    se_mtx_m=hcat(se_mtx_m,se_m)
    legend[1,i]="nu = $(nu_vec[i])" #the degree of freedom
 end

 len_tsamp=20 #20 generations
 tsamp=range(1,T+1,length=len_tsamp)
 tsamp=round.(Int,tsamp)

 return tsamp, legend, avn_mtx_g, se_mtx_g, avn_mtx_m, se_mtx_m
 end


##x:mig, y: asymptotic eta or t, series: nu

function onedim_xdispersal_yfix_nuseries(N,nDemes,twoNmu,s0,mig_vec,M,T,Gauss, nu_vec)
    len_mig=length(mig_vec)
    len_nu=length(nu_vec)
    #y: time to fixation
    t_avn_mtx=zeros(len_mig,len_nu)
    t_se_mtx=zeros(len_mig,len_nu)
    eta_avn_mtx=zeros(len_mig,len_nu)
    eta_se_mtx=zeros(len_mig,len_nu)
    fix_M_mtx=zeros(len_mig,len_nu)

     @threads for i in 1:len_nu
         for j in 1:len_mig
            t_avn, t_se, eta_avn, eta_se, fix_M=
            num_origin_1d_fix(N,nDemes,twoNmu,s0,mig_vec[j],M,T,Gauss, nu_vec[i])
            t_avn_mtx[j,i]=t_avn
            t_se_mtx[j,i]=t_se
            eta_avn_mtx[j,i]=eta_avn
            eta_se_mtx[j,i]=eta_se
            fix_M_mtx[j,i]=fix_M
        end#for j
    end#for i

    legend=Array{String,2}(undef,1,len_nu)
    for i in 1:len_nu
        legend[1,i] = "nu = $(nu_vec[i])"
    end
    fix_M_mtx=fix_M_mtx/M #fixation probability
    return legend, t_avn_mtx, t_se_mtx, eta_avn_mtx, eta_se_mtx,fix_M_mtx
end
