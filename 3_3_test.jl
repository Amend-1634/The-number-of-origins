pwd()
cd("D:\\Julia\\func")
using Profile
using BenchmarkTools
# include("2_1allele_var_mtx.jl")
include("2_2gaussian_multinomial.jl")
include("3_1_1dfrequency1.jl")
# @time gaussian_1d_localmig(1e6,10,0.05,5e-5,200,3, 1)
 #for the block before migration
# @time D, migLR, migrantCount=gaussian_1d_localmig(1e6,10,0.05,5e-5,200,3, 1)


##why deme population size change so much if tracking migrants? migration coefficient?
#total nubmer of case corresp onding to migration/utatio coefficient:
include("2_2gaussian_multinomial.jl")
include("3_1_1track_cons.jl")
Random.seed!(1234)
@time D= gaussian_1d_localmig(5e3,1,0.1,5e-5,200,5, 1)
using Plots

Plots.plot(D[1]')
using StatsPlots, Plots
violin(D[1],leg=false)

violin(["Series 1" "Series 2" "Series 3" "Series 4"], y, leg = true)



@time D, migrantCount, mutantCount = gaussian_1d_localmig(5e3,1,0.1,5e-5,200,5, 1)
                                   # gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss)

sum(migrantCount)
sum(mutantCount)

@time D, migrantCount, mutantCount = gaussian_1d_localmig(5e6,1,0.1,5e-5,200,5, 1) #8 seconds
                                   # gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss)
sum(migrantCount)
sum(mutantCount)

#plot deme population size along with time (how much fluctuate)
include("3_1_1track_cons_NDeme.jl")
@time D, NDeme_mtx = gaussian_1d_localmig(5e7,1,0.1,5e-5,2000,10, 1)
                    #gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss)
sum(NDeme_mtx[5,:])#0, why
plot(NDeme_mtx[:,2:2001]',linewidth=1,alpha=0.8,dpi=300,legend=:outertopright)
#same as 2001
xlabel!("generations")
ylabel!("deme population size (originally is 10^7)")
savefig("..\\photo\\track_cons_NDeme2")
plot(D[5]')#acceptable
plot(D[4]')
plot(D[3]')


##Compare track and no track
 #a. allele frequency #b. num of origins

#a.allele frequency
include("3_1_1track_cons.jl")
Random.seed!(1234)
@time D= gaussian_1d_localmig(5e3,10,0.1,5e-5,200,5, 1)
plot(D[2]',yscale=:log10,ylims=(1e-3,1),legend=false,dpi=300)
savefig("3_1_1track_cons.jpg")

include("3_1_2notrack_cons.jl")
Random.seed!(1234)
@time D= gaussian_1d_localmig(5e3,10,0.1,5e-5,200,5, 1) #8 seconds
        #gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss)#
#ok after ind=>frequency but return D in a wrong way-now solved scaling
 #within or outside of for loop

plot(D[2]',yscale=:log10,ylims=(1e-3,1),legend=false,dpi=300)
savefig("3_1_2notrack_cons.jpg")
#conclusion: differ a lot regardless of same sampling seed

#b.how about num of origins?
N=1e9;nDemes=5;twoNmu=1;s0=0.1;mig=0.05;M=5;T=200;Gauss=0

N=1e5;nDemes=5;twoNmu=1;s0=0.1;mig=0.05;M=5;T=200;Gauss=0
include("3_1_1track_cons.jl")
include("3_2_1d_tsamp.jl")
@time tsamp, eta_avn, eta_se, X_avn, X_se = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)

#13s, 8s much faster
scatter(tsamp, eta_avn, yerror=eta_se,lab="tracking mig",legend=:left)

plot(tsamp,X_avn)

N=1e5;nDemes=5;twoNmu=1;s0=0.1;mig=0.05;M=5;T=200;Gauss=0
include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
@time avn, se, tsamp = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3]
tsamp
avn
#18s,27s
scatter!(tsamp, avn, yerror=se,lab=" not tracking mig")
title!("N=1e5")
#debug

##"3_1_2notrack_cons.jl"

# Random.seed!(1234)
using BenchmarkTools
include("3_1_2notrack_cons.jl");println("dense") #42.243 ms
@btime D= gaussian_1d_localmig(5e3,10,0.1,5e-5,50,5, 1)

include("3_1_2notrack_cons_sparse.jl");println("sparse")
@time D= gaussian_1d_localmig(5e3,10,0.1,5e-5,50,5, 1) # 0.287747s
#7 times of time--3/4 space
#will be only used for large dataset
#use sparse matrix and avoid @threads to avoid Outofmemory() error

#used in fixsamp instead of tsamp
typeof(D[1])
typeof(sparse(rand(2,3)))

D[2]
extrema(D[2])
plot(D[2]',yscale=:log10,ylims=(1e-5,1),legend=false)
plot(D[2]')
#finally function well

include("4_1_1nonlocal_frequency_sparse.jl")
D, X_D= gaussian_1d_nonlocalmig(1e5,10,0.1,5e-5,100,5,1,1)
      # gaussian_1d_nonlocalmig(N,twoNmu,s0,mig,T,nDemes, Gauss,nu)#N




##"3_1_3notrack_osci.jl"

a=collect(1:50)
b=0.5(Nmax+Nmin) .+ 0.5(Nmax-Nmin)*sin.(2pi*a/12)
#so the lowest point won't reach 0
b=b/1e8
plot!(b,color=:red)

include("3_1_3notrack_osci.jl")
N=1e7; phi=10
@time D= gaussian_1d_localmig(N, phi, "H",       1, 0.1,5e-7,2000,5, 1) #132.452806 s
#function gaussian_1d_localmig(N, phi, har_geo, twoNmu,s0,mig,T,nDemes,Gauss)#N
D[1]

plot(D[5]',alpha=0.8,legend=false,dpi=300,yscale=:log10,ylims=(1e-7,1))
#the mutant frequency is extremly low?
#mig lower(mig=5e-7)=>bit more but same level
#longer time (2000)

plot!(b,color=:red)
savefig("..\\photo\\osci_period")

xlabel!("generations")
ylabel!("allele frequency or red total N*10^(-8)")
title!("allele frequency of oscillating total N")
savefig("..\\photo\\osci_freq")


##tsamp

N=1e7;
N=1e9; nDemes=10;twoNmu=1;s0=0.1;mig=5e-5;M=10;T=2000
 # N=1e8;

include("3_1_1track_cons.jl")
include("3_2_1d_tsamp.jl")
Gauss=1
@time avn, se, tsamp = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3] #40.9s
plot(tsamp, avn, yerror=se,lab="Gaussian",dpi=300)

Gauss=0
@time avn, se, tsamp = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3] #48s

N=1e9; nDemes=10;twoNmu=1;s0=0.1;mig=5e-5;M=10;T=2000;Gauss=1
include("3_1_1track_cons.jl")
include("3_2_1d_tsamp.jl")
using Profile
@profiler avn, se, tsamp = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3] #40.9s


plot!(tsamp, avn, yerror=se,lab="Multinomial",legend=:inside)
xlabel!("generations")
ylabel!("the number of origins")
# savefig("..\\photo\\1d_N1e7_Gaussian_multinomial_vary") #differ
savefig("..\\photo\\1d_N1e9_Gaus_multino_vary")#differ
#Conclusion: for 1d model, Gaussian may not substitute multinomial distribution
 #for 1e7,8,9 so later on I will use Multinomial distribution
-----------
include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
Gauss=1
@time avn, se, tsamp = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3] #40.9s
plot(tsamp, avn, yerror=se,lab="Gaussian",dpi=300)

Gauss=0
@time avn, se, tsamp = num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3] #48s
plot!(tsamp, avn, yerror=se,lab="Multinomial",legend=:inside)
xlabel!("generations")
ylabel!("the number of origins")
# savefig("..\\photo\\1d_N1e7_Gaussian_multinomial_vary") #differ
savefig("..\\photo\\1d_N1e8_Gaus_multino_vary")#differ

----------
nDemes=5;twoNmu=1;s0=0.1;mig=5e-5;M=3;T=200
Nmax=5.5e7
Nmin=5.5e6
include("3_2_1d_tsamp_osci.jl")
Gauss=0
har_geo="H"
@time avn, se, tsamp = num_origin_1d_t_osci(Nmax, Nmin, har_geo, nDemes,twoNmu,s0,mig,M,T,Gauss)[1:3]
#157s
plot(tsamp, avn, yerror=se,lab="Multinomial",dpi=300)
xlabel!("generations")
ylabel!("the number of origins")
savefig("..\\photo\\tsamp_osci_test")

##Why the number of origins in time-series is half of Joejyn's simulation?
#


##fixsamp
include("3_1_2notrack_cons.jl")
include("3_3_1d_fixsamp.jl")
using Profile
@profiler t_avn, t_se, eta_avn, eta_se, fix_M = num_origin_1d_fix(1e3,5,50,0.1,0.05,3,200,1)

include("3_1_2notrack_cons.jl")
include("3_3_1d_fixsamp.jl")
using Profile
@profiler t_avn, t_se, eta_avn, eta_se, fix_M = num_origin_1d_fix(1e3,5,50,0.1,0.05,3,200,1)

@time t_avn, t_se, eta_avn, eta_se, fix_M = num_origin_1d_fix(1e3,5,50,0.1,0.05,3,200,1)
                                # num_origin_1d_fix(N,nDemes,twoNmu,s0,mig,M,T,Gauss)
#checked, function well


##Large nDemes (100)
# include("3_1_1track_cons.jl")
include("3_1_2notrack_cons.jl") #Bhavin preferred
#Compare between different deme nubmers:
#nDemes changing:
 #1.should the NDeme be the constant (//I preferred)?
    #as constant N will make the NDeme in accordance with higher nDemes
    #more senesitive to genetic drift
    #as spatially expand-not influence local density
 #2.should population-scaled parameter (twoNmu, twoNs, twoNmig) be the constant
  #in comparison? or mu s mig themself (//I preferred)?
@time D1= gaussian_1d_localmig(1e9,100,0.1,5e-7,2000,100, 0) #14.17hr
         #gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss
@time D2= gaussian_1d_localmig(5e7,5,0.1,5e-7,2000,5, 0) #8s @Clifford
@time D3= gaussian_1d_localmig(1e10,1000,0.1,5e-7,2000,1000, 0) #
 #assume that twoNmu for nDemes=1 NDeme=1e7 is 1


using JLD
pwd()
readdir()
save("nDemes_100_D1.jld","D1",D1) #am
save("nDemes_5_D2.jld","D2",D2) #am2
save("nDemes_1000_D3.jld","D3",D3) #am2


readdir()
using JLD
@load "notracking_100nDemes.jld"


"notracking_NDeme1e7.jld"
"tracking_100nDemes.jld"
"tracking_100nDemes_1e7deme.jld"


using JLD
@load "nDemes_100_D1.jld"
@load "nDemes_5_D2.jld"


using Plots

# the mutation rate of spatial one seems to be less than non-spatial

#populatio size of each deme too small so the effect of genetic drift is too significant?
plot(D[1]',yscale=:log10,ylims=(1e-5,1),legend=false,dpi=300)
plot(D1[2]',yscale=:log10,ylims=(1e-7,1),legend=false,dpi=300)
plot(D1[3]',yscale=:log10,ylims=(1e-7,1),legend=false,dpi=300)

plot(D2[4]',yscale=:log10,ylims=(1e-7,1),legend=false,dpi=300)
plot(D2[5]',yscale=:log10,ylims=(1e-7,1),legend=false,dpi=300)

# why quick mutation and reach fixation during 500~1500 generations
# Is the dispersal  rate set too large? (range of dispersal rate)


##
N_vec=[5e3 5e2]

N=5e3
s0_vec=[0.5 0.1]
twoNmu_vec=[0.1 1]
nDemes=5
mig=0.005
Gauss=1
M=2
T=100

include("3_1_2notrack_cons.jl")
include("3_2_1d_tsamp.jl")
include("3_4_figjoejyn_series_threads.jl")
tsamp, legend, avn_mtx, se_mtx=s0_twoNmu_N_series(nDemes,N,twoNmu_vec,s0_vec, mig, M, T, Gauss)
plot(tsamp,avn_mtx,yerror=se_mtx)
legend
