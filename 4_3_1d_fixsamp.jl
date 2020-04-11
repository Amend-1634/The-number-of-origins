using Statistics
using Multibreak
# using Base.Threads #@threads
include("4_1_1nonlocal_frequency.jl")

#use dense matrix for tsamp (quicker)
#use sparse matrix (save memory) for fixsamp=>larger dataset (as longer time
     #employed to ensure the fixation reached)

#sampling data(time, number of origins) when reaching fixation
function num_origin_1d_fix(N,nDemes,twoNmu,s0,mig,M,T,Gauss,nu)

    #Output
        # t_avn: average (from replicates) time to fixation
        # t_se:  se
        # eta_avn: average number of origins
        # eta_se:  se
        # fix_M: fixation reached among replicates (check this in case
         # sampling bias due to sampling upper limit of time)

    #Input
    # N: population size of metapopulation
    # nDemes: the number of demes
    # twoNmu: the population-scaled mutation rate
    # s0 = twoNs/2/NDeme/nDemes #selection coefficient, s_ to avoid name crash with s in Julia
    # mig = twoNmig/2/NDeme/nDemes # migration rate
    # M: replicate number
    # T: number of generation
    # Gauss: 1 (using Gaussian approximation) or 0 (multinomial distribution)
    # nu: 𝝂 is the degree of freedom of the student t distribution deciding migration of individuals

    Nsamp = 1000 #total sample size of meta-population
    t_fix_vec=[] #vector of time to faxiation
    eta_fix_vec=[] #vector of asymptotic number of origins
    fix_M=0 #how many cases reaching fixation in replicates

    @multibreak for m in 1:M ##Replicates
        # println("replicate#######",m)  #Showing the progress of the simulation

        D = gaussian_1d_nonlocalmig(N,twoNmu,s0,mig,T,nDemes,Gauss,nu)[1]

        #Detect the data when reaching fixation (wild type frequency=0 in
         #all demes)
         for t in 1:T+1 #generation
            for i in 1:nDemes
                 n_i=D[i]
                 xwt=n_i[end,t]
                 if xwt==0
                     if i!=nDemes
                         continue
                     elseif i==nDemes
                         # println("reached fixation#####  t  ",t)
                         fix_M+=1 #+1 for counting fixation cases when reachign fixation
                         push!(t_fix_vec,t)

                         #the number of origins at fixation time point
                         n_mu=size(D[1])[1]-1 #total mutant haplotype amount(-1 for wild type)
                         x=zeros(n_mu,1)
                         for i in 1:nDemes
                             x_=D[i][1:end-1,t]
                             x = x + x_
                         end
                         eta_fix=count(a->a>0,x) #count the origins number
                         push!(eta_fix_vec, eta_fix)

                         break; break; continue # to start next replication
                     end
                 else
                     break; continue #break nDemes loop and skip this generation
                 end
                 # println(" m ", m, " t ", t," i ", i ) #ensure it follows break and jumpt through this end
            end #for i
            # println("m ", m, " t ", t)
        end #for t
        # println("m ", m)

     end #for m
     t_avn=mean(t_fix_vec) #average time to fixation among replicates
     t_se=std(t_fix_vec)/sqrt(M) #se

     eta_avn=mean(eta_fix_vec) #mean asymptotic eta
     eta_se=std(eta_fix_vec)/sqrt(M) #se

     return t_avn, t_se, eta_avn, eta_se, fix_M #output as one number
end #function num_origin_1d_fix()
