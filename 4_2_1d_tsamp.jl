using Distributions
using Random
using Base.Threads #@threads
include("4_1_1nonlocal_frequency.jl")

#sampling of num of origins at sampled time points
function num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss,nu)
 #Output
    #avn: average number of origins (eta) from replicates
    #se: standard error of eta from replicates
    #tsamp: sampled generation vector
    #medn: median of eta among replicates
    #lqrn: 1/4 quantile of eta among replicates
    #uqrn: 3/4 quantile of eta among replicates

    #X_avn: average (from replicates) total mutant frequency(X) from replicates
    #X_se: se (from replicates) of X from replicates
    #X_medn: median of X among replicates
    #X_lqrn: 1/4 quantile of X among replicates
    #X_uqrn: 3/4 quantile of X among replicates

 #Input
    # N: population size of metapopulation
    # twoNmu: the population-scaled mutation rate
    # s0 = twoNs/2/NDeme/nDemes #selection coefficient, s_ to avoid name crash with s in Julia
    # mig = twoNmig/2/NDeme/nDemes # migration rate
    # T: number of generation
    # nDemes: the number of demes
    # Gauss: 1 (using Gaussian approximation) or 0 (multinomial distribution)
    # nu: 𝝂 is the degree of freedom of the student t distribution deciding migration of individuals

    ## Setting the sampling generations (selective sampling points)
    t=collect(1:T+1)
    len_tsamp=20 #20 generations
    tsamp=range(1,T+1,length=len_tsamp)
    tsamp=round.(Int,tsamp) #Arbitrarily set to sample every 20 generations

    ns = 1000 #total sample size of meta-population freq_samp = Array{Array{Float64,2},2}(undef,len_tsamp,M)
    freq_samp_total = Array{Array{Float64,2},2}(undef,len_tsamp,M)
    freq_samp = Array{Array{Float64,2},2}(undef,len_tsamp,M)
    eta = zeros(len_tsamp,M) #Preallocating eta matrix
    X = zeros(len_tsamp,M) #Preallocating matrix for total frequency of mutants

    for m in 1:M ##Replicates
        println("replicate####### ",m) #Showing the progress of the simulation
        D, X_D = gaussian_1d_nonlocalmig(N,twoNmu,s0,mig,T,nDemes,Gauss,nu)
            #call 4_1_1nonlocal_frequency.jl to provide the allele frequency array

        ## Selecting ns of individuals
        for k in 1:len_tsamp
            println("generation##",k)
            indt=findfirst(isequal(tsamp[k]),t)

            #Randomly select ns individuals across all demes
            p = 1/nDemes*ones(nDemes) #vector of the probability of each deme
             #being sampled. Each deme has same chance of being picked
            K = rand(Multinomial(ns,p)) #Draw ns individuals,
             #where K is array holding number of individuals sampled from each deme
             #p has to be a column vector

            #Preallocating freq_sampkj
            nHaplotypes = size(D[1])[1] #nHaplotypes=number of rows in any one of the D arrays, since each element of the cell array should have same dimensions
            freq_sampkj = zeros(nDemes,nHaplotypes)

            #Sampling the frequency of K(j) number of haplotypes for each sampling generation
            for j in 1:nDemes
                p=vec(D[j][:,indt])
                 #p: frequncy vector in j-th deme in indt-th time point
                if !isempty(p[p.<0])
                    println("tsamp_p<0",p[p.<0])
                end

                if sum(p)!=1
                    println("tsamp_sum(p)", sum(p))
                end
                freq_sampkj[j,:] = rand(Multinomial(K[j],p))
            end

            freq_samp[k,m] = freq_sampkj #Putting freq_sampkj into an array for each tsamp and M

            nankm = zeros(len_tsamp,M) #Prellocating nankm
            if 1 in isnan.(freq_samp[k,m]) #any NaN found in the 2d array of kth generation and m-th replicate
                println("effing issues!")
                nankm[k,m] = NaN #??
            end

            #Finding the number of independent origins for each tsamp and M
            freq_sumkj=sum(freq_sampkj,dims=1) #2d #deme rows => 1 row of all demes

            X_ = freq_sumkj[1:end-1] #total frequency of all mutants of all demes
            eta[k,m] = count(a->a>0, X_) #eta: Number of individual origins

            X[k,m]=sum(X_)/ns #mean total frequency of a deme among demes
            freq_samp_total[k,m] = freq_sumkj #Sum over demes to give average total frequency of each haplotype. Any haplotype>0 indicates that it is present in at least one deme

        end #for k in 1:len_tsamp
    end #for m in 1:M

    ## Determining plotting variables
    #Preallocating vectors
    # len_tsamp=20 #20 time points (tsamp) at the outermost module in function module
    avn  = zeros(len_tsamp)
    se = zeros(len_tsamp)
    medn = zeros(len_tsamp)
    lqrn = zeros(len_tsamp)
    uqrn = zeros(len_tsamp)

    X_avn  = zeros(len_tsamp)#average total frequency of mutants among demes:X_
    X_se = zeros(len_tsamp)
    X_medn = zeros(len_tsamp)
    X_lqrn = zeros(len_tsamp)
    X_uqrn = zeros(len_tsamp)

    #Filling in vectors
    for k in 1:len_tsamp
          avn[k] = mean(eta[k,:])#mean
          se[k] = std(eta[k,:])/sqrt(M) #stadndard error
          medn[k] = median(eta[k,:])#median
          lqrn[k] = quantile(eta[k,:],0.25)#1/4
          uqrn[k]= quantile(eta[k,:],0.75)#3/4

          X_avn[k] = mean(X[k,:])
          X_se[k] = std(X[k,:])/sqrt(M)
          X_medn[k] = median(X[k,:])
          X_lqrn[k] = quantile(X[k,:],0.25)
          X_uqrn[k]= quantile(X[k,:],0.75)

    end

    return avn, se, tsamp, medn, lqrn, uqrn, X_avn, X_se, X_medn, X_lqrn, X_uqrn
    #all output are vectors
 end #function num_origin_1d_t()
