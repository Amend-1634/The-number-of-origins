#3_1_which.jl not included as it will be specified in test.jl according to setting
 #such as constant or oscillating population size
using Distributions
using Random
include("3_1_2notrack_cons.jl")

#1d local migration model

#sampling of num of origins along the time points
function num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss)

    # N: total population size of all demes
    # nDemes: deme number
    # twoNmu - population-scaled mutation rate (2Nμ)
    # s0: selection strength
    # mig: migration rate
    # M - number of replicates
    # T - genertions
    #Gauss: 0 (multinomial), 1 (using gaussian approximation)

    ## Setting the sampling generations (visually we need 20 or 30 poiints)
    t=collect(1:T+1)
    len_tsamp=30 #30 generations
    tsamp=range(1,T+1,length=len_tsamp)
    tsamp=round.(Int,tsamp)

    ns = 1000 #total sample size of meta-population freq_samp = Array{Array{Float64,2},2}(undef,len_tsamp,M)

    #Creat generic matrix
    #row: sampled time point; column: replicates
    freq_samp = Array{Array{Float64,2},2}(undef,len_tsamp,M)
    eta = zeros(len_tsamp,M) #eta matrix, eta(η) as number of origins
    X = zeros(len_tsamp,M) #matrix of total frequency of mutants


    for m in 1:M ##Replicates
        println("replicate#######",m) #Showing the progress of the simulation

        D = gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes,Gauss)

        ## Selecting ns of individuals
        for k in 1:len_tsamp
            println("generation##",k)
            indt=findfirst(isequal(tsamp[k]),t)

            #Randomly select ns individuals across all demes
            p = 1/nDemes*ones(nDemes) #vector of the probability of each deme
             #being sampled. Each deme has same chance of being picked
            K = rand(Multinomial(ns,p)) #Draw ns individuals, where K is array
             #holding number of individuals sampled from each deme
             #p has to be a column vector

            #Preallocating freq_sampkj
            nHaplotypes = size(D[1])[1] #nHaplotypes=number of rows in any one
             #of the D arrays, since each element of the cell array should have
             #same dimensions
            freq_sampkj = zeros(nDemes,nHaplotypes)

            #Sampling the frequency of K(j) number of haplotypes for each
             #sampling generation
            for j in 1:nDemes
                p2=vec(D[j][:,indt]) #indt: sampled time point
                                     #p2: probability vector of all haplotype at
                                      #this sampled time point
                freq_sampkj[j,:] = rand(Multinomial(K[j],p2))
                 #sampled individuals of each haplotype in corresponding deme
                 #row: deme number; col: haplotype
            end

            freq_samp[k,m] = freq_sampkj
            #each matrix 'freq_sampkj'--an element of the cell array 'freq_samp'
             #freq_samp:     row k: sampled time point, column m: replicate

            nankm = zeros(len_tsamp,M) #Prellocating nankm
            if 1 in isnan.(freq_samp[k,m]) #any NaN found in the 2d array in
             # kth generation and m-th replicate
                println("effing issues!")
                nankm[k,m] = NaN #??
            end

            #Sum of sampled individuals for each haplotype(column) across all
             #demes(row, now num of row=1)
            freq_sumkj=sum(freq_sampkj,dims=1) #freq_sumkj is still 2d array
            X_=freq_sumkj[1:end-1] #total mutants individuals across all demes
            eta[k,m] = count(a->a>0, X_) #eta: Number of individual origins
             #object in eta: the number of origins in k (time point) in m (replicate)

            X[k,m]=sum(X_)/ns #mean total mutant frequency among demes

        end #for k in 1:len_tsamp
    end #for m = 1:M


    ## Determining plotting variables
    #Preallocating vectors
    # len_tsamp=20 #20 time points (tsamp) at the outermost module in function module
    eta_avn  = zeros(len_tsamp)
    eta_se = zeros(len_tsamp)
    # medn = zeros(len_tsamp)
    # lqrn = zeros(len_tsamp)
    # uqrn = zeros(len_tsamp)

    X_avn  = zeros(len_tsamp)#average total frequency of mutants among demes:X_
    X_se = zeros(len_tsamp)


    #Filling in vectors
    for k in 1:len_tsamp
          eta_avn[k] = mean(eta[k,:])
          eta_se[k] = std(eta[k,:])/sqrt(M)
          #commented out cuz not outputed
          # medn[k] = median(eta[k,:])
          # lqrn[k] = quantile(eta[k,:],0.25)
          # uqrn[k]= quantile(eta[k,:],0.75)

          X_avn[k] = mean(X[k,:])
          X_se[k] = std(X[k,:])/sqrt(M)
    end


    return tsamp, eta_avn, eta_se, X_avn, X_se #vectors
 end #function num_origin_1d_t()
