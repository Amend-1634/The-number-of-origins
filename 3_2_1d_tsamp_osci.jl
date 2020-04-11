# include("3_1_1dfrequency3.jl")
using Distributions
using Random
include("3_1_3notrack_osci.jl")

#1d local migration model
#oscillating population size
#the only difference from 3_2_1d_tsamp.jl(constant population size) is
 #the function(oscillating population size) called in line48

#sampling of num of origins along the time points
function num_origin_1d_t_osci(Nmax, Nmin, har_geo, nDemes,twoNmu,s0,mig,M,T,Gauss)

    # Nmax: maximum population size in oscillation
    # Nmin: minimum population size in oscillation

    #har_geo: "H" as harmonically constrained mean
             #"G" as geometrically constrained mean
    # nDemes: deme number
    # twoNmu - population-scaled mutation rate (2Nμ)
    # s0: selection coefficient
    # mig: migration coefficient
    # M - number of replicates
    # T: generations
    # Gauss: 0 (multinomial), 1 (using gaussian approximation)

    ## Setting the sampling generations (since you dont want to sample every generation)
    t=collect(1:T+1)
    len_tsamp=20 #20 generations
    tsamp=range(1,T+1,length=len_tsamp)
    tsamp=round.(Int,tsamp) #Arbitrarily set to sample every 20 generations
    # T = maximum(tsamp) #In case T was not divisible by 20 #unnecessary

    ## Setting some variables
    # NDeme = N/nDemes #Population size of each deme=total population size/number of deme
    # N: total population size of all demes
    # NDemes: population size of each demes

    Nsamp = 1000 #total sample size of meta-population
    freq_samp_total = Array{Array{Float64,2},2}(undef,len_tsamp,M)
    freq_samp = Array{Array{Float64,2},2}(undef,len_tsamp,M)
    eta = zeros(len_tsamp,M) #Preallocating eta

    for m in 1:M ##Replicates

        println("replicate#######",m) #Showing the progress of the simulation

        ## Generating the frequency array of each deme
        D = gaussian_1d_localmig(Nmax, Nmin, har_geo, twoNmu,s0,mig,T,nDemes,Gauss)

        ## Selecting Nsamp of individuals
        for k in 1:len_tsamp
            println("generation##",k)
            # indt = t == tsamp(k);
            indt=findfirst(isequal(tsamp[k]),t)
            println("tsamp[k] ",tsamp[k], " t ", size(t))
            println("indt ", indt)
            #Randomly select Nsamp individuals across all demes
            p = 1/nDemes*ones(nDemes) #vector of the probability of each deme being sampled. Each deme has same chance of being picked
            K = rand(Multinomial(Nsamp,p)) #Draw Nsamp individuals, where K is array holding number of individuals sampled from each deme
             #p has too be a column vector

            #not necessary: as k-dimensional integer vector from Multinomial() that sums to n.
            # #Check that the total number sampled across all demes = Nsamp
            # if sum(K)!=Nsamp
            #     println("Issues!") #     println(["sum(k) = ",round(sum(k))])
            #     println(["Nsamp = ",round(Nsamp)])
            # end

            #Preallocating freq_sampkj
            nHaplotypes = size(D[1])[1] #nHaplotypes=number of rows in any one of the D arrays, since each element of the cell array should have same dimensions
            freq_sampkj = zeros(nDemes,nHaplotypes)

            #Sampling the frequency of K(j) number of haplotypes for each sampling generation
            for j in 1:nDemes
                p=vec(D[j][:,indt])#indt: sampled time point
                # println("p_demes ",D[j][:,indt],
                        # " size ",size(p)," sum(p) ",sum(p))
                #undef
                println("p ",size(p))
                println("sum(p) ",sum(p))
                freq_sampkj[j,:] = rand(Multinomial(K[j],p))
            end

            #Checking that there are no issues with sampling. Each cell should
             #be either 0 (no mutant) or some frequency (mutant present)

            freq_samp[k,m] = freq_sampkj #Putting freq_sampkj into an array for each tsamp and M

            nankm = zeros(len_tsamp,M) #Prellocating nankm
            # println(" size(freq_samp[k,m])  ",size(freq_samp[k,m])) #2d
            if 1 in isnan.(freq_samp[k,m]) #any NaN found in the 2d array of kth generation and m-th replicate
                println("effing issues!")
                nankm[k,m] = NaN #??
            end

            #Finding the number of independent origins for each tsamp and M
            # freq_samp_total = cell(len_tsamp,M); #Preallocating freq_samp_total
            freq_sumkj=sum(freq_sampkj,dims=1) #2d
            println("freq_sumkj", size(freq_sumkj)) #(1,49)
            # println(" eta ", freq_samp_total[k,m]) #all undef, not reported as NaN
              #undef object cannot be indexed before defined
            # println(" eta[k,m] ", freq_samp_total[k,m])

            #indexing problem?????????
            freq_samp_total[k,m] = freq_sumkj #Sum over demes to give total frequency of each haplotype. Any haplotype>0 indicates that it is present in at least one deme

            indsamp = findall(a->a>0,freq_samp_total[k,m][1:end-1]) #Finding the mutant haplotypes that did not die out (non-zero haplotypes)
            eta[k,m] = length(indsamp) #eta: Number of individual origins

        end #for k in 1:len_tsamp
    end #for m = 1:M


    ## Determining histogram variables
    #Preallocating vectors
    # len_tsamp=20 #20 generations at the outermost module in function module
    println(len_tsamp)
    avn  = zeros(len_tsamp)
    se = zeros(len_tsamp)
    medn = zeros(len_tsamp)
    lqrn = zeros(len_tsamp)
    uqrn = zeros(len_tsamp)

    #Filling in vectors
    for k in 1:len_tsamp
          avn[k] = mean(eta[k,:])
          se[k] = std(eta[k,:])/sqrt(M)
          medn[k] = median(eta[k,:])
          # lqrn[k] = prctile(eta[k,:],25);#lower percentiles of a sample.
          # uqrn[k]= prctile(eta[k,:],75);#upper percentiles of a sample.
          lqrn[k] = quantile(eta[k,:],0.25)
          uqrn[k]= quantile(eta[k,:],0.75)
    end

    # #Plot
    # plot(layout = 2)
    # scatter(tsamp,avn, markersize  = 12, linewidth = 1, yerror=stdn,subplot = 1)
    # scatter(tsamp,medn, markersize  = 12, linewidth = 1, yerror=stdn,subplot = 2)
    #? BoundsError: attempt to access 1-element Array{Plots.Subplot,1} at index [2]
    # # errorbar(tsamp,avn,stdn,'s:','MarkerSize',12,'LineWidth',1,
    # # 'MarkerFaceColor','k');
    # title!(["N=$N, 2Ns=$twoNs, 2N*mu=$twoNmu, K=$K"])
    # ylabel!("Mean number of originations")
    # xlabel!("Time (Generations)")
    # xlims!((0,1.1*maximum(tsamp)))
    # save("1d_simu_gaussian.jld")


    return avn, se, tsamp, medn, lqrn, uqrn #vectors
 end #function num_origin_1d_t()
