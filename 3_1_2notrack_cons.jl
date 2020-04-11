using PoissonRandom
using LazyArrays
using Multibreak
using Base.Threads
using Plots
include("2_2gaussian_multinomial.jl")

#1d local migration model
#scaling without traking the migrant individuals
 # =>constant deme size

# function gaussian_1d_localmig(NDeme,twoNmu,twoNs,twoNm,T,nDemes, Gauss)

function gaussian_1d_localmig(N,twoNmu,s0,mig,T,nDemes, Gauss)#N

    #Output
     #D: cell array containing frequency array in each deme -
      #rows of this array correspond to independent origins/haplotypes,
       #with common ordering across all demes (same num of rows)

    #Input:
    #N: population size of metapopulation
    #twoNmu: population-scaled mutation rate (2Nμ)
    #s0: initial selection coefficient
    #mig: migration rate
    #T: number of generation
    #nDemes: the number of demes
    #Gauss: 1 (using Gaussian approximation) or 0 (multinomial distribution)


    NDeme=round(Int,N/nDemes) #initial value as integer
    # NDeme: population size of each deme

    ## Population Parameters
    # twoNmu, twoNs & twoNm: the population scaled by the total population.
    mu = twoNmu/2/NDeme/nDemes # mutation coefficient

    #commented out as s0 and mig directly inputed
    # s0 = twoNs/2/NDeme/nDemes #selection coefficient
    # mig = twoNm/2/NDeme/nDemes # dispersal rate

    ## Preallocating dataframes
    migCell=Array{Array{Float64,1},1}(undef,nDemes) # stores how many migrants per deme
    migLR=Array{Array{Float64,2},1}(undef,nDemes) #left right migration amount
    mutantCount = zeros(nDemes,T+1) # counts number of de novo mutants generated
     #per deme for each time
    migrantCount = zeros(nDemes,T+1) # counts number of migrants generated per deme
     #for each time

    ## Creating the initial frequency array, with each array representing one deme
    D=Array{Array{Float64,2},1}(undef,nDemes)#deme cell array - each with their
     #own frequency array
    x = zeros(1,T+1)
    x[end,1] = 1 #Initial frequency array of the wt (no mutants present yet)

    for i in 1:nDemes # populate deme cell array with initial wt-frequency of 1
        D[i] = x
    end

 for t in 2:T+1 #Generations #T+1
    # println("generation######################",t)
     for i in 1:nDemes
        # println("deme################", i)

        ##Selection
        ## Filling up the frequency array
        x=D[i]
        ##mutation, selection, genetic drift
        xwt = x[end,t-1] #wild type frequency at t-1
        X = x[1:end-1,t-1] # mutant frequencies at t-1
        ss=xwt*s0 #selection modification
        #as selection coefficient are assumed to be same for all mutants
          #ss can be simplified to this

        xx=X+ss*X #adjusted probability from frequency to be sampled in
         #genetic drift

        #correction
        xx[xx.<0].=0
        xxK = 1-sum(xx) #Frequency of wild type (wt is the Kth allele)
        if xxK < 0 # in case wild type frequency is negative
            xxK = 0 # set wt frequency to 0
            xx = xx/sum(xx) #scale all mutant frequencies to sum to 1
        end

        ## Genetic drift （Gaussian approximation of multinomial sampling）
        if Gauss == 1 #true
            z=gaussian_multinomial(NDeme,[xx;xxK]) #z in frequency
            n_=round.(Int,z*NDeme) #n in individuals, rounded to be integer

        elseif Gauss == 0 #false
            n_ = rand(Multinomial(Int(NDeme),[xx;xxK])) #Multinomial random numbers
        end

        x[:,t] = n_ #n outputed in individuals in the column of this generation
        n=x #represent matrix with only using columnn in individuals
        D[i]=n

        # Mutation
        nwt = n[end,t] #wild type individuals in t generation
        m = pois_rand(mu*nwt) #number of de novo mutants generated
        mutantCount[i,t] = m #stores number of mutants generated into the array
         #i deme, t generation

        if m != 0
            #correction as cannnot generate more mutant than wild type individuals
            if  m > nwt
                m=nwt
            end

            #the matrix created for the added freuqnecy (1/N) in this generation
             #for new mutants
            xN = zeros(m,T+1) #rows of new mutant alleles
            for f in 1:nDemes #add xN to all demes
                x = D[f]
                x=ApplyArray(vcat,xN,x) #new mutant piled on top
                if f==i #frequency of the deme that mutant emerged as 1/NDeme
                    x[1:m,t] .=1 # 1 individual arising per mutant seeding into above array
                    x[end,t] = x[end,t] - m
                end
                D[f]=x #store back
            end
        end

    end # for i in 1:nDemes


    ##Migration
    ## Determining the number of migrants leaving the demes
    for i in 1:nDemes
        n_i = D[i] #extract from cell array, multi-step=>better computational efficiency
        n_B = n_i[:,t] #frequency column in t generation, i-th deme
        MIG = mig*n_B #expected number of migrants individuals #column vec
        migLV = pois_rand.(MIG)
            #vector of total migrants of each haplotyepes outside this deme

        for g in 1:length(n_B) # poisson may generate>there is=>negative frequency
            if migLV[g] > n_B[g] # looks for negative values and sets to zero
                migLV[g] = Int(n_B[g]) #both float vector
            end
        end

        n_lv = n_B - migLV #after leaving; both column vector
        n_i[:,t] = n_lv # storing population after migrants leave in
        D[i] = n_i # stores "frequency" array back in to deme cell array

        migCell[i] = migLV # leaving array! not left-right array (migLR)
        migrantCount[i,t] = sum(migLV) #sum of migrants per deme
    end

        ## Determining if migrants move left or right
    for i in 1:nDemes
        c = size(migCell[i])[1] #c=column number, registers size of migrant
        #array to preallocate size of migLR_i array
        migLR_i = Int.(zeros(c,2)) #column 1 amount going left, 2 going right
        tot_mig = migCell[i] #total: left+right
        tot_mig = Int.(tot_mig)
        if i == 1 #first deme
            migLR_i[:,2] = tot_mig # left end, only toward right
            migLR[i] = migLR_i
        elseif i == nDemes #last deme
            migLR_i[:,1] = tot_mig # right end, only toward left
            migLR[i] = migLR_i
        else # middle demes
            migLR_i[:,1] = rand.(Binomial.(tot_mig,0.5))#toward left
                            #two elemetwise dot is needed
            migLR_i[:,2] = tot_mig - migLR_i[:,1] #toward right
            migLR[i] = migLR_i
        end
    end
#migLR:    column1           column2
#leaving:  toward left       toward right
#entering: fr Right          fr Left

    ## Accounting for migrants entering in the "frequency" array, D
    for i in 1:nDemes
        if i == 1 # for first deme
            n_i = D[i]
            n_B = n_i[:,t] ## extract n before accounting for entering migrants
            n_frR = migLR[i+1] # n entering array !ROW VECTOR
            n_frR = n_frR[:,1] #only migration from right, columnn vector

            n_aft = n_B + n_frR #add entering cases
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array

            elseif i == nDemes # for last deme
            n_i = D[i]
            n_B = n_i[:,t]
            n_frL = migLR[i-1] # n entering array !ROW VECTOR
            n_frL = n_frL[:,2] #only migration from left, columnn vector
            n_aft = n_B + n_frL # adds entering array into existing array
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array
        else # for the middle demes
            n_i = D[i] # n before entering, COLUMN VECTOR
            n_B = n_i[:,t]
            # println("migLR_size",size(migLR))
            n_frR = migLR[i+1] # n entering array from right deme !ROW VECTOR
            n_frR = n_frR[:,1] # n entering (left migration therefore
            n_frL = migLR[i-1] # n entering array from left deme !ROW VECTOR
            n_frL = n_frL[:,2] # n entering (right migration therefore

            n_aft = n_B + n_frR + n_frL # adds entering array into existing array
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array
        end
    end

    #convert individual to frequency
    for i in 1:nDemes
       #not tracking migration individuals
        n_B=D[i][:,t]#access matrix in individuals
        x_aft=n_B/sum(n_B) #transform into frequency matrix
        D[i][:,t]=x_aft #store frequency matrix back to D cell array
     end


  end #for t = 2:T+1

  # Setting all frequencies to sum up to 1 after migration (convert
   #from population into frequency)
  for i in 1:nDemes
      x_i = D[i] # access frequency matrix of i-th deme
      for t in 1:T+1
          x_B = x_i[:,t]

          if sum(x_B)!=1
              x_B=round.(x_B, digits=10)
              #could be due to rounding in line 210
              x_aft = x_B/sum(x_B) #still with this problem
              #pls check Matlab codes in file num2strexp.m in 1D folder

              x_i[:,t]=x_aft #store frequency vector back

          end
      end
      D[i] = x_i # stores back into cell array
  end

  return D
 end #function Gaussian_1d_localmig()
