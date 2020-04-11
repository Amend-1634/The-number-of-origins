using PoissonRandom
using LazyArrays
using Multibreak
using Base.Threads
using Plots
include("2_2gaussian_multinomial.jl")


#compared to 3_1_2notrack_cons.jl, only the population size differ--oscillating
 #periodically in this function


# function gaussian_1d_localmig(NDeme,twoNmu,twoNs,twoNm,T,nDemes, Gauss)
function gaussian_1d_localmig(N, phi, har_geo, twoNmu,s0,mig,T,nDemes,Gauss)#N

    #Input:
    #N: meta-population size
    #phi: peak-to-trough ratio
    #har_geo: "H" as harmonically constrained mean
             #"G" as geometrically constrained mean
    #twoNmu: total mutation cases in this meta-population in each generation
    #s0: selection strength
    #mig: dispersal rate
    #T: time (generations)
    #nDemes: the number of demes
    # Gauss: 1 (using Gaussian approximation) or 0 (multinomial distribution)

    #Output
    #D: cell array containing frequency array in each deme -
    #rows of this array correspond to independent origins/haplotypes,
    #with common ordering across all demes (same num of rows)



    #Nmax, Nmin: maximum and minimum population size
    if har_geo=="H" #harmonically constrained mean
        Nmax=N/2*(phi+1)
        Nmin=N/2*(phi+1)/phi
    elseif har_geo=="G" #geometrically constrained mean
        Nmax=sqrt(phi)*N
        Nmin=N/sqrt(phi)
    end


    # twoNmu, twoNs & twoNm: the population scaled by the total population.
        #to enable input s0 and mig instead of population-scaled parameters
        twoNs=s0*N*2
        twoNmig=mig*N*2
        #constant population parameters later will be scaled by oscillating population size

    ## Preallocating dataframes
    migCell=Array{Array{Float64,1},1}(undef,nDemes) # stores how many migrants per deme
                               #one column is enough
    migLR=Array{Array{Float64,2},1}(undef,nDemes) #left right migration amount
    mutantCount = zeros(nDemes,T+1) # counts number of de novo mutants generated
     #per deme for each time
    migrantCount = zeros(nDemes,T+1) # counts number of migrants generated per deme
     #for each time
    ## Creating the initial frequency array, with each array representing one deme
    D=Array{Array{Float64,2},1}(undef,nDemes)#deme cell array - each with their own frequency array
    x = zeros(1,T+1)
    x[end,1] = 1 #Initial frequency array of the wt (no mutants present yet)

    for i in 1:nDemes # populate deme cell array with initial wt-frequency of 1 @ T= 0
        D[i] = x
    end

    for t in 2:T+1 #Generations #T+1
    # println("generation############",t)
    prd=12 #a period of 12 generations

    #Calculating oscillating population size in each generation
    N=0.5(Nmax+Nmin)+0.5(Nmax-Nmin)*sin(2pi*t/prd) #oscillating population size
    # println("oscillating NDeme ", NDeme)

    ## Population Parameters
    mu = twoNmu/2/N# mutation coefficient
    s0 = twoNs/2/N #selection coefficient, s_ to avoid name crash with s in Julia
    mig = twoNmig/2/N # dispersal rate

    for i in 1:nDemes
        # println("deme################", i)

        ##mutation, selection, genetic drift
        ## Filling up the frequency array
        x=D[i]

        ##Selection
        xwt = x[end,t-1] #wild type frequency at t-1
        X = x[1:end-1,t-1] # mutant frequencies at t-1

        # ss=s_[1:end-1] .- s_'*x[:,t-1] #ss: adjustment parameter for probability
                                         #simplified as below
        #(K-1, 1)            (1,K)   (K,1)
        #(K-1, 1)                1(constant now)
        #(x,y) bracket above is indicating matrix with x row and y columns

        ss=xwt*s0 #selection modification

        xx=X+ss*X

        xx[xx.<0].=0
        xxK = 1-sum(xx) #Frequency of wild type (wt is the Kth allele)

        if xxK < 0 # in case wild type frequency is negative
            xxK = 0 # set wt frequency to 0
            xx = xx/sum(xx) #resets all mutant frequencies to sum to 1
        end

        ## Genetic drift （Gaussian approximation of multinomial sampling）
        if Gauss == 1 #true
            z=gaussian_multinomial(NDeme,[xx;xxK]) #z in frequency
            n_=round.(Int,z*NDeme)#n in individuals, in integer individuals
        elseif Gauss == 0 #false
            NDeme = round(Int, NDeme)
            n_ = rand(Multinomial(NDeme,[xx;xxK])) #Multinomial random numbers
        end
        n=x #represent matrix with only using columnn in individuals
        n[:,t] = n_ #n outputed in individuals in the column of this generation
        D[i]=n

        # Mutation
        nwt = n[end,t] #wild type individuals in t generation
        m = pois_rand(mu*nwt) #number of de novo mutants generated
        mutantCount[i,t] = m #stores number of mutants generated into the array
         #i deme, t generation

        #selection, genetic drift then mutation(order make no difference as new
         #mutant not involved in drift in this generation)
        if m != 0
            if  m > nwt #correction of new mutant allele to be <= wild type individuals
                m=nwt
            end

            xN = zeros(m,T+1) #rows of new mutant alleles
            for f in 1:nDemes #add xN to all demes
                x = D[f]
                x=ApplyArray(vcat,xN,x) #new mutant piled on top (differ from former insertation)
                if f==i #frequency of the deme that mutant emerged as 1/NDeme
                    x[1:m,t] .=1 # 1 individual arising per mutant seeding into above array
                    x[end,t] = x[end,t] - m
                end
                D[f]=x
            end
        end

    end # for i in 1:nDemes


    ##Migration
    ## Determining the number of migrants leaving the demes
    for i in 1:nDemes
        n_i = D[i] #extract from cell array, multi-step=>better computational efficiency
        n_B = n_i[:,t] #frequency column in t generation, i-th deme
        MIG = mig*n_B #expected number of migrants individuals #column vec
        migLV = pois_rand.(MIG) #actual number of migrants
        #size: (allele_amount,)

                          #Before leaving in migration
                          #column vector
        for g in 1:length(n_B) # poisson may generate>there is=>negative frequency
            if migLV[g] > n_B[g] # looks for negative values and sets to zero
                migLV[g] = Int(n_B[g]) #both float vector
            end
        end

        n_lv = n_B - migLV #after leaving; both column vector
        n_i[:,t] = n_lv # storing population after migrants leave in
        D[i] = n_i # stores "frequency" array back in to deme cell array

        #recall "migCell[i] = poissrnd(MIG)' #actual number of migrants"--line111
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
            n_i = D[i] # n before accounting for entering migrants
            n_B = n_i[:,t]
            n_frR = migLR[i+1] # n entering array !ROW VECTOR
            n_frR = n_frR[:,1] #only migration from right, columnn vector
            #column vector

            n_aft = n_B + n_frR
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
            n_frR = migLR[i+1] # n entering array from right deme !ROW VECTOR
            n_frR = n_frR[:,1] # n entering (only left migration therefore                #row 1)!COLUMN VECTOR
            n_frL = migLR[i-1] # n entering array from left deme !ROW VECTOR
            n_frL = n_frL[:,2] # n entering (only right migration therefore
                #row 2)!COLUMN VECTOR

            n_aft = n_B + n_frR + n_frL # adds entering array into existing array
            n_i[:,t]=n_aft
            D[i] = n_i # stores back into deme cell array
        end
    end

    #convert individual to frequency
    for i in 1:nDemes
       #not tracking migration individuals
        n_i=D[i]
        n_B=n_i[:,t]#access matrix in individuals
        x_aft=n_B/sum(n_B) #transform into frequency matrix
        n_i[:,t]=x_aft #contain this column as in frequency
                        #lazy in renaming--extra computational time
        D[i]=n_i #store frequency matrix back to D cell array
     end


  end #for t = 2:T+1

  ## Setting all frequencies to sum up to 1 after migration (convert
  #from population into frequency)
  # for i in 1:nDemes
  #     x_i = D[i] # access frequency matrix of i-th deme
  #     for t in 1:T+1
  #         x_B = x_i[:,t]
  #
  #         if sum(x_B)!=1
  #             # println("need zoom at end, ", sum(x_B))
  #             println("x_B", x_B)
  #
  #             #could be due to rounding in line 266
  #             x_aft = x_B/sum(x_B)
  #             x_i[:,t]=x_aft #store frequency vector back
  #             println("x_after_final scaling",x_aft)
  #
  #         end
  #     end
  #     D[i] = x_i # stores back into cell array
  # end
 #the final scaling: why always print out the sum in level of e102 or e46...?
 #problem of sum(n_B)?
  #could be rounding issue
  #pls check Matlab codes in file num2strexp.m in 1D folder

  return D
 end #function Gaussian_1d_localmig()
