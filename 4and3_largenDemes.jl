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
