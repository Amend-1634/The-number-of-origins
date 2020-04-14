N=1e5;nDemes=5;twoNmu=10;s0=0.1;mig=0.5;M=3;T=1000;Gauss=1;nu=1

include("4_1_1nonlocal_frequency.jl")
D, X_D, migrantCount=gaussian_1d_nonlocalmig(N,twoNmu,s0,mig,T,nDemes, Gauss,nu)
plot(D[1]')
plot!(X_D[1]')

N=1e5;nDemes=5;twoNmu=10;s0=0.1;mig=0.5;M=3;T=1000;Gauss=0;nu=1
pwd()
cd("D:\\Julia\\func")
include("4_2_1d_tsamp.jl")
avn, se, tsamp, medn, lqrn, uqrn, X_avn, X_se, X_medn, X_lqrn, X_uqrn=num_origin_1d_t(N,nDemes,twoNmu,s0,mig,M,T,Gauss,nu)
using Plots
plot(tsamp,avn)


N=1e5;nDemes=5;twoNmu=10;s0=0.1;mig=0.5;M=3;T=1000;Gauss=0;nu=1
include("4_3_1d_fixsamp.jl")
t_avn, t_se, eta_avn, eta_se, fix_M=num_origin_1d_fix(N,nDemes,twoNmu,s0,mig,M,
  T,Gauss,nu)
eta_avn
t_avn
