using Plots
using Distributions #MvNormal
using StatsBase #isposdef
using LinearAlgebra #eigen, Hermitian
include("2_1allele_var_mtx.jl")

function gaussian_multinomial(N,xx)
  #N: population size for the whole panmictic population or each deme in meta-population
  #xx: frequecy vector with all allele
  # sum(xx)=1?
  beta=3/2
  #All frequencies less than beta threshold are zero - this is empirical
  #?do you mean cannot survive afterward? or mathmetically zero (<1/N)?
  # println("na in xx_input", findall(isnan,xx))

  xx[xx.<1/beta/N].=0
  #Reduce dimensionality by removing mean frequencies=0 in xx
  indnot0 = findall(!isequal(0),xx)
  # println("indnot0 ",indnot0)
  z = zeros(size(xx))#float

  if length(indnot0)==1
      z[indnot0]=xx[indnot0] #contains allele frequency >1/beta/N  and =0
   elseif !isempty(indnot0) #include length(indnot0)==1??

    #if there are n non-zero elements take only the n-1 elements to form xx,
     #on an effective simplex of size n
    indnot0_end = indnot0[end]
    indnot0 = indnot0[1:end-1]
    #why the last not-0 one element in z so special?

    xnz = xx[indnot0]
    #Draw random number from multivariate normal distirbution with
     #variance matrix of allele frequencies given by B
    # println("xnz ",xnz)
    # println("N ",N)

    B = allele_var_mtx(xnz)/N #AlleleVarianceMatrixWF(xnz,numel(xnz);

    # println("isposdef before transform_",isposdef(B))
    # println("issymmetric before transform_",issymmetric(B))#false #so not Hermitian

#transform block to make cov matrix positive definite
##trial1: cannot output Hermitian, so try another one
     # using StatsBase #eigen()
     # println("typeof(B)",typeof(B))
     # println("size(B)",size(B))
     # println("na in B", findall(isnan,B))
     # println("inf in B", findall(isinf,B))

     D, V=eigen(B)
      #eigenvalues in D(n,) and eigenvectors in V(n,n)
     # ##test if any eigenvalues (D) <0 if yes then err=false
     if isempty(D[D.<0]) == false  #err=false, not positive semidefinite
     #  #if err==false #not positive semi-definite(or definite?)

         D[D.<0].=0 #set to be semidefinite
          # using LinearAlgebra
         D=Diagonal(D)# Construct a matrix from the diagonal of A.
          #(n,) => (n,n)
         B = V*D*inv(V)
          #(n,n)*(n,n)*(n,n) ouput (n,n) covariance matrix

         # err=isposdef(B) ;#
         # while err==false ;#
         D, V=eigen(B)
         possemi=isempty(D[D.<0])#whether positive semidefinite matrix
         while possemi == false #if not

               nx = size(D)[1]
               D=Diagonal(D)
               D = D + 1e-10*Diagonal(randn(nx,nx))#normal-distributed
               D[D.<0].=0 #only give semidefinite
               # D[D.<1e-7].=1e-7 #definite (Matlab answer)
               B = V*D*inv(V) #transform into positive semidefinite
               D, V=eigen(B)

               possemi=isempty(D[D.<0]) #return to this while loop until not satisfied
          end
     end

     if issymmetric(B)==false
      B=1/2*(B+B')#make it Hermitian (diagonal-symmetric)
     end
     # println("isposdef after transform_",isposdef(B))
     # println("issymmetric after transform_",issymmetric(B))#false #so not Hermitian

##trail2: of transformation into positive definite
    #output issymmetric covariance matrix as well

    # D, V=eigen(B)
    # posdef=isempty(D[D.<0])#whether positive semidefinite
    # while posdef == false  #err=false, not positive semidefinite
    #   p=minimum(D[D.>0]) #lowest positive eigenvalues
    #   D_ndx=findall(a->a<=0,D)#I amended it to <=, "negative" originally
    #   s=sum(D[D_ndx])*2
    #   t=s^2*100+1
    #   for i in D_ndx
    #     D[i]=p*(s-D[i])^2/t
    #    end
    #   D=Diagonal(D)
    #   B = V*D*inv(V)
    #   D, V=eigen(B)
    #   posdef=isempty(D[D.<0])
    #   println("trans_while")
    # end
    # #main problem: positive definite but not issymmetric any more
    # println("issymmetric ",issymmetric(B))#false #so not Hermitian
    # println("posdef", posdef)
    # # println(B)
##
    z[indnot0] = rand(MvNormal(xnz,B)) #required to be positive definite, and Hermitian
    #output as frequency
    # ERROR: PosDefException: matrix is not positive definite;
     # Cholesky factorization failed.

    #B as covariance matrix,xnz as mean vector

    # indnegzero = findall(z<1/beta/N,z);
    # z[indnegzero] .= 0;
    # println("z[indnot0] ",z[indnot0])
    z[z.<1/beta/N].=0 #simpler

    zz = 1-sum(z)

    if zz < 1/beta/N
        zz = 0
    end

    z[indnot0_end] = zz
    z = z/sum(z)#zooming to sum=1

   end#if length(indnot0)==1
  return z #frequency
 end#function Gaussian_multinomial()
