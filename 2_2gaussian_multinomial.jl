using Plots
using Distributions #MvNormal
using StatsBase #isposdef
using LinearAlgebra #eigen, Hermitian
include("2_1allele_var_mtx.jl")

function gaussian_multinomial(N,xx)
  #N: population size for the whole panmictic population or each deme in meta-population
  #xx: frequecy vector with all allele
  # sum(xx)=1

  beta=3/2
  #All frequencies less than beta threshold are zero - this is empirical
  #?do you mean cannot survive afterward? or mathmetically zero (<1/N)?

  xx[xx.<1/beta/N].=0
  #Reduce dimensionality by removing mean frequencies=0 in xx
  indnot0 = findall(!isequal(0),xx)
  z = zeros(size(xx))#float

  if length(indnot0)==1
      z[indnot0]=xx[indnot0] #contains allele frequency >1/beta/N  and =0
   elseif !isempty(indnot0)

    #if there are n non-zero elements take only the n-1 elements to form xx,
     #on an effective simplex of size n
    indnot0_end = indnot0[end]
    indnot0 = indnot0[1:end-1]
    #why the last not-0 one element in z so special?

    xnz = xx[indnot0]
    #Draw random number from multivariate normal distirbution with
     #variance matrix of allele frequencies given by B

    B = allele_var_mtx(xnz)/N

    # println("isposdef before transform_",isposdef(B))
    # println("issymmetric before transform_",issymmetric(B))#false #so not Hermitian

#transform block to make cov matrix positive definite
     # using StatsBase #eigen()
     D, V=eigen(B)
     #eigenvalues in D(n,) and eigenvectors in V(n,n)

     D=collect(Diagonal(D))# Construct a matrix from the diagonal of D.
                           #collect() make (n,)Diagonal array to (n,n)normal array

      #eigenvalues in D(n,) and eigenvectors in V(n,n)
     # ##test if any eigenvalues (D) <0 if yes then err=false
     # if isempty(D[D.<0]) == false  #err=false, not positive semidefinite
     if isposdef(D)==false #if D is not positive definite matrix
     #  #if err==false #not positive semi-definite(or definite?)

         # D[D.<0].=0 #set to be semidefinite
         D[D.<1e-7].=1e-7 #definite (Matlab defines)

          # using LinearAlgebra
         D=Diagonal(D)# Construct a matrix from the diagonal of D.
          #(n,) => (n,n)
         B = V*D*inv(V)
          #(n,n)*(n,n)*(n,n) ouput (n,n) covariance matrix
          #(n,n) denote a matrix with n rows and n columns in comment

         D, V=eigen(B)
         # println("D")
         # println(D)
         D=Diagonal(D)
         posdef=isposdef(D)
         #possemi = isempty(D[D.<=0])#whether positive semidefinite(>=0) matrix
                                  #in Matlab codes only require the positive
                                  #semidefinite but due to MvNormal() input
                                  #requiremet, here we set positive definite

         while posdef == false   #if not positive definite matrix

               nx = size(D)[1]
               D=collect(Diagonal(D))
               D = D + 1e-10*Diagonal(randn(nx,nx))
                #randn():normally-distributed sampling
               println("type_size_D",typeof(D),"  ", size(D))
               # D[D.<0].=0 #only give semidefinite(>=0)
               D[D.<1e-7].=1e-7 #definite (Matlab defines)
               B = V*D*inv(V) #transform into positive semidefinite
               D, V=eigen(B)

               posdef=isposdef(D)
               #iterate within this while loop until the condition not satisfied
          end
     end

     #transfomration into symmetric matrix
     if issymmetric(B)==false
      B=1/2*(B+B')#make it Hermitian (diagonal-symmetric)
     end

     # println("isposdef after transform_",isposdef(B))
     # println("issymmetric after transform_",issymmetric(B))#false #so not Hermitian

    z[indnot0] = rand(MvNormal(xnz,B))
    #B as covariance matrix, xnz as mean vector
    #input B required to be positive definite, and diagonal-symmetric (Hermitian)
    #output frequency

    z[z.<1/beta/N].=0 #all element of z lower than 1/beta/N set to 0

    zz = 1-sum(z)

    if zz < 1/beta/N
        zz = 0
    end

    z[indnot0_end] = zz
    z = z/sum(z)#scaling to sum up to 1

   end#if length(indnot0)==1
  return z #frequency
 end#function Gaussian_multinomial()
