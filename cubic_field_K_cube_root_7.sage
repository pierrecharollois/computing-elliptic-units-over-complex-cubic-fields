### This file is loaded by 'BCG_cubic.sage'. It contains the data for a complex cubic field obtained by adjoining to QQ a cubic root of 7. It loads the following GLOBAL VARIABLES:

###  y, pol_K, yRR, yCC, K, K_basis_SAGE, K_basis_SAGE_sigmaC, ff, aa, bb, q_ff, N_aa, aa_comp, vareps, varepsRR, varepsCC, or_matrix, Eps
### TO DO: list of zeta derivatives of conductor ff.

###### ------ NUMBER FIELD DATA -------  #######

y = PolynomialRing(QQ, 'y').gen()
pol_K = y^3 - 7

# The real and complex roots of pol_K
yRR = pol_K.roots(ring=CC, multiplicities=False)[0]
yRR = real(yRR)
yCC = pol_K.roots(ring=CC, multiplicities=False)[1] # This is the complex root of pol_K with positive imaginary part. We denote the associated embedding by sigmaC


# The number field K and two integral bases, provided by SAGE and PARI

K.<y> = NumberField(pol_K) # This is one of the main examples in our paper
#K_basis = K.pari_zk()  # This is the basis for the ring of integers that PARI gives.
#K_basis_sigmaC = [**, **,**] # This is the embedding of the PARI basis under sigmaC
K_basis_SAGE = K.integral_basis() # This is the basis that SAGE gives: it is [1,y,y^2].
K_basis_SAGE_sigmaC = [1, yCC, yCC^2] # This is the embedding the SAGE basis under sigmaC




# Ideals in K and some of the associated invariants appearing in the paper.

ff = K.ideal(3).factor()[0][0] # ff is the conductor ideal
aa = K.ideal(5).factor()[1][0] # aa is the smoothing ideal
bb = K.ideal(2).factor()[1][0] # bb is the ideal bb
q_ff = 3 # Recall that q is the smallest positive integer in ff, which in this case is 3.
N_aa = aa.norm() # We use N_aa to denote the norm of aa, which in this case is 5.
aa_comp = N_aa * aa^(-1) # This ideal will be useful to compute admissible elements.

# The fundamental unit and its three by three matrix

vareps = 2-y # This is the fundamental unit. Note that vareps-1 belongs to ff
varepsRR = 2-yRR # The image of the fundamental unit in RR. Note in lands in (0,1).
varepsCC = 2-yCC


# This matrix has purely imaginary determinant with NEGATIVE imginary part. Note that if we switch to the other complex embedding we swap the orientation!
or_matrix = Matrix([[1,1,1],[yRR,yCC,conjugate(yCC)],[yRR^2,yCC^2,conjugate(yCC)^2]])



Eps = Matrix([vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[0])),vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[1])),vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[2]))])

# IMPORTANT NOTE: the convention in this file is that vectors representing elements of K are ROW VECTORS. That is: (1,0,0) represents K_basis_SAGE[0]=1, (0,1,0) represents K_basis_SAGE[1]=y and (0,0,1) represents K_basis_SAGE[2]=y^2. This means that the map x mapsto vareps x (i.e. multiplication by vareps) corresponds to the map v mapsto v * Eps. I have tested this in examples.
# The only time when we use column vectors is when we represent functionals on K (as a \Q-vector space or as a Z-module), for example the complex embedding sigmaC or the vectors appearing in the wedge W, called a and b. This is also why when we act on a by the unit vareps we write the transpose.

print("Loaded 'cubic_field_K_cube_root_7.sage'. Contains the following variables: ")
print()
print('K:', K)
print('K has discriminant ', K.disc(), '=', factor(K.disc()))
print('K has class number ', K.class_number())
print()
print('pol_K = defining polynomial = ' , pol_K)
print('K_basis_SAGE = basis of K = ', K_basis_SAGE)
print()
print('yRR = unique real root of pol_K  = ' , numerical_approx(yRR),'...')
print('yCC = unique complex root of pol_K with positive imaginary part ')
print('    = ' , numerical_approx(real(yCC)),'... + I * ',numerical_approx(imag(yCC)),'...')
print('sigmaC = complex embedding of K sending y to yCC')
print('K_basis_SAGE_sigmaC = [1, yCC, yCC^2]')
print('K_basis_SAGE_sigmaR = [1, yRR, yRR^2]')
print()
print('ff = ', ff , ' of norm ' , ff.norm())
print('q_ff = smallest positive integer in ff =' , q_ff)
print('The narrow ray class group of conductor ff has order 6.'  )
print()
print('aa = ', aa )
print('N_aa = norm of aa = ' , aa.norm())
print('aa_comp = N_aa * aa^(-1)')
print('bb = ', bb , ' of norm ' , bb.norm() )
print('(bb is a generator of the narrow ray class group modulo ff.)' )
print()
print('vareps = generator of totally positive units in 1+ff = ', vareps)
print('varepsRR = real image of vareps = ', numerical_approx(varepsRR))
print('varepsCC = image of vareps under sigmaC ')
print('or_matrix = Matrix with entries the real and complex images of K_basis_SAGE')
print('Eps =  matrix of the fundamental unit vareps in the basis K_basis_SAGE = ')
print(Eps)







##### -----------------------------------------------------    ######