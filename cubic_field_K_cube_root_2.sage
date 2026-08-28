### This file is loaded by 'BCG_cubic.sage'. It contains the data for the complex cubic field Q(beta) with beta^3 = 2. This is the example of Section 5.3 of the paper (following Ren-Sczech). It loads the following GLOBAL VARIABLES:

###  y, pol_K, yRR, yCC, K, K_basis_SAGE, K_basis_SAGE_sigmaC, K_basis_SAGE_sigmaR, ff, aa, bb, q_ff, N_aa, aa_comp, vareps, varepsRR, varepsCC, or_matrix, Eps

###### ------ NUMBER FIELD DATA -------  #######

y = PolynomialRing(QQ, 'y').gen()
pol_K = y^3 - 2

# The real and complex roots of pol_K
yRR = pol_K.roots(ring=CC, multiplicities=False)[0]
yRR = real(yRR)
yCC = pol_K.roots(ring=CC, multiplicities=False)[1] # This is the complex root of pol_K with positive imaginary part. We denote the associated embedding by sigmaC


# The number field K and two integral bases, provided by SAGE and PARI

K.<y> = NumberField(pol_K)
K_basis_SAGE = K.integral_basis()
K_basis_SAGE_sigmaC = [1, yCC, yCC^2]
K_basis_SAGE_sigmaR = [1, yRR, yRR^2]



# Ideals in K and some of the associated invariants appearing in the paper.
# Section 5.3: f = (3), narrow ray class field of relative degree 6, aa of norm 5,
# and b = (beta) generates the ray class group.

ff = K.ideal(3)                       # ff is the conductor ideal, f = (3)
aa = [P for P,e in K.ideal(5).factor() if P.norm()==5][0]  # aa is the smoothing ideal, the prime of norm 5
bb = K.ideal(y)                       # bb = (beta) generates the ray class group
q_ff = 3                              # q is the smallest positive integer in ff
N_aa = aa.norm()                      # norm of aa, which is 5
aa_comp = N_aa * aa^(-1)              # useful to compute admissible elements

# The fundamental unit and its three by three matrix.
# The paper states that (beta-1)^3 is a generator of O^{+,x}_f. We use it directly.
vareps = (y-1)^3                      # generator of totally positive units in 1+ff
varepsRR = (yRR-1)^3
varepsCC = (yCC-1)^3


# This matrix has purely imaginary determinant. Note that switching the complex
# embedding swaps the orientation.
or_matrix = Matrix([[1,1,1],[yRR,yCC,conjugate(yCC)],[yRR^2,yCC^2,conjugate(yCC)^2]])


Eps = Matrix([vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[0])),vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[1])),vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[2]))])

# IMPORTANT NOTE: vectors representing elements of K are ROW VECTORS. Multiplication
# by vareps corresponds to v mapsto v * Eps.

print("Loaded 'cubic_field_K_cube_root_2.sage'. Contains the following variables: ")
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
print()
print('ff = ', ff , ' of norm ' , ff.norm())
print('q_ff = smallest positive integer in ff =' , q_ff)
print()
print('aa = ', aa , ' of norm ', aa.norm())
print('N_aa = norm of aa = ' , aa.norm())
print('aa_comp = N_aa * aa^(-1)')
print('bb = ', bb , ' of norm ' , bb.norm())
print('(bb = (beta) generates the ray class group modulo ff.)')
print()
print('vareps = generator of totally positive units in 1+ff = (beta-1)^3 = ', vareps)
print('varepsRR = real image of vareps = ', numerical_approx(varepsRR))
print('  (should lie in (0,1):', bool(0 < varepsRR < 1), ')')
print('  (vareps - 1 in ff:', (vareps - 1) in ff, ')')
print('or_matrix = Matrix with entries the real and complex images of K_basis_SAGE')
print('Eps =  matrix of the fundamental unit vareps in the basis K_basis_SAGE = ')
print(Eps)


##### -----------------------------------------------------    ######
