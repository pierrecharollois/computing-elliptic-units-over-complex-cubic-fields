### This file is loaded after 'BCG_cubic.sage'. It contains the data for the complex cubic field Q(beta) with beta^3 = 2,
### of discriminant -108, with conductor f = (5). It loads the following GLOBAL VARIABLES:

###  y, pol_K, yRR, yCC, K, K_basis_SAGE, K_basis_SAGE_sigmaC, K_basis_SAGE_sigmaR, ff, aa, bb, q_ff, N_aa, aa_comp, vareps, varepsRR, varepsCC, or_matrix, Eps

###### ------ NUMBER FIELD DATA -------  #######

y = PolynomialRing(QQ, 'y').gen()
pol_K = y^3 - 2

# The real root and the complex root of positive imaginary part.
_rts = pol_K.roots(ring=CC, multiplicities=False)
yRR = [r for r in _rts if not (r.imag() > 0) and not (r.imag() < 0)][0].real()
yCC = [r for r in _rts if r.imag() > 0][0]


# The number field K and the integral basis.
K.<y> = NumberField(pol_K)
K_basis_SAGE = K.integral_basis()
K_basis_SAGE_sigmaC = [1, yCC, yCC^2]
K_basis_SAGE_sigmaR = [1, yRR, yRR^2]



# Ideals and invariants. f = (5): the order of beta-1 modulo 5 is 8 and no power of beta-1 is -1
# modulo 5, so that Cl^+_K(5) = Z/12 while Cl_K(5) = Z/6. The class of bb = (beta) generates
# Cl^+_K(5), and aa = (11, beta+4) is the smallest prime coprime to 5 that can serve, 7, 13 and 19
# being inert.

ff = K.ideal(5)                       # conductor, f = (5)
q_ff = 5                              # smallest positive integer in ff

aa = K.ideal(11, y+4)                 # smoothing ideal, prime of norm 11
N_aa = aa.norm()                      # 11
aa_comp = N_aa * aa^(-1)

bb = K.ideal(y)                       # bb = (beta), of norm 2; its class generates Cl^+_K(f) = Z/12

# The unit: (beta-1)^8 is the smallest power of the fundamental unit beta-1 lying in 1 + ff.
vareps = (y-1)^8                      # = 1 + 100*beta - 80*beta^2
varepsRR = (yRR-1)^8
varepsCC = (yCC-1)^8


# Orientation matrix: real and complex images of the basis. Purely imaginary det.
or_matrix = Matrix([[1,1,1],[yRR,yCC,conjugate(yCC)],[yRR^2,yCC^2,conjugate(yCC)^2]])


Eps = Matrix([vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[0])),vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[1])),vector(K.ring_of_integers().coordinates(vareps*K_basis_SAGE[2]))])

# NOTE: vectors representing elements of K are ROW VECTORS. Multiplication
# by vareps corresponds to v mapsto v * Eps.

# ---- consistency checks (printed below) ----
chk_disc  = (K.disc() == -108)
chk_class = (K.class_number() == 1)
chk_norm  = (N_aa == 11 and ZZ(N_aa).is_prime())
chk_unit0 = bool(0 < varepsRR < 1)
chk_unitf = ((vareps - 1) in ff)
chk_eps   = (Eps.det().abs() == 1)

print("Loaded 'cubic_field_K_cube_root_2_f5.sage' (disc -108, f = (5)).")
print()
print('K:', K)
print('K has discriminant ', K.disc(), '=', factor(K.disc()))
print('K has class number ', K.class_number())
print()
print('pol_K = ' , pol_K)
print('K_basis_SAGE = ', K_basis_SAGE)
print()
print('yRR = ' , numerical_approx(yRR),'...')
print('yCC = ' , numerical_approx(real(yCC)),'... + I * ',numerical_approx(imag(yCC)),'...')
print('sigmaC = complex embedding of K sending y to yCC')
print()
print('ff = ', ff , ' of norm ' , ff.norm())
print('q_ff =' , q_ff)
print('aa = ', aa , ' of norm ', aa.norm())
print('N_aa = ' , N_aa)
print('bb = ', bb , ' of norm ' , bb.norm(), ' (class generates Cl^+_K(f) = Z/12)')
print()
print('vareps = (beta-1)^8 = ', vareps)
print('varepsRR = ', RealField(53)(varepsRR))
print('Eps =')
print(Eps)
print()
print('--- consistency checks ---')
print('disc = -108           :', chk_disc)
print('class number = 1      :', chk_class)
print('N(aa) = 11 prime      :', chk_norm)
print('sigma_R(vareps) in (0,1):', chk_unit0)
print('vareps - 1 in ff      :', chk_unitf)
print('det(Eps) = +-1        :', chk_eps)
if not (chk_disc and chk_class and chk_norm and chk_unit0 and chk_unitf and chk_eps):
    print('*** WARNING: a consistency check FAILED, do not trust results ***')
print('--------------------------')


##### -----------------------------------------------------    ######
