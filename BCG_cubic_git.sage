

####### GLOBAL VARIABLES ######

# prec: Number of digits of precision
# CC: the complex field with prec digits of precision.
# Pi:  the number pi in the field CC.
# N_Ell_Gamma: Number of terms used to compute Elliptic Gamma

prec = 3400
N_Ell_Gamma = 10000 # Number of terms used to compute Elliptic Gamma
show_prec = 40      # bits used when PRINTING tau, sigma and the
                    # final value.  40 bits is about 12 digits
                    # It affects the display only, never the computation.

CC=ComplexBallField(prec)
RBF=RealBallField(prec)   #  needed by the certified tail bound below
Pi=CC(pi,0)

#------------------------------------------------------------------------------




###### ------ LOAD NUMBER FIELD DATA -------  #######

#load('cubic_field_K_disc_-588.sage')
#load ('cubic_field_K_cube_root_3.sage')
#load ('cubic_field_K_cube_root_7.sage')
#load ('cubic_field_K_cube_root_7_conj_sigmaC.sage')
#load ('cubic_field_K_cube_root_2.sage')
#load('cubic_field_K_disc_-23.sage')
#load('cubic_field_K_disc_-411.sage')

##### -----------------------------------------------------    ######


def check_hypotheses():
#  Checks the hypotheses that the field data files set by hand, and prints a line for
# each one that fails. Provides only a warning.  

   print()
   if N_aa != 2 and ZZ(N_aa).is_prime():
      print('N(aa) =', N_aa, 'is an odd prime, as (2.9) asks.')
   else:
      print('PROBLEM: N(aa) =', N_aa, 'should be an odd prime, see (2.9).')
   if gcd(ZZ(N_aa), ZZ(ff.norm())) == 1:
      print('aa is coprime to ff, as (2.9) asks.')
   else:
      print('PROBLEM: aa should be coprime to ff, see (2.9).')
   if 0 < varepsRR < 1:
      print('sigmaR(vareps) =', RealField(20)(varepsRR), 'lies in (0,1), as section 2.2 asks.')
   else:
      print('PROBLEM: sigmaR(vareps) =', RealField(20)(varepsRR), 'should lie in (0,1).')
      print('         Try the inverse unit.')
   if (vareps - 1) in ff:
      print('vareps lies in 1 + ff, as section 2.2 asks.')
   else:
      print('PROBLEM: vareps should lie in 1 + ff, see section 2.2.')
   if q_ff > 0 and q_ff in ff:
      print('q_ff =', q_ff, 'is a positive integer of ff, as section 2.2 asks.')
   else:
      print('PROBLEM: q_ff =', q_ff, 'should be a positive integer of ff.')
   print()


check_hypotheses()



### ROUTINES FOR WORKING WITH NUMBER FIELDS ###


def ideal_basis(I):
# This function computes an integral basis of the underlying Z-module to an ideal I of K. The output is a 3 by 3 matrix. The entries in each row are the coordinates of a basis of I with respect to the basis K_basis_SAGE 

   I0 = I.integral_basis()[0]
   I1 = I.integral_basis()[1]
   I2 = I.integral_basis()[2]
   v0 = K.ring_of_integers().coordinates(I0)
   v1 = K.ring_of_integers().coordinates(I1)
   v2 = K.ring_of_integers().coordinates(I2)
   M = matrix([v0, v1, v2])
   return  M
   
   
def sigmaR(v):
# This encodes the real embedding of K. The input v is a row vector whose entries are the coefficients of an element of K in the basis K_basis.
   s = v[0]*K_basis_SAGE_sigmaR[0] + v[1]*K_basis_SAGE_sigmaR[1]+v[2]*K_basis_SAGE_sigmaR[2]
   return s
   
def sigmaC(v):
# This encodes our fixed complex embedding. The input v is a row vector whose entries are the coefficients of an element of K in the basis K_basis.
   s = v[0]*K_basis_SAGE_sigmaC[0] + v[1]*K_basis_SAGE_sigmaC[1]+v[2]*K_basis_SAGE_sigmaC[2]
   return s
   


###### ROUTINES PRODUCING A LIST OF ADMISSIBLE ELEMENTS    ######## 
   
def insert_top_row(row,M):
    return matrix([row]+M.rows()[0:])   

   
def admissible_elements(N,L):

# this function generates a list of admissible elements in the lattice L wrt the ideal aa and the parameter q.
# The two conditions of section 2.2 are: lambda/q = 1_K mod L, and lambda/N(aa)
# generates aa^(-1)L/L. The second is tested below as primitivity in both L and
# aa_comp*L. A weaker test, lambda in aa_comp*L but not in N(aa)*L, only says
# that lambda/N(aa) is nonzero in the quotient. The two agree exactly when
# N(aa) is prime, since aa^(-1)L/L is then cyclic of prime order and every
# nonzero element generates it; hypothesis (2.9) of the paper requires N(aa) to
# be an odd prime, so this is the case here. In short, we are looking for the vectors lambda
# in q*1_K + qL that are primitive both in L and in aa_compL, a condition equivalent to lambda
# being admissible.

   qL = q_ff * L
   aa_compL = aa_comp * L
   M_L = ideal_basis(L)
   M_qL = ideal_basis(qL)
   M_aa_compL = ideal_basis(aa_compL)

   lambda0=vector([q_ff,0,0])
   list=Matrix([0,0,0])
   
   for j in [-N..N]:
      for k in [-N..N]:
         for l in [-N..N]:            
            candidate=lambda0 + j * M_qL[0] + k * M_qL[1] + l * M_qL[2]
            v=candidate*M_aa_compL^(-1)
            w=candidate*M_L^(-1)
            if (v[0].is_integer()) and (v[1].is_integer()) and (v[2].is_integer()) and (w[0].is_integer()) and (w[1].is_integer()) and (w[2].is_integer()):
               if (gcd(gcd(v[0],v[1]),v[2])==1) and (gcd(gcd(w[0],w[1]),w[2])==1):
                  list=insert_top_row(candidate,list)
   nrows = list.dimensions()[0]
   
   # Remove the last row ( meaningless [0,0,0])
   return list[0:nrows-1]

   
######## --------------------------------------------------------- ##########





########  ROUTINES TO HANDLE ADMISSIBLE VECTORS AND WEDGES ########


def primitive(lamb,L):
# This takes a row vector with rational entries (corresponding to the coordinates of an element of K in the basis K_basis_SAGE), and rescales by a unique positive rational number to make sure it belongs to the fractional ideal L.
   M=ideal_basis(L)
   coord_L_lamb = lamb*M^(-1)
   rho = abs(lcm([coord_L_lamb[0].denominator(),coord_L_lamb[1].denominator(),coord_L_lamb[2].denominator()]))
   r = rho*coord_L_lamb
   coord_L_lamb_int = vector([r[0].numerator(),r[1].numerator(),r[2].numerator()])
   sigma =abs(gcd(gcd(coord_L_lamb_int[0],coord_L_lamb_int[1]),coord_L_lamb_int[2]))
   coord_L_lamb_prim = sigma^(-1)*coord_L_lamb_int
   lamb_prim=coord_L_lamb_prim*M
   return lamb_prim
   

def compute_a_b_unoriented_wedge(lamb,L):
# This takes a row vector lamb and returns two column vectors a and b=vareps*a with H(a) cap H(b) the line spanned by lamb. These are primitive in Hom(L,Z). They are also not necessarily properly oriented with respect to sigmaC.
   M=Matrix([vector(lamb),vector(lamb*Eps^(-1))])
   A = kernel(M.transpose())
   a = A.basis()[0]
   b = a*(Eps^(-1)).transpose()
   
   # Now we rescale a and b to make them primitive in Hom(L,Z)
   L_basis = ideal_basis(L)
   vect_a = vector([a.dot_product(vector(L_basis[0])), a.dot_product(vector(L_basis[1])), a.dot_product(vector(L_basis[2]))])
   vect_b = vector([b.dot_product(vector(L_basis[0])), b.dot_product(vector(L_basis[1])), b.dot_product(vector(L_basis[2]))])
   rho_a = abs(lcm([vect_a[0].denominator(),vect_a[1].denominator(),vect_a[2].denominator()]))
   rho_b = abs(lcm([vect_b[0].denominator(),vect_b[1].denominator(),vect_b[2].denominator()]))
   r_a = rho_a * vect_a
   r_b = rho_b * vect_b
   vect_a_int = vector([r_a[0].numerator(),r_a[1].numerator(),r_a[2].numerator()])
   vect_b_int = vector([r_b[0].numerator(),r_b[1].numerator(),r_b[2].numerator()])
   sigma_a =abs(gcd(gcd(vect_a_int[0],vect_a_int[1]),vect_a_int[2]))
   sigma_b =abs(gcd(gcd(vect_b_int[0],vect_b_int[1]),vect_b_int[2]))
   a_prim = rho_a * sigma_a^(-1) * a
   b_prim = rho_b * sigma_b^(-1) * b
   
   
   W_unoriented = Matrix([a_prim,b_prim]).transpose()
   
   return W_unoriented
   
   
def compute_a_b_wedge(lamb,L):
# This takes a row vector lamb and returns two column vectors a and b=vareps*a with H(a) cap H(b) the line spanned by lamb. These are primitive in Hom(L,Z). They are also properly oriented with respect to sigmaC (we achieve this by calling orient_a_b_wedge below).
   M=Matrix([vector(lamb),vector(lamb*Eps^(-1))])
   A = kernel(M.transpose())
   a = A.basis()[0]
   b = a*(Eps^(-1)).transpose()
  
   
   # Now we rescale a and b to make them primitive in Hom(L,Z)
   L_basis = ideal_basis(L)
   vect_a = vector([a.dot_product(vector(L_basis[0])), a.dot_product(vector(L_basis[1])), a.dot_product(vector(L_basis[2]))])
   vect_b = vector([b.dot_product(vector(L_basis[0])), b.dot_product(vector(L_basis[1])), b.dot_product(vector(L_basis[2]))])
   rho_a = abs(lcm([vect_a[0].denominator(),vect_a[1].denominator(),vect_a[2].denominator()]))
   rho_b = abs(lcm([vect_b[0].denominator(),vect_b[1].denominator(),vect_b[2].denominator()]))
   r_a = rho_a * vect_a
   r_b = rho_b * vect_b
   vect_a_int = vector([r_a[0].numerator(),r_a[1].numerator(),r_a[2].numerator()])
   vect_b_int = vector([r_b[0].numerator(),r_b[1].numerator(),r_b[2].numerator()])
   sigma_a =abs(gcd(gcd(vect_a_int[0],vect_a_int[1]),vect_a_int[2]))
   sigma_b =abs(gcd(gcd(vect_b_int[0],vect_b_int[1]),vect_b_int[2]))
   a_prim = rho_a * sigma_a^(-1) * a
   b_prim = rho_b * sigma_b^(-1) * b
   
   
   W_unoriented = Matrix([a_prim,b_prim]).transpose()
   W = orient_a_b_wedge(W_unoriented,L)
   
   return W
   
def is_oriented_a_b_wedge_check(W,adm):
# This returns True if the wedge W=a wedge b is positively oriented wrt to adm and False otherwise.
   W_e = insert_top_row(adm,W.transpose())
   R = (det(W_e) * (-sign(imag(det(or_matrix)))) > 0)
   return R
   
def compute_gamma_wedge(W,L):
# This computes gamma given a wedge a,b and lattice L
   v = kernel(W).basis()[0] # This returns some non-zero vector in H(a) intersect H(b)
   v = primitive(vector(v),L) # Now we rescale the vector by a positive rational number so v is primitive in L
   W_e = insert_top_row(v,W.transpose())
   if (det(W_e) * (sign(-imag(det(or_matrix)))) < 0):
      v=-v   
   return v   
   
def orient_a_b_wedge(W,L):
# This changes the signs of the vectors a and b in the wedge W to orient them properly with respect to sigmaC.
# INPUTS: W is a matrix whose columns are the vectors a and b (remember, elements of K are row vectors and so a and b, which are functionals on K, are column vectors.
# L is a fractional ideal
# OUTPUT: the same matrix W or -W.

  # First we obtain a and b and we make sure they are integral vectors
   a = (W.columns())[0]
   b = (W.columns())[1]
   rho_a = lcm([a[0].denominator(),a[1].denominator(),a[2].denominator()])
   rho_b = lcm([b[0].denominator(),b[1].denominator(),b[2].denominator()])
   a = rho_a * a
   b = rho_b * b 
   a  = vector([a[0].numerator(),a[1].numerator(),a[2].numerator()])
   b  = vector([b[0].numerator(),b[1].numerator(),b[2].numerator()])  # Now we at least have vectors with integral entries
   
   # Next we check whether we have to change sign of a (and b) to ensure they are properly oriented.
   H_a = kernel(Matrix(a).transpose())
   H_b = kernel(Matrix(b).transpose())
   v0 = H_a.basis()[0]
   v1 = H_a.basis()[1]
   v0p = v0
   v1p = v1
   M = Matrix([v0, v1,a])
   if (det(M)*(sign(-imag(det(or_matrix)))) < 0):
      v0p = v1
      v1p = v0 # Now v0p and v1p form a positively oriented basis of H(a)
   # Now we use the notion of oriented appearing in Lemma 3.2 of FelderDuke.
   z = imag(sigmaC(v0p)*conjugate(sigmaC(v1p)))
   if (z<0):
      a=-a
      b=-b

   W=Matrix([a,b]).transpose()   
   return W
   


def compute_alpha_beta(W,L):
# Computes alpha and beta for the wedge W, normalised by alpha(b) = beta(a) = 0 and
# alpha(a) = beta(b) = s, where s is the width, given by det(a,b,.) = s*gamma. Working in
# coordinates relative to a basis of L, so that L and Hom(L,Z) are both Z^3, the integer s
# is the gcd of the cross product of a and b, and any c in Hom(L,Z) with gamma(c) = 1 works :
# the rows alpha, beta, s*gamma are the comatrix of the rows a, b, c.
# INPUTS: a 2x3 matrix W whose columns are the vectors a and b forming the wedge, and a fractional ideal L.
# OUTPUT: a 3x2 matrix whose first row is alpha and second row is beta.

   # First we make sure that a and b are properly oriented with respect to sigmaC
   W = orient_a_b_wedge(W,L)
   a = (W.columns())[0]
   b = (W.columns())[1]

   # Coordinates of a and b in the basis of L. They are integral since a and b are
   # primitive in Hom(L,Z).
   M = ideal_basis(L)
   a_L = vector(ZZ,[a.dot_product(vector(M[0])), a.dot_product(vector(M[1])), a.dot_product(vector(M[2]))])
   b_L = vector(ZZ,[b.dot_product(vector(M[0])), b.dot_product(vector(M[1])), b.dot_product(vector(M[2]))])

   v = a_L.cross_product(b_L)
   s = gcd(gcd(v[0],v[1]),v[2])
   gamma_L = vector(ZZ,[v[0]/s, v[1]/s, v[2]/s])

   # A third element c of Hom(L,Z), with gamma(c) = 1, by Bezout: gamma is primitive in L.
   d, p, q = xgcd(gamma_L[0], gamma_L[1])
   d, u, w = xgcd(d, gamma_L[2])
   c = vector(ZZ,[u*p, u*q, w])

   alpha = b_L.cross_product(c) * M
   beta = c.cross_product(a_L) * M

   return Matrix([alpha, beta])
   
   
   
   
def wedge_width(adm,L):
# This computes the width 's' of the wedge associated with an admissible vector adm, that is the
# index in L of the sublattice generated by alpha, beta and gamma, which is also the size of F/Zgamma.
# With the alpha and beta returned by compute_alpha_beta, which are the minimal ones, this index
# is also the modulus mod(a,b) of Felder's Duke paper, defined there by det(a,b,.) = s*gamma. 

   W = compute_a_b_wedge(adm,L) # The columns of W are the vectors a and b of the wedge corresponding to adm.
   gamma = compute_gamma_wedge(W,L)
   [alpha, beta] = compute_alpha_beta(W,L) # These are alpha and beta for our wedge
   M = Matrix([alpha, beta, gamma])
   width = ZZ(abs(det(M)/det(ideal_basis(L))))
   return width
   

   
   
  
def find_small_width_admissible_vectors(N,L,max_width):
# This returns a list of admissible vectors such that the associated wedge has width at most max_width, by looking amongst the first (2N+1)^3 choices.
   all_adm_list = admissible_elements(N,L) # Generate a list of all admissible elements in a certain cube of size more or less (2N)^3
   n_adm = all_adm_list.dimensions()[0]
   list = Matrix([0,0,0])
   for j in range(1,n_adm+1):
      adm = all_adm_list[j-1]
      w = wedge_width(adm,L)
      if (w <= max_width):
         list = insert_top_row(adm,list)
         
   nrows = list.dimensions()[0]
   # Remove the last row (the meaningless [0,0,0])
   return list[0:nrows-1]
         
def compute_set_F(alpha, beta,gamma, L):
# Returns a set of representatives of F/Zgamma
   M = Matrix([alpha, beta, gamma])
   w = ZZ(abs(det(M)/det(ideal_basis(L))))
   list = Matrix([0,0,0])
   for j in range(0,w):
      for k in range(0,w):
         for l in range(0,w):
            lam = j/w * alpha + k/w * beta + l/w * gamma
            if (lam in L):
               list = insert_top_row(lam, list)
            
   nrows = list.dimensions()[0]
   # Remove the last row (the meaningless [0,0,0])
   return list[0:nrows-1]
   
   
def orientation_check(alpha,beta, gamma):
# Checks that (alpha,beta,gamma) is a positively oriented basis. The remaining
# checks, the values alpha(a), alpha(b), beta(a), beta(b), are printed by
# display_wedge_info.
   M = Matrix([alpha, beta, gamma])
   s = imag(det(M * or_matrix))
   if (s<0):
      print('(alpha,beta,gamma) is positively oriented basis')
   else:
      print('(alpha,beta,gamma) is negatively oriented basis (ERROR!)')


#------------------------------------------------------------------------------



##########     FUNCTIONS    #################

# EllGamma(z,tau,sigma): The basic elliptic gamma function, truncated using N_Ell_Gamma.


##### Provides bound for the tail of the series. With a = pi|Im tau|, b = pi|Im sigma|,
##### c = pi|Im(2z-tau-sigma)| and rho = a+b-c > 0, the terms decay like exp(-rho j)/j and
##### |sum_{j>=N} T_j| <= 4 exp(-rho N) / ( N (1-exp(-rho))(1-exp(-2Na))(1-exp(-2Nb)) ).

def decay_rate(z,tau,sigma):
# The terms of the series used by EllGamma decay like exp(-rho*j)/j. This rho is positive
# exactly in the region of convergence (16), and decides how many terms are needed:
# about d*ln(10)/rho of them for d decimals.
   a = (RBF(pi)*abs(CC(tau).imag())).lower()
   b = (RBF(pi)*abs(CC(sigma).imag())).lower()
   c = (RBF(pi)*abs(CC(2*z-tau-sigma).imag())).upper()
   return a + b - c


def tail_bound(z,tau,sigma,N):
# Upper bound for what truncating the series at N terms leaves out:
#    |sum_{j>=N} T_j| <= 4 exp(-rho N)/( N (1-exp(-rho))(1-exp(-2Na))(1-exp(-2Nb)) ),
# obtained from |sin(x+iy)| <= cosh(y) on the numerator and |sin(x+iy)| >= |sinh(y)| on each
# factor of the denominator. EllGamma adds it to the partial sum, so that the ball it returns
# encloses Gamma_ell itself and not merely the partial sum.
   a = (RBF(pi)*abs(CC(tau).imag())).lower()
   b = (RBF(pi)*abs(CC(sigma).imag())).lower()
   rho = decay_rate(z,tau,sigma)
   if not (a > 0 and b > 0 and rho > 0):
      raise ValueError('outside the region of convergence')
   Nb = RBF(N)
   E = 4*(-RBF(rho)*Nb).exp()
   E = E/(Nb*(1-(-RBF(rho)).exp())*(1-(-2*Nb*RBF(a)).exp())*(1-(-2*Nb*RBF(b)).exp()))
   return E.upper()


def certified_digits(x):
# Number of certified significant decimals of a ball, that is -log10(radius/|centre|).
   r = RealField(40)(max(x.real().rad(), x.imag().rad()))
   m = RealField(40)(CC(x).abs().lower())
   if (m <= 0) or (r <= 0):
      return 0
   return floor(-log(r/m)/log(RealField(40)(10)))


def EllGamma(z,tau,sigma):
# Uses the formula (15) in Felder-Varchenko to write Gamma_ell as exp of an infinite sum.

   E = tail_bound(z,tau,sigma,N_Ell_Gamma)
   s = 0
   for j in range(1,N_Ell_Gamma):
      s=s+(sin(Pi*j*(2*z-tau-sigma)))/(j*sin(Pi*j*tau)*sin(Pi*j*sigma))
   s = CC(s).add_error(E)
   gamma=exp(-I*s/2)
   return gamma


def EllGamma_aa_h(w,h,L):
# This is the smoothed elliptic Gamma function. The products below are the
# expression of Gamma_{a,b} (equation (10) of Felder's Duke paper) given by
# Prop. 3.5 of Felder-Henriques-Rossi-Zhu, i.e. formula (2.7) of the paper.
   global tau, sigma   
   adm = h * q_ff
   W_un = compute_a_b_unoriented_wedge(adm,L)
   W = compute_a_b_wedge(adm,L)
   
   gamma = compute_gamma_wedge(W,L)
   [alpha, beta] = compute_alpha_beta(W,L)
   F = compute_set_F(alpha,beta,gamma,L)
   s = F.dimensions()[0] # The cardinality of the set F
   tau = sigmaC(alpha)/sigmaC(gamma)
   sigma = sigmaC(beta)/sigmaC(gamma)
   
   display_wedge_info(h,W,W_un,alpha,beta,gamma,tau,sigma,F,s)

   # First the denominator
   denominator = 1
   for j in range(0,s):
      delta=F[j]
      denominator = denominator * EllGamma((w+sigmaC(delta))/sigmaC(gamma),tau,sigma)^(N_aa)
   
   numerator = 1
   for j in range(0,s):
      delta=F[j]
      numerator = numerator * EllGamma(N_aa*(w+sigmaC(delta))/sigmaC(gamma),N_aa*tau,N_aa*sigma)
   R = numerator/denominator   
  
   
   #  check the orientation.
   if is_oriented_a_b_wedge_check(W,adm):
      print("Wedge positively oriented with respect to h")
   else:
      print("Wedge negatively oriented with respect to h")

   print('certified decimals = ', certified_digits(R))

   return R


def in_K(z):
# The element  z of K recovered by lindep on 1, yCC and yCC^2.
   c = pari([1, yCC.mid(), (yCC^2).mid(), z.mid()]).lindep(100)
   return -(ZZ(c[0]) + ZZ(c[1])*K.gen() + ZZ(c[2])*K.gen()^2)/ZZ(c[3])


def display_wedge_info(h,W,W_un,alpha,beta,gamma,tau,sigma,F,s):
# This prints on screen all the info regarding the wedge, tau, sigma, etc, associated with a choice of h.
 
   if (W_un == W):
      print('Did not have to change signs of a and b to satisfy Lemma 2.1 (ii).')
   else:
      print('Had to change signs of a and b to satisfy Lemma 2.1 (ii).')
   print('Evaluating at h = adm/'+str(q_ff)+' =', h, ' with adm =', h*q_ff)
   print('The wedge for ', h , ' is ')
   a=W.columns()[0]
   b=W.columns()[1]
   print('a = ',a)
   print('b = ',b)
   print('alpha =',alpha, ', alpha(b) = ',alpha*b, 'alpha(a) = ', alpha*a)
   print('beta =',beta, ', beta(a) = ',beta*a, 'beta(b) = ', beta*b)
   print('gamma =',gamma)
   orientation_check(alpha,beta,gamma)
   print('tau = ', ComplexField(show_prec)(tau.mid()))
   print('sigma = ', ComplexField(show_prec)(sigma.mid()))
   print('The set F has ',s,'elements. They are')
   print(F)
# The triples at which the elliptic Gamma is evaluated, as exact elements of K: the
# denominator uses (z,tau,sigma) and the numerator (N_aa*z, N_aa*tau, N_aa*sigma).
   for j in range(0,s):
      z = (sigmaC(h) + sigmaC(F[j]))/sigmaC(gamma)
      print('(z, tau, sigma) in K = ', [in_K(z), in_K(tau), in_K(sigma)])

   ######   
   
   
def EllGammaunit_aa_h(h,L):
# This is the number we are looking after 
   return EllGamma_aa_h(sigmaC(h),h,L)

   

def scan_admissible(N,L,max_width,D=1000):
# Prints the five admissible vectors of L whose wedge has width at most max_width and the
# smallest ratio s/rho, together with the number of terms N_Ell_Gamma each of them needs for D
# decimals, and returns the first of them. One evaluation costs 2s calls to EllGamma of about
# D*ln(10)/rho terms each, so it is s/rho, and not s alone, that has to be made small.

   A = find_small_width_admissible_vectors(N,L,max_width)
   n = A.dimensions()[0]
   if (n == 0):
      print('No admissible vector of width at most', max_width, 'in the box of radius', N, '.')
      print('Raise max_width, or raise N if admissible_elements(N,L) is empty too.')
      return 0

   # For each admissible vector we record its width s, its decay rate rho and the ratio s/rho.
   width = []
   rate = []
   ratio = []
   for j in range(0,n):
      adm = A[j]
      W = compute_a_b_wedge(adm,L)
      gamma = compute_gamma_wedge(W,L)
      [alpha, beta] = compute_alpha_beta(W,L)
      F = compute_set_F(alpha,beta,gamma,L)
      s = F.dimensions()[0]
      tau = sigmaC(alpha)/sigmaC(gamma)
      sigma = sigmaC(beta)/sigmaC(gamma)
      w = sigmaC(vector(adm)/q_ff)

      # The decay rate depends on the point z at which EllGamma is called, and the slowest of
      # the s points of F is the one that decides. The numerator is evaluated at
      # (N_aa*z, N_aa*tau, N_aa*sigma), whose rate is N_aa times the rate computed here, so the
      # smallest rate always comes from the denominator.
      rho = decay_rate((w + sigmaC(F[0]))/sigmaC(gamma),tau,sigma)
      for k in range(1,s):
         r = decay_rate((w + sigmaC(F[k]))/sigmaC(gamma),tau,sigma)
         if (r < rho):
            rho = r

      width.append(s)
      rate.append(rho)
      ratio.append(s/rho)

   # We bring the five smallest ratios to the front. 
   order = list(range(0,n))
   for i in range(0,min(5,n)):
      for k in range(i+1,n):
         if (ratio[order[k]] < ratio[order[i]]):
            temp = order[i]
            order[i] = order[k]
            order[k] = temp

   print('the five best admissible vectors, with the number of terms they need for', D, 'decimals:')
   print('s/rho      s     rho         N_Ell_Gamma  adm')
   for i in range(0,min(5,n)):
      j = order[i]
      N_needed = ceil(D*log(10)/rate[j])
      print(str(RealField(20)(ratio[j])).ljust(10),
            str(width[j]).ljust(5),
            str(RealField(20)(rate[j])).ljust(11),
            str(N_needed).ljust(12),
            A[j])

   return A[order[0]]




##########   ALGEBRAIC RECOGNITION    ##########

# Two routines. algdep_Q looks for a relation over Q, algdep_K for one over the ring of
# integers of K, which is the form the polynomials of the paper take.
# In both cases the number of trusted digits must be supplied: it is the value printed as
# 'certified decimals' by EllGamma_aa_h, that is certified_digits(u), up to a small margin.



def algdep_Q(u, degree, digits):
# Looks for the minimal polynomial over Q of u, of the given degree, and returns it in ZZ[x].
# This is the algdep of Sage

   z = ComplexField(prec)(u.mid())
   return algdep(z, degree, known_digits = digits)


def algdep_K(u, degree, digits, beta=None):
# Looks for a relation sum_{i<3, j<=degree} c_{ij} beta^i u^j = 0 with integer c_{ij}, and
# returns it as a monic polynomial in K[x] as  in the paper. The relation has
# to be assembled in K[x].

   if beta is None:
      beta = yCC
   z = ComplexField(prec)(u.mid())
   bC = ComplexField(prec)(beta)
   v = []
   for j in range(0,degree+1):
      for i in range(0,3):
         v.append(bC^i * z^j)
   c = pari(v).lindep(digits)
   coeffs = []
   for j in range(0,degree+1):
      coeffs.append(ZZ(c[3*j]) + ZZ(c[3*j+1])*K.gen() + ZZ(c[3*j+2])*K.gen()^2)
   P = PolynomialRing(K,'x')(coeffs)
   if P == 0:
      return P
   return P/P.leading_coefficient()

