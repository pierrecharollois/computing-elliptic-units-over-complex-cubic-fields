#### TO DO: 1) commnent every function, including INPUT and OUTPUT, and the code inside. 2) Write a README file explaining how to use all this. 


####### GLOBAL VARIABLES ######

# prec: Number of digits of precision
# CC: the complex field with prec digits of precision.
# Pi:  the number pi in the field CC.
# N_Ell_Gamma: Number of terms used to compute Elliptic Gamma
# N_theta0: Number of terms used to compute theta0 (not strictly needed, only want to check the vanishing of the boundary term in the proof of the main theorem).

prec = 1000
N_Ell_Gamma = 10000 # Number of terms used to compute Elliptic Gamma
N_theta0 = 500
CC=ComplexBallField(prec)
RBF=RealBallField(prec)   # ADDED in v5: needed by the certified tail bound below
Pi=CC(pi,0)

#------------------------------------------------------------------------------




###### ------ LOAD NUMBER FIELD DATA -------  #######

#load('cubic_field_K_disc_-588.sage')
#load ('cubic_field_K_cube_root_3.sage')
load ('cubic_field_K_cube_root_7.sage')
#load ('cubic_field_K_cube_root_7_conj_sigmaC.sage')
#load ('cubic_field_K_cube_root_2.sage')
#load('cubic_field_K_disc_-23.sage')
#load('cubic_field_K_disc_-411.sage')

##### -----------------------------------------------------    ######


##### SOME ROUTINES FOR WORKING WITH NUMBER FIELDS #####


#def convert_bases_matrix:
# TO DO: This returns a 3 by 3 matrix that gives the change of basis from K_basis to K_basis_SAGE


def ideal_basis(I):
# This function computes an integral basis of the underlying Z-module to an ideal I of K. The output is a 3 by 3 matrix. The entries in each row are the coordinates of a basis of I with respect to the basis K_basis_SAGE (to do: return also the matrix with respect to K_basis)
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
   

######## -------------------------------------------------  #########


###### ROUTINES THAT PRODUCE A LIST OF ADMISSIBLE ELEMENTS    ######## 
   
def insert_top_row(row,M):
    return matrix([row]+M.rows()[0:])   

   
def admissible_elements(N,L):
# This function generates a list of admissible elements in the lattice L wrt the ideal aa and the parameter q.  The idea is that we are looking for vectors lambd in q* 1_K + qL that are also primitive in L and in aa_compL, a condition that is equivalent to lambd being admissible. 
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
            #print candidate
            #list=insert_top_row(candidate,list)
            v=candidate*M_aa_compL^(-1)
            w=candidate*M_L^(-1)
            if (v[0].is_integer()) and (v[1].is_integer()) and (v[2].is_integer()) and (w[0].is_integer()) and (w[1].is_integer()) and (w[2].is_integer()):
               if (gcd(gcd(v[0],v[1]),v[2])==1) and (gcd(gcd(w[0],w[1]),w[2])==1):
                  list=insert_top_row(candidate,list)
   nrows = list.dimensions()[0]
   
   # Remove the last row (the meaningless [0,0,0])
   return list[0:nrows-1]

   
def weakly_admissible_elements(N,L):
# This function generates a list of weakly admissible elements in the lattice L wrt the ideal aa and the parameter q.  The idea is that we are looking for vectors lambd in q* 1_K + qL that are also in aacomp_L but not in N(\mathfrak{a})L (equivalently: lambd belongs to q*1_K +qL and lambd/N_aa generates aa^(-1)L/L).

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
            if (candidate in aa_compL) and not(candidate in N_aa * L):
                  list=insert_top_row(candidate,list)
   nrows = list.dimensions()[0]
   
   # Remove the last row (the meaningless [0,0,0])
   return list[0:nrows-1]




######## --------------------------------------------------------- ##########





########  ROUTINES TO HANDLE ADMISSIBLE VECTORS AND WEDGES ########


def primitive(lamb,L):
# This takes a row vector with rational entries (corresponding to the coordinates of an element of K in the basis K_basis_SAGE), and rescales by a unique positive rational number to make sure it belongs to the fractional ideal L.
   M=ideal_basis(L)
   coord_L_lamb = lamb*M^(-1)
   #print coord_L_lamb
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
   #print lamb
   A = kernel(M.transpose())
   a = A.basis()[0]
   b = a*(Eps^(-1)).transpose()
   #print a
   #print b
   
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
   #print lamb
   A = kernel(M.transpose())
   a = A.basis()[0]
   b = a*(Eps^(-1)).transpose()
   #print a
   #print b
   
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
   #print v
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
   #print z
   if (z<0):
      a=-a
      b=-b

   W=Matrix([a,b]).transpose()   
   return W
   

   

def basis_L_mod_Zgamma(adm,L):
#TO DO: This computes a pair of vectors in the fractional ideal L that give an isom Z^2 simeq L/Zgamma
   return adm
   
def compute_alpha_beta(W,L):
# Compute a pair of vectors alpha and beta for the wedge W. NOTE: alpha and beta are not unique. Other choices may lead to smaller width!
# INPUTS: a 2x3 matrix W whose columns are the vectors a and b forming the wedge, and a fractional ideal L.
# OUTPUT: a 3x2 matrix whose first row is alpha and second row is beta.

   # First we make sure that a and b are properly oriented with respect to sigmaC
   W = orient_a_b_wedge(W,L)
   
   # Next we compute first the hyperplanes H(a) and H(b) forming the wedge, and bases (v0,v1) for H(a) and (w0,w1) for H(b).
   a = (W.columns())[0]
   b = (W.columns())[1]
   H_a = kernel(Matrix(a).transpose())
   H_b = kernel(Matrix(b).transpose())
   v0 = H_a.basis()[0]
   v1 = H_a.basis()[1]
   w0 = H_b.basis()[0]
   w1 = H_b.basis()[1]
   
   # Next we compute gamma for this wedge.
   gamma = compute_gamma_wedge(W,L)
   
   
  # Now we can choose alpha among w0, -w0, w1, -w1. We only need to choose the one that is not proportional to gamma (equivalently the one that does not lie on H(b)) and that gives a positive result when evaluating at a. First we find alpha.
   alpha = w0
   if (a.dot_product(w0) == 0):
      alpha = w1 # Now we are sure that alpha at least satisfies alpha(b)=0 and alpha(a) neq 0
   if (alpha.dot_product(a) < 0):
      alpha = -alpha # Now we know that alpha lies on H(b) and satisfies alpha(a)>0.

   #Now the same with beta   
   beta = v0
   if (b.dot_product(v0) == 0):
      beta = v1 # Now we are sure that beta at least satisfies beta(a)=0 and beta(b) neq 0
   if (beta.dot_product(b) < 0):
      beta = -beta # Now we know that alpha lies on H(b) and satisfies alpha(a)>0.
   
   M = Matrix([primitive(vector(alpha),L),primitive(vector(beta),L)])
   
   return M
   
   
   
   
   
def wedge_width(adm,L):
# This computes the width of the wedge associated with an admissible vector adm and our choice of data (alpha, beta). Careful! This is different from the original thing we called the width, which Felder in the Duke paper calls the modulus. What this function computes is simply the index in L of the sublattice generated by alpha, beta and gamma. This will be the size of F/Zgamma.

   W = compute_a_b_wedge(adm,L) # The columns of W are the vectors a and b of the wedge corresponding to adm.
   gamma = compute_gamma_wedge(W,L)
   [alpha, beta] = compute_alpha_beta(W,L) # These are alpha and beta for our wedge
   M = Matrix([alpha, beta, gamma])
   width = ZZ(abs(det(M)/det(ideal_basis(L))))
   return width
   

#def wedge_modulus(adm,L):
# To do: compute the modulus of the wedge in the sense of Felder et al. Not necessary.

   
   
  
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
         
def find_small_width_weakly_admissible_vectors(N,L,max_width):
# This returns a list of weakly admissible vectors such that the associated wedge has width at most max_width, by looking amongst the first (2N+1)^3 choices.
   all_w_adm_list = weakly_admissible_elements(N,L) # Generate a list of all weakly admissible elements in a certain cube of size more or less (2N)^3
   n_adm = all_w_adm_list.dimensions()[0]
   list = Matrix([0,0,0])
   for j in range(1,n_adm+1):
      adm = all_w_adm_list[j-1]
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
   #w = wedge_width(gamma,L)
   #print w
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
# TO DO: print all the checks that everything is well oriented
   M = Matrix([alpha, beta, gamma])
   s = imag(det(M * or_matrix))
   if (s<0):
      print('(alpha,beta,gamma) is positively oriented basis')
   else:
      print('(alpha,beta,gamma) is negatively oriented basis (ERROR!)')


#------------------------------------------------------------------------------



##########     FUNCTIONS    #################

# EllGamma(z,tau,sigma): The basic elliptic gamma function, truncated using N_Ell_Gamma.
# EllGamma_ab(w,x,L,a,b): The elliptic Gamma function for the lattice L and the wedge defined by a and b. 
# theta0 and theta0_smoothed: these are not really necessary to compute our numbers (they are not called by any of the other functions in this file) but they are easy so I include them here.
# theta and theta_smoothed: again not necessary but I include them to test the Kronecker limit formula that we write down in our paper.


def theta0(z,tau):
# We evaluate theta0 using the infinite product definition. This is probably very slow but that's fine.
   p=1-exp(2*Pi*I*z)
   for j in range(1,N_theta0):
      p = p*(1-exp(2*Pi*I*(j*tau+z)))*(1-exp(2*Pi*I*(j*tau-z)))
   return p
      
      
def theta0_smoothed(z,tau,N):
# This is the smoothed theta0 function
   numerator = theta0(z,tau)^N
   denominator = theta0(N*z,N*tau)
   R = numerator/denominator
   return R
         

def theta(z,tau):
# This is the theta function as it appears in Colmez's Bourbaki expos\'e ``La conjecture de Birch and Swinnerton--Dyer p-adique'', middle of p. 274. The relation with theta0 is theta = -q_tau^(1/12)*q_z^(-1/2)*theta0.
  f = -exp(2*Pi*I*tau/12)*exp(-Pi*I*z)
  t = f* theta0(z,tau)
  return t
  
def theta_smoothed(z,tau,N):
# This is the smoothed theta function. Note that for N odd it is the same as the smoothed theta0!
   numerator = theta(z,tau)^N
   denominator = theta(N*z,N*tau)
   R = numerator/denominator
   return R
  

   


      

##### ADDED in v5: certification of the truncation error #####
# The three functions below, together with the global statistics, are the only
# additions of v5 to the numerical part. They bound the tail of the series in
# EllGamma so that the returned ball is a rigorous enclosure of Gamma_ell and
# not merely of its partial sum. The mathematics is unchanged.
#
# Bound. With T_j the term of index j, a = pi*|Im(tau)|, b = pi*|Im(sigma)|,
# c = pi*|Im(2z-tau-sigma)| and rho = a+b-c > 0, one has for every N >= 1
#   |sum_{j>=N} T_j| <= 4*exp(-rho*N)
#                        / ( N*(1-exp(-rho))*(1-exp(-2*N*a))*(1-exp(-2*N*b)) ),
# from |sin(x+iy)| <= cosh(y) applied to the numerator, |sin(x+iy)| >= |sinh(y)|
# applied to each denominator factor, cosh(t) <= exp(t) and
# sinh(t) = (exp(t)/2)*(1-exp(-2t)), then summing the geometric series.

rho_min = None       # smallest decay rate encountered in the current evaluation
tail_max = None      # largest tail bound used in the current evaluation
n_calls = 0          # number of EllGamma calls in the current evaluation

def reset_stats():
   global rho_min, tail_max, n_calls
   rho_min = None
   tail_max = None
   n_calls = 0

def decay_data(z,tau,sigma):
# Returns (rho, a, b) with rho, a, b LOWER bounds of the true quantities. The
# tail bound decreases in each of rho, a, b, so lower bounds keep it valid.
   a = (RBF(pi)*abs(CC(tau).imag())).lower()
   b = (RBF(pi)*abs(CC(sigma).imag())).lower()
   c = (RBF(pi)*abs(CC(2*z-tau-sigma).imag())).upper()
   if not (a > 0) or not (b > 0):
      raise ValueError("Im(tau) or Im(sigma) is not separated from zero: a = %s, b = %s" % (a,b))
   rho = a + b - c
   if not (rho > 0):
      raise ValueError("outside the region of convergence: rho = %s" % rho)
   return rho, a, b

def tail_bound(N,rho,a,b):
# Rigorous upper bound for |sum_{j>=N} T_j|, as a RealNumber.
   Nb = RBF(N); rb = RBF(rho); ab = RBF(a); bb = RBF(b)
   E = 4*(-rb*Nb).exp()/(Nb*(1-(-rb).exp())*(1-(-2*Nb*ab).exp())*(1-(-2*Nb*bb).exp()))
   return E.upper()

def certified_digits(x):
# Certified significant decimal digits of a ball x, i.e. floor(-log10(rad/|centre|)).
# Returns 0 when the ball contains zero.
   try:
      r = RealField(53)(x.rad())
   except (AttributeError, TypeError):
      r = RealField(53)(max(x.real().rad(), x.imag().rad()))
   m = RealField(53)(CC(x).abs().lower())
   if m <= 0 or r <= 0:
      return Integer(0)
   return max(Integer(floor(-log(r/m)/log(RealField(53)(10)))), Integer(0))


def EllGamma(z,tau,sigma):
# We use the formula (15) in Felder-Varchenko to write Gamma_ell as exp of an infinite sum. The region of convergence is specified in (16).
# CHANGED in v5: the convergence condition (16) is now checked rather than
# assumed, and the truncation error is added to the partial sum before
# exponentiating. The summation loop itself is unchanged.
   global rho_min, tail_max, n_calls
   rho, a, b = decay_data(z,tau,sigma)
   E = tail_bound(N_Ell_Gamma,rho,a,b)
   s = 0
   for j in range(1,N_Ell_Gamma):
      s=s+(sin(Pi*j*(2*z-tau-sigma)))/(j*sin(Pi*j*tau)*sin(Pi*j*sigma))
   s = CC(s).add_error(E)
   if rho_min is None or rho < rho_min:
      rho_min = rho
   if tail_max is None or E > tail_max:
      tail_max = E
   n_calls = n_calls + 1
   gamma=exp(-I*s/2)
   return gamma
   
  
 
   
   
#def EllGamma_a_b(w,W,L):
# This computes the elliptic gamma function Gamma_{a,b} in Felder's Duke paper, equation (10), using the formula in Prop. 3.5. It is obsolete, I don't use it anymore. The problem is it works, but it is very slow because every time it gets called it computes F, etc. , which is very slow because to compute the smoothed elliptic zeta we are computing the same things multiple times.
#   #W = compute_a_b_wedge(adm,L)
#   gamma = compute_gamma_wedge(W,L)
#   [alpha, beta] = compute_alpha_beta(W,L)
#   F = compute_set_F(alpha,beta,gamma,L)
#   s = F.dimensions()[0] # The cardinality of the set F
#   tau = sigmaC(alpha)/sigmaC(gamma)
#   sigma = sigmaC(beta)/sigmaC(gamma)
#   p=1
#   for j in range(0,s):
#      delta=F[j]
#      p = p* EllGamma((w+sigmaC(delta))/sigmaC(gamma),tau,sigma)
#   return p
      
def EllGamma_aa_h(w,h,L):
# This is the smoothed elliptic Gamma function
   reset_stats()   # ADDED in v5
   adm = h * q_ff
   W_un = compute_a_b_unoriented_wedge(adm,L)
   W = compute_a_b_wedge(adm,L)
   #check_weak_admissibility(h,L)
   
   gamma = compute_gamma_wedge(W,L)
   [alpha, beta] = compute_alpha_beta(W,L)
   F = compute_set_F(alpha,beta,gamma,L)
   s = F.dimensions()[0] # The cardinality of the set F
   tau = sigmaC(alpha)/sigmaC(gamma)
   sigma = sigmaC(beta)/sigmaC(gamma)
   
   display_wedge_info(h,W,W_un,alpha,beta,gamma,tau,sigma,F,s)

   # The commented code below is obsolete.
   # First the denominator
   #denominator = EllGamma_a_b(w,W,L)^(N_aa)
   # Now we compute the numerator using the distribution relation
   #numerator = 1
   #for k in range(0,N_aa):
   #   numerator = numerator * EllGamma_a_b(w+k*sigmaC(gamma)/N_aa,W,L)
   #R = numerator/denominator
   
   # First the denominator
   denominator = 1
   for j in range(0,s):
      delta=F[j]
      denominator = denominator * EllGamma((w+sigmaC(delta))/sigmaC(gamma),tau,sigma)^(N_aa)
   
   # Now we compute the numerator using the distribution relation
   # CHANGED in v5: the product over k of Gamma(z+k/N_aa, tau, sigma) is contracted
   # into the single value Gamma(N_aa*z, N_aa*tau, N_aa*sigma), reducing the number
   # of EllGamma calls from s*(1+N_aa) to 2*s. Same value, faster convergence.
   numerator = 1
   for j in range(0,s):
      delta=F[j]
      numerator = numerator * EllGamma(N_aa*(w+sigmaC(delta))/sigmaC(gamma),N_aa*tau,N_aa*sigma)
   R = numerator/denominator   
  
   
   # We also check the orientation.
   if is_oriented_a_b_wedge_check(W,adm):
      print("Wedge positively oriented with respect to h")
   else:
      print("Wedge negatively oriented with respect to h")

   # ADDED in v5: certification report
   print()
   print('--- certification ---')
   print('N_Ell_Gamma          = ', N_Ell_Gamma, 'terms,', n_calls, 'calls')
   print('smallest rho         = ', RealField(30)(rho_min))
   print('largest tail bound   = ', RealField(15)(tail_max))
   print('CERTIFIED DECIMALS   = ', certified_digits(R))
   print('---------------------')

   return R
      
   
#def display_wedge_info(h,L):
## This prints on screen all the info regarding the wedge, tau, sigma, etc, associated with a choice of h.
#   adm = h * q_ff
#   W_unoriented = compute_a_b_unoriented_wedge(adm,L)
#   W = compute_a_b_wedge(adm,L)
#   if (W_unoriented == W):
#      print 'Did not have to change signs of a and b to satisfy (ii).'
#   else:
#      print 'Had to change signs of a and b to satisfy (ii).'
#   print 'Evaluating at h =',h,'.'
#   #check_weak_admissibility(h,L)
#   print 'The wedge for ', h , ' is '
#   a=W.columns()[0]
#   b=W.columns()[1]
#   print 'a = ',a
#   print 'b = ',b
#
#   
#   gamma = compute_gamma_wedge(W,L)
#   [alpha, beta] = compute_alpha_beta(W,L)
#   F = compute_set_F(alpha,beta,gamma,L)
#   s = F.dimensions()[0] # The cardinality of the set F
#   tau = sigmaC(alpha)/sigmaC(gamma)
#   sigma = sigmaC(beta)/sigmaC(gamma)
#      
#   print 'alpha =',alpha, ', alpha(b) = ',alpha*b, 'alpha(a) = ', alpha*a
#   print 'beta =',beta, ', beta(a) = ',beta*a, 'beta(b) = ', beta*b
#   print 'gamma =',gamma
#   orientation_check(alpha,beta,gamma)
#   print 'tau = ', tau
#   print 'sigma = ', sigma
#   print 'The set F has ',s,'elements. They are'
#   print F
#   ######   
   
   
def display_wedge_info(h,W,W_un,alpha,beta,gamma,tau,sigma,F,s):
# This prints on screen all the info regarding the wedge, tau, sigma, etc, associated with a choice of h.
 
   if (W_un == W):
      print('Did not have to change signs of a and b to satisfy (ii).')
   else:
      print('Had to change signs of a and b to satisfy (ii).')
   print('Evaluating at h =',h,'.')
   print('The wedge for ', h , ' is ')
   a=W.columns()[0]
   b=W.columns()[1]
   print('a = ',a)
   print('b = ',b)
   print('alpha =',alpha, ', alpha(b) = ',alpha*b, 'alpha(a) = ', alpha*a)
   print('beta =',beta, ', beta(a) = ',beta*a, 'beta(b) = ', beta*b)
   print('gamma =',gamma)
   orientation_check(alpha,beta,gamma)
   print('tau = ', tau)
   print('sigma = ', sigma)
   print('The set F has ',s,'elements. They are')
   print(F)
   ######   
   
   
def EllGammaunit_aa_h(h,L):
# This is the number.
   #display_wedge_info(h,L)
   return EllGamma_aa_h(sigmaC(h),h,L)

   
   



def rho_of_admissible(adm,L):
# ADDED in v5. Returns the decay rate rho of the slowest series involved, without
# running the evaluation. Use it to compare admissible vectors of equal width
# before committing: the truncation error decreases like exp(-rho*N_Ell_Gamma).
   W = compute_a_b_wedge(adm,L)
   gamma = compute_gamma_wedge(W,L)
   [alpha, beta] = compute_alpha_beta(W,L)
   F = compute_set_F(alpha,beta,gamma,L)
   s = F.dimensions()[0]
   tau = sigmaC(alpha)/sigmaC(gamma)
   sigma = sigmaC(beta)/sigmaC(gamma)
   w = sigmaC(vector(adm)/q_ff)
   r = None
   for j in range(0,s):
      delta = F[j]
      z0 = (w+sigmaC(delta))/sigmaC(gamma)
      for (zz,tt,ss) in [(z0,tau,sigma),(N_aa*z0,N_aa*tau,N_aa*sigma)]:
         rr = decay_data(zz,tt,ss)[0]
         if r is None or rr < r:
            r = rr
   return r
