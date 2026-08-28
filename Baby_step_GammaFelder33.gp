\\ ---------------------------------------------------------------------------
\\ Baby_step_GammaFelder33.gp
\\
\\ Evaluation of the Felder--Varchenko elliptic Gamma function
\\      Gamma(z; tau, sigma) = prod_{j,k >= 0}
\\        (1 - e^{2 i Pi ((j+1)tau + (k+1)sigma - z)}) / (1 - e^{2 i Pi (j tau + k sigma + z)})
\\ by SL_3(Z) modular moves ("baby steps"), theta cocycle factors and
\\ truncated Lambert series.  Companion code for BCG, arXiv:2311.04110.
\\
\\ References: Felder--Varchenko, "The elliptic gamma function and SL(3,Z)xZ^3",
\\ Adv. Math. 156 (2000), formulas (15)-(16); Felder--Henriques--Rossi--Zhu,
\\ Prop. 3.5.  The three-term modular relation used in split() is
\\   Gamma(z/sigma; tau/sigma, -1/sigma)
\\     = e^{i Pi Q(z,tau,sigma)} Gamma((z-sigma)/tau; -1/tau, -sigma/tau)
\\       * Gamma(z; tau, sigma),
\\ with Q the cubic polynomial of testexpfactor() below.
\\
\\ Public interface:
\\   lessnaifGamma(z, TAU, SIGMA, BBd, BDterm)   main evaluator
\\   THETA0(z, TAU, BOUND)                        normalized Jacobi theta
\\   lessnaifEuler(x, q, r, BBd, BDterm)          double product prod_{j,k>=0}(1 - x q^j r^k)
\\   basicGamma(z, TAU, SIGMA, BD)                naive oracle (slow, reference)
\\
\\ Accuracy contract: lessnaifGamma aims at absolute error O(10^-BDterm) on
\\ log Gamma; take BDterm = target digits + 10 and BBd >= 8*BDterm.  This is
\\ an empirical contract, NOT a certified bound (see the Sage chain for that).
\\ lessnaifEuler prints a warning (once) if a sum exits through the BBd cap
\\ before reaching the 10^-BDterm threshold: raise BBd in that case.
\\
\\ Notes and limits: z is reduced mod 1 on entry (Gamma is 1-periodic in z);
\\ real z (the strip boundary) is handled exactly.  The theta series uses
\\ BDspecial = 40 terms after reduction to Im(tau) > 0.865, which is exact
\\ below roughly 3800 decimal digits (term j ~ |q|^(j^2), |q| < 0.00437);
\\ raise BDspecial beyond that.
\\ ---------------------------------------------------------------------------

\\ ------------------------------- globals -----------------------------------

BDspecial   = 40;      \\ theta terms after reduction; enough below ~3800 digits
SPLITBOUND  = 0.8;     \\ below this Im, normalize (tau,sigma) through split()
EULERWARN   = 0;       \\ set to 1 after the first truncation warning

\\ ---------------------------------------------------------------------------
\\ basicGamma: naive evaluation from the log series
\\   log Gamma = -(I/2) sum_{j>=1} sin(Pi j (2z-tau-sigma)) / (j sin(Pi j tau) sin(Pi j sigma))
\\ Converges iff |Im(2z-tau-sigma)| < Im(tau)+Im(sigma), at geometric rate
\\ rho = Pi (Im tau + Im sigma - |Im(2z-tau-sigma)|).  Slow; kept as an
\\ independent naive reference.

basicGamma(z, TAU, SIGMA, BD)={
local(S,termnum, den1,den2);
if(imag(2*z-TAU-SIGMA) < abs(imag(TAU))+abs(imag(SIGMA)), ,
   print("basicGamma: argument outside the convergence region"));
S=0;
for(j=1, BD,
  termnum=sin(Pi*j*(2*z-TAU-SIGMA));
  den1=sin(Pi*j*TAU);
  den2=sin(Pi*j*SIGMA);
  S=S+termnum/(j*den1*den2));
return(exp(-I/2*S))
}

\\ ---------------------------------------------------------------------------
\\ newTHETAinit: the theta sum  -sum_j e^{i Pi (tau (j+1/2)^2 + 2 (j+1/2)(z+1/2))}
\\ after three normalizations: Re(tau) reduced mod 1 (factor e^{n i Pi/4}),
\\ tau -> -1/tau when Im(tau) <= 0.865 (Jacobi transformation), and z pushed
\\ up by m*tau so that the series over |j| <= BDspecial suffices.
\\ BOUND is kept in the interface for compatibility; the truncation level is
\\ the global BDspecial (see header).

newTHETAinit(z,TAU, BOUND)={
local(qq, SUM, TERM, n,m,x1,x2,TT, TEMP, TERM1, TERM2, ZEROterm, FACTOR, NEWZ);
qq=exp(2*I*Pi*TAU);
if(abs(qq)>1, print("newTHETAinit: |q| > 1, wrong half-plane"); return("STOP"));
n=round(real(TAU));
if(n==0, , return(newTHETAinit(z,TAU-n, BOUND)*exp(n*I*Pi/4)));
if(imag(TAU)>0.865, ,
  [x1,x2]=polroots(x^2-TAU/I);
  if(real(x1)>0, TT=x1, TT=x2);
  TEMP=newTHETAinit(z/TAU,-1/TAU, BOUND)*exp(-I*Pi*z*z/TAU)/(-I*TT); return(TEMP));
m=ceil((imag(TAU-z)/imag(TAU)));
if(m==0, , return(newTHETAinit(z+m*TAU,TAU, BOUND)*(-1)^m*exp(I*Pi*TAU*m^2+2*I*Pi*m*z)));
FACTOR=(-1)^m;
NEWZ=z+m*TAU;
SUM=0;
for(j=1, BDspecial,
  TERM1=(j+1/2)*(TAU*(j+1/2)+(2*NEWZ+1));
  TERM2=(-j+1/2)*(TAU*(-j+1/2)+(2*NEWZ+1));
  SUM=SUM+exp(I*Pi*(TERM1+m*(m*TAU+2*z)))+exp(I*Pi*(TERM2+m*(m*TAU+2*z))));
ZEROterm=exp(I*Pi/2*(TAU/2+2*NEWZ+1+2*m*(m*TAU+2*z)));
if(abs(ZEROterm)>10, print("newTHETAinit: suspicious zero term ", ZEROterm));
return((-SUM-ZEROterm)*FACTOR)
}

\\ ---------------------------------------------------------------------------
\\ THETA0: the normalized theta of Felder, theta_0(z,tau) =
\\   -i e^{i Pi (z - tau/4)} theta_sum / (q;q)_infty ,  extended to Im(tau)<0
\\ by  theta_0(z,tau) = 1/theta_0(z-tau,-tau).  It is 1-periodic in z and
\\ satisfies the standard modular relation under tau -> -1/tau.

THETA0(z, TAU, BOUND)={
local(ETA, temp ,vv);
vv=imag(TAU);
if(vv<0, return(1/THETA0(z-TAU, -TAU, BOUND)));
ETA=eta(TAU, 0);      \\ PARI's eta(,0) is the plain (q;q)_infty product
temp=-I*exp(I*Pi*(z-TAU/4))* newTHETAinit(z, TAU, BOUND)/ETA;
return(temp)
}

\\ ---------------------------------------------------------------------------
\\ lessnaifEuler: with x, q, r of modulus < 1, evaluates the double product
\\   E(x; q, r) = prod_{j,k >= 0} (1 - x q^j r^k)
\\ through its truncated Lambert series
\\   log E(x) = - sum_{n>=1} x^n / ( n (1-q^n)(1-r^n) )
\\ when |x| < 1 - 1e-7, and otherwise through the q-Pochhammer factorization
\\   E(x) = (1-x) prod_{k>=1} (1-x q^k)(1-x r^k) * E(x q r),
\\ which lowers |x| to |x q r| and recurses.
\\ Stops when the current term drops below 10^-BDterm; BBd caps the number
\\ of terms and a diagnostic is printed (once) if the cap is the reason the
\\ loop ends -- callers should then raise BBd.

lessnaifEuler(xxx,q,r, BBd, BDterm)={
local(PRO1,localx, S, term, tterm, TERM1, TERM2);
PRO1=1;S=0;
localx=xxx;
if(abs(q)>1, print("lessnaifEuler: |q| > 1"));
if(abs(r)>1, print("lessnaifEuler: |r| > 1"));
if(abs(localx)<1-10^(-7),
  for(jj=1, BBd,
    tterm=localx^jj/(jj*(1-q^jj)*(1-r^jj));
    S=S+tterm;
    if(abs(tterm)<10^(-BDterm), return(exp(-S))));
  if(!EULERWARN, EULERWARN=1;
     print("lessnaifEuler WARNING: sum truncated by BBd=",BBd,
           " before reaching 10^-",BDterm,"; raise BBd"));
  return(exp(-S)));
PRO1=lessnaifEuler(xxx*q*r,q,r, BBd, BDterm);
PRO1=PRO1*(1-localx);
for(k=1, BBd,
  TERM1= localx*q^k; TERM2=localx*r^k;
  if(max(abs(TERM1),abs(TERM2))<10^(-BDterm), return(PRO1),
     PRO1=PRO1*(1-TERM1)*(1-TERM2)));
if(!EULERWARN, EULERWARN=1;
   print("lessnaifEuler WARNING: product truncated by BBd=",BBd,
         " before reaching 10^-",BDterm,"; raise BBd"));
return(PRO1)
}

\\ ---------------------------------------------------------------------------
\\ testexpfactor: the cubic polynomial Q(z,tau,sigma) in the exponential
\\ factor e^{i Pi Q} of the SL_3 three-term relation (see header).

testexpfactor(z,TAU,SIGMA)={
local(Q);
Q=z^3/(3*TAU*SIGMA)-(TAU+SIGMA-1)/(2*TAU*SIGMA)*z^2
  +(TAU^2+SIGMA^2+3*TAU*SIGMA-3*TAU-3*SIGMA+1)*z/(6*TAU*SIGMA)
  +1/12*(TAU+SIGMA-1)*(1/TAU+1/SIGMA-1);
return(Q)
}

\\ ---------------------------------------------------------------------------
\\ split: when both imaginary parts are small (< SPLITBOUND), rewrite Gamma
\\ through the three-term relation; the two Gammas on the right have larger
\\ imaginary parts (-1/tau, -1/sigma directions).  Real parts of tau, sigma
\\ are reduced mod 1 first (exact for Gamma).

split(z, toto, SSIG, BOUND, BDterm)={
local(TAU, SIGMA, BB, CC);
TAU=toto-round(real(toto));
SIGMA=SSIG-round(real(SSIG));
BB=lessnaifGamma((z-SIGMA)/TAU, -1/TAU, -SIGMA/TAU, BOUND, BDterm);
CC=lessnaifGamma(z/SIGMA, TAU/SIGMA, -1/SIGMA, BOUND, BDterm);
return(CC*exp(-I*Pi*testexpfactor(z,TAU,SIGMA))/(BB))
}

\\ ---------------------------------------------------------------------------
\\ lessnaifGamma: main evaluator.  Route of the moves, in order:
\\   0. reduce z mod 1 (exact);
\\   1. sign normalizations Im(tau) > 0, Im(sigma) > 0 (inversion formulas);
\\   2. ensure Im(sigma) >= Im(tau) (Gamma is symmetric in tau, sigma);
\\   3. if Im(sigma) < SPLITBOUND, normalize through split();
\\   4. translate z by sigma into the strip |Im z| <= Im(sigma)/2, one theta
\\      cocycle factor theta_0(., tau) per step;
\\   5. reflect z -> sigma - z into the upper half strip (Gamma reflection);
\\   6. recenter z by M1*tau near (sigma +- tau)/2, collecting the product
\\      of theta_0(z + k tau, sigma);
\\   7. finish with the two lessnaifEuler evaluations U/V at the
\\      recentered point (Gamma = E(qr/x; q,r) / E(x; q,r)).

lessnaifGamma(z,TAU, SIGMA, BBd, BDterm)={
local(xx,r,q, U, V,a,b,c, DELT, DELT1, zim,xxm,qm,rm, aminus, aplus,
      test5, test6, objectifinf, objectifsup, goal, M1, ww, FACTEUR,
      TAUtemp, SIGtemp, TESTIFYz1, RATIO1, RATIO2, RATIOnew1, RATIOnew2);
z = z - round(real(z));                            \\ exact: Gamma is 1-periodic
if(imag(TAU)<0,
  return(1/lessnaifGamma(z-TAU, -TAU, SIGMA, BBd, BDterm)));
if(imag(SIGMA)<0,
  return(1/lessnaifGamma(z-SIGMA, TAU, -SIGMA, BBd, BDterm)));
a=imag(SIGMA);
zim=imag(z);
b=imag(TAU);
TAUtemp=TAU; SIGtemp=SIGMA;
if(a<=0, print("lessnaifGamma: Im(sigma) <= 0 after normalization"));
if(b<=0, print("lessnaifGamma: Im(tau) <= 0 after normalization"));
if(a<b,
  return(lessnaifGamma(z, SIGtemp, TAUtemp, BBd, BDterm)));
if(a<SPLITBOUND, return(split(z, SIGtemp, TAUtemp, BBd, BDterm)));
if(zim<=-a/2,
  return( lessnaifGamma(z+SIGtemp,TAU, SIGMA, BBd, BDterm)/THETA0(z, TAU, BBd)));
if(zim>a/2,
  DELT1=THETA0(z-SIGtemp, TAU, BBd);
  return( lessnaifGamma(z-SIGtemp, TAU, SIGMA, BBd, BDterm)*DELT1));
\\ here -a/2 < Im(z) <= a/2
TESTIFYz1=(zim<=a/2)*(-a/2<zim);
if(TESTIFYz1==0, print("lessnaifGamma: z escaped the sigma strip"));
if(imag(z)<0,
  RATIO1=lessnaifGamma(SIGMA-z,TAU, SIGMA, BBd, BDterm);
  RATIO2=THETA0(SIGMA-z, SIGMA, BBd);
  return(1/(RATIO1*RATIO2)));
\\ recenter z by multiples of tau near (sigma - tau)/2 .. (sigma + tau)/2
objectifinf= 1/2*(imag(SIGMA)-imag(TAU));
objectifsup= 1/2*(imag(SIGMA)+imag(TAU));
goal=imag(SIGMA/2-TAU/2-z)/imag(TAU);
M1= ceil(goal);
test5=(objectifinf<imag(z+M1*TAU));
test6=(imag(z+M1*TAU)<objectifsup);
if(test5*test6, ,
  RATIOnew1=lessnaifGamma(z-TAU,TAU, SIGMA, BBd, BDterm);
  RATIOnew2=THETA0(z-TAU, SIGMA, BBd);
  return(RATIOnew1*RATIOnew2));
FACTEUR=1;
if(M1>0, for(k=0, M1-1, FACTEUR=FACTEUR*THETA0(z+k*TAU, SIGMA, BBd)), M1=0);
ww= z+M1*TAU;
xxm=exp(2*I*Pi*ww);
qm=exp(2*I*Pi*TAU);
rm=exp(2*I*Pi*SIGMA);
U=lessnaifEuler(qm*rm/xxm,qm,rm, BBd, BDterm);
V=lessnaifEuler(xxm,qm,rm, BBd, BDterm);
return(U/(V*FACTEUR))
}
