# computing-elliptic-units-over-complex-cubic-fields

## Description

This project contains the code used to run the numerical examples of the article [Bergeron, Charollois and García, *Elliptic units for complex cubic fields*](https://arxiv.org/abs/2311.04110).

Let K be a complex cubic field, f an ideal of its ring of integers, and K(f) the narrow ray class field of conductor f. Given a smoothing ideal a and an admissible base point h in K, the code evaluates the elliptic gamma value `Gamma_{a,h}(f b^-1)`, a finite product of values of the Felder–Varchenko elliptic gamma function Γ(z, τ, σ) at points of K³. The main conjecture predicts that this complex number is the image of a unit `u_a(f b^-1)` of K(f), that it does not depend on h, and that `u_a(f b^-1)` is the image of `u_a(f)` under the Artin symbol of b. Letting b run through the narrow ray class group classes, the resulting complex numbers are matched, to 1000 decimal digits, against the roots of a polynomial whose splitting field is K(f).
 
The evaluation rests on the summation formula of Felder–Varchenko for log Γ(z, τ, σ), whose convergence is slow when τ and σ are close to the real axis.
This is the main bottleneck, the other being the length 's' of the finite product of the ellitic Gamma values.  These are the two reasons the  implementation looks for   τ and σ away from the real axis, and wedges of small width 's'  before evaluating.

Two independent implementations are provided: a SageMath chain, with a certified truncation bound, and a GP/PARI 2.15.2 chain, faster but working to an empirical accuracy contract. The field K = Q(∛2) is run in both, and the two agree.

## Results

See [tables](./Tables_of_Examples.pdf) for the commands used, their outputs to 1000 digits, and the polynomials these are recognised as roots of.

The examples covered are the following.

| file | K = Q(β) | disc(K) | h(K) | conductor f | smoothing a | `Cl^+(f)` | paper |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `cubic_field_K_disc_-588.sage` | `β^3 - β^2 + 5β + 1 = 0` | −588 = −2²·3·7² | 3 | (3, β+1) | norm 5 | Z/6 | §2.4 |
| `cubic_field_K_cube_root_7.sage` | `β^3 = 7` | −1323 = −3³·7² | 3 | the ideal of norm 3 | norm 5 | Z/6 | §5.2 |
| `cubic_field_K_cube_root_3.sage` | `β^3 = 3` | −243 = −3⁵ | 1 | the ideal of norm 3 | norm 5, then norm 2 | Z/2 | §5.1 |
| `cubic_field_K_cube_root_2.sage` | `β^3 = 2` | −108 = −2²·3³ | 1 | (3) | norm 5 | Z/6, gen. (β) | §5.3 |
| `setup108.gp` | `β^3 = 2` | −108 | 1 | (3) | norm 5 | Z/6, gen. (β) | §5.3 |
| `setup23c.gp` | `β^3 - β + 1 = 0` | −23 | 1 | (5) | norm 7 | Z/4, gen. (2) | §5.4 |

The first one is Dasgupta's original example, the complex cubic field of smallest discriminant with class number 3; the second is the example of the introduction of the paper. The field of discriminant −108 (Ren–Sczech §4.1) is treated by both implementations, which is how they are cross-checked. The field of minimal discriminant −23 (Ren–Sczech §4.2) is treated in GP/PARI only.

## Instructions

### SageMath

Download `BCG_cubic_git.sage` together with the data file of the example you want to run, and load them in this order from a SageMath session:

```
sage: load('BCG_cubic_git.sage')
sage: load('cubic_field_K_cube_root_2.sage')
```

The second file sets the global variables describing the example: the field `K`, the conductor `ff` and the smallest positive integer `q_ff` it contains, the smoothing ideal `aa`, an ideal `bb` whose class generates the narrow ray class group, the generator `vareps` of the totally positive units in 1+f and its matrix `Eps`, the roots `yRR` and `yCC` defining the real and complex embeddings, and the orientation matrix `or_matrix`. It prints them, with a few consistency checks.

One then chooses an admissible λ for `L = f b^(-k)` and evaluates typically as :

```
sage: L = ff*bb^(-1)
sage: A = find_small_width_admissible_vectors(8, L, 30)
sage: u = EllGammaunit_aa_h(A[0]/q_ff, L)
sage: certified_digits(u)
sage: algdep(u, 18)
```

`admissible_elements(N, L)` returns all admissible elements in the search box of radius `N`, and `find_small_width_admissible_vectors(N, L, max_width)` keeps those whose wedge has small width, which is what governs the cost of one evaluation. Different rows of `A` give different h and should return the same value, one of the checks the conjecture predicts.

`algdep` recognises `u` over Q. The polynomial over K printed in the paper is recovered by factoring that polynomial over K, or directly by a `lindep` of the numbers `β^i u^j`. In the `β^3 = 3` example the coefficients are already rational and `algdep(u, 6)` suffices.

The Galois conjugates are obtained by letting b vary:

```
sage: for k in range(6):
....:     L = ff * bb^(-k)
....:     A = find_small_width_admissible_vectors(8, L, 30)
....:     print(EllGammaunit_aa_h(A[0]/q_ff, L))
```

The working precision and the truncation bounds are set among the global variables at the top of `BCG_cubic_git.sage`, and `certified_digits` returns the number of decimal digits the truncation error bound guarantees. Reaching 1000 certified digits takes a while, the run time being dominated by the width of the wedge and by the decay rate ρ of the series.

### GP/PARI

`Baby_step_GammaFelder33.gp` evaluates the elliptic gamma function itself, through SL₃(Z) moves and theta cocycle factors; its public interface is `lessnaifGamma(z, TAU, SIGMA, BBd, BDterm)`. `setup108.gp` and `setup23c.gp` each set up one example and run all of its narrow ray classes. The three files must sit in the same directory, and the setup scripts are meant to be run directly:

```
$ gp -q setup108.gp        # six classes, disc -108, BCG §5.3
$ gp -q setup23c.gp        # four classes, disc  -23, BCG §5.4
```

Precision is set by `TARGETDIGITS` at the top of each setup file, 25 by default; raise it to 1000 to reproduce the tables. Each run prints the wedge data and the invariant of every class, locates that invariant among the roots of the polynomial of the paper (the sextic of §5.3, the quartic of §5.4), and writes it to `invariant_disc-108_class<k>.txt`. `runclass(k)` runs a single class.

At 1000 digits the six classes of disc −108 take a few minutes. The wedges of disc −23 are far wider, and two of its classes run for hours. The remaining comments, on the marking normalization `FACTOR = -rhosign`, on the accuracy contract of the evaluator and on the `parisizemax` pitfall, are in the headers of the files.

## Contributors

Nicolas Bergeron, Pierre Charollois and Luis E. García.

## License

GPL-3.0, see [LICENSE](./LICENSE). If you use this code, please cite the article and this repository, see [CITATION.cff](./CITATION.cff).
