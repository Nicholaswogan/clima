# Dry radiative-convective equilibrium: smooth fluxes and convective adjustment

This note accompanies [`smooth_rce_toy.py`](smooth_rce_toy.py). The script is a
self-contained experiment for comparing four ways to calculate a dry,
one-dimensional radiative-convective equilibrium (RCE):

1. a smooth, finite-rate convective flux with automatic continuation in its
   strength;
2. conventional energy-conserving convective adjustment after explicit radiative
   steps; and
3. exact adjustment coupled to radiation in a constrained backward-Euler
   timestep; and
4. the same exact adjustment written as a projected steady residual and solved
   with pseudo-transient continuation (PTC).

All four formulations approach the same dry RCE, but they are no longer equally
attractive. If convection is intended to be instantaneous, weighted PAVA is
cleaner than a smooth finite-$K$ flux: it has no convective-strength parameter,
conserves thermal energy to roundoff, and enforces dry stability exactly. The
smooth flux remains useful only when a genuinely finite-rate, differentiable
convective closure is desired.

Projected PTC is the leading steady-state method from these experiments. It solves
the exact adjustment fixed point without physically waiting for slow radiative or
ocean modes, and it used far fewer full radiative-transfer (RT) evaluations than
explicit relaxation. It is a leading candidate rather than a universal winner:
real opacity feedbacks, difficult active-set changes, moist adjustment, and
production-scale Jacobian approximations still need to be tested.

For time-accurate evolution, the recommendation is conditional. Explicit
radiative stepping followed by PAVA decisively wins the present toy benchmark
because its stable timestep is already comparable to the timestep required for
accuracy. Constrained backward Euler is mathematically sound but currently much
more expensive. Higher order, analytical projection derivatives, and Jacobian
reuse could narrow that gap substantially, and implicit integration may become
preferable in atmospheres where fast radiative modes make explicit stability much
more restrictive than physical accuracy.

The toy deliberately excludes chemistry, condensation, latent heating, clouds,
scattering, and temperature-dependent opacity. It tests numerical formulations,
not a complete habitable-planet climate model.

## Current method recommendations

| Objective | Preferred method | Current interpretation |
|---|---|---|
| Dry steady RCE | Projected PTC on the PAVA natural residual | Best steady solver tested; bypasses slow physical relaxation |
| Time-accurate evolution without strong stability restrictions | Explicit radiation followed by PAVA | Simplest method and clear winner in the fixed-time toy benchmark |
| Time-accurate evolution with strong radiative stiffness | Constrained implicit adjustment, after optimization | Promising for extreme problems, but not yet competitive in this implementation |
| A deliberately finite-rate differentiable convection model | Smooth finite-$K$ flux | Valid closure, but $K$ is physical/modeling input rather than a route to exact adjustment |

Weighted PAVA is the preferred dry instantaneous-convection operator in the first
three rows. The numerical method used to couple it to radiation should be selected
according to whether the objective is a steady root or a time-resolved trajectory.
None of these dry results establishes a complete moist-adjustment method.

## Running the experiment

For continuity with the original experiment, the script default remains the
smooth hybrid-continuation calculation; this is not the current recommendation
for exact dry convective neutrality:

```bash
python smooth_rce_toy.py
```

The three exact-adjustment calculations are

```bash
python smooth_rce_toy.py --continuation-mode projected
python smooth_rce_toy.py --continuation-mode projected-explicit
python smooth_rce_toy.py --continuation-mode projected-implicit
```

For terminal-only runs, add `--no-plot`. The available modes are:

| Mode | Convection | Integration or steady solver |
|---|---|---|
| `staged` | Smooth finite-$K$ flux | Fully converge each $K$ stage with PTC |
| `coupled` | Smooth finite-$K$ flux | Increase $K$ during one PTC solve |
| `hybrid` | Smooth finite-$K$ flux | Guard transient lapse rates and increase $K$ near equilibrium |
| `projected-explicit` | Exact dry adjustment | Forward-Euler radiation, then adjustment |
| `projected-implicit` | Exact dry adjustment | Constrained backward Euler with step doubling |
| `projected` | Exact dry adjustment | PTC on a projected natural residual |

An apples-to-apples transient benchmark integrates the exact-adjustment methods
from the same initial state to the same physical end time:

```bash
python smooth_rce_toy.py --time-benchmark --no-plot
```

It compares fixed-step explicit runs and error-controlled constrained backward
Euler against a timestep-refined explicit reference. The benchmark end time and
reference timestep are controlled by `--benchmark-end-time` and
`--benchmark-reference-dt`.

Important controls are:

```text
--nlev N
--k-conv-initial K
--k-conv-max KMAX
--superadiabatic-tolerance S
--superadiabatic-guard S_GUARD
--projection-time ALPHA
--explicit-dt-max DTMAX
--explicit-max-temperature-step DTEMP
--implicit-dt-initial DT
--implicit-dt-max DTMAX
--implicit-temperature-tolerance DTEMP
--implicit-flux-tolerance DFLUX
--implicit-max-steps NSTEPS
--time-benchmark
--benchmark-end-time TEND
--benchmark-reference-dt DTREF
```

`--k-conv` remains an alias for `--k-conv-initial`.

## Shared finite-volume climate model

Pressure increases downward and net flux is positive upward. The dry atmospheric
energy equation in pressure coordinates is

$$
\frac{c_p}{g}\frac{\partial T}{\partial t}
=\frac{\partial F^{\mathrm{tot}}}{\partial p},
\qquad
F^{\mathrm{tot}}=F^{\mathrm{rad}}+F^{\mathrm{conv}}.
$$

RT is a nonlocal diagnostic operator,

$$
F^{\mathrm{rad}}=\mathcal R[\mathbf T,T_s],
$$

so the continuous model is more accurately an integro-differential equation than
a local diffusion PDE. After vertical discretization it is a nonlinear ODE system.

For pressure cell $i$,

$$
C_i=\frac{c_p\Delta p_i}{g},
\qquad
C_i\frac{dT_i}{dt}
=F^{\mathrm{tot}}_{i+1/2}-F^{\mathrm{tot}}_{i-1/2}.
$$

One shared flux at each interface makes interior energy exchange conservative. The
surface is prognostic as well:

$$
C_s\frac{dT_s}{dt}=-F^{\mathrm{tot}}_{N+1/2}.
$$

The default $C_s=4\times10^8$ J m$^{-2}$ K$^{-1}$ is approximately the heat
capacity of a 100 m ocean. It introduces a slow mode that is useful for testing
stiff solvers.

### Gray radiative transfer

Longwave optical depth is prescribed at cell edges as

$$
\tau(p)=\tau_s
\left(\frac{p-p_{\mathrm{top}}}{p_s-p_{\mathrm{top}}}\right)^n,
$$

with default values $\tau_s=4$ and $n=2$. For layer optical depth $\Delta\tau_i$,

$$
t_i=\exp(-D\Delta\tau_i), \qquad D=1.66.
$$

The pure-absorption, isothermal-layer recurrences are

$$
F^\uparrow_{i-1/2}
=t_iF^\uparrow_{i+1/2}+(1-t_i)\sigma T_i^4,
$$

$$
F^\downarrow_{i+1/2}
=t_iF^\downarrow_{i-1/2}+(1-t_i)\sigma T_i^4.
$$

The surface emits $\sigma T_s^4$ and the downward longwave flux at the top is zero.
The default shortwave flux passes through the atmosphere and deposits 240 W
m$^{-2}$ at the surface. Thus

$$
F^{\mathrm{rad}}
=F^\uparrow_{\mathrm{LW}}-F^\downarrow_{\mathrm{LW}}
-F^\downarrow_{\mathrm{SW}}.
$$

This intentionally crude RT scheme supplies a smooth, nonlocal temperature-to-flux
map with correct column-energy bookkeeping.

### Dry stability

For a dry ideal gas,

$$
\nabla=\frac{\partial\ln T}{\partial\ln p},
\qquad
\nabla_{\mathrm{ad}}=\frac{R}{c_p}.
$$

At an internal interface the discrete superadiabaticity is

$$
s_{i+1/2}
=\frac{\ln T_{i+1}-\ln T_i}{\ln p_{i+1}-\ln p_i}
-\nabla_{\mathrm{ad}}.
$$

Dry stability requires $s\le0$. The default $c_p=1506$ J kg$^{-1}$ K$^{-1}$ gives
$\nabla_{\mathrm{ad}}\simeq0.1906$. This deliberately shallow dry adiabat creates a
substantial convective region; it should not be interpreted as modern Earth's
physical dry heat capacity.

## Formulation A: a smooth finite-rate convective flux

The smooth model closes the temperature ODE with

$$
F^{\mathrm{conv}}_{i+1/2}
=\rho_{i+1/2}c_pK_{\mathrm{conv}}\frac{g}{R}A_\epsilon(s_{i+1/2}),
$$

where $K_{\mathrm{conv}}$ is an eddy diffusivity and $A_\epsilon$ is a one-sided
$C^1$ Huber activation:

$$
A_\epsilon(s)=
\begin{cases}
0, & s\le0,\\[3pt]
s^2/(2\epsilon), & 0<s<\epsilon,\\[3pt]
s-\epsilon/2, & s\ge\epsilon.
\end{cases}
$$

The boundary convective fluxes are zero, so this closure only redistributes column
energy. It is an engineered potential-temperature diffusion, not mixing-length
theory and not a prediction of convective velocity.

The exact zero for $s\le0$ matters. An earlier softplus activation had
$\operatorname{softplus}_\epsilon(0)=\epsilon\ln2>0$; increasing $K$ then produced
an increasingly large artificial flux even at neutral stability.

### Why $K_{\mathrm{conv}}$ cannot be chosen once and forgotten

In the approximately linear active region,

$$
F^{\mathrm{conv}}
\simeq \rho c_pK_{\mathrm{conv}}\frac{g}{R}s,
\qquad
s\simeq
\frac{F^{\mathrm{conv}}R}{\rho c_pK_{\mathrm{conv}}g}.
$$

A finite flux carried by finite $K$ therefore requires finite
superadiabaticity. The profile approaches an exact adiabat only as $K\to\infty$,
while the ODE becomes increasingly stiff. Here $K$ is best interpreted as a
penalty strength selected to meet a lapse-rate tolerance, not necessarily as a
physical eddy diffusivity.

### PTC for the smooth ODE

Writing the temperature tendency as $\mathbf f(\mathbf T)$, PTC solves

$$
\left(\frac{I}{\Delta\tau}-J_f\right)\delta\mathbf T
=\mathbf f(\mathbf T),
\qquad
\mathbf T\leftarrow\mathbf T+\delta\mathbf T.
$$

Small $\Delta\tau$ resembles a damped implicit time step; large $\Delta\tau$
approaches Newton's method. Rejected trials reduce $\Delta\tau$ while reusing the
same Jacobian. Accepted trials adjust it using residual progress. The trajectory is
pseudo-time and should not be interpreted as physical climate evolution.

The toy currently uses a dense central-difference Jacobian. The steady convergence
test is based on maximum flux imbalance, not temperature tendency, because the
layer and ocean heat capacities span a wide range.

### Three $K$-continuation strategies

The user specifies a target

$$
s_{\max}=\max_i(s_i)\le s_{\mathrm{target}}
$$

rather than manually selecting the final $K$.

- **Staged:** converge RCE at fixed $K$, estimate the next $K$ from
  $s_{\max}\propto K^{-1}$, and warm-start the next stage. This is simple and
  provides a useful reference, but its preliminary radiative equilibrium can be
  extremely superadiabatic.
- **Coupled:** increase $K$ during one PTC integration after sufficient residual
  progress, and reduce the pseudo-timestep when $K$ grows. This can save accepted
  steps but is more controller-sensitive.
- **Hybrid:** begin with weak convection, increase $K$ if an accepted state reaches
  a transient superadiabatic guard, and also increase it whenever the current
  fixed-$K$ equilibrium remains above the final target. A harder guard rejects an
  excessive candidate before accepting it. This avoids both a separate radiative
  equilibrium and unrestricted superadiabatic transients.

The default hybrid guard is

$$
s_{\mathrm{guard}}
=\max(10s_{\mathrm{target}},\,0.05\nabla_{\mathrm{ad}}).
$$

Setting the guard very large reduces the hybrid behavior approximately to staged
equilibrium-triggered continuation, except that it still starts from the initial
profile instead of first calculating radiative equilibrium.

## Formulation B: exact energy-conserving convective adjustment

Dry potential temperature can be written

$$
\theta_i=\frac{T_i}{\Pi_i},
\qquad
\Pi_i=\left(\frac{p_i}{p_s}\right)^{\nabla_{\mathrm{ad}}}.
$$

Because the grid is ordered from the top toward increasing pressure, dry stability
requires

$$
\theta_1\ge\theta_2\ge\cdots\ge\theta_N.
$$

The adjustment operator $P$ is a weighted pool-adjacent-violators algorithm
(PAVA). It finds the nearest stable sequence in weighted potential-temperature
space:

$$
P(\theta^*)=
\arg\min_{\theta_1\ge\cdots\ge\theta_N}
\frac12\sum_i w_i(\theta_i-\theta_i^*)^2,
\qquad
w_i=C_i\Pi_i.
$$

Whenever neighboring blocks violate the ordering, they are pooled and assigned

$$
\bar\theta_B
=\frac{\sum_{i\in B}C_i\Pi_i\theta_i^*}
{\sum_{i\in B}C_i\Pi_i}.
$$

Since $T_i=\Pi_i\theta_i$, this exactly conserves thermal energy in every pooled
block:

$$
\sum_{i\in B}C_iT_i^{\mathrm{adjusted}}
=\sum_{i\in B}C_iT_i^*.
$$

PAVA is $O(N)$, deterministic, and removes the need to choose
$K_{\mathrm{conv}}$. It formalizes conventional dry convective adjustment as a
projection/optimization problem rather than an informal profile repair.

At steady state this corresponds to the complementarity conditions

$$
F^{\mathrm{conv}}\ge0,
\qquad
\nabla_{\mathrm{ad}}-\nabla\ge0,
\qquad
F^{\mathrm{conv}}(\nabla_{\mathrm{ad}}-\nabla)=0.
$$

The exact-adjustment problem is therefore naturally a variational inequality or
complementarity problem. It need not be represented as a general-purpose DAE: the
projection eliminates the convective-flux multipliers from the state solve.

### Diagnosing the convective flux

PAVA determines temperatures but does not explicitly evolve a convective flux. At
equilibrium the flux can be reconstructed inside each neutral mixed block. If the
block begins at edge $a$, then at an internal edge $j$

$$
F^{\mathrm{conv}}_j=F^{\mathrm{rad}}_a-F^{\mathrm{rad}}_j.
$$

The total flux is consequently constant through the block. Convective flux is zero
at the block boundaries, and the remaining block-integrated radiative imbalance is
placed in its final cell. At a converged RCE that final imbalance also vanishes.

## Explicit radiative steps followed by adjustment

The conventional exo_k-like split algorithm is

$$
\mathbf T^*=\mathbf T^n+\Delta t\,\mathbf f_{\mathrm{rad}}(\mathbf T^n),
\qquad
\mathbf T^{n+1}=P(\mathbf T^*).
$$

Only atmospheric temperatures are projected; the ocean temperature takes the same
explicit radiative step but is not part of PAVA. Each accepted step needs one full
RT evaluation and one inexpensive adjustment.

The toy grows one global timestep up to `--explicit-dt-max`. A trial is retried
with a smaller timestep if the post-adjustment temperature change exceeds
`--explicit-max-temperature-step` or leaves broad physical temperature bounds.
Limiting the post-adjustment change is less restrictive than limiting raw
layer-by-layer radiative heating, because adjustment legitimately redistributes
that heating through a convective block. Individual layer increments are not
clipped, so the projected forward-Euler map is preserved.

The explicit solver requires both a small diagnosed flux imbalance and a small
projected change per unit time. The current map-tendency threshold is $10^{-13}$ K
s$^{-1}$. A preliminary $10^{-10}$ K s$^{-1}$ threshold stopped with the uppermost
temperature still 0.27 K from the PTC result, even though the flux criterion had
passed. The optically thin top cell has so little radiative leverage that apparently
small flux and tendency residuals can still permit a visible temperature error.

This method is attractive because it is transparent, conservative, and has no
nonlinear linear solve. Its main weakness is the explicit stability limit. A large
timestep can oscillate indefinitely even while remaining bounded, and a very slow
optically thin or ocean mode can require many RT calls after the bulk profile looks
converged.

## Constrained backward Euler

The implicit exact-adjustment timestep is

$$
\mathbf T^{n+1}
=P\!\left[\mathbf T^n+\Delta t\,
\mathbf f_{\mathrm{rad}}(\mathbf T^{n+1})\right].
$$

It is computed by solving the nonsmooth residual

$$
\mathbf H(\mathbf Y)
=\mathbf Y-P\!\left[\mathbf T^n+\Delta t\,
\mathbf f_{\mathrm{rad}}(\mathbf Y)\right]=0.
$$

This is fundamentally different from applying PAVA after an arbitrary implicit
radiative or Newton correction: the radiation operator is evaluated at the new
constrained state inside the equation defining the timestep. The surface is
implicit but unprojected. The atmospheric PAVA step conserves the energy of the
implicit radiative trial state.

The toy solves $\mathbf H=0$ with damped Newton iterations and a one-sided
finite-difference Jacobian through both RT and PAVA. The projection is piecewise
linear, so active-set changes make the residual semismooth. Invalid trials or a
failed nonlinear solve reject the physical timestep.

Time adaptivity uses one full step and two half steps. Their maximum temperature
difference estimates the local error; the two-half-step solution is accepted when

$$
\left\|\mathbf T_{\Delta t}
-\mathbf T_{\Delta t/2,\Delta t/2}\right\|_\infty
\le \epsilon_T.
$$

This costs three nonlinear solves per attempted macro-step, but cleanly separates
temporal accuracy from nonlinear convergence. Small-step tests show that the
backward-Euler versus forward-Euler difference grows by approximately a factor of
four when $\Delta t$ doubles, as expected for the $O(\Delta t^2)$ local difference
between two first-order methods.

The method is time-accurate for the ideal instantaneous-adjustment evolution
defined by this weighted projection. That does not make the transient a unique
model of real three-dimensional convection; the projection and its weighting are
part of the physical idealization.

## PTC on the projected natural residual

Define the projected map for any positive scaling time $\alpha$:

$$
M_\alpha(\mathbf T)
=P\!\left[\mathbf T+\alpha\mathbf f_{\mathrm{rad}}(\mathbf T)\right].
$$

The surface component is left unprojected. The natural residual is

$$
\mathbf G_\alpha(\mathbf T)
=\frac{M_\alpha(\mathbf T)-\mathbf T}{\alpha}.
$$

A root of $\mathbf G_\alpha=0$ is a fixed point of energy-conserving convective
adjustment. For this convex dry constraint, the root represents the same
complementarity equilibrium for any valid positive $\alpha$; $\alpha$ changes
residual scaling and active-set behavior, not the desired physical solution.

The projected solver applies PTC directly to $\mathbf G$:

$$
\left(\frac{I}{\Delta\tau}-J_G\right)\delta\mathbf T
=\mathbf G(\mathbf T).
$$

Every accepted iterate is projected back into the stable set. The toy computes a
one-sided finite-difference Jacobian through both RT and PAVA. The residual is
continuous and piecewise differentiable; active-set changes make it semismooth
rather than globally smooth. Candidate acceptance is therefore relaxed modestly
when the mixed-block structure changes.

The PTC pseudo-timestep $\Delta\tau$ and the projection time $\alpha$ have distinct
roles. The former globalizes the nonlinear solve. The latter scales the natural
residual and influences which convective constraints become active in the
radiative trial state. They should not be increased together as though they were
the same continuation parameter.

### Automatic selection of the projection time

The code now chooses one $\alpha$ from the initially projected profile and holds
it fixed for the complete PTC solve. The default rule limits the largest relative
unconstrained radiative trial change:

$$
\alpha
=\frac{\eta}{\displaystyle\max_i
\left(
\frac{|f_{\mathrm{rad},i}|}{\max(|T_i|,T_{\mathrm{floor}})}
\right)},
$$

with $\eta=5\times10^{-3}$ and $T_{\mathrm{floor}}=100$ K. Thus, the initial
radiative trial changes no component by more than approximately $0.5\%$ in the
scaled norm. A separate cooling bound keeps the trial state conservatively above
50 K. The command-line option `--projection-relative-change` changes $\eta$;
providing `--projection-time` bypasses automatic selection.

It is important to freeze $\alpha$ while constructing and using a Jacobian.
Changing it every accepted PTC step would change the nonlinear residual at the
same time that the pseudo-timestep and active set are changing. A production
solver could change $\alpha$ only on a controlled restart, invalidating the old
Jacobian, and should freeze it completely near convergence.

For the default 40-layer problem the rule selects
$\alpha=3.0212\times10^4$ s. It converges in 31 accepted steps with no rejections,
the same performance as the former fixed $10^4$ s value. The two final profiles
differ by only $4.4\times10^{-7}$ K. Targets from $10^{-3}$ through
$5\times10^{-2}$ all required 31 steps; the deliberately small $10^{-4}$ target
required 36 steps. Grid tests with 8, 20, and 80 layers and additional thin-gray,
thick-gray, high-forcing, and low-$c_p$ cases also converged without manual
selection.

After convergence the code reevaluates the same state at $\alpha/10$, $\alpha$,
and $10\alpha$. In the default case all three checks have the same convective
block, maximum flux imbalance $6.60\times10^{-8}$ W m$^{-2}$, and maximum natural
residual $8.37\times10^{-14}$ K s$^{-1}$. This is a useful production diagnostic:
for an exact dry convex projection, the converged complementarity state should not
depend on a positive scalar $\alpha$.

The tests also expose an important limitation. An artificial case with
$c_p=4000$ J kg$^{-1}$ K$^{-1}$ fails with the default automatically selected
$\alpha\simeq8.0\times10^4$ s under the toy's simple active-set acceptance rule,
whereas a $2\%$ target selects $\alpha\simeq3.2\times10^5$ s and converges in 32
steps. Dependence is nonmonotonic: making $\alpha$ smaller does not necessarily
repair convergence. The relative-change rule is therefore a good initial scaling,
not yet a production-grade globalization method. A robust implementation should
combine it with a semismooth line search or trust region and permit a controlled
restart with a larger or smaller $\alpha$ after genuine active-set stagnation.

### Longer-term option: metric preconditioning

A single scalar $\alpha$ asks every layer and the surface to use the same scaling
time. That may be inadequate for a real column spanning many orders of magnitude
in pressure, heat capacity, and radiative timescale. The mathematically clean
generalization is a positive-definite scaling operator $H$:

$$
\mathbf G_H(\mathbf T)
=H^{-1}\left\{
P_H[\mathbf T+H\mathbf f_{\mathrm{rad}}(\mathbf T)]-\mathbf T
\right\},
$$

where $P_H$ is the convective projection in the corresponding $H^{-1}$ metric.
The diagonal entries of $H$ could represent layer-dependent temperature or
radiative-timescale scaling, preventing one fast, low-mass upper layer from
forcing an unhelpfully small scalar $\alpha$ for the entire atmosphere.

The matching metric is essential. Simply using unrelated layerwise values
$\alpha_i$ inside the trial state while retaining the original energy-weighted
PAVA projection can change the fixed-point/complementarity problem and therefore
the climate solution. Metric preconditioning must be derived together with the
projection and its conservation constraints, then checked for fixed-point
equivalence, idempotence, and energy conservation.

It is not worth adding this complexity to the toy yet. Scalar automatic selection
has a broad successful plateau in the default and most perturbed cases, while the
high-$c_p$ failure is at least partly caused by the deliberately simple semismooth
globalization. The next priority should be a production-quality active-set-aware
line search or trust region and tests with real RT. Metric preconditioning becomes
worth considering if those tests show that a few fast layers repeatedly dictate
$\alpha$, make the projected Jacobian poorly conditioned, or cause scalar-alpha
restarts across otherwise reasonable atmospheric types.

### Why an implicit radiative step followed by projection is not enough

We also tested a tempting split method: calculate a linearly implicit or Newton
radiative correction and then project the corrected profile. It did not converge
to the RCE complementarity solution. In the large-pseudo-timestep limit it projects
a radiative-Newton target; that map is not
$P[\mathbf T+\alpha\mathbf f_{\mathrm{rad}}(\mathbf T)]$ and has the wrong fixed
points.

Thus, implicit acceleration must be applied to the projected residual itself (or
to an equivalent active-set/complementarity system). Projection is not generally a
valid post-processing operation on an arbitrary Newton correction.

## Numerical results

All values below use the current 40-layer defaults unless noted otherwise.

### Baseline and smooth-$K$ behavior

Pure radiative equilibrium has a surface temperature of 359.63 K and maximum
superadiabaticity of about 0.166. It is a valid radiative solution but a poor
intermediate climate profile.

The default hybrid smooth solve starts at $K=10$ m$^2$ s$^{-1}$ and obtains:

| Quantity | Result |
|---|---:|
| Accepted PTC steps | 66 |
| Rejected attempts | 12 |
| Guard-triggered $K$ increases | 14 |
| Equilibrium-triggered increases | 1 |
| Final $K$ | 3117 m$^2$ s$^{-1}$ |
| Final maximum superadiabaticity | $8.93\times10^{-4}$ |
| Maximum upward convective flux | 119.51 W m$^{-2}$ |
| Surface temperature | 337.481 K |

Tightening the staged smooth penalty demonstrates convergence toward exact
adjustment:

| Smooth target | Final $K$ [m$^2$ s$^{-1}$] | Surface temperature [K] | Maximum $s$ |
|---:|---:|---:|---:|
| $10^{-3}$ | 3,247 | 337.476744 | $8.58\times10^{-4}$ |
| $10^{-4}$ | 48,711 | 337.373074 | $8.08\times10^{-5}$ |
| $10^{-5}$ | 3,711,912 | 337.362643 | $8.56\times10^{-6}$ |
| Exact projection | -- | 337.361368 | $3.66\times10^{-15}$ |

The escalating $K$ and solve cost illustrate why an exact constraint becomes
attractive when a nearly neutral profile is required.

### Exact-adjustment methods

The explicit and projected-PTC methods converge to the same surface temperature,
and their full temperature profiles agree within $4.5\times10^{-5}$ K under their
current stopping criteria. Constrained backward Euler agrees in the dynamically
important lower column but, with its looser default physical-integration stopping
tolerance, stops earlier along the slow optically thin upper-atmosphere mode.

| Method | Accepted steps | Full RT evaluations | Jacobian evaluations | Surface temperature [K] |
|---|---:|---:|---:|---:|
| Explicit adjustment, $\Delta t_{\max}=10^6$ s | 23,971 | 23,972 | 0 | 337.361368 |
| Constrained backward Euler, $\epsilon_T=0.02$ K | 548 | 111,705 | 2,598 | 337.361368 |
| Projected PTC | 31 | 1,335 | 31 | 337.361368 |

The PTC count includes every RT call used by its one-sided finite-difference
Jacobians and candidate residuals. Despite rebuilding a dense Jacobian at every
accepted step, it uses about 18 times fewer full RT evaluations than the
conservative explicit configuration.

The constrained backward-Euler run conserves each PAVA projection to about
$8\times10^{-16}$ relative error and remains stable while its physical timestep
grows far beyond the explicit stability boundary. It is nevertheless the most
expensive calculation here because step doubling requires three nonlinear solves
per attempted macro-step and every finite-difference Newton Jacobian uses a full
set of RT calls. Jacobian reuse or an analytical PAVA derivative would reduce this
cost substantially.

Its default $10^{-5}$ W m$^{-2}$ stopping tolerance leaves the very optically thin
top cell about 0.7 K from the more tightly converged PTC profile, even though the
surface, convective region, TOA balance, and surface temperature agree. Tightening
`--implicit-flux-tolerance` continues the physical integration toward the same
upper-atmosphere solution, but exposes the genuinely long radiative timescale.
This is an important distinction: an implicit time integrator removes stability
limits but does not remove slow physical modes when time accuracy is retained.

The explicit timestep is consequential:

| Maximum explicit timestep | Outcome |
|---:|---|
| $5\times10^5$ s | Converged in about 47,800 RT calls |
| $10^6$ s | Converged in about 24,000 RT calls |
| $2\times10^6$ s | Converged in about 12,100 RT calls |
| $3\times10^6$ s | Persistent oscillation; did not converge |

This is the expected explicit tradeoff: a conservative default is robust but slow,
while the fastest stable timestep lies close to a problem-dependent stability
boundary. The final slow convergence came mainly from an optically thin upper-layer
mode that remained after the surface and convective region appeared equilibrated.

### Fixed-time apples-to-apples comparison

The `--time-benchmark` calculation removes equilibrium stopping criteria from the
comparison. Every candidate starts from the same projected profile and ends at
$t=10^8$ s. The reference is projected forward Euler with $\Delta t=1.25\times
10^4$ s. Comparing it with $\Delta t=2.5\times10^4$ s changes the maximum
atmospheric temperature by only $2.22\times10^{-3}$ K, comfortably below the
candidate errors.

| Method and control | Steps | Full RT evaluations | Maximum atmospheric error [K] | Atmospheric RMS error [K] | Surface error [K] |
|---|---:|---:|---:|---:|---:|
| Explicit, $\Delta t=2.5\times10^5$ s | 400 | 401 | 0.0423 | 0.00481 | 0.00339 |
| Explicit, $\Delta t=5\times10^5$ s | 200 | 201 | 0.0872 | 0.00944 | 0.00585 |
| Explicit, $\Delta t=10^6$ s | 100 | 101 | 0.178 | 0.0178 | 0.00688 |
| Explicit, $\Delta t=2\times10^6$ s | 50 | 51 | 0.379 | 0.0610 | 0.0590 |
| Implicit, $\epsilon_T=0.2$ K | 49 | 15,653 | 0.510 | 0.118 | 0.135 |
| Implicit, $\epsilon_T=0.05$ K | 87 | 28,279 | 0.277 | 0.0640 | 0.0734 |
| Implicit, $\epsilon_T=0.02$ K | 131 | 38,325 | 0.175 | 0.0401 | 0.0459 |

At nearly identical maximum atmospheric error, explicit $\Delta t=10^6$ s and
implicit $\epsilon_T=0.02$ K take 100 and 131 accepted macro-steps, respectively.
The current implicit code is vastly more expensive because each macro-step uses a
full step plus two half steps, each nonlinear iteration rebuilds a dense
finite-difference Jacobian, and rejected attempts repeat that work.

More importantly, the comparable step counts show that this particular toy
trajectory does not yet exhibit a large separation between the explicit stability
limit and the timestep required for accuracy. An optimized BDF2 or TR-BDF2 method
could still become competitive by reducing temporal error per step, using the
analytical block derivative of PAVA, and reusing radiative Jacobians. This test does
not support replacing the simple explicit reference yet; it supplies the baseline
against which those optimizations must demonstrate an advantage.

### Conservation and resolution

Random unstable profiles at 8, 20, 40, and 80 layers were projected to machine
precision. The largest measured relative thermal-energy error was approximately
$6\times10^{-16}$.

Projected PTC gave:

| Layers | Accepted steps | Rejections | Surface temperature [K] | Maximum flux imbalance [W m$^{-2}$] |
|---:|---:|---:|---:|---:|
| 20 | 31 | 0 | 334.2974 | $4.78\times10^{-8}$ |
| 40 | 31 | 0 | 337.3614 | $6.65\times10^{-8}$ |
| 80 | 31 | 0 | 338.3749 | $7.19\times10^{-8}$ |
| 120 | 31 | 0 | 338.5821 | $7.34\times10^{-8}$ |

The nonlinear iteration count is essentially grid-independent in this test.
Surface temperature still changes with resolution, especially from 20 to 40
layers, so this table demonstrates solver behavior rather than complete spatial
convergence.

## What the comparison establishes

1. A single conservative temperature equation can combine RT with a smooth
   convective flux and can be solved robustly with PTC.
2. Smooth finite-rate convection is an ODE closure, but a finite $K$ necessarily
   leaves finite superadiabaticity.
3. Standard dry convective adjustment can be formalized as a conservative weighted
   projection rather than an ad hoc profile correction.
4. Explicit adjustment, constrained backward Euler, and projected PTC approach the
   same dry constrained equilibrium when converged tightly.
5. Explicit adjustment is an excellent correctness baseline but can be limited by
   stability and slow thermal modes.
6. Constrained backward Euler removes the explicit stability restriction while
   preserving a first-order physical-time trajectory, but does not accelerate
   genuinely slow radiative modes.
7. In the present fixed-time benchmark, explicit stability is not substantially
   more restrictive than accuracy, so the unoptimized implicit method has no
   efficiency advantage.
8. PTC on the projected residual is substantially more RT-efficient for steady
   RCE in this toy.
9. An implicit radiative correction followed by projection is not equivalent to
   solving the projected steady problem.
10. Scaling $\alpha$ from a target relative radiative trial change removes manual
    selection in the tested baseline and most perturbed cases, but does not replace
    semismooth globalization; the high-$c_p$ counterexample demonstrates that its
    convergence behavior can be nonmonotonic.
11. Layer-dependent metric preconditioning is a plausible response to extreme
    timescale separation, but must be derived with a matching projection metric to
    preserve the intended complementarity equilibrium.
12. Neither the smooth nor exact dry formulation in this script yet solves moist
   convection.

## Moving to production radiative transfer

The most useful next experiment is a dry, fixed-composition implementation using
the production RT operator. Projected PTC should be treated as the leading steady
solver, while explicit adjustment should be retained as both the reference
trajectory and the initial time-accurate implementation. Constrained implicit
integration should be pursued when a fixed-time benchmark demonstrates that the
explicit stability limit is substantially smaller than the timestep required for
accuracy. The smooth hybrid method remains useful only when a finite-rate
differentiable convective closure is itself desired.

### When implicit time integration may become preferable

The gray toy does not exhibit a large separation between its explicit stability
limit and its accuracy-limited timestep. That result should not be generalized to
all one-dimensional atmospheres. A dry gas giant spanning roughly $10^{-6}$ to
$10^3$ bar, for example, combines nine orders of magnitude in pressure with
temperatures ranging from a few hundred kelvin aloft to perhaps 2000 K at depth.

On a logarithmic pressure grid, layer heat capacity scales approximately as

$$
C_i=\frac{c_p\Delta p_i}{g}\propto p_i,
$$

while local thermal-emission derivatives scale as $4\sigma T^3$. Optical depth,
radiative diffusion, weakly coupled upper layers, and a massive deep convective
reservoir can broaden the timescale spectrum further. Exact PAVA adjustment removes
the artificial finite-$K$ convective timescale, but it does not remove this
radiative separation.

Such a problem may have an explicit timestep fixed by fast upper or photospheric
modes while the scientifically relevant deep evolution occurs over years or
longer. An optimized BDF2 or TR-BDF2 constrained method could then be strongly
favored. Pressure range alone does not prove stiffness, because very optically thin
layers may also be weakly coupled. The decision should be based on the same
fixed-end-time accuracy-versus-work benchmark used here, supplemented by the
spectrum or conditioning of the projected radiative Jacobian.

### Always converge the full physical residual

At every accepted nonlinear state:

1. reconstruct density and any temperature-dependent geometry;
2. update opacities and other RT properties;
3. compute complete shortwave and longwave interface fluxes;
4. construct either the smooth convective flux or exact adjustment map;
5. form conservative layer and surface energy imbalances; and
6. test convergence using scaled flux residuals as well as temperature/map
   residuals.

Approximate Jacobians may change convergence, but they do not bias the final state
if every accepted residual uses the full physics and is converged tightly.

### Jacobian strategy for smooth convection

A practical approximation is

$$
J\approx J_{\mathrm{rad}}^{\mathrm{frozen\ opacity}}
+J_{\mathrm{conv}}^{\mathrm{analytic}}.
$$

The radiative block can use one-sided temperature differences through the fast
`radiate` path while holding opacity and shortwave properties fixed. The complete
opacity dependence remains in accepted-state residuals. The analytic convective
block is sparse because an interface flux depends only on its two adjacent
temperatures; it captures the stiffest modes introduced by large $K$.

This work is most valuable for PTC, where a new Jacobian would otherwise be formed
at every pseudo-step. The existing HYBRJ path reuses Jacobians and is mostly
residual-dominated, so an analytic Jacobian is less likely to transform its runtime.

### Jacobian strategy for projected PTC

The adjustment itself is $O(N)$ and cheap. The expensive part is differentiating
the projected RT map. A first implementation can:

- finite-difference temperature through the PAVA map;
- use frozen opacity for the perturbed `radiate` calls;
- evaluate unperturbed and accepted candidate residuals with complete opacity;
- lag the projected radiative Jacobian while residual reduction remains adequate;
- rebuild after active-set changes accompanied by poor progress, repeated rejected
  trials, or large temperature changes; and
- reuse the same factorization while only the PTC pseudo-timestep changes.

Because PAVA is piecewise linear, its derivative is simple for a fixed block
partition: perturbations within a mixed block are replaced by their weighted block
average. A later implementation could combine that analytic projection derivative
with the approximate radiative Jacobian. Active-set boundaries still require a
semismooth globalization strategy.

The toy currently rebuilds every projected Jacobian. Its 1,335-RT count is therefore
a conservative baseline; production-style lagging could reduce it substantially.

### Scaling and stopping criteria

Temperature tendencies alone are poorly scaled because atmospheric cell heat
capacities vary with pressure and the ocean heat capacity is enormous. Production
convergence tests should include:

- maximum and normed flux imbalance;
- TOA and surface energy balance;
- the projected-map residual for exact adjustment;
- lapse-rate feasibility; and
- scaled temperature corrections.

A trust-region or line-search merit function should account for all of these rather
than relying on an unweighted tendency norm.

For projected PTC, production diagnostics should additionally record the selected
$\alpha$, the scaled radiative trial displacement, active-set changes, and the
natural residual evaluated at nearby scalar values. Agreement at $\alpha/10$,
$\alpha$, and $10\alpha$ is a direct check that numerical projection details have
not introduced an unwanted scaling dependence.

## Coupling steady RCE to chemistry and clouds

The projected residual need not contain every physical variable as an explicit
Newton unknown. The appropriate coupling depends on whether another process is a
cheap diagnostic closure, an expensive steady subproblem, or a genuinely
time-dependent state equation. In all cases, approximate derivatives may lag
expensive physics, but final convergence must be tested with all coupled physics
updated consistently.

### Equilibrium chemistry: eliminate it inside the residual

If chemical equilibrium is unique and inexpensive, define composition as a
diagnostic function of temperature, pressure, and elemental abundances,

$$
\mathbf q_{\mathrm{eq}}=Q(\mathbf T,p,\mathbf b),
$$

and construct the reduced radiative tendency

$$
\mathbf f_{\mathrm{rad}}(\mathbf T)
=R[\mathbf T,Q(\mathbf T),\kappa(\mathbf T,Q(\mathbf T))].
$$

Chemical equilibrium should then be solved during every complete base-state and
candidate residual evaluation. Temperature remains the exposed projected-PTC
unknown, while composition is eliminated by the inner equilibrium solve. This
lets candidate acceptance see feedbacks such as
$T\rightarrow q\rightarrow\kappa\rightarrow R$ without greatly enlarging the
nonlinear system.

The approximate Jacobian need not include all of that response initially. A
practical inexact-Newton implementation can freeze equilibrium composition and
opacity in the perturbed `radiate` calls, while recomputing both in every complete
accepted-state and candidate residual. Rebuild the approximation after large
temperature changes, relevant chemical transitions, poor residual reduction, or
active-set changes accompanied by stagnation. Complete equilibrium chemistry and
opacity must always be restored for the final convergence test.

A completely separate outer equilibrium-chemistry iteration is useful for initial
implementation and debugging, but it is a weaker final architecture. Solving RCE
to tight tolerance with obsolete composition can overshoot strongly coupled
opacity transitions and lead to oscillation or heavy under-relaxation.

### Steady photochemistry: an inexact outer block solve

A full steady photochemical model is different. Solving reaction--transport
kinetics during every projected residual or every Jacobian perturbation would
usually dominate the calculation. The natural first production architecture is a
block nonlinear iteration for

$$
G_T(\mathbf T,\mathbf q)=0,
\qquad
G_q(\mathbf T,\mathbf q)=0,
$$

where $G_T$ is the projected RCE residual and $G_q$ is the steady
reaction--transport residual. One outer cycle is

1. hold composition fixed and take several projected-PTC climate steps;
2. hold the updated climate and atmospheric structure fixed and solve
   photochemistry to steady state, warm-started from the previous composition;
3. update mean molecular weight, heat capacity, hydrostatic structure, opacity,
   photolysis rates, and other shared quantities; and
4. repeat until both blocks and their shared state stop changing.

Neither block should be tightly converged early in the outer iteration. It is
wasteful to solve RCE accurately with composition that will soon change, and
similarly wasteful to solve photochemistry accurately at a climate state far from
the coupled solution. Early cycles can use a fixed small number of climate PTC
steps and looser photochemical tolerances; both tolerances should tighten as the
coupled residual falls.

Photochemistry should be updated adaptively rather than necessarily after every
climate step. Useful triggers include:

- a sufficiently large accumulated temperature change;
- a material change in convective topology or tropopause location;
- a large change in UV or thermal optical depth;
- a prescribed reduction in the climate residual;
- stagnation of projected PTC with frozen composition; or
- a maximum number of accepted climate steps since the previous update.

Warm starts are essential. Previous abundances, timestep history, Jacobian
sparsity, symbolic factorizations, and unchanged photolysis intermediates should be
reused where possible. Strong climate--chemistry feedback may require safeguarded
under-relaxation, preferably in logarithmic abundance variables with explicit
positivity and elemental-conservation checks. Anderson acceleration is a natural
next step for the expensive outer fixed point; limited-memory Broyden or a block
Newton method should be considered only if simpler acceleration is inadequate.

A converged outer fixed point is a valid coupled steady solution if the subproblems
are solved consistently and are locally unique. The final verification should use
a handshake cycle: converge photochemistry at the latest climate, reevaluate and
converge climate with that composition, rerun photochemistry, and verify that the
climate residual, chemical residual, temperature, abundances, and opacities all
remain within tolerance.

An advanced future alternative is a reduced coupled Jacobian. If
$F_q(T,q)=0$, then

$$
\frac{d\mathbf q}{d\mathbf T}
=-\left(\frac{\partial F_q}{\partial\mathbf q}\right)^{-1}
\frac{\partial F_q}{\partial\mathbf T}.
$$

A photochemical solver already constructs the stiff chemical Jacobian
$\partial F_q/\partial q$. Reusing its factorization for sensitivity solves could
capture chemical feedback in the climate Jacobian without rerunning a complete
steady photochemical integration for every temperature perturbation. This
Schur-complement approach is promising but substantially more complex than the
inexact outer block iteration.

### Clouds: choose the coupling from the cloud model

Clouds range from inexpensive diagnostic closures to prognostic microphysical
systems, so there is no single correct update location.

For a deterministic diagnostic cloud model,

$$
\mathbf c=C(\mathbf T,p,\mathbf q,K_{zz},\ldots),
$$

the eventual reduced steady residual can evaluate cloud properties during every
complete candidate residual, just like equilibrium chemistry. Cloud opacity may
still be frozen in approximate Jacobian calls. This is attractive only when the
diagnostic map is unique, inexpensive, and sufficiently continuous.

Many 1-D cloud models instead contain an internal steady closure involving cloud
base, particle size, sedimentation, mixing, supersaturation, rainout, and elemental
depletion. They may also depend on a convective flux diagnosed from the climate
solution. The safer initial coupling is then a damped outer cloud iteration:

1. hold cloud properties fixed and take several projected-PTC steps;
2. update the cloud solution;
3. limit abrupt changes in cloud optical depth or use continuation from clear to
   cloudy opacity; and
4. repeat until climate, cloud structure, and opacity are all stationary.

Cloud bases appearing, disappearing, or crossing grid cells introduce
nonsmoothness and possibly hysteresis beyond the clean convex nonsmoothness of dry
PAVA. Useful safeguards include cloud-opacity continuation, limits on accepted
optical-depth changes, lagged cloud updates, trust-region globalization, and
Jacobian rebuilds after cloud-topology changes. Anderson acceleration may help a
well-defined outer cloud fixed point, but it does not cure a nonunique or
hysteretic closure.

If clouds affect only opacity, they can remain part of the radiative closure. If
condensation changes latent heating, molecular weight, condensable inventory, or
the critical lapse rate, that thermodynamic effect belongs inside the moist
convective projection. An outer opacity update cannot substitute for a
conservative $P_{\mathrm{moist}}$.

For prognostic nucleation, growth, evaporation, sedimentation, or aerosol
transport, cloud variables are genuine unknowns. Their steady equations should
eventually be included in a block-coupled nonlinear solve, although a split outer
iteration remains a useful first implementation and preconditioner.

### What carries over to time evolution

The block-organization ideas carry over, but projected PTC itself is a
steady-state solver and its pseudo-trajectory is not time accurate. A physical
time integration must distinguish diagnostic assumptions from prognostic
processes:

- equilibrium chemistry may be reevaluated each physical timestep if
  instantaneous equilibrium is the intended model;
- finite-rate photochemistry must be advanced over the same physical interval,
  using operator splitting, subcycling, or a coupled implicit integrator rather
  than a steady photochemical solve at every step;
- diagnostic clouds may be updated each step, while prognostic clouds must evolve
  through their physical microphysical and transport equations; and
- instantaneous convection can use explicit adjustment or constrained implicit
  integration, but projected PTC should be reserved for finding the steady state.

For time-accurate splitting, the ordering and update frequency create temporal
error and must be checked by timestep refinement; Strang splitting may be useful
when individual subsolvers permit it. For a steady solve, update scheduling affects
efficiency and globalization but not the desired fixed point, provided every block
is ultimately converged against the latest values of all the others.

## Extension to moist atmospheres

Nothing in the toy yet establishes a self-consistent moist lapse rate. The dry
PAVA algorithm is unusually simple because $T=\Pi\theta$, $c_p$ is fixed, and dry
thermal energy is linear in the adjusted variable.

A moist adjustment operator would need to define, at minimum:

- the moist stability criterion;
- which quantities are mixed or conserved, such as moist enthalpy and total water;
- phase equilibrium and condensate retention or rainout;
- latent heating during block adjustment; and
- consistent surface/ocean volatile exchange.

The block problem would generally be nonlinear rather than a weighted average of
potential temperature. Nevertheless, the architecture can survive: if a robust,
energy-conserving moist adjustment operator $P_{\mathrm{moist}}$ can be defined,
then explicit adjustment and the natural residual

$$
G(\mathbf x)
=\frac{P_{\mathrm{moist}}[\mathbf x+\alpha f(\mathbf x)]-\mathbf x}{\alpha}
$$

remain available, as does constrained backward Euler with the same moist
projection inside each implicit timestep. The state $\mathbf x$ would likely
include composition or total water as well as temperature. Whether either implicit
method remains efficient will depend on the uniqueness, smoothness, and cost of
that moist block solve.

The smooth-flux route can also use a state-dependent critical moist gradient, but
latent heating cannot be represented consistently by changing the lapse-rate
threshold alone. Condensation, vapor transport, condensate treatment, and energy
conservation must be coupled to the temperature equations.

## Recommended validation sequence

1. Reproduce dry radiative equilibrium with production RT and verify column,
   surface, and TOA energy conservation.
2. Implement dry PAVA adjustment and test random-profile stability and block-energy
   conservation independently of RT.
3. Run conservative explicit adjustment over a timestep sweep; use it as the
   reference dry RCE.
4. Implement projected PTC with a complete finite-difference Jacobian and verify
   agreement with explicit adjustment.
5. Compare projected PTC directly with the existing active-set Newton RCE solver
   over difficult initial profiles and changing convective-zone structures.
6. Substitute frozen-opacity `radiate` derivatives and confirm that converged full
   residuals and profiles are unchanged.
7. Add Jacobian lagging, refresh logic, and active-set-aware globalization; compare
   work in full opacity updates, fast `radiate` calls, and factorizations.
8. For every intended time-dependent application, compare explicit adjustment at
   a fixed end time with a timestep-refined reference. Implement optimized
   constrained implicit integration only if stability is materially more
   restrictive than accuracy.
9. Repeat over vertical resolution, initial profiles, optical depths, stellar flux,
   surface heat capacity, and deep gas-giant pressure ranges. Retain the smooth
   hybrid path only as a finite-rate comparison.
10. Design and unit-test a conservative moist block adjustment before coupling it
    to explicit, constrained-implicit, or projected-PTC evolution.
11. Add equilibrium chemistry inside complete projected residual evaluations while
    initially freezing its composition and opacity response in approximate
    Jacobians.
12. Couple expensive steady photochemistry through an inexact, warm-started outer
    block iteration and test safeguarded relaxation and Anderson acceleration.
13. Introduce clouds first through a safeguarded outer continuation; move a cloud
    closure inside the complete residual only after its uniqueness, continuity,
    conservation, and cost are understood.

## Present conclusion

The main modeling decision is whether convection is intended to have a resolved
finite timescale. If it is, the smooth conservative flux is a viable differentiable
closure and hybrid $K$ continuation is its most useful controller. If convection
is instead assumed to adjust instantaneously, finite-$K$ evolution is an indirect
penalty approximation: it requires a strength parameter, leaves residual
superadiabaticity, and becomes stiff as exact neutrality is approached. Weighted
PAVA is then the cleaner dry-convection operator.

For steady dry RCE, projected PTC is the leading method from this study. It solves
the PAVA fixed point directly and avoids waiting for slow physical modes. In the
default toy it used roughly 18 times fewer full RT evaluations than conservative
explicit relaxation even though every projected Jacobian was rebuilt with forward
differences. A production implementation should combine it with frozen-opacity
radiative derivatives, analytical PAVA block derivatives, Jacobian reuse, and
active-set-aware globalization. Automatic scalar-$\alpha$ selection is a useful
default, but the high-$c_p$ experiment shows that restart logic or better
globalization is still required. Its superiority is well supported by the toy but
must still be tested against the existing active-set Newton solver under realistic
opacity, composition, and convective-zone changes. Metric preconditioning should
remain a second-stage optimization unless real RT demonstrates that scalar scaling
is a persistent limitation.

For coupled steady problems, the same reduced-space principle should be used where
it is economical. Fast equilibrium chemistry and simple diagnostic clouds can be
eliminated inside complete climate residual evaluations while remaining frozen in
approximate Jacobians. Expensive steady photochemistry and complex cloud closures
should begin as inexact, warm-started outer blocks, with all component residuals
and shared opacities checked together at final convergence. This retains the
steady-state efficiency of projected PTC without multiplying every radiative
Jacobian call by another expensive nonlinear solve.

For time-accurate evolution, explicit radiation followed by PAVA is the current
recommendation. In the fixed-time benchmark it reached essentially the same
accuracy as constrained backward Euler in 100 versus 131 accepted macro-steps, and
its one-RT-call steps made it vastly cheaper. This demonstrates that the gray toy
is not strongly stability-limited over the tested trajectory.

The constrained implicit method should not be discarded. Its present cost is
dominated by an intentionally conservative first implementation: first-order
backward Euler, three nonlinear solves for step doubling, complete finite-difference
Jacobians, and no reuse. BDF2 or TR-BDF2, an analytical generalized derivative of
PAVA, lagged radiative Jacobians, and factorization reuse could reduce its cost by
orders of magnitude. It becomes the more plausible time-accurate choice when a
problem—potentially a deep gas giant—has fast radiative modes that constrain
explicit stability far more strongly than the desired temporal accuracy. The jury
therefore remains open for production atmospheres, and the fixed-time benchmark is
the appropriate decision test.

The practical hierarchy is:

1. use projected PTC as the leading dry steady-state candidate;
2. use explicit adjustment as the time-accurate baseline and correctness
   reference;
3. optimize constrained implicit integration only where measured stiffness or
   robustness requirements justify it; and
4. retain smooth finite-$K$ convection when finite-rate convection is the intended
   closure rather than an approximation to instantaneous adjustment.

Moist convection is a separate unresolved problem. These recommendations can be
extended only after a unique, conservative moist adjustment map has specified
water conservation, latent heating, phase equilibrium, condensate retention or
rainout, and surface volatile exchange.
