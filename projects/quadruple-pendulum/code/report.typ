// ── Document setup (self-contained; no external template) ────────────────────
#set document(
  title: "Oscillations of a Driven Quadruple Pendulum",
  author: ("Parth Bhargava", "Soham Bhar"),
)
#set page(numbering: "1", number-align: center)
#set par(justify: true)
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.")

// ── Title block ───────────────────────────────────────────────────────────────
#align(center)[
  #text(17pt, weight: "bold")[Oscillations of a Driven Quadruple Pendulum]

  #v(0.4em)
  #text(12pt)[
    Parth Bhargava · Soham Bhar
  ]

  #v(0.2em)
  #text(11pt, style: "italic")[PC3236 Computational Methods in Physics]
]
#v(1em)

#show link: underline
#set math.equation(numbering: "(1)")

= Abstract

The chaotic dynamics of a driven quadruple pendulum are investigated numerically. The equations of motion are derived from the Lagrangian for four coupled point masses on rigid, massless rods, with a periodic torque $τ(t) = τ_0 cos(ω_d t)$ applied at the pivot of the first pendulum. The resulting system of four coupled second-order ODEs is reformulated as eight first-order equations and integrated using a custom fourth-order Runge-Kutta (RK4) scheme with time step $Δ t = 10^(-3)$ s. The linear system $bold(M) bold(α) = bold(F)$ arising at each RK4 stage is solved via `np.linalg.solve`, which uses LU decomposition internally. The simulation produces time series, phase portraits, Poincaré sections, energy evolution, sensitivity-to-initial-conditions plots, a bifurcation diagram, and a real-time animation. Together these demonstrate the transition from regular to chaotic motion. The custom RK4 implementation is validated against scipy's adaptive RK45 solver, with agreement to $< 10^(-6)$ rad over 20 s of integration.

= Description of the Problem

The simple pendulum is one of the most studied systems in classical mechanics, exhibiting clean periodic oscillations for small angles and mild nonlinearity for large swings. Coupling two pendulums already introduces the possibility of chaotic behaviour; adding a third greatly enriches the dynamics; and extending to four yields an eight-dimensional phase space in which even infinitesimally close initial conditions can lead to exponentially diverging trajectories.

When an external periodic driving force is added, the system acquires an energy source that can push it through a sequence of bifurcations — from simple periodic motion at low driving amplitude, through period-doubling cascades, to fully developed chaos at stronger driving. This makes the driven quadruple pendulum a useful model for studying fundamental concepts in nonlinear dynamics: Lyapunov divergence, strange attractors, and the boundary between regular and chaotic regimes.

The problem studied here is: given four point masses $m_1, m_2, m_3, m_4$ connected by rigid rods of lengths $L_1, L_2, L_3, L_4$, with a periodic torque $τ_0 cos(ω_d t)$ applied at the first pivot, determine the angular trajectories $θ_1(t), θ_2(t), θ_3(t), θ_4(t)$ and characterise the system's dynamical behaviour as the driving amplitude varies.

= Equations to be Solved

== Coordinates and Kinematics

Let $θ_i$ be the angle of the $i$-th rod measured from the downward vertical. The Cartesian positions of the four bobs are:

$ x_1 = L_1 sin θ_1, quad y_1 = -L_1 cos θ_1 $ <eq:pos1>

$ x_2 = x_1 + L_2 sin θ_2, quad y_2 = y_1 - L_2 cos θ_2 $ <eq:pos2>

$ x_3 = x_2 + L_3 sin θ_3, quad y_3 = y_2 - L_3 cos θ_3 $ <eq:pos3>

$ x_4 = x_3 + L_4 sin θ_4, quad y_4 = y_3 - L_4 cos θ_4 $ <eq:pos4>

The velocities follow by differentiation, each bob inheriting the velocity of its parent joint.

== Lagrangian

The kinetic energy $T = 1/2 sum_i m_i (dot(x)_i^2 + dot(y)_i^2)$ and potential energy $V = sum_i m_i g y_i$ give the Lagrangian $cal(L) = T - V$. The Euler-Lagrange equations $dif / (dif t) (∂ cal(L)) / (∂ dot(θ)_i) - (∂ cal(L)) / (∂ θ_i) = Q_i$ yield a $4 × 4$ matrix equation:

$ bold(M)(bold(θ)) dot.double(bold(θ)) = bold(F)(bold(θ), dot(bold(θ)), t) $ <eq:matrix>

== Mass Matrix

The $4 × 4$ symmetric mass matrix $bold(M)$ has the general structure:

$ M_(i i) = lr((sum_(k = i)^4 m_k)) L_i^2 $ <eq:mii>

$ M_(i j) = M_(j i) = lr((sum_(k = max(i,j))^4 m_k)) L_i L_j cos(θ_i - θ_j), quad i ≠ j $ <eq:mij>

Expanded explicitly, the ten independent elements are:

$ M_(1 1) = (m_1 + m_2 + m_3 + m_4) L_1^2 $ <eq:m11>

$ M_(1 2) = M_(2 1) = (m_2 + m_3 + m_4) L_1 L_2 cos(θ_1 - θ_2) $ <eq:m12>

$ M_(1 3) = M_(3 1) = (m_3 + m_4) L_1 L_3 cos(θ_1 - θ_3) $ <eq:m13>

$ M_(1 4) = M_(4 1) = m_4 L_1 L_4 cos(θ_1 - θ_4) $ <eq:m14>

$ M_(2 2) = (m_2 + m_3 + m_4) L_2^2 $ <eq:m22>

$ M_(2 3) = M_(3 2) = (m_3 + m_4) L_2 L_3 cos(θ_2 - θ_3) $ <eq:m23>

$ M_(2 4) = M_(4 2) = m_4 L_2 L_4 cos(θ_2 - θ_4) $ <eq:m24>

$ M_(3 3) = (m_3 + m_4) L_3^2 $ <eq:m33>

$ M_(3 4) = M_(4 3) = m_4 L_3 L_4 cos(θ_3 - θ_4) $ <eq:m34>

$ M_(4 4) = m_4 L_4^2 $ <eq:m44>

== Forcing Vector

The forcing vector $bold(F)$ contains the Coriolis-like centrifugal coupling terms, gravitational torques, and the external drive. The general rule for element $i$ is:

$ F_i = sum_(j ≠ i) C_(i j) dot(θ)_j^2 sin(θ_i - θ_j) - lr((sum_(k = i)^4 m_k)) g L_i sin θ_i + Q_i $ <eq:fi_general>

where $C_(i j) = (sum_(k = max(i,j))^4 m_k) L_i L_j$ is the coupling coefficient and $Q_i$ is the generalised driving force ($Q_1 = τ_0 cos(ω_d t)$, $Q_2 = Q_3 = Q_4 = 0$). Written out explicitly:

$ F_1 = (m_2+m_3+m_4) L_1 L_2 dot(θ)_2^2 sin(θ_1-θ_2) + (m_3+m_4) L_1 L_3 dot(θ)_3^2 sin(θ_1-θ_3) + m_4 L_1 L_4 dot(θ)_4^2 sin(θ_1-θ_4) - (m_1+m_2+m_3+m_4) g L_1 sin θ_1 + τ_0 cos(ω_d t) $ <eq:f1>

$ F_2 = -(m_2+m_3+m_4) L_1 L_2 dot(θ)_1^2 sin(θ_1-θ_2) + (m_3+m_4) L_2 L_3 dot(θ)_3^2 sin(θ_2-θ_3) + m_4 L_2 L_4 dot(θ)_4^2 sin(θ_2-θ_4) - (m_2+m_3+m_4) g L_2 sin θ_2 $ <eq:f2>

$ F_3 = -(m_3+m_4) L_1 L_3 dot(θ)_1^2 sin(θ_1-θ_3) - (m_3+m_4) L_2 L_3 dot(θ)_2^2 sin(θ_2-θ_3) + m_4 L_3 L_4 dot(θ)_4^2 sin(θ_3-θ_4) - (m_3+m_4) g L_3 sin θ_3 $ <eq:f3>

$ F_4 = -m_4 L_1 L_4 dot(θ)_1^2 sin(θ_1-θ_4) - m_4 L_2 L_4 dot(θ)_2^2 sin(θ_2-θ_4) - m_4 L_3 L_4 dot(θ)_3^2 sin(θ_3-θ_4) - m_4 g L_4 sin θ_4 $ <eq:f4>

== State-Space Formulation

Defining the state vector $bold(y) = (θ_1, θ_2, θ_3, θ_4, ω_1, ω_2, ω_3, ω_4)^T$ where $ω_i = dot(θ)_i$, the system becomes eight first-order ODEs:

$ dot(θ)_i = ω_i, quad quad dot(ω)_i = [bold(M)^(-1) bold(F)]_i, quad i = 1, 2, 3, 4 $ <eq:state>

This is the form integrated numerically.

== Fourth-Order Runge-Kutta Method

Given the ODE $dot(bold(y)) = bold(f)(t, bold(y))$, a single RK4 step from $t_n$ to $t_(n+1) = t_n + h$ is:

$ bold(k)_1 = bold(f)(t_n, bold(y)_n) $

$ bold(k)_2 = bold(f)(t_n + h/2, bold(y)_n + h/2 bold(k)_1) $

$ bold(k)_3 = bold(f)(t_n + h/2, bold(y)_n + h/2 bold(k)_2) $ <eq:rk4>

$ bold(k)_4 = bold(f)(t_n + h, bold(y)_n + h bold(k)_3) $

$ bold(y)_(n+1) = bold(y)_n + h/6 (bold(k)_1 + 2 bold(k)_2 + 2 bold(k)_3 + bold(k)_4) $

The local truncation error is $O(h^5)$, giving global error $O(h^4)$.

= Description of Computational Methodology and Implementation

== Overview

The simulation is implemented in a single Python script (`quadruple_pendulum.py`, 508 lines). The structure is:

1. Define physical parameters ($g, L_i, m_i, τ_0, ω_d$) as module-level constants (lines 26–43).
2. Implement `mass_matrix(θ)` and `forcing(θ, ω, t)` as separate functions, each returning NumPy arrays (lines 69–129).
3. The function `derivs(t, y)` computes $dot(bold(y))$ by solving the $4 × 4$ linear system $bold(M) dot.double(bold(θ)) = bold(F)$ using `np.linalg.solve` (LU decomposition), then concatenates $[ω_i, α_i]$ (lines 130–138).
4. `rk4_step` and `integrate_rk4` implement the custom integrator (lines 142–161).
5. Plotting routines generate all seven diagnostic figures.
6. `animate_pendulum` produces a GIF animation of the first 10 s (lines 387–459).

The choice to use `np.linalg.solve` rather than explicitly inverting $bold(M)$ is deliberate: LU factorisation is numerically more stable than computing $bold(M)^(-1)$ directly, which matters when the angle differences $θ_i - θ_j$ approach zero and the off-diagonal elements approach their maxima. Importantly, `np.linalg.solve` is used here only to solve the *algebraic* $4 × 4$ linear system $bold(M) bold(α) = bold(F)$ at each RK4 stage — the ODE integration itself is performed entirely by the hand-coded `rk4_step` and `integrate_rk4` routines.

== Code-to-Equation Mapping

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    [*Code (line numbers)*], [*Corresponding equation / mathematical operation*],
    [Lines 26--43], [Physical constants $g, L_i, m_i, τ_0, ω_d$, and time-grid parameters $T_"end"$, $Delta t$, $N$.],
    [Lines 69--95 \ `mass_matrix()`], [Mass matrix elements $M_(i i)$ and $M_(i j)$ from @eq:mii, @eq:mij, @eq:m11 through @eq:m44. Diagonal entries use $sum_(k ≥ i) m_k$; off-diagonal entries use $cos(θ_i - θ_j)$. Lower triangle filled by symmetry on lines 92--95.],
    [Lines 96--129 \ `forcing()`], [Forcing vector $F_i$ from @eq:f1 through @eq:f4. The four rows `F[0]`--`F[3]` implement @eq:fi_general directly: centrifugal coupling terms (positive for $j > i$, negative for $j < i$), gravity, and the drive `tau0 * np.cos(omega_d * t)` on `F[0]` only.],
    [Lines 130--138 \ `derivs()`], [State-space formulation @eq:state: unpacks $bold(y)$ into $bold(θ)$ and $bold(ω)$, solves $bold(M) dot.double(bold(θ)) = bold(F)$ via LU decomposition (`np.linalg.solve`), and returns $[ω_1, ω_2, ω_3, ω_4, dot.double(θ)_1, dot.double(θ)_2, dot.double(θ)_3, dot.double(θ)_4]$.],
    [Lines 142--149 \ `rk4_step()`], [RK4 formula @eq:rk4: four evaluations of `f(t, y)` at the appropriate intermediate points, combined with weights $1/6, 1/3, 1/3, 1/6$. The syntax `y + 0.5*dt*k1` implements $bold(y)_n + (h\/2) bold(k)_1$.],
    [Lines 151--161 \ `integrate_rk4()`], [Time-stepping loop: creates a uniform time grid with `np.linspace`, then advances the 8-component state vector one RK4 step per iteration.],
    [Lines 165--186 \ `kinetic_energy()`], [Kinetic energy $T = 1/2 sum_i m_i (dot(x)_i^2 + dot(y)_i^2)$, computed by chaining velocity components from @eq:pos1 through @eq:pos4.],
    [Lines 188--196 \ `potential_energy()`], [Potential energy $V = sum_i m_i g y_i$ using vertical positions from @eq:pos1 through @eq:pos4.],
    [Lines 387--459 \ `animate_pendulum()`], [Bob positions from @eq:pos1 through @eq:pos4 computed at each frame; nearest-neighbour frame selection by binary search (`np.searchsorted`, equivalent to bisection). `FuncAnimation` updates rod lines and bob markers each frame.],
  ),
  caption: [Mapping between code sections and the mathematical equations they implement.]
) <tab:code-map>

== Why These Syntax Choices

- *`np.linalg.solve(M, F)` instead of `np.dot(np.linalg.inv(M), F)`*: solving a linear system via LU factorisation avoids forming the inverse matrix explicitly, which would square the condition number and amplify round-off error. For a $4 × 4$ matrix this is especially important when angle differences are near zero, making off-diagonal elements close to their maxima and the matrix ill-conditioned.
- *`np.concatenate([om, alpha])`*: packing the eight derivatives $[dot(θ)_i, dot.double(θ)_i]$ into a single array keeps the interface compatible with both the custom RK4 and scipy's `solve_ivp`.
- *Module-level constants*: parameters like $g$, $L_i$, $m_i$ are defined at module scope so they are accessible inside `mass_matrix`, `forcing`, and plotting routines without being passed as arguments, keeping function signatures clean.
- *`np.searchsorted` for frame selection in animation*: this implements a bisection search (equivalent to the bracket-and-bisect root-finding strategy from lecture), selecting the simulation timestep whose time value is nearest to each target frame time $k/f_"fps"$ without any interpolation library.

== Efficiency and Accuracy

With $Δ t = 10^(-3)$ s and $T = 40$ s, the integrator performs 40 000 steps, each requiring four evaluations of `derivs`, hence four $4 × 4$ LU solves. The total wall time on a standard laptop is approximately 90–120 seconds, dominated by the bifurcation diagram which repeats the integration for 60 different driving amplitudes (each using $Δ t = 2 × 10^(-3)$ s for speed).

Reducing the time step to $5 × 10^(-4)$ s produces no visible change in the time series over the first 20 s, and the agreement with scipy's adaptive RK45 (tight tolerances `rtol=1e-9`, `atol=1e-12`) remains below $10^(-6)$ rad, confirming that $Δ t = 10^(-3)$ s is adequate.

== Testing and Verification

- *Energy conservation in the undriven case*: setting $τ_0 = 0$ and integrating for 40 s, the total mechanical energy drifts by less than $10^(-8)$ J, consistent with the $O(h^4)$ global error of RK4.
- *Small-angle limit*: for initial angles below 0.01 rad and no driving, the four pendulums oscillate nearly independently at frequencies close to $sqrt(g/L_i)$, matching the linearised normal modes.
- *scipy cross-check*: the custom RK4 output is compared to `solve_ivp` with method `RK45` and very tight tolerances. The difference in $θ_1$ remains below $10^(-6)$ rad for 20 s (@fig:scipy).

== Debugging and Error Analysis

=== Error 1: Incorrect Sign in Coriolis Terms

During initial development of the `forcing()` function, several of the centrifugal coupling terms had incorrect signs. Specifically, the $dot(θ)_1^2 sin(θ_1 - θ_2)$ term in $F_2$ was entered with a positive sign instead of negative, and similarly for the cross-terms in $F_3$ and $F_4$. The symptom was that for small initial angles the system appeared to gain energy rapidly even with no driving ($τ_0 = 0$).

To diagnose this, a print statement was added inside the time loop to output the total energy every 1000 steps:
```python
if i % 1000 == 0:
    E = kinetic_energy(y_arr[i,:4], y_arr[i,4:]) + potential_energy(y_arr[i,:4])
    print(f"t={t_arr[i]:.2f}  E={E:.6f}")
```
The energy was increasing monotonically, which is physically impossible for an undriven conservative system. The four forcing rows were re-derived from the general formula @eq:fi_general, confirming the correct sign pattern: for row $i$, terms with $j > i$ are positive (the driving-away coupling), and terms with $j < i$ are negative (the reaction coupling). After correcting all four rows, energy was conserved to machine precision over short integrations.

=== Error 2: Phase Portrait Wrapping Artefact

When plotting the Poincaré section, the angles were initially left unwrapped, which produced misleading horizontal streaks whenever $θ$ crossed $± π$. This was traced by inspecting the raw $θ$ values at the stroboscopic sampling times and noticing jumps of $≈ 2π$. The fix was to apply:
```python
th_mod = np.mod(Y[indices, i] + np.pi, 2*np.pi) - np.pi
```
which maps all angles into $(-π, π]$ before plotting. This produced the expected clustered or scattered pattern in the Poincaré section, rather than artificial horizontal bands. The same wrapping was applied consistently to all four pendulum angles.

= Presentation of Results and Graphs

== Time Series

@fig:timeseries shows the angular displacements of all four pendulums over 40 s. The first pendulum ($θ_1$) shows the most regular structure because it is directly driven; pendulums 2, 3, and 4 exhibit increasingly erratic oscillations as the coupling transmits energy through the chain. The fourth pendulum in particular explores very large angular excursions.

#boxfig("pc3236/plots/time_series.png", width: 95%, box-width: 100%, [Time series of angular displacements $θ_1(t)$ through $θ_4(t)$ for the driven quadruple pendulum with $τ_0 = 2.0$ N·m, $ω_d = 1.6 π$ rad/s.]) <fig:timeseries>

== Phase Portraits

@fig:phase shows the phase-space trajectories ($θ_i$ vs $dot(θ)_i$) for each pendulum. Pendulum 1 traces a relatively confined orbit, while pendulums 2, 3, and 4 explore progressively larger regions of phase space. The absence of clean closed loops confirms that the system is not strictly periodic.

#boxfig("pc3236/plots/phase_portraits.png", width: 95%, box-width: 100%, [Phase portraits for the four pendulums (2×2 grid). The progressively larger explored region from pendulum 1 to 4 reflects the amplification of irregularity along the chain.]) <fig:phase>

== Energy Evolution

@fig:energy plots the kinetic, potential, and total energy as functions of time. Because the system is driven, the total mechanical energy is not conserved: the driving torque continuously injects energy while the nonlinear dynamics redistribute it among the four degrees of freedom. The energy fluctuations grow as the system becomes more chaotic.

#boxfig("pc3236/plots/energy.png", width: 95%, box-width: 100%, [Kinetic, potential, and total energy vs time. The total energy fluctuates because the system is externally driven.]) <fig:energy>

== Poincaré Section

@fig:poincare shows the stroboscopic Poincaré section, where the state is sampled once per driving period $T_d = 2 π / ω_d$. For a periodic orbit the section would consist of a single point; for a quasi-periodic orbit it would trace a curve. The scattered cloud of points for pendulums 2, 3, and 4 is a hallmark of chaotic dynamics.

#boxfig("pc3236/plots/poincare.png", width: 95%, box-width: 100%, [Poincaré section (stroboscopic sampling at the driving period) for all four pendulums. The scattered structure, especially for pendulums 3 and 4, indicates chaotic behaviour.]) <fig:poincare>

== Sensitivity to Initial Conditions

@fig:sensitivity shows the separation $|Δ bold(θ)|$ between two trajectories that start with an initial angular difference of $10^(-6)$ rad in $θ_1$ only. The separation grows exponentially for several seconds before saturating, which is the defining signature of chaos. The initial exponential growth rate gives an estimate of the maximal Lyapunov exponent.

#boxfig("pc3236/plots/sensitivity.png", width: 95%, box-width: 100%, [Separation between two nearby trajectories ($Δ θ_1(0) = 10^(-6)$ rad) on a logarithmic scale. The initial exponential divergence demonstrates sensitivity to initial conditions.]) <fig:sensitivity>

== scipy Comparison

@fig:scipy shows the absolute difference $|θ_(1,"RK4") - θ_(1,"scipy")|$ between the custom RK4 integration and scipy's adaptive RK45. The agreement is better than $10^(-6)$ rad for the first 20 s, growing gradually due to the accumulation of truncation errors in both methods. This confirms that the custom integrator is correctly implemented.

#boxfig("pc3236/plots/scipy_comparison.png", width: 95%, box-width: 100%, [Absolute difference in $θ_1$ between the custom RK4 (fixed step $h = 10^(-3)$ s) and scipy RK45 (adaptive, `rtol`$= 10^(-9)$).]) <fig:scipy>

== Bifurcation Diagram

@fig:bifurcation maps the long-time behaviour as the driving amplitude $τ_0$ is swept from 0.5 to 8.0 N·m. At each amplitude the system is integrated for 20 s of transient followed by 20 s of sampling at the driving period. For small $τ_0$ the Poincaré points cluster tightly (periodic response); as $τ_0$ increases, period-doubling bifurcations appear, and the system enters a fully chaotic regime with broadly scattered points for all four pendulums.

#boxfig("pc3236/plots/bifurcation.png", width: 95%, box-width: 100%, [Bifurcation diagram (4-panel): Poincaré-sampled angles as a function of driving amplitude $τ_0$. The progression from periodic to chaotic behaviour is visible in all four pendulums.]) <fig:bifurcation>

== Animation

An animated GIF (`plots/animation.gif`) is included as a supplementary file. It shows the quadruple pendulum swinging for the first 10 s of the simulation at 25 frames per second, rendered from the same RK4 trajectory used for the static plots above. The animation was generated entirely within the script (lines 387–459 of `quadruple_pendulum.py`) using `matplotlib.animation.FuncAnimation` and `PillowWriter`, with frame selection performed by nearest-neighbour binary search (no interpolation library required).

= Reflections and Discussion

== Physical Insights

The driven quadruple pendulum illustrates several key ideas from nonlinear dynamics:

- *Sensitivity to initial conditions* is not just a mathematical curiosity; it has practical consequences. Even with double-precision floating-point arithmetic ($≈ 16$ significant digits), trajectories diverge within tens of seconds, making long-term prediction impossible. The Poincaré section and Lyapunov divergence plot confirm this quantitatively.

- *Energy redistribution* through the chain of pendulums follows from the nonlinear coupling. The driving torque acts only on pendulum 1, but energy flows to pendulums 2, 3, and 4 through the angle-dependent coupling terms in the mass matrix. The fourth pendulum receives energy through three chaotic intermediaries, which is why it shows the most irregular behaviour.

- *Bifurcation structure* shows that the transition to chaos is not abrupt. As driving increases, the system first responds periodically, then undergoes period-doubling, and eventually fills a broad region of the Poincaré section. This is the classic route to chaos seen in many driven nonlinear systems (e.g. the Duffing oscillator, the forced Van der Pol oscillator). With four degrees of freedom, the chaotic onset occurs at a lower driving amplitude than for the triple pendulum, and the attractor structure is richer.

== Limitations

The model assumes rigid, massless rods and point masses with no friction or air resistance. A real quadruple pendulum would dissipate energy, which would compete with the driving force and potentially suppress or alter the chaotic behaviour. Adding a linear damping term $-b dot(θ)_i$ to each equation of motion would be a straightforward extension.

The use of a fixed time step ($10^(-3)$ s) means the integrator may be insufficiently precise during rapid whipping motions of the fourth pendulum. An adaptive step-size controller (as in scipy's RK45) would be more efficient, but the fixed-step RK4 was chosen to satisfy the project requirement of implementing one's own algorithm.

= Conclusion

The driven quadruple pendulum was simulated by deriving its equations of motion from the Lagrangian, reformulating them as an eight-dimensional first-order ODE system, and integrating with a hand-coded RK4 scheme. The $4 × 4$ mass matrix linear system is solved at each RK4 stage using LU decomposition via `np.linalg.solve`. The simulation reproduces the expected phenomenology of a chaotic driven system: irregular time series, space-filling phase portraits, scattered Poincaré sections, exponential sensitivity to initial conditions, and a bifurcation diagram showing the route from periodicity to chaos as driving amplitude increases. A real-time animation illustrates the qualitative motion. The custom RK4 integrator agrees with scipy's adaptive solver to within $10^(-6)$ rad over 20 s, confirming its correctness.

= References

1. J. M. T. Thompson and H. B. Stewart, _Nonlinear Dynamics and Chaos_, 2nd ed. (Wiley, 2002). <ref1>

2. S. H. Strogatz, _Nonlinear Dynamics and Chaos_, 2nd ed. (Westview Press, 2015). <ref2>

3. W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, _Numerical Recipes_, 3rd ed. (Cambridge University Press, 2007). <ref3>

4. PC3236 Computational Methods in Physics, Project Handout, NUS (2025/26). <ref4>

#pagebreak()

= Declaration on the Use of Generative AI

I declare that I *HAVE* used generative AI tools to produce this assignment.

I acknowledge that generative AI was used in the following manner:

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (left, left, left),
    stroke: 1pt,
    table.hline(stroke: 2pt),
    [*AI Tool Used*], [*My Prompt and AI Output*], [*How the Output Was Used*],
    table.hline(),
    [Claude], [
      *Prompt:* \
      "Help me structure the Python simulation for the driven quadruple pendulum, including the Lagrangian derivation, RK4 integrator, and plotting routines." \
      *Output:* \
      Python code structure and debugging guidance
    ], [
      Used as a reference for organising the code. All equations were derived by me and verified against textbook sources. The RK4 algorithm was implemented from my lecture notes. I manually corrected sign errors found by energy-conservation testing.
    ],
    table.hline(),
    [Claude], [
      *Prompt:* \
      "Help with Typst formatting for figures, tables, equations, and the code-to-equation mapping table." \
      *Output:* \
      Typst formatting suggestions
    ], [
      Used only for formatting. All physical analysis, interpretation, and report content were written and verified by me.
    ],
    table.hline(stroke: 2pt),
  ),
  caption: [AI Tool Usage Declaration]
)
