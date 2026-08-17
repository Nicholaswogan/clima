#!/usr/bin/env python3
"""Toy smooth radiative-convective equilibrium model.

This script tests a temperature-only atmospheric evolution equation,

    C_i dT_i/dt = F_total[i + 1/2] - F_total[i - 1/2],

on a pressure-coordinate finite-volume grid. The model contains:

* fixed-composition, pure-absorption gray longwave radiative transfer;
* a prescribed downward shortwave flux with optional gray attenuation;
* a one-sided, continuously differentiable dry-convective flux;
* a finite-difference Jacobian and pseudo-transient continuation (PTC);
* automatic continuation in convective diffusivity until a requested lapse-rate
  tolerance is reached; and
* an alternative exact, energy-conserving dry convective adjustment, evolved
  explicitly, with constrained backward Euler, or as a projected residual with
  PTC.

It intentionally has no chemistry, condensation, clouds, or mixing-length theory.
The gray optical depths are prescribed functions of pressure and do not depend on
temperature.

Run:

    python smooth_rce_toy.py
    python smooth_rce_toy.py --no-plot

Only NumPy is required for a non-plotting run. Matplotlib is used when available.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import numpy as np


SIGMA = 5.670374419e-8  # Stefan-Boltzmann constant [W m-2 K-4]


@dataclass(frozen=True)
class Model:
    """Parameters and pressure grid for the toy atmosphere."""

    nlev: int = 40
    p_top: float = 100.0  # Pa
    p_surf: float = 1.0e5  # Pa
    t_surf: float = 288.0  # K; initial surface temperature
    surface_heat_capacity: float = 4.0e8  # J m-2 K-1; roughly a 100 m ocean
    gravity: float = 9.81  # m s-2
    # cp: float = 1004.0  # J kg-1 K-1
    cp: float = 1004.0*1.5 # J kg-1 K-1
    gas_constant: float = 287.0  # J kg-1 K-1
    tau_lw_surf: float = 4.0
    lw_pressure_exponent: float = 2.0
    tau_sw_surf: float = 0.0  # default: SW passes through atmosphere to surface
    sw_pressure_exponent: float = 2.0
    stellar_flux: float = 240.0  # globally averaged incident SW [W m-2]
    diffusivity_factor: float = 1.66
    k_conv: float = 10.0  # initial convective eddy diffusivity [m2 s-1]
    smooth_width: float = 5.0e-5  # smoothing width in log lapse rate

    def __post_init__(self) -> None:
        if self.nlev < 3:
            raise ValueError("nlev must be at least 3")
        if not 0.0 < self.p_top < self.p_surf:
            raise ValueError("Require 0 < p_top < p_surf")
        if self.cp <= self.gas_constant:
            raise ValueError("cp must exceed the specific gas constant")

    @property
    def p_edge(self) -> np.ndarray:
        """Cell-edge pressures ordered from TOA to surface."""

        return np.geomspace(self.p_top, self.p_surf, self.nlev + 1)

    @property
    def p_cell(self) -> np.ndarray:
        """Logarithmic cell-center pressures ordered from top to bottom."""

        edge = self.p_edge
        return np.sqrt(edge[:-1] * edge[1:])

    @property
    def atmospheric_heat_capacity(self) -> np.ndarray:
        """Layer heat capacities per unit horizontal area [J m-2 K-1]."""

        return self.cp * np.diff(self.p_edge) / self.gravity

    @property
    def heat_capacity(self) -> np.ndarray:
        """Heat capacities for atmospheric cells followed by the surface."""

        return np.concatenate(
            (self.atmospheric_heat_capacity, [self.surface_heat_capacity])
        )

    @property
    def nabla_ad(self) -> float:
        """Dry adiabatic logarithmic temperature gradient."""
        return self.gas_constant / self.cp


def one_sided_huber(x: np.ndarray, width: float) -> np.ndarray:
    """One-sided C1 approximation to ``max(x, 0)``.

    Unlike a softplus, this activation is exactly zero for stable and neutral
    interfaces. Consequently its neutral-state flux does not grow when the
    convective diffusivity is increased during continuation.
    """

    if width <= 0.0:
        return np.maximum(x, 0.0)

    result = np.zeros_like(x)
    transition = (x > 0.0) & (x < width)
    unstable = x >= width
    result[transition] = x[transition] ** 2 / (2.0 * width)
    result[unstable] = x[unstable] - 0.5 * width
    return result


def gray_radiation(
    model: Model, temperature: np.ndarray, surface_temperature: float
) -> dict[str, np.ndarray]:
    """Compute gray LW and prescribed SW fluxes at cell edges.

    Longwave transfer uses an isothermal-layer, pure-absorption two-stream
    recurrence. The diffusivity factor converts vertical optical depth into the
    effective slant optical depth of the hemispheric flux approximation.

    All returned net fluxes are positive upward.
    """

    p_edge = model.p_edge

    # A p^2-like optical-depth profile crudely represents pressure broadening and
    # deliberately makes the deep radiative equilibrium superadiabatic. With
    # tau proportional only to p, the classic deep gray radiative gradient tends
    # toward 1/4, which is slightly shallower than Earth's dry adiabat.
    pressure_coordinate = (
        (p_edge - model.p_top) / (model.p_surf - model.p_top)
    )
    tau_lw_edge = (
        model.tau_lw_surf
        * pressure_coordinate**model.lw_pressure_exponent
    )
    delta_tau = np.diff(tau_lw_edge)
    transmission = np.exp(-model.diffusivity_factor * delta_tau)
    source = SIGMA * temperature**4

    up_lw = np.zeros(model.nlev + 1)
    down_lw = np.zeros(model.nlev + 1)

    # Integrate surface thermal emission upward.
    up_lw[-1] = SIGMA * surface_temperature**4
    for i in range(model.nlev - 1, -1, -1):
        tr = transmission[i]
        up_lw[i] = tr * up_lw[i + 1] + (1.0 - tr) * source[i]

    # Integrate atmospheric back-radiation downward from a dark upper boundary.
    down_lw[0] = 0.0
    for i in range(model.nlev):
        tr = transmission[i]
        down_lw[i + 1] = tr * down_lw[i] + (1.0 - tr) * source[i]

    net_lw = up_lw - down_lw

    # Prescribed downward shortwave flux. It is independent of temperature.
    tau_sw_edge = model.tau_sw_surf * (
        pressure_coordinate**model.sw_pressure_exponent
    )
    down_sw = model.stellar_flux * np.exp(-tau_sw_edge)
    net_sw = -down_sw

    return {
        "up_lw": up_lw,
        "down_lw": down_lw,
        "net_lw": net_lw,
        "down_sw": down_sw,
        "net_sw": net_sw,
        "net_rad": net_lw + net_sw,
    }


def convective_flux(
    model: Model, temperature: np.ndarray, k_conv: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Compute smooth upward convective flux at cell edges.

    The internal-interface closure is

        F_conv = rho cp K (g/R) huber(nabla - nabla_ad).

    Convective flux is zero at the top and bottom boundaries. The returned
    superadiabaticity has one entry per internal interface.
    """

    if k_conv is None:
        k_conv = model.k_conv

    p = model.p_cell
    nabla = np.diff(np.log(temperature)) / np.diff(np.log(p))
    superadiabaticity = nabla - model.nabla_ad

    # Logarithmic interface means are natural on the logarithmic pressure grid.
    p_int = np.sqrt(p[:-1] * p[1:])
    t_int = np.sqrt(temperature[:-1] * temperature[1:])
    rho_int = p_int / (model.gas_constant * t_int)

    activation = one_sided_huber(superadiabaticity, model.smooth_width)
    flux = np.zeros(model.nlev + 1)
    flux[1:-1] = (
        rho_int
        * model.cp
        * k_conv
        * model.gravity
        / model.gas_constant
        * activation
    )
    return flux, superadiabaticity


def evaluate(
    model: Model, state: np.ndarray, k_conv: float | None = None
) -> dict[str, np.ndarray]:
    """Evaluate fluxes, finite-volume imbalance, and temperature tendency."""

    if state.shape != (model.nlev + 1,):
        raise ValueError(f"Expected state shape {(model.nlev + 1,)}")
    if np.any(~np.isfinite(state)) or np.any(state <= 0.0):
        raise ValueError("Temperatures must be positive and finite")

    temperature = state[:-1]
    surface_temperature = float(state[-1])
    radiation = gray_radiation(model, temperature, surface_temperature)
    flux_conv, superadiabaticity = convective_flux(model, temperature, k_conv)
    flux_total = radiation["net_rad"] + flux_conv

    # Lower-edge upward flux minus upper-edge upward flux heats a cell.
    atmospheric_imbalance = flux_total[1:] - flux_total[:-1]
    # The surface gains the energy represented by a downward net bottom flux.
    surface_imbalance = -flux_total[-1]
    imbalance = np.concatenate((atmospheric_imbalance, [surface_imbalance]))
    tendency = imbalance / model.heat_capacity

    return {
        **radiation,
        "conv": flux_conv,
        "total": flux_total,
        "imbalance": imbalance,
        "tendency": tendency,
        "superadiabaticity": superadiabaticity,
    }


def finite_difference_jacobian(
    model: Model,
    temperature: np.ndarray,
    k_conv: float,
    relative_step: float = 1.0e-5,
    absolute_step: float = 1.0e-3,
) -> np.ndarray:
    """Central-difference Jacobian of dT/dt with respect to temperature."""

    n = temperature.size
    jac = np.empty((n, n))
    for j in range(n):
        step = max(absolute_step, relative_step * abs(temperature[j]))
        forward = temperature.copy()
        backward = temperature.copy()
        forward[j] += step
        backward[j] -= step
        f_plus = evaluate(model, forward, k_conv)["tendency"]
        f_minus = evaluate(model, backward, k_conv)["tendency"]
        jac[:, j] = (f_plus - f_minus) / (2.0 * step)
    return jac


def dry_convective_projection(
    model: Model,
    temperature: np.ndarray,
) -> tuple[np.ndarray, list[tuple[int, int]], float]:
    """Project atmospheric temperatures onto the dry-stable set.

    This is the weighted pool-adjacent-violators algorithm used by standard dry
    convective adjustment. Potential temperature must be non-increasing from the
    top of this pressure-ordered grid toward the surface. The weights make each
    pooled block conserve its dry thermal energy exactly for constant cp.

    Returns the adjusted temperatures, inclusive bounds of mixed blocks, and the
    absolute column-energy conservation error in J/m2.
    """

    if temperature.shape != (model.nlev,):
        raise ValueError(f"Expected atmospheric shape {(model.nlev,)}")
    if np.any(~np.isfinite(temperature)) or np.any(temperature <= 0.0):
        raise ValueError("Temperatures must be positive and finite")

    exner = (model.p_cell / model.p_surf) ** model.nabla_ad
    theta = temperature / exner
    weights = model.atmospheric_heat_capacity * exner

    # Each entry is [start, end, summed weight, weighted-mean theta].
    blocks: list[list[float | int]] = []
    for i in range(model.nlev):
        blocks.append([i, i, float(weights[i]), float(theta[i])])
        while len(blocks) >= 2 and blocks[-2][3] < blocks[-1][3]:
            right = blocks.pop()
            left = blocks.pop()
            combined_weight = float(left[2]) + float(right[2])
            combined_theta = (
                float(left[2]) * float(left[3])
                + float(right[2]) * float(right[3])
            ) / combined_weight
            blocks.append(
                [int(left[0]), int(right[1]), combined_weight, combined_theta]
            )

    adjusted_theta = np.empty(model.nlev)
    mixed_blocks: list[tuple[int, int]] = []
    for start, end, _, block_theta in blocks:
        i_start = int(start)
        i_end = int(end)
        adjusted_theta[i_start : i_end + 1] = float(block_theta)
        if i_end > i_start:
            mixed_blocks.append((i_start, i_end))

    adjusted = adjusted_theta * exner
    energy_before = float(
        np.dot(model.atmospheric_heat_capacity, temperature)
    )
    energy_after = float(np.dot(model.atmospheric_heat_capacity, adjusted))
    return adjusted, mixed_blocks, abs(energy_after - energy_before)


def evaluate_projected_state(
    model: Model,
    state: np.ndarray,
    mixed_blocks: list[tuple[int, int]],
) -> dict[str, np.ndarray]:
    """Evaluate a projected state and diagnose its implied convective flux.

    Within each neutral mixed block, convection cancels variations of radiative
    flux relative to the block's upper boundary. The residual at the lower cell of
    the block is consequently the net radiative energy input to the whole block.
    It vanishes at steady radiative-convective equilibrium.
    """

    temperature = state[:-1]
    surface_temperature = float(state[-1])
    radiation = gray_radiation(model, temperature, surface_temperature)
    flux_conv = np.zeros(model.nlev + 1)

    for start, end in mixed_blocks:
        upper_edge = start
        for edge in range(start + 1, end + 1):
            flux_conv[edge] = (
                radiation["net_rad"][upper_edge]
                - radiation["net_rad"][edge]
            )

    flux_total = radiation["net_rad"] + flux_conv
    atmospheric_imbalance = flux_total[1:] - flux_total[:-1]
    surface_imbalance = -flux_total[-1]
    imbalance = np.concatenate((atmospheric_imbalance, [surface_imbalance]))
    tendency = imbalance / model.heat_capacity
    _, superadiabaticity = convective_flux(model, temperature, 0.0)

    return {
        **radiation,
        "conv": flux_conv,
        "total": flux_total,
        "imbalance": imbalance,
        "tendency": tendency,
        "superadiabaticity": superadiabaticity,
    }


def _projected_residual_only(
    model: Model,
    state: np.ndarray,
    projection_time: float,
) -> tuple[np.ndarray, list[tuple[int, int]], float]:
    """Evaluate the projected residual without extra flux diagnostics."""

    if projection_time <= 0.0:
        raise ValueError("projection_time must be positive")

    radiative = evaluate(model, state, 0.0)
    trial = state + projection_time * radiative["tendency"]
    if np.any(~np.isfinite(trial)) or np.any(trial <= 0.0):
        raise ValueError("Projected residual produced an invalid trial state")

    adjusted, mixed_blocks, energy_error = dry_convective_projection(
        model, trial[:-1]
    )
    mapped = trial.copy()
    mapped[:-1] = adjusted
    residual = (mapped - state) / projection_time
    return residual, mixed_blocks, energy_error


def projected_natural_residual(
    model: Model,
    state: np.ndarray,
    projection_time: float,
) -> tuple[
    np.ndarray,
    dict[str, np.ndarray],
    list[tuple[int, int]],
    float,
]:
    """Evaluate a natural residual and diagnostics for constrained equilibrium.

    A root of

        G(T) = [P(T + alpha f_rad(T)) - T] / alpha

    is a fixed point of an energy-conserving convective adjustment applied to a
    radiative trial step. Unlike appending a projection to a Newton correction,
    the roots of G represent the constrained steady problem for any positive
    projection time alpha.
    """

    residual, mixed_blocks, energy_error = _projected_residual_only(
        model, state, projection_time
    )
    diagnostics = evaluate_projected_state(model, state, mixed_blocks)
    return residual, diagnostics, mixed_blocks, energy_error


def forward_difference_projected_jacobian(
    model: Model,
    state: np.ndarray,
    projection_time: float,
    base_residual: np.ndarray,
    relative_step: float = 1.0e-5,
    absolute_step: float = 1.0e-3,
) -> np.ndarray:
    """Forward-difference Jacobian of the projected natural residual."""

    n = state.size
    jac = np.empty((n, n))
    for j in range(n):
        step = max(absolute_step, relative_step * abs(state[j]))
        perturbed = state.copy()
        perturbed[j] += step
        residual, _, _ = _projected_residual_only(
            model, perturbed, projection_time
        )
        jac[:, j] = (residual - base_residual) / step
    return jac


@dataclass
class PTCResult:
    temperature: np.ndarray
    diagnostics: dict[str, np.ndarray]
    steps: int
    jacobian_evaluations: int
    rejected_steps: int
    converged: bool


@dataclass
class ContinuationResult:
    """Result and aggregate work from continuation in convective diffusivity."""

    result: PTCResult
    k_conv: float
    stages: int
    total_steps: int
    total_jacobian_evaluations: int
    total_rejected_steps: int
    converged: bool
    max_superadiabaticity_encountered: float = 0.0
    guard_triggered_updates: int = 0
    equilibrium_triggered_updates: int = 0


@dataclass
class ProjectedSolveResult:
    """Result and work diagnostics for a projected convective solve."""

    result: PTCResult
    projection_evaluations: int
    max_projection_energy_error: float
    max_projection_relative_energy_error: float
    radiative_evaluations: int = 0
    nonlinear_iterations: int = 0
    simulated_time: float = 0.0
    final_timestep: float = 0.0
    max_estimated_local_error: float = 0.0


def solve_ptc(
    model: Model,
    initial_temperature: np.ndarray,
    k_conv: float,
    *,
    flux_tolerance: float = 1.0e-5,
    max_steps: int = 250,
    dt_initial: float = 100.0,
    dt_max: float = 1.0e12,
    growth: float = 1.5,
    verbose: bool = True,
) -> PTCResult:
    """Drive the atmospheric ODE to a steady state with dense PTC."""

    temperature = np.asarray(initial_temperature, dtype=float).copy()
    if temperature.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    dt = dt_initial
    rejected = 0
    njev = 0

    diagnostics = evaluate(model, temperature, k_conv)
    tendency_norm = np.linalg.norm(diagnostics["tendency"])

    if verbose:
        print(
            f"{'step':>5s} {'dt [s]':>12s} {'max |dF| [W/m2]':>18s} "
            f"{'max superad.':>14s}"
        )

    for step_number in range(max_steps + 1):
        max_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))
        max_superadiabaticity = float(
            max(0.0, np.max(diagnostics["superadiabaticity"]))
        )

        if verbose and (step_number < 8 or step_number % 10 == 0):
            print(
                f"{step_number:5d} {dt:12.4e} {max_flux_error:18.6e} "
                f"{max_superadiabaticity:14.6e}"
            )

        if max_flux_error < flux_tolerance:
            if verbose:
                print(
                    f"Converged in {step_number} accepted PTC steps "
                    f"({rejected} rejected)."
                )
            return PTCResult(
                temperature,
                diagnostics,
                step_number,
                njev,
                rejected,
                True,
            )

        if step_number == max_steps:
            break

        jac = finite_difference_jacobian(model, temperature, k_conv)
        njev += 1

        # Retry the same state with progressively smaller pseudo-timesteps.
        accepted = False
        for _ in range(15):
            matrix = np.eye(model.nlev + 1) / dt - jac
            try:
                delta_temperature = np.linalg.solve(
                    matrix, diagnostics["tendency"]
                )
            except np.linalg.LinAlgError:
                dt *= 0.25
                rejected += 1
                continue

            candidate = temperature + delta_temperature
            if (
                np.any(~np.isfinite(candidate))
                or np.min(candidate) < 50.0
                or np.max(candidate) > 1000.0
            ):
                dt *= 0.25
                rejected += 1
                continue

            candidate_diagnostics = evaluate(model, candidate, k_conv)
            candidate_norm = np.linalg.norm(candidate_diagnostics["tendency"])

            # A modest increase is allowed because PTC is not a line search.
            if candidate_norm > 1.25 * tendency_norm:
                dt *= 0.5
                rejected += 1
                continue

            accepted = True
            break

        if not accepted:
            raise RuntimeError("PTC failed to find an acceptable pseudo-step")

        old_norm = tendency_norm
        temperature = candidate
        diagnostics = candidate_diagnostics
        tendency_norm = candidate_norm

        # Switched-evolution-relaxation-style pseudo-timestep update.
        if tendency_norm > 0.0:
            ratio = old_norm / tendency_norm
            dt_factor = np.clip(growth * ratio, 0.5, 5.0)
        else:
            dt_factor = 5.0
        dt = min(dt_max, dt * dt_factor)

    if verbose:
        print(
            f"Did not converge in {max_steps} accepted PTC steps; "
            f"max flux imbalance is "
            f"{np.max(np.abs(diagnostics['imbalance'])):.6e} W/m2."
        )
    return PTCResult(
        temperature, diagnostics, max_steps, njev, rejected, False
    )


def solve_projected_explicit(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    flux_tolerance: float = 1.0e-5,
    tendency_tolerance: float = 1.0e-13,
    max_steps: int = 200_000,
    dt_initial: float = 100.0,
    dt_max: float = 1.0e6,
    growth: float = 1.5,
    max_temperature_step: float = 1.0,
    verbose: bool = True,
) -> ProjectedSolveResult:
    """Explicit radiative evolution followed by exact convective adjustment.

    A single global timestep is limited so that the largest temperature change
    after adjustment does not exceed ``max_temperature_step``. This retains the
    projected forward-Euler map instead of clipping individual layer updates.
    """

    if max_temperature_step <= 0.0:
        raise ValueError("max_temperature_step must be positive")
    if tendency_tolerance <= 0.0:
        raise ValueError("tendency_tolerance must be positive")
    if not 0.0 < dt_initial <= dt_max:
        raise ValueError("Require 0 < dt_initial <= dt_max")

    state = np.asarray(initial_temperature, dtype=float).copy()
    if state.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    adjusted, mixed_blocks, initial_energy_error = dry_convective_projection(
        model, state[:-1]
    )
    state[:-1] = adjusted
    diagnostics = evaluate_projected_state(model, state, mixed_blocks)
    radiative_evaluations = 1
    projection_evaluations = 1
    max_energy_error = initial_energy_error
    column_energy = abs(
        float(np.dot(model.atmospheric_heat_capacity, state[:-1]))
    )
    max_relative_energy_error = initial_energy_error / max(
        column_energy, np.finfo(float).tiny
    )
    dt = dt_initial
    rejected = 0
    last_map_residual = np.inf

    if verbose:
        print(
            f"{'step':>7s} {'dt [s]':>12s} {'mixed blocks':>13s} "
            f"{'max |dT/dt| [K/s]':>20s} {'max |dF| [W/m2]':>18s}"
        )

    for step_number in range(max_steps + 1):
        max_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))
        if verbose and (step_number < 8 or step_number % 1000 == 0):
            print(
                f"{step_number:7d} {dt:12.4e} {len(mixed_blocks):13d} "
                f"{last_map_residual:20.6e} {max_flux_error:18.6e}"
            )

        if (
            max_flux_error < flux_tolerance
            and last_map_residual < tendency_tolerance
        ):
            if verbose:
                print(
                    f"Converged in {step_number} explicit projected steps "
                    f"({rejected} rejected)."
                )
            result = PTCResult(
                state, diagnostics, step_number, 0, rejected, True
            )
            return ProjectedSolveResult(
                result,
                projection_evaluations,
                max_energy_error,
                max_relative_energy_error,
                radiative_evaluations,
            )

        if step_number == max_steps:
            break

        net_rad = diagnostics["net_rad"]
        radiative_imbalance = np.concatenate(
            (net_rad[1:] - net_rad[:-1], [-net_rad[-1]])
        )
        radiative_tendency = radiative_imbalance / model.heat_capacity
        accepted = False
        for _ in range(15):
            trial = state + dt * radiative_tendency
            if (
                np.any(~np.isfinite(trial))
                or np.min(trial) < 50.0
                or np.max(trial) > 1000.0
            ):
                dt *= 0.25
                rejected += 1
                continue

            adjusted, candidate_blocks, energy_error = dry_convective_projection(
                model, trial[:-1]
            )
            projection_evaluations += 1
            candidate = trial.copy()
            candidate[:-1] = adjusted
            max_temperature_change = float(
                np.max(np.abs(candidate - state))
            )
            if max_temperature_change > max_temperature_step:
                dt *= max(
                    0.1,
                    0.9 * max_temperature_step / max_temperature_change,
                )
                rejected += 1
                continue
            accepted = True
            break

        if not accepted:
            raise RuntimeError(
                "Explicit projected solve failed to find a valid timestep"
            )

        last_map_residual = float(np.max(np.abs((candidate - state) / dt)))
        state = candidate
        mixed_blocks = candidate_blocks
        diagnostics = evaluate_projected_state(model, state, mixed_blocks)
        radiative_evaluations += 1

        candidate_energy = abs(
            float(np.dot(model.atmospheric_heat_capacity, state[:-1]))
        )
        max_energy_error = max(max_energy_error, energy_error)
        max_relative_energy_error = max(
            max_relative_energy_error,
            energy_error
            / max(candidate_energy, np.finfo(float).tiny),
        )
        if max_temperature_change > 0.0:
            controller_growth = np.clip(
                0.9 * max_temperature_step / max_temperature_change,
                0.5,
                growth,
            )
        else:
            controller_growth = growth
        dt = min(dt_max, dt * float(controller_growth))

    if verbose:
        print(
            f"Explicit projected solve did not converge in {max_steps} steps; "
            f"max flux imbalance is "
            f"{np.max(np.abs(diagnostics['imbalance'])):.6e} W/m2."
        )
    result = PTCResult(state, diagnostics, max_steps, 0, rejected, False)
    return ProjectedSolveResult(
        result,
        projection_evaluations,
        max_energy_error,
        max_relative_energy_error,
        radiative_evaluations,
    )


@dataclass
class _ImplicitStepResult:
    """Internal result from one constrained backward-Euler solve."""

    state: np.ndarray
    mixed_blocks: list[tuple[int, int]]
    converged: bool
    nonlinear_iterations: int
    jacobian_evaluations: int
    radiative_evaluations: int
    projection_evaluations: int
    max_projection_energy_error: float


def _constrained_backward_euler_residual(
    model: Model,
    old_state: np.ndarray,
    new_state: np.ndarray,
    dt: float,
) -> tuple[np.ndarray, list[tuple[int, int]], float]:
    """Residual of one exact-adjustment backward-Euler timestep.

    The root satisfies

        T_new = P[T_old + dt f_rad(T_new)],

    where P acts on atmospheric temperatures and leaves the surface unchanged.
    Radiation is evaluated at the new state, while the dry adjustment exactly
    conserves the atmospheric energy of the implicit radiative trial state.
    """

    if dt <= 0.0:
        raise ValueError("dt must be positive")
    radiative = evaluate(model, new_state, 0.0)
    trial = old_state + dt * radiative["tendency"]
    if np.any(~np.isfinite(trial)) or np.any(trial <= 0.0):
        raise ValueError("Implicit radiative trial state is invalid")

    adjusted, mixed_blocks, energy_error = dry_convective_projection(
        model, trial[:-1]
    )
    mapped = trial.copy()
    mapped[:-1] = adjusted
    return new_state - mapped, mixed_blocks, energy_error


def _solve_constrained_backward_euler_step(
    model: Model,
    old_state: np.ndarray,
    dt: float,
    *,
    initial_guess: np.ndarray | None = None,
    nonlinear_tolerance: float = 1.0e-8,
    max_iterations: int = 18,
) -> _ImplicitStepResult:
    """Solve one nonsmooth constrained backward-Euler step with damped Newton."""

    n = old_state.size
    if initial_guess is None:
        new_state = old_state.copy()
    else:
        new_state = np.asarray(initial_guess, dtype=float).copy()
    if new_state.shape != old_state.shape:
        raise ValueError("Initial guess has the wrong shape")

    radiative_evaluations = 0
    projection_evaluations = 0
    jacobian_evaluations = 0
    max_energy_error = 0.0
    mixed_blocks: list[tuple[int, int]] = []

    def residual_at(
        state: np.ndarray,
    ) -> tuple[np.ndarray, list[tuple[int, int]], float] | None:
        nonlocal radiative_evaluations, projection_evaluations, max_energy_error
        if (
            np.any(~np.isfinite(state))
            or np.min(state) < 50.0
            or np.max(state) > 1000.0
        ):
            return None
        try:
            residual, blocks, energy_error = (
                _constrained_backward_euler_residual(
                    model, old_state, state, dt
                )
            )
        except ValueError:
            # Radiation was evaluated before the radiative trial was rejected.
            radiative_evaluations += 1
            return None
        radiative_evaluations += 1
        projection_evaluations += 1
        max_energy_error = max(max_energy_error, energy_error)
        return residual, blocks, energy_error

    base = residual_at(new_state)
    if base is None:
        return _ImplicitStepResult(
            new_state,
            mixed_blocks,
            False,
            0,
            0,
            radiative_evaluations,
            projection_evaluations,
            max_energy_error,
        )
    residual, mixed_blocks, _ = base

    for iteration in range(max_iterations + 1):
        residual_norm = float(np.max(np.abs(residual)))
        if residual_norm <= nonlinear_tolerance:
            return _ImplicitStepResult(
                new_state,
                mixed_blocks,
                True,
                iteration,
                jacobian_evaluations,
                radiative_evaluations,
                projection_evaluations,
                max_energy_error,
            )
        if iteration == max_iterations:
            break

        jacobian = np.empty((n, n))
        jacobian_ok = True
        for j in range(n):
            difference_step = max(1.0e-3, 1.0e-5 * abs(new_state[j]))
            perturbed = new_state.copy()
            perturbed[j] += difference_step
            perturbed_result = residual_at(perturbed)
            if perturbed_result is None:
                perturbed[j] = new_state[j] - difference_step
                perturbed_result = residual_at(perturbed)
                if perturbed_result is None:
                    jacobian_ok = False
                    break
                jacobian[:, j] = (
                    residual - perturbed_result[0]
                ) / difference_step
            else:
                jacobian[:, j] = (
                    perturbed_result[0] - residual
                ) / difference_step
        jacobian_evaluations += 1
        if not jacobian_ok:
            break

        try:
            correction = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError:
            break

        accepted = False
        damping = 1.0
        for _ in range(16):
            candidate = new_state + damping * correction
            candidate_result = residual_at(candidate)
            if candidate_result is not None:
                candidate_residual, candidate_blocks, _ = candidate_result
                candidate_norm = float(
                    np.max(np.abs(candidate_residual))
                )
                # A weak Armijo condition works across PAVA active-set changes.
                if candidate_norm < (1.0 - 1.0e-4 * damping) * residual_norm:
                    new_state = candidate
                    residual = candidate_residual
                    mixed_blocks = candidate_blocks
                    accepted = True
                    break
            damping *= 0.5
        if not accepted:
            break

    return _ImplicitStepResult(
        new_state,
        mixed_blocks,
        False,
        min(iteration + 1, max_iterations),
        jacobian_evaluations,
        radiative_evaluations,
        projection_evaluations,
        max_energy_error,
    )


def solve_projected_implicit(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    flux_tolerance: float = 1.0e-5,
    tendency_tolerance: float = 1.0e-13,
    max_steps: int = 1000,
    dt_initial: float = 100.0,
    dt_max: float = 1.0e12,
    temperature_tolerance: float = 2.0e-2,
    nonlinear_tolerance: float = 1.0e-8,
    end_time: float | None = None,
    verbose: bool = True,
) -> ProjectedSolveResult:
    """Integrate exact dry adjustment with constrained backward Euler.

    A full step and two half steps estimate the local temporal error. The two
    half-step solution is accepted, providing a conservative, feasible physical
    trajectory and separating timestep accuracy from nonlinear convergence.
    """

    if temperature_tolerance <= 0.0:
        raise ValueError("temperature_tolerance must be positive")
    if flux_tolerance <= 0.0:
        raise ValueError("flux_tolerance must be positive")
    if tendency_tolerance <= 0.0:
        raise ValueError("tendency_tolerance must be positive")
    if nonlinear_tolerance <= 0.0:
        raise ValueError("nonlinear_tolerance must be positive")
    if max_steps < 1:
        raise ValueError("max_steps must be positive")
    if not 0.0 < dt_initial <= dt_max:
        raise ValueError("Require 0 < dt_initial <= dt_max")
    if end_time is not None and end_time <= 0.0:
        raise ValueError("end_time must be positive")

    state = np.asarray(initial_temperature, dtype=float).copy()
    if state.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    adjusted, mixed_blocks, initial_energy_error = dry_convective_projection(
        model, state[:-1]
    )
    state[:-1] = adjusted
    diagnostics = evaluate_projected_state(model, state, mixed_blocks)
    radiative_evaluations = 1
    projection_evaluations = 1
    jacobian_evaluations = 0
    nonlinear_iterations = 0
    rejected = 0
    max_energy_error = initial_energy_error
    max_relative_energy_error = 0.0
    max_estimated_error = 0.0
    last_estimated_error = 0.0
    simulated_time = 0.0
    dt = dt_initial
    last_rate = np.inf

    if verbose:
        print(
            f"{'step':>6s} {'time [s]':>12s} {'dt [s]':>12s} "
            f"{'mixed blocks':>13s} {'error [K]':>12s} "
            f"{'max |dF| [W/m2]':>18s}"
        )

    for step_number in range(max_steps + 1):
        max_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))
        if verbose and (step_number < 8 or step_number % 10 == 0):
            print(
                f"{step_number:6d} {simulated_time:12.4e} {dt:12.4e} "
                f"{len(mixed_blocks):13d} {last_estimated_error:12.4e} "
                f"{max_flux_error:18.6e}"
            )

        reached_end_time = end_time is not None and simulated_time >= (
            end_time
            - 16.0 * np.finfo(float).eps * max(1.0, abs(end_time))
        )
        if reached_end_time:
            if verbose:
                print(
                    f"Reached t = {simulated_time:.6e} s in {step_number} "
                    f"constrained backward-Euler steps ({rejected} rejected)."
                )
            result = PTCResult(
                state,
                diagnostics,
                step_number,
                jacobian_evaluations,
                rejected,
                True,
            )
            return ProjectedSolveResult(
                result,
                projection_evaluations,
                max_energy_error,
                max_relative_energy_error,
                radiative_evaluations,
                nonlinear_iterations,
                simulated_time,
                dt,
                max_estimated_error,
            )

        # A flux tolerance implies a looser tendency tolerance in layers with
        # small heat capacity. Do not demand a temporal stopping criterion that
        # is inconsistent with the requested finite-volume flux balance.
        effective_tendency_tolerance = max(
            tendency_tolerance,
            flux_tolerance / float(np.min(model.heat_capacity)),
        )
        if end_time is None and (
            max_flux_error < flux_tolerance
            and last_rate < effective_tendency_tolerance
        ):
            if verbose:
                print(
                    f"Converged in {step_number} constrained backward-Euler "
                    f"steps ({rejected} rejected)."
                )
            result = PTCResult(
                state,
                diagnostics,
                step_number,
                jacobian_evaluations,
                rejected,
                True,
            )
            return ProjectedSolveResult(
                result,
                projection_evaluations,
                max_energy_error,
                max_relative_energy_error,
                radiative_evaluations,
                nonlinear_iterations,
                simulated_time,
                dt,
                max_estimated_error,
            )
        if step_number == max_steps:
            break

        if end_time is not None:
            dt = min(dt, end_time - simulated_time)

        accepted = False
        for _ in range(20):
            full = _solve_constrained_backward_euler_step(
                model,
                state,
                dt,
                initial_guess=state,
                nonlinear_tolerance=nonlinear_tolerance,
            )
            half_1 = _solve_constrained_backward_euler_step(
                model,
                state,
                0.5 * dt,
                initial_guess=state,
                nonlinear_tolerance=nonlinear_tolerance,
            )
            for solve in (full, half_1):
                radiative_evaluations += solve.radiative_evaluations
                projection_evaluations += solve.projection_evaluations
                jacobian_evaluations += solve.jacobian_evaluations
                nonlinear_iterations += solve.nonlinear_iterations
                max_energy_error = max(
                    max_energy_error, solve.max_projection_energy_error
                )
            if not full.converged or not half_1.converged:
                dt *= 0.25
                rejected += 1
                continue

            half_2 = _solve_constrained_backward_euler_step(
                model,
                half_1.state,
                0.5 * dt,
                initial_guess=full.state,
                nonlinear_tolerance=nonlinear_tolerance,
            )
            radiative_evaluations += half_2.radiative_evaluations
            projection_evaluations += half_2.projection_evaluations
            jacobian_evaluations += half_2.jacobian_evaluations
            nonlinear_iterations += half_2.nonlinear_iterations
            max_energy_error = max(
                max_energy_error, half_2.max_projection_energy_error
            )
            if not half_2.converged:
                dt *= 0.25
                rejected += 1
                continue

            estimated_error = float(
                np.max(np.abs(half_2.state - full.state))
            )
            if estimated_error > temperature_tolerance:
                factor = max(
                    0.1,
                    0.9 * np.sqrt(temperature_tolerance / estimated_error),
                )
                dt *= factor
                rejected += 1
                continue
            accepted = True
            break

        if not accepted:
            raise RuntimeError(
                "Constrained backward Euler failed to find an acceptable step"
            )

        old_state = state
        state = half_2.state
        mixed_blocks = half_2.mixed_blocks
        diagnostics = evaluate_projected_state(model, state, mixed_blocks)
        radiative_evaluations += 1
        simulated_time += dt
        last_rate = float(np.max(np.abs(state - old_state))) / dt
        last_estimated_error = estimated_error
        max_estimated_error = max(max_estimated_error, estimated_error)

        column_energy = abs(
            float(np.dot(model.atmospheric_heat_capacity, state[:-1]))
        )
        max_relative_energy_error = max(
            max_relative_energy_error,
            max_energy_error / max(column_energy, np.finfo(float).tiny),
        )

        if estimated_error > 0.0:
            factor = float(
                np.clip(
                    0.9 * np.sqrt(temperature_tolerance / estimated_error),
                    0.5,
                    2.0,
                )
            )
        else:
            factor = 2.0
        dt = min(dt_max, dt * factor)

    if verbose:
        print(
            f"Constrained backward Euler did not converge in {max_steps} "
            f"steps; max flux imbalance is {np.max(np.abs(diagnostics['imbalance'])):.6e} W/m2."
        )
    result = PTCResult(
        state,
        diagnostics,
        max_steps,
        jacobian_evaluations,
        rejected,
        False,
    )
    return ProjectedSolveResult(
        result,
        projection_evaluations,
        max_energy_error,
        max_relative_energy_error,
        radiative_evaluations,
        nonlinear_iterations,
        simulated_time,
        dt,
        max_estimated_error,
    )


def integrate_projected_explicit_fixed_time(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    end_time: float,
    dt: float,
) -> ProjectedSolveResult:
    """Integrate projected forward Euler to an exact physical end time."""

    if end_time <= 0.0 or dt <= 0.0:
        raise ValueError("end_time and dt must be positive")
    state = np.asarray(initial_temperature, dtype=float).copy()
    if state.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    adjusted, mixed_blocks, initial_energy_error = dry_convective_projection(
        model, state[:-1]
    )
    state[:-1] = adjusted
    max_energy_error = initial_energy_error
    max_relative_energy_error = 0.0
    projection_evaluations = 1
    radiative_evaluations = 0
    simulated_time = 0.0
    steps = 0

    while simulated_time < end_time:
        step_dt = min(dt, end_time - simulated_time)
        radiative = evaluate(model, state, 0.0)
        radiative_evaluations += 1
        trial = state + step_dt * radiative["tendency"]
        if (
            np.any(~np.isfinite(trial))
            or np.min(trial) < 50.0
            or np.max(trial) > 1000.0
        ):
            raise RuntimeError(
                f"Projected explicit integration became invalid at "
                f"t={simulated_time:.6e} s with dt={step_dt:.6e} s"
            )
        adjusted, mixed_blocks, energy_error = dry_convective_projection(
            model, trial[:-1]
        )
        projection_evaluations += 1
        trial[:-1] = adjusted
        state = trial
        simulated_time += step_dt
        steps += 1
        max_energy_error = max(max_energy_error, energy_error)
        column_energy = abs(
            float(np.dot(model.atmospheric_heat_capacity, state[:-1]))
        )
        max_relative_energy_error = max(
            max_relative_energy_error,
            energy_error / max(column_energy, np.finfo(float).tiny),
        )

    diagnostics = evaluate_projected_state(model, state, mixed_blocks)
    radiative_evaluations += 1
    result = PTCResult(state, diagnostics, steps, 0, 0, True)
    return ProjectedSolveResult(
        result,
        projection_evaluations,
        max_energy_error,
        max_relative_energy_error,
        radiative_evaluations,
        simulated_time=simulated_time,
        final_timestep=dt,
    )


@dataclass
class TimeBenchmarkRow:
    """Accuracy and work for one fixed-end-time projected integration."""

    method: str
    control: float
    steps: int
    radiative_evaluations: int
    wall_time: float
    max_atmospheric_error: float
    rms_atmospheric_error: float
    surface_error: float


def _time_benchmark_errors(
    model: Model,
    state: np.ndarray,
    reference: np.ndarray,
) -> tuple[float, float, float]:
    difference = state - reference
    weights = model.atmospheric_heat_capacity
    rms_atmospheric_error = float(
        np.sqrt(np.dot(weights, difference[:-1] ** 2) / np.sum(weights))
    )
    return (
        float(np.max(np.abs(difference[:-1]))),
        rms_atmospheric_error,
        float(abs(difference[-1])),
    )


def run_time_integration_benchmark(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    end_time: float = 1.0e8,
    reference_dt: float = 2.5e4,
    explicit_timesteps: tuple[float, ...] = (2.5e5, 5.0e5, 1.0e6, 2.0e6),
    implicit_tolerances: tuple[float, ...] = (2.0e-1, 5.0e-2, 2.0e-2),
) -> list[TimeBenchmarkRow]:
    """Compare projected explicit and implicit evolution at one physical time.

    The reference is projected forward Euler with half ``reference_dt``. A run
    at ``reference_dt`` reports the remaining reference discretization change.
    Every candidate starts from the same dry-stable initial profile and reaches
    exactly ``end_time``; no equilibrium convergence criterion is involved.
    """

    if end_time <= 0.0 or reference_dt <= 0.0:
        raise ValueError("Benchmark times must be positive")

    print(
        f"\nFixed-time projected-integration benchmark: nlev={model.nlev}, "
        f"t_end={end_time:.6e} s"
    )
    start = perf_counter()
    coarse_reference = integrate_projected_explicit_fixed_time(
        model,
        initial_temperature,
        end_time=end_time,
        dt=reference_dt,
    )
    coarse_wall = perf_counter() - start
    start = perf_counter()
    fine_reference = integrate_projected_explicit_fixed_time(
        model,
        initial_temperature,
        end_time=end_time,
        dt=0.5 * reference_dt,
    )
    fine_wall = perf_counter() - start
    reference_state = fine_reference.result.temperature
    reference_change = _time_benchmark_errors(
        model, coarse_reference.result.temperature, reference_state
    )
    print(
        "Reference refinement, dt -> dt/2: "
        f"max atmosphere={reference_change[0]:.3e} K, "
        f"atmosphere RMS={reference_change[1]:.3e} K, "
        f"surface={reference_change[2]:.3e} K"
    )
    print(
        f"Reference work: {coarse_reference.radiative_evaluations} + "
        f"{fine_reference.radiative_evaluations} RT evaluations, "
        f"{coarse_wall + fine_wall:.3f} s"
    )

    rows: list[TimeBenchmarkRow] = []
    for explicit_dt in explicit_timesteps:
        start = perf_counter()
        candidate = integrate_projected_explicit_fixed_time(
            model,
            initial_temperature,
            end_time=end_time,
            dt=explicit_dt,
        )
        wall_time = perf_counter() - start
        errors = _time_benchmark_errors(
            model, candidate.result.temperature, reference_state
        )
        rows.append(
            TimeBenchmarkRow(
                "explicit dt",
                explicit_dt,
                candidate.result.steps,
                candidate.radiative_evaluations,
                wall_time,
                *errors,
            )
        )

    for tolerance in implicit_tolerances:
        start = perf_counter()
        candidate = solve_projected_implicit(
            model,
            initial_temperature,
            max_steps=5000,
            dt_initial=min(100.0, end_time),
            dt_max=end_time,
            temperature_tolerance=tolerance,
            nonlinear_tolerance=min(1.0e-10, 1.0e-3 * tolerance),
            end_time=end_time,
            verbose=False,
        )
        wall_time = perf_counter() - start
        if not candidate.result.converged:
            raise RuntimeError(
                "Implicit benchmark failed to reach the requested end time"
            )
        errors = _time_benchmark_errors(
            model, candidate.result.temperature, reference_state
        )
        rows.append(
            TimeBenchmarkRow(
                "implicit tol",
                tolerance,
                candidate.result.steps,
                candidate.radiative_evaluations,
                wall_time,
                *errors,
            )
        )

    print(
        f"\n{'method':>14s} {'control':>12s} {'steps':>8s} {'RT evals':>10s} "
        f"{'wall [s]':>10s} {'max atm [K]':>13s} {'rms atm [K]':>13s} "
        f"{'surface [K]':>12s}"
    )
    for row in rows:
        print(
            f"{row.method:>14s} {row.control:12.4e} {row.steps:8d} "
            f"{row.radiative_evaluations:10d} {row.wall_time:10.3f} "
            f"{row.max_atmospheric_error:13.4e} "
            f"{row.rms_atmospheric_error:13.4e} {row.surface_error:12.4e}"
        )
    return rows


def solve_projected_ptc(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    projection_time: float = 1.0e4,
    flux_tolerance: float = 1.0e-5,
    max_steps: int = 300,
    dt_initial: float = 100.0,
    dt_max: float = 1.0e12,
    growth: float = 1.5,
    verbose: bool = True,
) -> ProjectedSolveResult:
    """Apply PTC to the natural residual of exact dry convective adjustment."""

    state = np.asarray(initial_temperature, dtype=float).copy()
    if state.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    # Begin in the feasible set without changing atmospheric thermal energy.
    adjusted, _, initial_energy_error = dry_convective_projection(
        model, state[:-1]
    )
    state[:-1] = adjusted
    residual, diagnostics, mixed_blocks, energy_error = (
        projected_natural_residual(model, state, projection_time)
    )
    residual_norm = np.linalg.norm(residual)

    column_energy = abs(
        float(np.dot(model.atmospheric_heat_capacity, state[:-1]))
    )
    max_energy_error = max(initial_energy_error, energy_error)
    max_relative_energy_error = max_energy_error / max(
        column_energy, np.finfo(float).tiny
    )
    projection_evaluations = 2
    radiative_evaluations = 2
    dt = dt_initial
    rejected = 0
    njev = 0

    if verbose:
        print(
            f"{'step':>5s} {'dt [s]':>12s} {'mixed blocks':>13s} "
            f"{'max |G| [K/s]':>16s} {'max |dF| [W/m2]':>18s}"
        )

    for step_number in range(max_steps + 1):
        max_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))
        max_residual = float(np.max(np.abs(residual)))

        if verbose and (step_number < 8 or step_number % 10 == 0):
            print(
                f"{step_number:5d} {dt:12.4e} {len(mixed_blocks):13d} "
                f"{max_residual:16.6e} {max_flux_error:18.6e}"
            )

        if max_flux_error < flux_tolerance and max_residual < 1.0e-10:
            if verbose:
                print(
                    f"Converged in {step_number} accepted projected PTC steps "
                    f"({rejected} rejected)."
                )
            result = PTCResult(
                state,
                diagnostics,
                step_number,
                njev,
                rejected,
                True,
            )
            return ProjectedSolveResult(
                result,
                projection_evaluations,
                max_energy_error,
                max_relative_energy_error,
                radiative_evaluations,
            )

        if step_number == max_steps:
            break

        jac = forward_difference_projected_jacobian(
            model,
            state,
            projection_time,
            residual,
        )
        njev += 1
        projection_evaluations += model.nlev + 1
        radiative_evaluations += model.nlev + 1

        accepted = False
        for _ in range(15):
            matrix = np.eye(model.nlev + 1) / dt - jac
            try:
                delta_temperature = np.linalg.solve(matrix, residual)
            except np.linalg.LinAlgError:
                dt *= 0.25
                rejected += 1
                continue

            candidate = state + delta_temperature
            if (
                np.any(~np.isfinite(candidate))
                or np.min(candidate) < 50.0
                or np.max(candidate) > 1000.0
            ):
                dt *= 0.25
                rejected += 1
                continue

            # Keep every accepted iterate in the dry-stable feasible set.
            adjusted, _, candidate_projection_error = dry_convective_projection(
                model, candidate[:-1]
            )
            projection_evaluations += 1
            candidate[:-1] = adjusted

            try:
                (
                    candidate_residual,
                    candidate_diagnostics,
                    candidate_blocks,
                    candidate_mapping_error,
                ) = projected_natural_residual(
                    model,
                    candidate,
                    projection_time,
                )
                radiative_evaluations += 2
            except ValueError:
                # The residual RT was evaluated before the invalid trial was
                # detected, but flux diagnostics were not.
                radiative_evaluations += 1
                dt *= 0.25
                rejected += 1
                continue
            projection_evaluations += 1
            candidate_norm = np.linalg.norm(candidate_residual)

            candidate_energy = abs(
                float(
                    np.dot(
                        model.atmospheric_heat_capacity,
                        candidate[:-1],
                    )
                )
            )
            candidate_energy_error = max(
                candidate_projection_error,
                candidate_mapping_error,
            )
            max_energy_error = max(max_energy_error, candidate_energy_error)
            max_relative_energy_error = max(
                max_relative_energy_error,
                candidate_energy_error
                / max(candidate_energy, np.finfo(float).tiny),
            )

            # The natural residual is only piecewise smooth. Allow a temporary
            # increase when the convective active set changes; PTC damping still
            # limits invalid or very large temperature steps.
            active_set_changed = candidate_blocks != mixed_blocks
            acceptance_factor = 10.0 if active_set_changed else 1.25
            if candidate_norm > acceptance_factor * residual_norm:
                dt *= 0.5
                rejected += 1
                continue

            accepted = True
            break

        if not accepted:
            raise RuntimeError(
                "Projected PTC failed to find an acceptable pseudo-step"
            )

        old_norm = residual_norm
        state = candidate
        residual = candidate_residual
        diagnostics = candidate_diagnostics
        mixed_blocks = candidate_blocks
        residual_norm = candidate_norm

        if residual_norm > 0.0:
            ratio = old_norm / residual_norm
            dt_factor = np.clip(growth * ratio, 0.5, 5.0)
        else:
            dt_factor = 5.0
        dt = min(dt_max, dt * dt_factor)

    if verbose:
        print(
            f"Projected PTC did not converge in {max_steps} steps; "
            f"max flux imbalance is "
            f"{np.max(np.abs(diagnostics['imbalance'])):.6e} W/m2."
        )
    result = PTCResult(
        state,
        diagnostics,
        max_steps,
        njev,
        rejected,
        False,
    )
    return ProjectedSolveResult(
        result,
        projection_evaluations,
        max_energy_error,
        max_relative_energy_error,
        radiative_evaluations,
    )


def solve_convective_continuation(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    superadiabatic_tolerance: float = 1.0e-3,
    k_initial: float | None = None,
    k_max: float = 1.0e7,
    max_stages: int = 12,
    safety_factor: float = 1.2,
    minimum_growth: float = 1.5,
    maximum_growth: float = 100.0,
    verbose: bool = True,
) -> ContinuationResult:
    """Increase convective diffusivity until the lapse-rate target is met.

    In an actively convecting region, the residual superadiabaticity scales
    approximately as 1/K. The measured lapse-rate error therefore supplies a
    useful update for K, while each stage is warm-started from the preceding
    equilibrium.
    """

    if superadiabatic_tolerance <= 0.0:
        raise ValueError("superadiabatic_tolerance must be positive")
    if k_initial is None:
        k_initial = model.k_conv
    if not 0.0 < k_initial <= k_max:
        raise ValueError("Require 0 < k_initial <= k_max")
    if max_stages < 1:
        raise ValueError("max_stages must be positive")
    if not 1.0 < minimum_growth <= maximum_growth:
        raise ValueError("Require 1 < minimum_growth <= maximum_growth")
    if safety_factor <= 1.0:
        raise ValueError("safety_factor must exceed one")

    temperature = np.asarray(initial_temperature, dtype=float).copy()
    k_conv = float(k_initial)
    total_steps = 0
    total_njev = 0
    total_rejected = 0

    for stage in range(1, max_stages + 1):
        if verbose:
            print(f"\nK-continuation stage {stage}: k_conv = {k_conv:.6g} m2/s")

        result = solve_ptc(
            model,
            temperature,
            k_conv,
            verbose=verbose,
        )
        total_steps += result.steps
        total_njev += result.jacobian_evaluations
        total_rejected += result.rejected_steps

        if not result.converged:
            return ContinuationResult(
                result,
                k_conv,
                stage,
                total_steps,
                total_njev,
                total_rejected,
                False,
            )

        max_superadiabaticity = float(
            max(0.0, np.max(result.diagnostics["superadiabaticity"]))
        )
        if verbose:
            print(
                "  continuation target: max superadiabaticity "
                f"{max_superadiabaticity:.6e} <= "
                f"{superadiabatic_tolerance:.6e}"
            )

        if max_superadiabaticity <= superadiabatic_tolerance:
            return ContinuationResult(
                result,
                k_conv,
                stage,
                total_steps,
                total_njev,
                total_rejected,
                True,
            )

        growth = np.clip(
            safety_factor
            * max_superadiabaticity
            / superadiabatic_tolerance,
            minimum_growth,
            maximum_growth,
        )
        new_k_conv = min(k_max, k_conv * float(growth))
        if new_k_conv <= k_conv:
            break

        temperature = result.temperature
        k_conv = new_k_conv

    return ContinuationResult(
        result,
        k_conv,
        stage,
        total_steps,
        total_njev,
        total_rejected,
        False,
    )


def solve_coupled_continuation(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    superadiabatic_tolerance: float = 1.0e-3,
    k_initial: float | None = None,
    k_max: float = 1.0e7,
    flux_tolerance: float = 1.0e-5,
    max_steps: int = 400,
    dt_initial: float = 100.0,
    dt_max: float = 1.0e12,
    ptc_growth: float = 1.5,
    residual_reduction_for_k_update: float = 0.1,
    minimum_steps_between_k_updates: int = 3,
    k_safety_factor: float = 1.2,
    k_minimum_growth: float = 1.5,
    k_maximum_growth: float = 3.0,
    k_dt_reduction_power: float = 1.0,
    verbose: bool = True,
) -> ContinuationResult:
    """Continue pseudo-time and convective diffusivity in one PTC solve.

    K is held fixed during every attempted PTC step and its retries. It is only
    increased after enough accepted progress at the current K. When K grows by a
    factor q, the pseudo-timestep is divided by q so that the dimensionless
    strength of the fastest diffusion-like modes does not jump abruptly.
    """

    if superadiabatic_tolerance <= 0.0:
        raise ValueError("superadiabatic_tolerance must be positive")
    if k_initial is None:
        k_initial = model.k_conv
    if not 0.0 < k_initial <= k_max:
        raise ValueError("Require 0 < k_initial <= k_max")
    if not 0.0 < residual_reduction_for_k_update < 1.0:
        raise ValueError("residual_reduction_for_k_update must be between 0 and 1")
    if minimum_steps_between_k_updates < 1:
        raise ValueError("minimum_steps_between_k_updates must be positive")
    if not 1.0 < k_minimum_growth <= k_maximum_growth:
        raise ValueError("Require 1 < k_minimum_growth <= k_maximum_growth")
    if k_dt_reduction_power < 1.0:
        raise ValueError("k_dt_reduction_power must be at least one")

    temperature = np.asarray(initial_temperature, dtype=float).copy()
    if temperature.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    k_conv = float(k_initial)
    dt = dt_initial
    rejected = 0
    njev = 0
    k_updates = 0
    steps_since_k_update = 0

    diagnostics = evaluate(model, temperature, k_conv)
    tendency_norm = np.linalg.norm(diagnostics["tendency"])
    reference_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))

    if verbose:
        print(
            f"{'step':>5s} {'dt [s]':>12s} {'k_conv [m2/s]':>15s} "
            f"{'max |dF| [W/m2]':>18s} {'max superad.':>14s}"
        )

    for step_number in range(max_steps + 1):
        max_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))
        max_superadiabaticity = float(
            max(0.0, np.max(diagnostics["superadiabaticity"]))
        )

        if verbose and (step_number < 8 or step_number % 10 == 0):
            print(
                f"{step_number:5d} {dt:12.4e} {k_conv:15.6g} "
                f"{max_flux_error:18.6e} {max_superadiabaticity:14.6e}"
            )

        if (
            max_flux_error < flux_tolerance
            and max_superadiabaticity <= superadiabatic_tolerance
        ):
            result = PTCResult(
                temperature,
                diagnostics,
                step_number,
                njev,
                rejected,
                True,
            )
            return ContinuationResult(
                result,
                k_conv,
                k_updates + 1,
                step_number,
                njev,
                rejected,
                True,
            )

        # Change K only at a fixed accepted state, never during a rejected-step
        # retry. Requiring residual reduction lets the profile respond to the
        # current K before its remaining superadiabaticity is used as an estimate.
        ready_for_k_update = (
            max_superadiabaticity > superadiabatic_tolerance
            and k_conv < k_max
            and steps_since_k_update >= minimum_steps_between_k_updates
            and (
                max_flux_error
                <= residual_reduction_for_k_update * reference_flux_error
                or max_flux_error < flux_tolerance
            )
        )
        if ready_for_k_update:
            estimated_growth = np.sqrt(
                k_safety_factor
                * max_superadiabaticity
                / superadiabatic_tolerance
            )
            k_growth = float(
                np.clip(
                    estimated_growth,
                    k_minimum_growth,
                    k_maximum_growth,
                )
            )
            new_k_conv = min(k_max, k_conv * k_growth)

            if new_k_conv > k_conv:
                actual_growth = new_k_conv / k_conv
                k_conv = new_k_conv
                dt = max(
                    np.finfo(float).tiny,
                    dt / actual_growth**k_dt_reduction_power,
                )
                k_updates += 1
                steps_since_k_update = 0

                # K changes the residual at fixed temperature. Re-evaluate it and
                # restart the progress reference before forming the next Jacobian.
                diagnostics = evaluate(model, temperature, k_conv)
                tendency_norm = np.linalg.norm(diagnostics["tendency"])
                reference_flux_error = float(
                    np.max(np.abs(diagnostics["imbalance"]))
                )

                if verbose:
                    print(
                        f"      increased k_conv by {actual_growth:.3g} to "
                        f"{k_conv:.6g} m2/s; reduced dt to {dt:.4e} s"
                    )

        if step_number == max_steps:
            break

        jac = finite_difference_jacobian(model, temperature, k_conv)
        njev += 1

        accepted = False
        for _ in range(15):
            matrix = np.eye(model.nlev + 1) / dt - jac
            try:
                delta_temperature = np.linalg.solve(
                    matrix, diagnostics["tendency"]
                )
            except np.linalg.LinAlgError:
                dt *= 0.25
                rejected += 1
                continue

            candidate = temperature + delta_temperature
            if (
                np.any(~np.isfinite(candidate))
                or np.min(candidate) < 50.0
                or np.max(candidate) > 1000.0
            ):
                dt *= 0.25
                rejected += 1
                continue

            candidate_diagnostics = evaluate(model, candidate, k_conv)
            candidate_norm = np.linalg.norm(candidate_diagnostics["tendency"])
            if candidate_norm > 1.25 * tendency_norm:
                dt *= 0.5
                rejected += 1
                continue

            accepted = True
            break

        if not accepted:
            raise RuntimeError("Coupled PTC failed to find an acceptable pseudo-step")

        old_norm = tendency_norm
        temperature = candidate
        diagnostics = candidate_diagnostics
        tendency_norm = candidate_norm
        steps_since_k_update += 1

        if tendency_norm > 0.0:
            ratio = old_norm / tendency_norm
            dt_factor = np.clip(ptc_growth * ratio, 0.5, 5.0)
        else:
            dt_factor = 5.0
        dt = min(dt_max, dt * dt_factor)

    result = PTCResult(
        temperature,
        diagnostics,
        max_steps,
        njev,
        rejected,
        False,
    )
    return ContinuationResult(
        result,
        k_conv,
        k_updates + 1,
        max_steps,
        njev,
        rejected,
        False,
    )


def solve_hybrid_continuation(
    model: Model,
    initial_temperature: np.ndarray,
    *,
    superadiabatic_tolerance: float = 1.0e-3,
    superadiabatic_guard: float | None = None,
    k_initial: float | None = None,
    k_max: float = 1.0e7,
    flux_tolerance: float = 1.0e-5,
    max_steps: int = 400,
    dt_initial: float = 100.0,
    dt_max: float = 1.0e12,
    ptc_growth: float = 1.5,
    hard_guard_factor: float = 2.0,
    minimum_steps_between_k_updates: int = 1,
    k_safety_factor: float = 1.2,
    transient_k_minimum_growth: float = 1.2,
    transient_k_maximum_growth: float = 3.0,
    equilibrium_k_minimum_growth: float = 1.5,
    equilibrium_k_maximum_growth: float = 20.0,
    k_dt_reduction_power: float = 1.5,
    verbose: bool = True,
) -> ContinuationResult:
    """PTC with equilibrium- and lapse-guard-triggered continuation in K.

    The equilibrium trigger increases K when the flux residual has converged but
    the requested lapse-rate tolerance has not. The transient trigger increases K
    before the evolving profile becomes excessively superadiabatic. Candidate
    states beyond a harder lapse-rate guard are rejected before evaluating RT.
    """

    if superadiabatic_tolerance <= 0.0:
        raise ValueError("superadiabatic_tolerance must be positive")
    if superadiabatic_guard is None:
        superadiabatic_guard = max(
            10.0 * superadiabatic_tolerance,
            0.05 * model.nabla_ad,
        )
    if superadiabatic_guard <= superadiabatic_tolerance:
        raise ValueError(
            "superadiabatic_guard must exceed superadiabatic_tolerance"
        )
    if hard_guard_factor <= 1.0:
        raise ValueError("hard_guard_factor must exceed one")
    if k_initial is None:
        k_initial = model.k_conv
    if not 0.0 < k_initial <= k_max:
        raise ValueError("Require 0 < k_initial <= k_max")
    if minimum_steps_between_k_updates < 1:
        raise ValueError("minimum_steps_between_k_updates must be positive")
    if not 1.0 < transient_k_minimum_growth <= transient_k_maximum_growth:
        raise ValueError("Invalid transient K growth limits")
    if not 1.0 < equilibrium_k_minimum_growth <= equilibrium_k_maximum_growth:
        raise ValueError("Invalid equilibrium K growth limits")
    if k_dt_reduction_power < 1.0:
        raise ValueError("k_dt_reduction_power must be at least one")

    temperature = np.asarray(initial_temperature, dtype=float).copy()
    if temperature.shape != (model.nlev + 1,):
        raise ValueError(f"Expected temperature shape {(model.nlev + 1,)}")

    k_conv = float(k_initial)
    dt = dt_initial
    hard_guard = hard_guard_factor * superadiabatic_guard
    rejected = 0
    njev = 0
    k_updates = 0
    equilibrium_k_updates = 0
    guard_k_updates = 0
    steps_since_k_update = minimum_steps_between_k_updates

    diagnostics = evaluate(model, temperature, k_conv)
    tendency_norm = np.linalg.norm(diagnostics["tendency"])
    max_superadiabaticity_encountered = float(
        max(0.0, np.max(diagnostics["superadiabaticity"]))
    )

    if verbose:
        print(
            f"  target superadiabaticity: {superadiabatic_tolerance:.6e}\n"
            f"  transient guard:          {superadiabatic_guard:.6e}\n"
            f"  hard candidate guard:     {hard_guard:.6e}"
        )
        print(
            f"{'step':>5s} {'dt [s]':>12s} {'k_conv [m2/s]':>15s} "
            f"{'max |dF| [W/m2]':>18s} {'max superad.':>14s}"
        )

    for step_number in range(max_steps + 1):
        max_flux_error = float(np.max(np.abs(diagnostics["imbalance"])))
        max_superadiabaticity = float(
            max(0.0, np.max(diagnostics["superadiabaticity"]))
        )

        if verbose and (step_number < 8 or step_number % 10 == 0):
            print(
                f"{step_number:5d} {dt:12.4e} {k_conv:15.6g} "
                f"{max_flux_error:18.6e} {max_superadiabaticity:14.6e}"
            )

        flux_converged = max_flux_error < flux_tolerance
        lapse_converged = (
            max_superadiabaticity <= superadiabatic_tolerance
        )
        if flux_converged and lapse_converged:
            if verbose:
                print(
                    f"Converged in {step_number} accepted PTC steps "
                    f"({rejected} rejected, {guard_k_updates} guard-triggered "
                    f"and {equilibrium_k_updates} equilibrium-triggered "
                    "K updates)."
                )
            result = PTCResult(
                temperature,
                diagnostics,
                step_number,
                njev,
                rejected,
                True,
            )
            return ContinuationResult(
                result,
                k_conv,
                k_updates + 1,
                step_number,
                njev,
                rejected,
                True,
                max_superadiabaticity_encountered,
                guard_k_updates,
                equilibrium_k_updates,
            )

        update_reason: str | None = None
        k_growth = 1.0
        if (
            flux_converged
            and not lapse_converged
            and steps_since_k_update >= minimum_steps_between_k_updates
        ):
            update_reason = "equilibrium"
            k_growth = float(
                np.clip(
                    k_safety_factor
                    * max_superadiabaticity
                    / superadiabatic_tolerance,
                    equilibrium_k_minimum_growth,
                    equilibrium_k_maximum_growth,
                )
            )
        elif (
            max_superadiabaticity > superadiabatic_guard
            and steps_since_k_update >= minimum_steps_between_k_updates
        ):
            update_reason = "guard"
            k_growth = float(
                np.clip(
                    np.sqrt(
                        k_safety_factor
                        * max_superadiabaticity
                        / superadiabatic_guard
                    ),
                    transient_k_minimum_growth,
                    transient_k_maximum_growth,
                )
            )

        if update_reason is not None:
            new_k_conv = min(k_max, k_conv * k_growth)
            if new_k_conv <= k_conv:
                if flux_converged:
                    break
            else:
                actual_growth = new_k_conv / k_conv
                k_conv = new_k_conv
                dt = max(
                    np.finfo(float).tiny,
                    dt / actual_growth**k_dt_reduction_power,
                )
                k_updates += 1
                steps_since_k_update = 0
                if update_reason == "guard":
                    guard_k_updates += 1
                else:
                    equilibrium_k_updates += 1

                # Changing K changes the flux residual but not the temperature or
                # lapse rate. Reevaluate before constructing the new Jacobian.
                diagnostics = evaluate(model, temperature, k_conv)
                tendency_norm = np.linalg.norm(diagnostics["tendency"])

                if verbose:
                    print(
                        f"      {update_reason}-triggered K update: "
                        f"x{actual_growth:.3g} -> {k_conv:.6g} m2/s; "
                        f"dt -> {dt:.4e} s"
                    )

        if step_number == max_steps:
            break

        jac = finite_difference_jacobian(model, temperature, k_conv)
        njev += 1

        accepted = False
        for _ in range(15):
            matrix = np.eye(model.nlev + 1) / dt - jac
            try:
                delta_temperature = np.linalg.solve(
                    matrix, diagnostics["tendency"]
                )
            except np.linalg.LinAlgError:
                dt *= 0.25
                rejected += 1
                continue

            candidate = temperature + delta_temperature
            if (
                np.any(~np.isfinite(candidate))
                or np.min(candidate) < 50.0
                or np.max(candidate) > 1000.0
            ):
                dt *= 0.25
                rejected += 1
                continue

            # Lapse rate is local and cheap. Enforce the hard guard before an
            # expensive production RT residual would be evaluated.
            candidate_lapse = (
                np.diff(np.log(candidate[:-1]))
                / np.diff(np.log(model.p_cell))
                - model.nabla_ad
            )
            if float(max(0.0, np.max(candidate_lapse))) > hard_guard:
                dt *= 0.5
                rejected += 1
                continue

            candidate_diagnostics = evaluate(model, candidate, k_conv)
            candidate_norm = np.linalg.norm(candidate_diagnostics["tendency"])
            if candidate_norm > 1.25 * tendency_norm:
                dt *= 0.5
                rejected += 1
                continue

            accepted = True
            break

        if not accepted:
            raise RuntimeError("Hybrid PTC failed to find an acceptable pseudo-step")

        old_norm = tendency_norm
        temperature = candidate
        diagnostics = candidate_diagnostics
        tendency_norm = candidate_norm
        max_superadiabaticity_encountered = max(
            max_superadiabaticity_encountered,
            float(max(0.0, np.max(diagnostics["superadiabaticity"]))),
        )
        steps_since_k_update += 1

        if tendency_norm > 0.0:
            ratio = old_norm / tendency_norm
            dt_factor = np.clip(ptc_growth * ratio, 0.5, 5.0)
        else:
            dt_factor = 5.0
        dt = min(dt_max, dt * dt_factor)

    if verbose:
        print(
            f"Hybrid continuation did not converge in {max_steps} steps; "
            f"max flux imbalance is "
            f"{np.max(np.abs(diagnostics['imbalance'])):.6e} W/m2."
        )
    result = PTCResult(
        temperature,
        diagnostics,
        min(step_number, max_steps),
        njev,
        rejected,
        False,
    )
    return ContinuationResult(
        result,
        k_conv,
        k_updates + 1,
        min(step_number, max_steps),
        njev,
        rejected,
        False,
        max_superadiabaticity_encountered,
        guard_k_updates,
        equilibrium_k_updates,
    )


def initial_profile(model: Model) -> np.ndarray:
    """A deliberately simple atmospheric profile and surface temperature."""

    # Isothermal upper atmosphere transitioning smoothly to the fixed surface.
    pressure_fraction = np.log(model.p_cell / model.p_top) / np.log(
        model.p_surf / model.p_top
    )
    atmosphere = 210.0 + (model.t_surf - 210.0) * pressure_fraction**0.65
    return np.concatenate((atmosphere, [model.t_surf]))


def print_summary(model: Model, label: str, result: PTCResult) -> None:
    """Print compact equilibrium diagnostics."""

    diag = result.diagnostics
    print(f"\n{label}")
    print(f"  converged:                  {result.converged}")
    print(f"  accepted solver steps:      {result.steps}")
    print(f"  rejected attempts:          {result.rejected_steps}")
    print(f"  Jacobian evaluations:       {result.jacobian_evaluations}")
    print(
        "  max flux imbalance:        "
        f"{np.max(np.abs(diag['imbalance'])):.6e} W/m2"
    )
    print(f"  OLR:                        {diag['up_lw'][0]:.6f} W/m2")
    print(f"  net TOA flux:               {diag['total'][0]:.6f} W/m2")
    print(
        "  max superadiabaticity:     "
        f"{max(0.0, np.max(diag['superadiabaticity'])):.6e}"
    )
    print(
        "  max upward convective flux:"
        f" {np.max(diag['conv']):.6f} W/m2"
    )
    print(
        "  temperature range:         "
        f"{np.min(result.temperature[:-1]):.3f}--"
        f"{np.max(result.temperature[:-1]):.3f} K"
    )
    print(f"  surface temperature:       {result.temperature[-1]:.3f} K")


def print_continuation_summary(continuation: ContinuationResult) -> None:
    """Print aggregate diagnostics for the outer K continuation."""

    print("\nConvective-diffusivity continuation")
    print(f"  converged:                  {continuation.converged}")
    print(f"  continuation stages:        {continuation.stages}")
    print(f"  selected k_conv:             {continuation.k_conv:.6g} m2/s")
    print(f"  total accepted PTC steps:    {continuation.total_steps}")
    print(f"  total rejected attempts:     {continuation.total_rejected_steps}")
    print(
        "  total Jacobian evaluations: "
        f"{continuation.total_jacobian_evaluations}"
    )
    if (
        continuation.guard_triggered_updates
        + continuation.equilibrium_triggered_updates
        > 0
    ):
        print(
            "  guard-triggered K updates: "
            f"{continuation.guard_triggered_updates}"
        )
        print(
            "  equilibrium K updates:     "
            f"{continuation.equilibrium_triggered_updates}"
        )
        print(
            "  largest accepted superad.: "
            f"{continuation.max_superadiabaticity_encountered:.6e}"
        )


def print_projected_summary(
    projected: ProjectedSolveResult,
    method: str,
) -> None:
    """Print work and conservation diagnostics for a projected solve."""

    print(f"\nProjected dry-convective solve ({method})")
    print(f"  converged:                  {projected.result.converged}")
    print(f"  accepted integration steps: {projected.result.steps}")
    print(f"  rejected attempts:          {projected.result.rejected_steps}")
    print(
        "  projected Jacobian evals:   "
        f"{projected.result.jacobian_evaluations}"
    )
    print(
        "  full radiative evaluations:"
        f" {projected.radiative_evaluations}"
    )
    print(
        "  projection evaluations:     "
        f"{projected.projection_evaluations}"
    )
    print(
        "  max projection energy error:"
        f" {projected.max_projection_energy_error:.6e} J/m2"
    )
    print(
        "  max relative energy error:  "
        f"{projected.max_projection_relative_energy_error:.6e}"
    )
    if projected.nonlinear_iterations > 0:
        print(
            "  nonlinear iterations:      "
            f"{projected.nonlinear_iterations}"
        )
        print(
            "  simulated physical time:   "
            f"{projected.simulated_time:.6e} s"
        )
        print(
            "  final physical timestep:   "
            f"{projected.final_timestep:.6e} s"
        )
        print(
            "  max estimated local error: "
            f"{projected.max_estimated_local_error:.6e} K"
        )


def make_plot(
    model: Model,
    initial: np.ndarray,
    radiative_result: PTCResult | None,
    convective_result: PTCResult,
    superadiabatic_tolerance: float,
    superadiabatic_guard: float | None,
    rce_label: str,
    output: Path,
) -> None:
    """Plot profiles and equilibrium flux diagnostics."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Matplotlib is required for plotting; use --no-plot instead"
        ) from exc

    p_bar = model.p_cell / 1.0e5
    p_edge_bar = model.p_edge / 1.0e5

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.2))

    ax = axes[0]
    ax.plot(initial[:-1], p_bar, color="0.65", linestyle="--", label="initial")
    if radiative_result is not None:
        ax.plot(
            radiative_result.temperature[:-1],
            p_bar,
            color="#3366aa",
            label="radiative only",
        )
    ax.plot(
        convective_result.temperature[:-1],
        p_bar,
        color="#bb3322",
        label=rce_label,
    )
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.grid(alpha=0.25)
    ax.legend()
    if radiative_result is not None:
        ax.scatter(
            [radiative_result.temperature[-1]],
            [model.p_surf / 1.0e5],
            color=["#3366aa"],
            marker="s",
            s=25,
            zorder=3,
        )
    ax.scatter(
        [convective_result.temperature[-1]],
        [model.p_surf / 1.0e5],
        color=["#bb3322"],
        marker="s",
        s=25,
        zorder=3,
    )

    ax = axes[1]
    diag = convective_result.diagnostics
    ax.plot(diag["net_lw"], p_edge_bar, label="net LW")
    ax.plot(diag["net_sw"], p_edge_bar, label="net SW")
    ax.plot(diag["net_rad"], p_edge_bar, label="net radiation")
    ax.plot(diag["conv"], p_edge_bar, label="convection")
    ax.plot(diag["total"], p_edge_bar, color="black", linewidth=2, label="total")
    ax.axvline(0.0, color="0.7", linewidth=0.8)
    ax.set_xlabel("Upward flux [W m$^{-2}$]")
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[2]
    p_int_bar = np.sqrt(model.p_cell[:-1] * model.p_cell[1:]) / 1.0e5
    if radiative_result is not None:
        ax.plot(
            radiative_result.diagnostics["superadiabaticity"],
            p_int_bar,
            color="#3366aa",
            label="radiative only",
        )
    ax.plot(
        convective_result.diagnostics["superadiabaticity"],
        p_int_bar,
        color="#bb3322",
        label=rce_label,
    )
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.axvline(
        superadiabatic_tolerance,
        color="#dd8800",
        linestyle="--",
        linewidth=1.3,
        label=r"final target $s_{\rm target}$",
    )
    if superadiabatic_guard is not None:
        ax.axvline(
            superadiabatic_guard,
            color="#8844aa",
            linestyle=":",
            linewidth=1.6,
            label=r"transient guard $s_{\rm guard}$",
        )
    ax.set_xlabel(r"$\nabla-\nabla_{\rm ad}$")
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    fig.suptitle("Toy gray radiative-convective equilibrium")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    print(f"\nSaved plot to {output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-plot", action="store_true", help="run without importing Matplotlib"
    )
    parser.add_argument(
        "--time-benchmark",
        action="store_true",
        help=(
            "compare explicit and constrained-implicit adjustment at one "
            "physical end time, then exit"
        ),
    )
    parser.add_argument(
        "--benchmark-end-time",
        type=float,
        default=1.0e8,
        help="physical end time for --time-benchmark, in seconds (default: 1e8)",
    )
    parser.add_argument(
        "--benchmark-reference-dt",
        type=float,
        default=2.5e4,
        help=(
            "coarser of two explicit reference timesteps for --time-benchmark, "
            "in seconds (default: 2.5e4)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("smooth_rce_toy.png"),
        help="output plot path (default: smooth_rce_toy.png)",
    )
    parser.add_argument(
        "--nlev", type=int, default=40, help="number of atmospheric layers"
    )
    parser.add_argument(
        "--k-conv",
        "--k-conv-initial",
        dest="k_conv_initial",
        type=float,
        default=10.0,
        help="initial convective eddy diffusivity in m2/s (default: 10)",
    )
    parser.add_argument(
        "--superadiabatic-tolerance",
        type=float,
        default=1.0e-3,
        help="maximum allowed nabla - nabla_ad (default: 1e-3)",
    )
    parser.add_argument(
        "--superadiabatic-guard",
        type=float,
        default=None,
        help=(
            "transient lapse-rate guard for hybrid continuation; default is "
            "max(10*tolerance, 0.05*nabla_ad)"
        ),
    )
    parser.add_argument(
        "--k-conv-max",
        type=float,
        default=1.0e7,
        help="maximum diffusivity allowed during continuation (default: 1e7)",
    )
    parser.add_argument(
        "--continuation-mode",
        choices=(
            "staged",
            "coupled",
            "hybrid",
            "projected",
            "projected-explicit",
            "projected-implicit",
        ),
        default="hybrid",
        help=(
            "choose staged, coupled, or hybrid smooth-K continuation, or exact "
            "dry adjustment via projected PTC, explicit integration, or "
            "constrained backward Euler "
            "(default: hybrid)"
        ),
    )
    parser.add_argument(
        "--projection-time",
        type=float,
        default=1.0e4,
        help=(
            "scaling time in the projected natural residual, in seconds "
            "(default: 1e4)"
        ),
    )
    parser.add_argument(
        "--explicit-dt-max",
        type=float,
        default=1.0e6,
        help=(
            "maximum timestep for projected-explicit integration, in seconds "
            "(default: 1e6)"
        ),
    )
    parser.add_argument(
        "--explicit-max-temperature-step",
        type=float,
        default=1.0,
        help=(
            "maximum adjusted temperature change per explicit step, in K "
            "(default: 1)"
        ),
    )
    parser.add_argument(
        "--implicit-dt-initial",
        type=float,
        default=100.0,
        help=(
            "initial physical timestep for constrained backward Euler, in "
            "seconds (default: 100)"
        ),
    )
    parser.add_argument(
        "--implicit-dt-max",
        type=float,
        default=1.0e12,
        help=(
            "maximum physical timestep for constrained backward Euler, in "
            "seconds (default: 1e12)"
        ),
    )
    parser.add_argument(
        "--implicit-temperature-tolerance",
        type=float,
        default=2.0e-2,
        help=(
            "step-doubling local temperature error tolerance for constrained "
            "backward Euler, in K (default: 2e-2)"
        ),
    )
    parser.add_argument(
        "--implicit-flux-tolerance",
        type=float,
        default=1.0e-5,
        help=(
            "steady maximum flux-imbalance tolerance for constrained backward "
            "Euler, in W/m2 (default: 1e-5)"
        ),
    )
    parser.add_argument(
        "--implicit-max-steps",
        type=int,
        default=1000,
        help=(
            "maximum accepted constrained backward-Euler steps "
            "(default: 1000)"
        ),
    )
    parser.add_argument(
        "--quiet", action="store_true", help="suppress per-step solver output"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    model = Model(nlev=args.nlev, k_conv=args.k_conv_initial)
    initial = initial_profile(model)

    if args.time_benchmark:
        run_time_integration_benchmark(
            model,
            initial,
            end_time=args.benchmark_end_time,
            reference_dt=args.benchmark_reference_dt,
        )
        return 0

    radiative_result: PTCResult | None = None
    if args.continuation_mode in ("staged", "coupled"):
        print("Solving radiative equilibrium...")
        radiative_result = solve_ptc(
            model, initial, 0.0, verbose=not args.quiet
        )
        print_summary(model, "Radiative equilibrium", radiative_result)
    else:
        print(
            "Skipping a separate radiative-equilibrium solve; the "
            f"{args.continuation_mode} solver starts from the initial profile."
        )

    projected_result: ProjectedSolveResult | None = None
    continuation: ContinuationResult | None = None
    if args.continuation_mode == "projected-explicit":
        print(
            "\nSolving radiative-convective equilibrium with explicit "
            "radiative steps and convective adjustment..."
        )
        projected_result = solve_projected_explicit(
            model,
            initial,
            dt_max=args.explicit_dt_max,
            max_temperature_step=args.explicit_max_temperature_step,
            verbose=not args.quiet,
        )
        convective_result = projected_result.result
        print_projected_summary(projected_result, "explicit")
        print_summary(
            model,
            "Explicit adjusted radiative-convective equilibrium",
            convective_result,
        )
    elif args.continuation_mode == "projected-implicit":
        print(
            "\nSolving radiative-convective evolution with constrained "
            "backward Euler..."
        )
        projected_result = solve_projected_implicit(
            model,
            initial,
            flux_tolerance=args.implicit_flux_tolerance,
            max_steps=args.implicit_max_steps,
            dt_initial=args.implicit_dt_initial,
            dt_max=args.implicit_dt_max,
            temperature_tolerance=args.implicit_temperature_tolerance,
            verbose=not args.quiet,
        )
        convective_result = projected_result.result
        print_projected_summary(projected_result, "constrained backward Euler")
        print_summary(
            model,
            "Implicit adjusted radiative-convective equilibrium",
            convective_result,
        )
    elif args.continuation_mode == "projected":
        print(
            "\nSolving radiative-convective equilibrium with projected PTC..."
        )
        projected_result = solve_projected_ptc(
            model,
            initial,
            projection_time=args.projection_time,
            verbose=not args.quiet,
        )
        convective_result = projected_result.result
        print_projected_summary(projected_result, "PTC")
        print_summary(
            model,
            "Projected radiative-convective equilibrium",
            convective_result,
        )
    elif args.continuation_mode == "hybrid":
        print(
            "\nSolving smooth radiative-convective equilibrium by "
            "hybrid K continuation..."
        )
        # The hybrid guard is intended to avoid first converging a potentially
        # extreme radiative equilibrium, so start it from the original profile.
        continuation = solve_hybrid_continuation(
            model,
            initial,
            superadiabatic_tolerance=args.superadiabatic_tolerance,
            superadiabatic_guard=args.superadiabatic_guard,
            k_initial=args.k_conv_initial,
            k_max=args.k_conv_max,
            verbose=not args.quiet,
        )
        convective_result = continuation.result
        print_continuation_summary(continuation)
        print_summary(
            model,
            "Smooth radiative-convective equilibrium",
            convective_result,
        )
    else:
        print(
            "\nSolving smooth radiative-convective equilibrium by "
            f"{args.continuation_mode} K continuation..."
        )
        continuation_solver = (
            solve_convective_continuation
            if args.continuation_mode == "staged"
            else solve_coupled_continuation
        )
        continuation = continuation_solver(
            model,
            radiative_result.temperature,
            superadiabatic_tolerance=args.superadiabatic_tolerance,
            k_initial=args.k_conv_initial,
            k_max=args.k_conv_max,
            verbose=not args.quiet,
        )
        convective_result = continuation.result
        print_continuation_summary(continuation)
        print_summary(
            model,
            "Smooth radiative-convective equilibrium",
            convective_result,
        )

    solve_converged = (
        projected_result.result.converged
        if projected_result is not None
        else continuation.converged
    )
    if not solve_converged:
        return 1
    if radiative_result is not None and not radiative_result.converged:
        return 1

    if not args.no_plot:
        if radiative_result is None:
            print("\nSolving radiative equilibrium for plotting only...")
            radiative_result = solve_ptc(
                model,
                initial,
                0.0,
                verbose=not args.quiet,
            )
            print_summary(model, "Radiative equilibrium", radiative_result)
            if not radiative_result.converged:
                return 1

        plot_guard = None
        if args.continuation_mode == "hybrid":
            plot_guard = args.superadiabatic_guard
            if plot_guard is None:
                plot_guard = max(
                    10.0 * args.superadiabatic_tolerance,
                    0.05 * model.nabla_ad,
                )
        make_plot(
            model,
            initial,
            radiative_result,
            convective_result,
            args.superadiabatic_tolerance,
            plot_guard,
            (
                "projected RCE"
                if args.continuation_mode
                in ("projected", "projected-explicit", "projected-implicit")
                else "smooth RCE"
            ),
            args.output,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
