#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin

"""FVQual against measured DeGray temperature profiles.

CE-QUAL-W2 ships `prf_dam.npt` with the DeGray example, and constituent 22 in
it is measured temperature at the dam: 24 casts through 1980, roughly 50 points
each at 1 m resolution. That is field data for the exact dataset the model was
run on, and it is a stronger reference than W2's own output would be — one
tells you whether FVQual reproduces another model, the other whether it
reproduces the reservoir.

Reports bias and RMSE by date and by depth band, and draws the comparison.

Usage: python3 compare_observed.py [results-dir] [prf_dam.npt]

NOTE: written while the build sandbox was unavailable and therefore never
executed. Run it before trusting anything it prints.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_degray import style, months, INK, MUTED

#! Constituent number carrying temperature in the DeGray prf file.
TEMPERATURE = 22

#! Mixed-layer criterion: first depth colder than the surface by this much.
MIXED_LAYER_DROP = 1.0


def read_prf(path, constituent=TEMPERATURE):
    """Parses W2's observed-profile format into {day: (depths, values)}.

    The format is a run of blocks, each a `CNO DAYO NRO` line followed by a
    `ZO C2O` header and then NRO depth/value pairs, five pairs to a line and
    wrapping onto as many lines as needed.
    """
    with open(path) as handle:
        lines = handle.readlines()

    profiles = {}
    i = 0

    while i < len(lines):
        if "CNO" not in lines[i] or "DAYO" not in lines[i]:
            i += 1
            continue

        # The values line follows, after any blank line.
        i += 1

        while i < len(lines) and not lines[i].strip():
            i += 1

        fields = lines[i].split()

        if len(fields) < 3:
            i += 1
            continue

        cno, day, count = int(fields[0]), float(fields[1]), int(fields[2])

        # Skip to past the ZO/C2O header.
        i += 1

        while i < len(lines) and "ZO" not in lines[i]:
            i += 1

        i += 1

        numbers = []

        while i < len(lines) and len(numbers) < 2 * count:
            stripped = lines[i].strip()

            if stripped:
                numbers.extend(float(f) for f in stripped.split())

            i += 1

        if cno != constituent or len(numbers) < 2 * count:
            continue

        pairs = np.array(numbers[:2 * count]).reshape(-1, 2)

        # Guard against an unsorted cast.
        order = np.argsort(pairs[:, 0])

        profiles[day] = (pairs[order, 0], pairs[order, 1])

    return profiles


def mixed_layer_depth(depth, temperature, drop=MIXED_LAYER_DROP):
    """Shallowest depth colder than the surface value by `drop`."""
    if len(depth) < 2:
        return np.nan

    cooler = np.where(temperature <= temperature[0] - drop)[0]

    if cooler.size == 0:
        # No crossing: the column is within `drop` of isothermal top to bottom,
        # so it has no mixed layer to report rather than one as deep as the
        # reservoir. Returning the bed depth here made a February profile look
        # like a 54 m mixed layer, which is true and useless.
        return np.nan

    first = cooler[0]

    if first == 0:
        return float(depth[0])

    # Linear interpolation onto the exact crossing.
    target = temperature[0] - drop
    hi, lo = temperature[first - 1], temperature[first]
    span = hi - lo

    weight = (hi - target) / span if span > 0 else 0.0

    return float(depth[first - 1] + weight * (depth[first] - depth[first - 1]))


def model_profile(profiles, day):
    """The modelled profile nearest `day`, surface first."""
    days = profiles["jday"].unique()
    nearest = days[np.abs(days - day).argmin()]

    slab = profiles[profiles.jday == nearest].sort_values("depth_m")

    return nearest, slab.depth_m.values, slab.temperature_C.values


def compare(profiles, observed):
    """Matches every observed cast to a modelled profile."""
    rows = []

    for day in sorted(observed):
        depth, measured = observed[day]

        matched, mdepth, modelled = model_profile(profiles, day)

        # Only compare where the model has water; the casts run slightly
        # deeper than the deepest layer centre.
        inside = (depth >= mdepth.min()) & (depth <= mdepth.max())

        if inside.sum() < 3:
            continue

        interpolated = np.interp(depth[inside], mdepth, modelled)

        error = interpolated - measured[inside]

        rows.append({
            "observed_day": day,
            "model_day": matched,
            "n": int(inside.sum()),
            "bias": float(error.mean()),
            "rmse": float(np.sqrt((error ** 2).mean())),
            "surface_observed": float(measured[0]),
            "surface_model": float(modelled[0]),
            "bottom_observed": float(measured[-1]),
            "bottom_model": float(modelled[-1]),
            "mld_observed": mixed_layer_depth(depth, measured),
            "mld_model": mixed_layer_depth(mdepth, modelled),
        })

    return pd.DataFrame(rows)


def depth_bands(profiles, observed, edges=(0, 5, 10, 20, 30, 60)):
    """Bias and RMSE pooled over every cast, by depth band."""
    collected = {}

    for day in sorted(observed):
        depth, measured = observed[day]
        _, mdepth, modelled = model_profile(profiles, day)

        inside = (depth >= mdepth.min()) & (depth <= mdepth.max())

        if inside.sum() < 3:
            continue

        error = np.interp(depth[inside], mdepth, modelled) - measured[inside]

        for z, e in zip(depth[inside], error):
            for lo, hi in zip(edges[:-1], edges[1:]):
                if lo <= z < hi:
                    collected.setdefault((lo, hi), []).append(e)
                    break

    rows = []

    for (lo, hi), errors in sorted(collected.items()):
        errors = np.array(errors)

        rows.append({
            "band": f"{lo}-{hi} m",
            "n": len(errors),
            "bias": errors.mean(),
            "rmse": np.sqrt((errors ** 2).mean()),
        })

    return pd.DataFrame(rows)


def figure_profiles(profiles, observed, path, count=9):
    """Small multiples: measured against modelled, through the year."""
    days = sorted(observed)
    chosen = [days[i] for i in np.linspace(0, len(days) - 1, count).astype(int)]

    fig, axes = plt.subplots(3, 3, figsize=(10.5, 10.5), sharey=True)

    for ax, day in zip(axes.flat, chosen):
        depth, measured = observed[day]
        matched, mdepth, modelled = model_profile(profiles, day)

        ax.plot(measured, depth, color="#1b1b1b", lw=1.6, label="measured")
        ax.plot(modelled, mdepth, color="#c0392b", lw=1.6, label="FVQual")

        ax.set_title(f"day {day:.0f}", fontsize=9.5)
        ax.set_xlim(4, 36)
        ax.grid(alpha=0.4)

    axes[0, 0].invert_yaxis()
    axes[0, 0].legend(loc="lower left", fontsize=8)

    for ax in axes[-1]:
        ax.set_xlabel("temperature (°C)")

    for ax in axes[:, 0]:
        ax.set_ylabel("depth (m)")

    fig.suptitle("DeGray at the dam — measured against modelled, 1980",
                 fontsize=12, fontweight="semibold")

    fig.savefig(path)
    plt.close(fig)


def figure_bias(profiles, observed, path):
    """Where and when the model is wrong, as a depth-time field."""
    days = sorted(observed)
    grid = np.arange(0.0, 50.0, 1.0)

    field = np.full((len(grid), len(days)), np.nan)

    for j, day in enumerate(days):
        depth, measured = observed[day]
        _, mdepth, modelled = model_profile(profiles, day)

        usable = (grid >= max(depth.min(), mdepth.min())) & \
                 (grid <= min(depth.max(), mdepth.max()))

        field[usable, j] = (np.interp(grid[usable], mdepth, modelled) -
                            np.interp(grid[usable], depth, measured))

    fig, ax = plt.subplots(figsize=(10.5, 4.6))

    limit = np.nanmax(np.abs(field))

    mesh = ax.pcolormesh(np.array(days), grid, field, cmap="RdBu_r",
                         vmin=-limit, vmax=limit, shading="nearest")

    ax.invert_yaxis()
    ax.set_ylabel("depth (m)")
    ax.set_title("Model minus measured (°C) — red is too warm, blue too cold")
    months(ax, min(days), max(days))

    bar = fig.colorbar(mesh, ax=ax, pad=0.015, aspect=28)
    bar.set_label("error (°C)")
    bar.outline.set_visible(False)

    fig.savefig(path)
    plt.close(fig)


def figure_mixed_layer(table, path):
    """The diagnosis in one panel: is the epilimnion the right thickness?"""
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(10.5, 6.0), sharex=True)

    top.plot(table.observed_day, table.surface_observed, "o-", color=INK,
             lw=1.4, ms=4, label="measured")
    top.plot(table.observed_day, table.surface_model, "o-", color="#c0392b",
             lw=1.4, ms=4, label="FVQual")
    top.plot(table.observed_day, table.bottom_observed, "o--", color=INK,
             lw=1.0, ms=3, alpha=0.6)
    top.plot(table.observed_day, table.bottom_model, "o--", color="#1f5f8b",
             lw=1.0, ms=3, alpha=0.8)

    top.set_ylabel("temperature (°C)")
    top.set_title("Surface (solid) and bottom (dashed)")
    top.legend(loc="upper left", fontsize=8.5)

    bottom.plot(table.observed_day, table.mld_observed, "o-", color=INK,
                lw=1.4, ms=4, label="measured")
    bottom.plot(table.observed_day, table.mld_model, "o-", color="#c0392b",
                lw=1.4, ms=4, label="FVQual")

    bottom.invert_yaxis()
    bottom.set_ylabel("mixed-layer depth (m)")
    bottom.set_title(f"Mixed layer, first depth {MIXED_LAYER_DROP:.0f} °C "
                     "below the surface")
    bottom.legend(loc="lower left", fontsize=8.5)

    months(bottom, table.observed_day.min(), table.observed_day.max())

    fig.savefig(path)
    plt.close(fig)


def main():
    where = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    prf = (sys.argv[2] if len(sys.argv) > 2 else
           "/sessions/festive-sharp-feynman/mnt/v5/Examples/DeGray/prf_dam.npt")

    if not os.path.exists(prf):
        print(f"no observed profiles at {prf}")
        return 1

    style()

    profiles = pd.read_csv(os.path.join(where, "degray_profiles.csv"))
    observed = read_prf(prf)

    print(f"{len(observed)} measured casts, "
          f"days {min(observed):.1f} to {max(observed):.1f}")

    table = compare(profiles, observed)

    if table.empty:
        print("no casts overlap the modelled period")
        return 1

    print("\nPer cast:")
    print(table[["observed_day", "n", "bias", "rmse", "surface_observed",
                 "surface_model", "mld_observed", "mld_model"]]
          .to_string(index=False, float_format=lambda v: f"{v:8.2f}"))

    print("\nBy depth band, pooled over every cast:")
    print(depth_bands(profiles, observed)
          .to_string(index=False, float_format=lambda v: f"{v:8.2f}"))

    pooled = np.sqrt((table.rmse ** 2 * table.n).sum() / table.n.sum())

    print(f"\nOverall: bias {np.average(table.bias, weights=table.n):+.2f} °C, "
          f"RMSE {pooled:.2f} °C over {int(table.n.sum())} measurements")

    table.to_csv(os.path.join(where, "degray_observed_comparison.csv"),
                 index=False)

    figure_profiles(profiles, observed,
                    os.path.join(where, "degray_vs_observed_profiles.png"))
    figure_bias(profiles, observed,
                os.path.join(where, "degray_vs_observed_bias.png"))
    figure_mixed_layer(table,
                       os.path.join(where, "degray_vs_observed_mixedlayer.png"))

    print(f"\nwrote the comparison table and three figures to {where}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
