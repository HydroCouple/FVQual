#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin
"""Figures for the DeGray Reservoir verification run.

Reads the three CSV files the DeGray application writes and draws what a
limnologist would ask to see: the depth-time temperature contour, the surface
temperature against the equilibrium temperature the weather implies, the heat
and volume budgets, and seasonal profiles.

Every number annotated on a figure is computed from the CSV here rather than
typed in, so a figure cannot drift away from the run it describes.

Usage: python3 plot_degray.py [directory]
"""

import sys
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# 1980 was a leap year, which is what DeGray's Julian days are counted in.
MONTH_STARTS = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336]
MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

INK = "#1b1b1b"
MUTED = "#6b6b6b"
GRID = "#d8d8d8"


def style():
    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 9.5,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlecolor": INK,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.alpha": 0.9,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def months(ax, lo, hi):
    ticks = [d for d in MONTH_STARTS if lo - 5 <= d <= hi + 5]
    ax.set_xticks(ticks)
    ax.set_xticklabels([MONTH_NAMES[MONTH_STARTS.index(t)] for t in ticks])
    ax.set_xlim(lo, hi)


def isotherms(daily, profiles, path):
    """Depth-time temperature contour at the dam, the defining figure."""
    days = np.sort(profiles.jday.unique())
    layers = np.sort(profiles.layer.unique())

    # The sigma grid moves with the surface, so depth is a function of both day
    # and layer. Building 2-D coordinates keeps that exact instead of pretending
    # a layer sits at a fixed depth.
    temp = profiles.pivot(index="layer", columns="jday", values="temperature_C").values
    depth = profiles.pivot(index="layer", columns="jday", values="depth_m").values
    jday = np.tile(days, (len(layers), 1))

    fig, ax = plt.subplots(figsize=(10.5, 4.6))

    levels = np.arange(5, 33, 1.0)
    fill = ax.contourf(jday, depth, temp, levels=levels, cmap="RdYlBu_r", extend="both")

    lines = ax.contour(jday, depth, temp, levels=np.arange(6, 32, 2.0),
                       colors="white", linewidths=0.5, alpha=0.55)
    ax.clabel(lines, fmt="%.0f", fontsize=7, colors="white")

    # The steepest-gradient tracker is only meaningful while there *is* a
    # thermocline. Spread 2 °C over 55 m and the mean gradient is 0.04 °C/m,
    # at which point the tracker is choosing between layer differences of the
    # same size as the transport's own round-off and jumps the full depth of the
    # column from one day to the next. Masked rather than drawn, because a line
    # that wanders is read as a signal.
    strat = daily.surface_temperature_C - daily.bottom_temperature_C
    tracked = daily.thermocline_depth_m.where(strat >= 2.0)

    ax.plot(daily.jday, tracked, color="black", lw=1.4,
            label="steepest gradient (while stratified)")

    ax.invert_yaxis()
    ax.set_ylabel("depth below surface (m)")
    ax.set_ylim(depth.max(), 0)
    ax.grid(False)
    months(ax, days.min(), days.max())

    bar = fig.colorbar(fill, ax=ax, pad=0.015, aspect=28)
    bar.set_label("temperature (°C)")
    bar.outline.set_visible(False)

    peak = daily.loc[daily.surface_temperature_C.idxmax()]

    ax.set_title("DeGray Reservoir at the dam — temperature with depth and time, 1980")
    ax.legend(loc="lower left", fontsize=8, labelcolor="black")

    ax.annotate(f"surface peaks at {peak.surface_temperature_C:.1f} °C",
                xy=(peak.jday, 1.0), xytext=(peak.jday - 62, 12.0),
                fontsize=8, color="black",
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8))

    fig.savefig(path)
    plt.close(fig)


def temperatures(daily, deepest, path):
    """Surface and bottom against the forcing that drives them."""
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(10.5, 6.4), sharex=True,
                                      gridspec_kw={"height_ratios": [3, 2]})

    # Sampled at the same instant each day — a shade after noon, because the run
    # starts at day 64.5 and writes every whole day after. So this is the
    # equilibrium temperature under peak sun, not a daily mean, which is why it
    # runs to 40 °C in July and below freezing on clear December nights. A
    # reservoir with this much thermal mass never gets near it.
    top.plot(daily.jday, daily.equilibrium_C, color="#c9a227", lw=0.8, alpha=0.7,
             label="equilibrium temperature (near solar noon)")
    top.plot(daily.jday, daily.surface_temperature_C, color="#c0392b", lw=1.8,
             label="surface")
    top.plot(daily.jday, daily.bottom_temperature_C, color="#1f5f8b", lw=1.8,
             label=f"bottom layer, {deepest:.0f} m down")

    top.set_ylabel("temperature (°C)")
    top.set_ylim(top.get_ylim()[0], daily.equilibrium_C.max() + 12.0)
    top.legend(loc="upper left", ncol=3, fontsize=8.5)
    top.yaxis.set_major_locator(MultipleLocator(10))
    top.set_title("Surface follows the weather; the hypolimnion does not")

    strat = daily.surface_temperature_C - daily.bottom_temperature_C
    peak = daily.jday[strat.idxmax()]

    top.axvline(peak, color=MUTED, lw=0.8, ls=":")
    top.annotate(f"{strat.max():.1f} °C apart", xy=(peak, 20), xytext=(peak + 6, 27),
                 fontsize=8, color=INK)

    bottom.fill_between(daily.jday, 0, strat, color="#8e44ad", alpha=0.18)
    bottom.plot(daily.jday, strat, color="#8e44ad", lw=1.4)
    bottom.set_ylabel("stratification (°C)")
    bottom.axhline(1.0, color=MUTED, lw=0.8, ls="--")
    bottom.annotate("1 °C — effectively isothermal", xy=(daily.jday.min() + 4, 1.6),
                    fontsize=7.5, color=MUTED)

    months(bottom, daily.jday.min(), daily.jday.max())

    fig.savefig(path)
    plt.close(fig)


def budgets(daily, path):
    """Heat and water, the two things that must add up."""
    fig, (heat, water) = plt.subplots(2, 1, figsize=(10.5, 5.4), sharex=True)

    heat.plot(daily.jday, daily.heat_content_J / 1e16, color="#c0392b", lw=1.6)
    heat.set_ylabel("heat content (10$^{16}$ J)")
    heat.set_title("Whole-reservoir heat and volume")

    change = 100.0 * (daily.volume_m3 / daily.volume_m3.iloc[0] - 1.0)

    water.plot(daily.jday, change, color="#1f5f8b", lw=1.6)
    water.axhline(0.0, color=MUTED, lw=0.8)
    water.set_ylabel("volume change (%)")
    water.annotate(f"net {change.iloc[-1]:+.2f}% over the run, inflow less release",
                   xy=(daily.jday.min() + 4, change.max() * 0.62),
                   fontsize=8, color=MUTED)

    months(water, daily.jday.min(), daily.jday.max())

    fig.savefig(path)
    plt.close(fig)


def profiles_figure(profiles, final, path):
    """Seasonal profiles, which is how a thermal survey reports a reservoir."""
    days = np.sort(profiles.jday.unique())

    wanted = [(days.min(), "early March — nearly isothermal", "#1f5f8b"),
              (150, "late May — stratifying", "#27ae60"),
              (200, "mid July — fully stratified", "#c0392b"),
              (280, "early October — thermocline deepening", "#e67e22"),
              (days.max(), "late December — turned over", "#6b6b6b")]

    fig, ax = plt.subplots(figsize=(5.6, 6.4))

    for target, label, colour in wanted:
        day = days[np.abs(days - target).argmin()]
        slab = profiles[profiles.jday == day].sort_values("depth_m")

        ax.plot(slab.temperature_C, slab.depth_m, color=colour, lw=1.8, label=label)

    ax.invert_yaxis()
    ax.set_xlabel("temperature (°C)")
    ax.set_ylabel("depth below surface (m)")
    ax.set_title("Profiles through the year")
    ax.legend(loc="lower right", fontsize=8)

    fig.savefig(path)
    plt.close(fig)


def storage_curve(bathymetry_path, bottom_elevation=66.35):
    """Stage-storage for DeGray, read straight from the W2 bathymetry file.

    Deliberately a second implementation, in a different language, from the raw
    file rather than from anything FVQual computed. Its only job is to disagree
    with the model if the model got the storage curve wrong.
    """
    lengths, thickness, widths = None, [], []

    with open(bathymetry_path) as handle:
        in_layers = False

        for line in handle:
            fields = [f.strip() for f in line.strip().split(",")]

            if not fields or not fields[0]:
                if fields and len(fields) > 2 and lengths is None:
                    continue

            label = fields[0]

            if label == "DLX":
                lengths = np.array([float(f) for f in fields[1:] if f])
            elif label == "LAYERH":
                in_layers = True
            elif in_layers and label:
                values = [f for f in fields if f]
                thickness.append(float(values[0]))
                widths.append([float(f) for f in values[1:1 + len(lengths)]])

    thickness = np.array(thickness)
    widths = np.array(widths)

    # K runs downward in the file, so the last row sits on the bed.
    area = widths * lengths[None, :]          # plan area of each cell (m²)
    cell = area * thickness[:, None]          # its volume (m³)

    bottoms = bottom_elevation + np.cumsum(thickness[::-1])[::-1] - thickness

    def storage(z):
        wet = np.clip((z - bottoms) / thickness, 0.0, 1.0)
        return float((cell * wet[:, None]).sum())

    def plan(z):
        wet = (z > bottoms) & (z <= bottoms + thickness)
        return float(area[wet].sum()) if wet.any() else 0.0

    return storage, plan, bottoms


def stage_storage(daily, bathymetry_path, path):
    """Does the model's water sit on the curve the raw bathymetry file implies?"""
    storage, plan, bottoms = storage_curve(bathymetry_path)

    levels = np.linspace(daily.surface_elevation_m.min() - 0.4,
                         daily.surface_elevation_m.max() + 0.4, 120)
    curve = np.array([storage(z) for z in levels])

    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.4))

    left.plot(curve / 1e8, levels, color=MUTED, lw=3.0, alpha=0.45,
              label="W2 bathymetry file, read independently")
    left.scatter(daily.volume_m3 / 1e8, daily.surface_elevation_m, s=7,
                 color="#1f5f8b", zorder=3, label="FVQual, every day of the run")

    left.set_xlabel("storage (10$^8$ m³)")
    left.set_ylabel("surface elevation (m)")
    left.set_title("Stage and storage agree")
    left.legend(loc="lower right", fontsize=8)

    modelled = np.array([storage(z) for z in daily.surface_elevation_m])
    error = 100.0 * (daily.volume_m3.values - modelled) / modelled

    right.plot(daily.jday, error, color="#c0392b", lw=1.2)
    right.axhline(0.0, color=MUTED, lw=0.8)
    right.set_ylabel("model storage less file storage (%)")
    right.set_title("Disagreement over the year")
    months(right, daily.jday.min(), daily.jday.max())

    right.annotate(f"worst {np.abs(error).max():.2f}%, mean {error.mean():+.2f}%",
                   xy=(0.03, 0.9), xycoords="axes fraction", fontsize=8, color=INK)

    fig.savefig(path)
    plt.close(fig)

    return error


def main():
    where = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    bathymetry = (sys.argv[2] if len(sys.argv) > 2 else
                  "/sessions/festive-sharp-feynman/mnt/v5/Examples/DeGray/InputFiles/bth1.csv")

    style()

    daily = pd.read_csv(os.path.join(where, "degray_daily.csv"))
    profiles = pd.read_csv(os.path.join(where, "degray_profiles.csv"))
    final = pd.read_csv(os.path.join(where, "degray_final_profile.csv"))

    isotherms(daily, profiles, os.path.join(where, "degray_isotherms.png"))
    temperatures(daily, profiles.depth_m.max(),
                 os.path.join(where, "degray_temperatures.png"))
    budgets(daily, os.path.join(where, "degray_budgets.png"))
    profiles_figure(profiles, final, os.path.join(where, "degray_profiles.png"))

    if os.path.exists(bathymetry):
        error = stage_storage(daily, bathymetry,
                              os.path.join(where, "degray_stage_storage.png"))
        print(f"stage-storage: worst {np.abs(error).max():.3f}%, "
              f"mean {error.mean():+.3f}%")
    else:
        print(f"skipped the stage-storage check: no bathymetry at {bathymetry}")

    print(f"wrote figures to {where}")


if __name__ == "__main__":
    main()
