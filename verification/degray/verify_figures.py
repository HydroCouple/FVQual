#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin

"""Checks that every claim the DeGray figures make is true of the CSV data.

A figure is an assertion about a run, and an assertion that nobody checks is a
caption. The annotated numbers in plot_degray.py are computed from the data
rather than typed, so those cannot drift; what *can* drift are the words — a
panel titled "the hypolimnion does not follow the weather" or a profile labelled
"fully stratified" is a claim in prose, and prose does not recompute itself.

So each one is restated here as an inequality over the CSV. Run after plotting.
Exits non-zero if a figure says something the run does not.
"""

import os
import sys

import numpy as np
import pandas as pd

FAILURES = []


def check(name, condition, detail):
    mark = "ok  " if condition else "FAIL"
    print(f"  [{mark}] {name}: {detail}")

    if not condition:
        FAILURES.append(name)


def main():
    where = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

    daily = pd.read_csv(os.path.join(where, "degray_daily.csv"))
    profiles = pd.read_csv(os.path.join(where, "degray_profiles.csv"))

    strat = daily.surface_temperature_C - daily.bottom_temperature_C
    days = np.sort(profiles.jday.unique())

    print("degray_isotherms.png")

    peak = daily.loc[daily.surface_temperature_C.idxmax()]
    check("annotated peak is the true maximum",
          np.isclose(peak.surface_temperature_C, daily.surface_temperature_C.max()),
          f"{peak.surface_temperature_C:.2f} °C on day {peak.jday:.1f}")

    # The contour levels run 5 to 32; anything outside would be silently clipped
    # into the "extend" arrows and read as a flat region.
    check("all temperatures fall inside the drawn contour range",
          profiles.temperature_C.min() > 4.0 and profiles.temperature_C.max() < 33.0,
          f"{profiles.temperature_C.min():.2f} to {profiles.temperature_C.max():.2f} °C")

    check("every day has a full profile",
          len(profiles) == len(days) * profiles.layer.nunique(),
          f"{len(days)} days x {profiles.layer.nunique()} layers = {len(profiles)} rows")

    print("\ndegray_temperatures.png")

    # "Surface follows the weather; the hypolimnion does not."
    surface_r = np.corrcoef(daily.surface_temperature_C, daily.equilibrium_C)[0, 1]
    bottom_r = np.corrcoef(daily.bottom_temperature_C, daily.equilibrium_C)[0, 1]

    check("surface tracks the equilibrium temperature more closely than the bottom",
          surface_r > 2.0 * bottom_r,
          f"r = {surface_r:+.2f} at the surface against {bottom_r:+.2f} at the bed")

    summer = daily[(daily.jday > 150) & (daily.jday < 260)]

    check("the hypolimnion stays cold all summer",
          summer.bottom_temperature_C.max() < 10.0,
          f"bottom peaks at {summer.bottom_temperature_C.max():.2f} °C "
          f"between days 150 and 260")

    check("the annotated maximum stratification is the true maximum",
          np.isclose(strat.max(), strat.max()),
          f"{strat.max():.2f} °C on day {daily.jday[strat.idxmax()]:.1f}")

    check("the run begins and ends effectively isothermal",
          strat.iloc[0] < 1.0 and strat.iloc[-1] < 1.0,
          f"{strat.iloc[0]:.2f} °C in March, {strat.iloc[-1]:.2f} °C in December")

    # The equilibrium curve is labelled "near solar noon", so the sample instant
    # had better be near noon and had better be the same each day.
    hour = (daily.jday % 1.0) * 24.0

    check("the daily sample really is taken near solar noon",
          hour.min() > 11.0 and hour.max() < 13.0 and hour.std() < 0.1,
          f"{hour.mean():.2f} h local, spread {hour.std():.3f} h")

    print("\ndegray_budgets.png")

    heat = daily.heat_content_J

    check("heat content peaks in summer, not at an end of the run",
          150 < daily.jday[heat.idxmax()] < 280,
          f"maximum on day {daily.jday[heat.idxmax()]:.0f}")

    change = 100.0 * (daily.volume_m3 / daily.volume_m3.iloc[0] - 1.0)

    check("the volume excursion is a drawdown and refill, not a drift",
          change.min() < -5.0 and abs(change.iloc[-1]) < 3.0,
          f"low {change.min():.2f}%, ends {change.iloc[-1]:+.2f}%")

    print("\ndegray_profiles.png")

    labels = [(days.min(), "nearly isothermal", lambda s, t: s < 1.5),
              (150, "stratifying", lambda s, t: 5.0 < s < 20.0),
              (200, "fully stratified", lambda s, t: s > 20.0),
              (280, "thermocline deepening", lambda s, t: t > 8.0),
              (days.max(), "turned over", lambda s, t: s < 1.5)]

    for target, label, test in labels:
        day = days[np.abs(days - target).argmin()]
        row = daily.iloc[np.abs(daily.jday.values - day).argmin()]

        spread = row.surface_temperature_C - row.bottom_temperature_C

        check(f"day {day:.0f} is '{label}'", test(spread, row.thermocline_depth_m),
              f"{spread:.2f} °C top to bottom, steepest gradient at "
              f"{row.thermocline_depth_m:.1f} m")

    # The figure claims the thermocline deepens through the season.
    tracked = daily[strat >= 2.0]
    midsummer = tracked[(tracked.jday > 170) & (tracked.jday < 210)]
    autumn = tracked[(tracked.jday > 270) & (tracked.jday < 310)]

    check("the thermocline is deeper in autumn than at midsummer",
          autumn.thermocline_depth_m.mean() > 2.0 * midsummer.thermocline_depth_m.mean(),
          f"{midsummer.thermocline_depth_m.mean():.1f} m in July, "
          f"{autumn.thermocline_depth_m.mean():.1f} m in October")

    print("\ndegray_stage_storage.png")

    # dV/dz recovered from the run must be the reservoir's plan area. This is
    # the mesh and the barotropic mode agreeing through the storage curve, and
    # neither of them is told the answer.
    dv = np.diff(daily.volume_m3)
    dz = np.diff(daily.surface_elevation_m)
    moved = np.abs(dz) > 1.0e-4

    area = np.median(dv[moved] / dz[moved])

    check("dV/dz recovers the reservoir's plan area",
          4.5e7 < area < 5.3e7,
          f"{area:.4g} m² against 5.09e7 m² at full pool")

    print()

    if FAILURES:
        print(f"{len(FAILURES)} figure claim(s) not supported by the data:")

        for name in FAILURES:
            print(f"  - {name}")

        return 1

    print("every figure claim checks out against the CSV data")
    return 0


if __name__ == "__main__":
    sys.exit(main())
