# Forecast verification study (research only, not used by the site)

Scores as-issued Open-Meteo forecasts (Previous Runs API, no look-ahead) against
SNOTEL (US) and JMA AMeDAS (Japan) daily observations, winters 2024-25 and 2025-26.

    python3 fetch.py          # US: SNOTEL + forecasts -> cache/
    python3 fetch.py japan    # Japan: AMeDAS HTML + forecasts -> cache/
    python3 analyze.py > results.txt

Truth = liquid-equivalent precipitation (SNOTEL PREC gauge increments; AMeDAS daily precip).
Snowfall is scored separately against SNOTEL depth gain / AMeDAS 降雪の深さ, both crude.

Caveats: SNOTEL gauges under-catch snow in wind (10-30%+), so model bias ~1.1-1.2 may be
"correct". AMeDAS stations are valley sites (Hakuba 703 m, Sekiyama 350 m), not on-hill.
Lead-N daily totals splice hours from runs issued ~N days earlier, not one single run.
