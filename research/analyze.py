"""Score as-issued Open-Meteo forecasts against SNOTEL / AMeDAS daily observations.

Forecast for local day D at lead N = sum of hourly precipitation_previous_dayN over D.
72h forecast starting D = day1(D) + day2(D+1) + day3(D+2), i.e. roughly one issue time.
"""
import glob, json, math, os, re, sys
from collections import defaultdict
from datetime import date, timedelta
from fetch import CACHE, SITES, JP_SITES, SEASONS, MODELS, JP_MODELS

IN = 25.4
ENS_US = ["gfs_seamless", "ecmwf_ifs025", "icon_seamless", "gem_seamless", "ncep_nbm_conus"]
ENS_JP = ["jma_seamless", "ecmwf_ifs025", "icon_seamless", "gem_seamless", "gfs_seamless"]


def d(s):
    return date.fromisoformat(s[:10])


def season_of(day):
    for i, (a, b) in enumerate(SEASONS):
        if d(a) <= day <= d(b):
            return i
    return None


def load_snotel(site):
    out = defaultdict(dict)  # day -> {precip, swe, snwd, tavg}
    for f in glob.glob(f"{CACHE}/obs_{site}_*.json"):
        j = json.load(open(f))
        el = {e["stationElement"]["elementCode"]: {d(v["date"]): v.get("value") for v in e["values"]}
              for e in j[0]["data"]}
        for day, p in el.get("PREC", {}).items():
            nxt = el["PREC"].get(day + timedelta(1))
            if p is not None and nxt is not None:
                out[day]["precip"] = max(0.0, (nxt - p) * IN)  # SNOTEL PREC is start-of-day cumulative
            w, w1 = el.get("WTEQ", {}).get(day), el.get("WTEQ", {}).get(day + timedelta(1))
            if w is not None and w1 is not None:
                out[day]["swe"] = max(0.0, (w1 - w) * IN)
            s, s1 = el.get("SNWD", {}).get(day), el.get("SNWD", {}).get(day + timedelta(1))
            if s is not None and s1 is not None:
                out[day]["newsnow"] = max(0.0, (s1 - s) * 2.54)  # cm, crude: depth gain
        for day, t in el.get("TAVG", {}).items():
            if t is not None:
                out[day]["tavg"] = (t - 32) / 1.8
    return out


def num(x):
    x = x.strip().rstrip(")]").strip()
    try:
        return float(x)
    except ValueError:
        return None


def load_amedas(site):
    out = defaultdict(dict)
    for f in glob.glob(f"{CACHE}/amedas_{site}_*.html"):
        y, m = map(int, re.search(r"(\d{4})-(\d\d)\.html", f).groups())
        t = open(f, encoding="utf-8").read()
        for r in re.findall(r'<tr class="mtx"[^>]*>(.*?)</tr>', t, re.S):
            c = [re.sub("<[^>]+>", "", x) for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
            if len(c) < 18 or not c[0].strip().isdigit():
                continue
            day = date(y, m, int(c[0]))
            for k, i in (("precip", 1), ("tavg", 4), ("newsnow", 16)):
                v = num(c[i])
                if v is not None:
                    out[day][k] = v
    return out


def load_fc(site, model):
    """-> {var_lead: {day: daily value}}; precip/snow summed, temp averaged; needs full 24h."""
    acc = defaultdict(lambda: defaultdict(list))
    elev = None
    for f in glob.glob(f"{CACHE}/fc_{site}_{model}_*.json"):
        j = json.load(open(f))
        if "hourly" not in j:
            continue
        elev = j.get("elevation")
        h = j["hourly"]
        for k, vals in h.items():
            if k == "time":
                continue
            base = k.replace(f"_{model}", "")
            for t, v in zip(h["time"], vals):
                acc[base][d(t)].append(v)
    out = {}
    for k, days in acc.items():
        agg = {}
        for day, vs in days.items():
            if len(vs) == 24 and all(v is not None for v in vs):
                agg[day] = sum(vs) / 24 if k.startswith("temperature") else sum(vs)
        out[k] = agg
    return out, elev


def series72(fc):
    p1, p2, p3 = (fc.get(f"precipitation_previous_day{n}", {}) for n in (1, 2, 3))
    out = {}
    for day in p1:
        a, b = day + timedelta(1), day + timedelta(2)
        if a in p2 and b in p3:
            out[day] = p1[day] + p2[a] + p3[b]
    return out


def obs72(obs, key="precip"):
    out = {}
    for day in obs:
        vs = [obs.get(day + timedelta(i), {}).get(key) for i in range(3)]
        if all(v is not None for v in vs):
            out[day] = sum(vs)
    return out


def scores(pairs, thr):
    """pairs: list of (fc, obs). Returns dict of metrics."""
    n = len(pairs)
    if n == 0:
        return None
    f = [p[0] for p in pairs]; o = [p[1] for p in pairs]
    mae = sum(abs(a - b) for a, b in pairs) / n
    bias = sum(f) / sum(o) if sum(o) else float("nan")
    mf, mo = sum(f) / n, sum(o) / n
    cov = sum((a - mf) * (b - mo) for a, b in pairs)
    sf = math.sqrt(sum((a - mf) ** 2 for a in f)); so = math.sqrt(sum((b - mo) ** 2 for b in o))
    r = cov / (sf * so) if sf and so else float("nan")
    hit = sum(1 for a, b in pairs if a >= thr and b >= thr)
    miss = sum(1 for a, b in pairs if a < thr and b >= thr)
    fa = sum(1 for a, b in pairs if a >= thr and b < thr)
    pod = hit / (hit + miss) if hit + miss else float("nan")
    prec = hit / (hit + fa) if hit + fa else float("nan")
    return dict(n=n, mae=mae, bias=bias, r=r, pod=pod, prec=prec, events=hit + miss)


def fmt(s):
    if not s:
        return "    (no data)"
    return (f"n={s['n']:4d} MAE={s['mae']:5.1f} bias={s['bias']:4.2f} r={s['r']:4.2f} "
            f"| ev={s['events']:3d} POD={s['pod']:4.2f} prec={s['prec']:4.2f}")


def mean_series(series_list):
    common = set.intersection(*[set(s) for s in series_list]) if series_list else set()
    return {k: sum(s[k] for s in series_list) / len(series_list) for k in common}


def run(sites, models, ens_members, load_obs, label):
    print(f"\n######## {label} ########")
    # thresholds in mm SWE: 24h >= 12.7 (~0.5", ~6-8" snow); 72h >= 38 (~1.5", ~18"+ snow)
    T24, T72 = 12.7, 38.0
    agg = defaultdict(list)       # (model, kind, season) -> pairs, pooled over sites
    temp_bias = defaultdict(list)
    big = defaultdict(list)       # model -> (site, day, fc, obs) for top obs 72h storms
    for site in sites:
        obs = load_obs(site)
        o72 = obs72(obs)
        fcs = {}
        for m in models:
            fc, elev = load_fc(site, m)
            if fc:
                fcs[m] = (fc, elev)
        members = [fcs[m][0] for m in ens_members if m in fcs]
        ens = {f"precipitation_previous_day{n}": mean_series([x.get(f"precipitation_previous_day{n}", {}) for x in members])
               for n in (1, 2, 3, 5)}
        fcs["ENSMEAN"] = (ens, None)
        top = sorted(o72.items(), key=lambda kv: -kv[1])
        # de-overlap storms: keep top windows at least 3 days apart
        storms = []
        for day, v in top:
            if all(abs((day - s).days) >= 3 for s, _ in storms):
                storms.append((day, v))
            if len(storms) == 8:
                break
        print(f"\n--- {site}: grid elev used by Open-Meteo = {fcs.get('best_match', (0, '?'))[1]} m")
        for m, (fc, _) in fcs.items():
            line = []
            for n in (1, 2, 3, 5):
                ser = fc.get(f"precipitation_previous_day{n}", {})
                pairs = [(ser[k], obs[k]["precip"]) for k in ser if "precip" in obs.get(k, {}) and season_of(k) is not None]
                for k in ser:
                    if "precip" in obs.get(k, {}) and season_of(k) is not None:
                        agg[(m, f"24h_L{n}", season_of(k))].append((ser[k], obs[k]["precip"]))
                s = scores(pairs, T24)
                line.append(f"L{n} MAE={s['mae']:4.1f} b={s['bias']:4.2f}" if s else f"L{n}   --")
            f72 = series72(fc)
            for k in f72:
                if k in o72 and season_of(k) is not None:
                    agg[(m, "72h", season_of(k))].append((f72[k], o72[k]))
            s72 = scores([(f72[k], o72[k]) for k in f72 if k in o72], T72)
            print(f"  {m:22s} " + " | ".join(line) + f" || 72h {fmt(s72)}")
            for day, v in storms:
                if day in f72:
                    big[m].append((site, day, f72[day], v))
            t = fc.get("temperature_2m_previous_day1", {})
            for k in t:
                if "tavg" in obs.get(k, {}):
                    temp_bias[m].append(t[k] - obs[k]["tavg"])
    # pooled tables per season (only models with data that season)
    for si, (a, b) in enumerate(SEASONS):
        print(f"\n=== POOLED over sites, season {a[:4]}-{b[2:4]} (mm SWE; events 24h>={T24}, 72h>={T72}) ===")
        for kind, thr in (("24h_L1", T24), ("24h_L3", T24), ("24h_L5", T24), ("72h", T72)):
            print(f" [{kind}]")
            rows = [(m, scores(agg[(m, kind, si)], thr)) for m in list(models) + ["ENSMEAN"]]
            for m, s in sorted([r for r in rows if r[1]], key=lambda r: r[1]["mae"]):
                print(f"   {m:22s} {fmt(s)}")
    # cross-season multiplicative bias correction per site, 72h
    print("\n=== D: per-site multiplicative bias correction (train on other season), 72h, pooled ===")
    for m in list(models) + ["ENSMEAN"]:
        raw, cor = [], []
        for site in sites:
            obs = load_obs(site); o72 = obs72(obs)
            fc = (load_fc(site, m)[0] if m != "ENSMEAN" else None)
            if m == "ENSMEAN":
                members = [load_fc(site, x)[0] for x in ens_members]
                members = [x for x in members if x]
                fc = {f"precipitation_previous_day{n}": mean_series([x.get(f"precipitation_previous_day{n}", {}) for x in members]) for n in (1, 2, 3)}
            f72 = series72(fc)
            by = defaultdict(list)
            for k in f72:
                if k in o72 and season_of(k) is not None:
                    by[season_of(k)].append((f72[k], o72[k]))
            if len(by) < 2:
                continue
            for si in by:
                train = by[1 - si]
                kfac = sum(o for _, o in train) / max(1e-6, sum(f for f, _ in train))
                raw += by[si]
                cor += [(f * kfac, o) for f, o in by[si]]
        if raw:
            print(f"   {m:22s} raw: {fmt(scores(raw, T72))}\n   {'':22s} cor: {fmt(scores(cor, T72))}")
    print("\n=== Large storms (top-8 observed 72h windows per site per both seasons), 72h fc/obs ===")
    for m in list(models) + ["ENSMEAN"]:
        if not big[m]:
            continue
        ratios = [f / o for _, _, f, o in big[m] if o > 0]
        ratios.sort()
        med = ratios[len(ratios) // 2]
        under = sum(1 for r in ratios if r < 0.5)
        mae = sum(abs(f - o) for _, _, f, o in big[m]) / len(big[m])
        print(f"   {m:22s} n={len(ratios):2d} median fc/obs={med:4.2f}  MAE={mae:5.1f}mm  share <50% of obs={under/len(ratios):4.2f}")
    print("\n=== Daily-mean 2m temp bias vs station (lead 1), degC ===")
    for m, v in temp_bias.items():
        if v:
            mae = sum(abs(x) for x in v) / len(v)
            print(f"   {m:22s} n={len(v):4d} bias={sum(v)/len(v):+5.2f} MAE={mae:4.2f}")


def snow_check(sites, models, load_obs, label):
    print(f"\n######## {label}: model snowfall (lead1, cm) vs observed new snow (cm) ########")
    for m in models:
        pairs = []
        for site in sites:
            obs = load_obs(site)
            fc, _ = load_fc(site, m)
            s = fc.get("snowfall_previous_day1", {})
            pairs += [(s[k], obs[k]["newsnow"]) for k in s if "newsnow" in obs.get(k, {}) and season_of(k) is not None]
        if pairs and sum(p[0] for p in pairs) > 0:
            print(f"   {m:22s} {fmt(scores(pairs, 15.0))}   (event: >=15cm/24h)")


if __name__ == "__main__":
    us = [s for s in SITES]
    run(us, MODELS, ENS_US, load_snotel, "US / SNOTEL (precip gauge)")
    snow_check(us, MODELS, load_snotel, "US (SNOTEL depth gain)")
    run(list(JP_SITES), JP_MODELS, ENS_JP, load_amedas, "JAPAN / AMeDAS valley stations")
    snow_check(list(JP_SITES), JP_MODELS, load_amedas, "JAPAN (AMeDAS 降雪の深さ)")
