"""Fetch as-issued Open-Meteo forecasts (Previous Runs API) + SNOTEL daily obs into ./cache."""
import json, os, sys, time, urllib.request, urllib.parse

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
os.makedirs(CACHE, exist_ok=True)

# name: (snotel triplet, lat, lon, fixed-offset tz matching SNOTEL local standard time)
SITES = {
    "palisades": ("784:CA:SNTL", 39.18986, -120.26576, "Etc/GMT+8"),
    "kirkwood_carsonpass": ("1067:CA:SNTL", 38.6927, -120.0022, "Etc/GMT+8"),
    "snowbird": ("766:UT:SNTL", 40.56914, -111.65852, "Etc/GMT+7"),
    "alta_atwater": ("1308:UT:SNTL", 40.59124, -111.63775, "Etc/GMT+7"),
    "jackson_phillips": ("689:WY:SNTL", 43.51687, -110.91258, "Etc/GMT+7"),
}
SEASONS = [("2024-11-01", "2025-04-30"), ("2025-11-01", "2026-04-30")]
MODELS = ["best_match", "gfs_seamless", "ecmwf_ifs025", "ecmwf_ifs", "ecmwf_aifs025_single",
          "ncep_nbm_conus", "icon_seamless", "gem_seamless", "gfs_hrrr"]
LEADS = [1, 2, 3, 5]


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            print("  retry", attempt, e, file=sys.stderr)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(url)


def fetch_forecast(site, model, start, end, sites=SITES):
    path = f"{CACHE}/fc_{site}_{model}_{start}.json"
    if os.path.exists(path):
        return
    _, lat, lon, tz = sites[site]
    hourly = [f"precipitation_previous_day{d}" for d in LEADS]
    hourly += [f"snowfall_previous_day{d}" for d in (1, 2, 3)]
    hourly += ["temperature_2m_previous_day1"]
    q = dict(latitude=lat, longitude=lon, hourly=",".join(hourly), models=model,
             start_date=start, end_date=end, timezone=tz)
    j = get("https://previous-runs-api.open-meteo.com/v1/forecast?" + urllib.parse.urlencode(q))
    json.dump(j, open(path, "w"))
    print("fc", site, model, start, "grid elev", j.get("elevation"))


def fetch_snotel(site, start, end):
    path = f"{CACHE}/obs_{site}_{start}.json"
    if os.path.exists(path):
        return
    trip = SITES[site][0]
    q = dict(stationTriplets=trip, elements="PREC,WTEQ,SNWD,TAVG", duration="DAILY",
             beginDate=start, endDate=end)
    j = get("https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data?" + urllib.parse.urlencode(q))
    json.dump(j, open(path, "w"))
    print("obs", site, start)


# JMA AMeDAS: (prec_no, block_no) as id; valley stations, well below the lifts
JP_SITES = {
    "hakuba_amedas": ((48, "0396"), 36.6983, 137.8617, "Asia/Tokyo"),  # 703 m
    "arai_sekiyama": ((54, "0998"), 36.9333, 138.2217, "Asia/Tokyo"),  # 350 m
}
JP_MODELS = ["best_match", "jma_seamless", "jma_msm", "ecmwf_ifs025", "ecmwf_ifs",
             "ecmwf_aifs025_single", "gfs_seamless", "icon_seamless", "gem_seamless"]


def fetch_amedas(site, year, month):
    path = f"{CACHE}/amedas_{site}_{year}-{month:02d}.html"
    if os.path.exists(path):
        return
    prec, block = JP_SITES[site][0]
    url = (f"https://www.data.jma.go.jp/stats/etrn/view/daily_a1.php?prec_no={prec}"
           f"&block_no={block}&year={year}&month={month}&day=&view=")
    with urllib.request.urlopen(url, timeout=60) as r:
        open(path, "wb").write(r.read())
    print("amedas", site, year, month)


if __name__ == "__main__" and sys.argv[1:] == ["japan"]:
    for site in JP_SITES:
        for start, end in SEASONS:
            y0 = int(start[:4])
            for y, m in [(y0, 11), (y0, 12), (y0 + 1, 1), (y0 + 1, 2), (y0 + 1, 3), (y0 + 1, 4)]:
                fetch_amedas(site, y, m)
                time.sleep(1)
            for m in JP_MODELS:
                fetch_forecast(site, m, start, end, JP_SITES)
                time.sleep(1)
elif __name__ == "__main__":
    for site in SITES:
        for start, end in SEASONS:
            fetch_snotel(site, start, end)
            for m in MODELS:
                fetch_forecast(site, m, start, end)
                time.sleep(1)
