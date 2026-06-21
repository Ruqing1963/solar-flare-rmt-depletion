# Data

All flare timings extracted from the science-quality **NOAA/NCEI GOES-R XRS
Flare Report** (one record per flare). Source directory:
`https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/multi/l2/data/xrsf-l2-flrpt_science/`
Yearly files: `sci_xrsf-l2-flrpt_geo_y2014_v1-0-0.csv`, `sci_xrsf-l2-flrpt_geo_y2017_v1-0-0.csv`.
Public domain.

## ar12673_mx_flares.csv
41 M- and X-class flares from AR 12673, 2017-09-03 to 2017-09-11. AR 12673 is the
only active region tagged in the flare report anywhere in September 2017 and
produced the two largest flares of Solar Cycle 24 (X9.3 on Sep 6, X8.2 on Sep 10).

**Columns:** time (peak, UTC), start_time, end_time, flare_class, xrsb_irrad (W/m²), active_region

## ar12192_mx_flares.csv
20 M/X flares from AR 12192 (Oct 2014) carrying *direct* active-region tags in the
flare report — a contamination-free single-AR sample. AR 12192 was the largest
sunspot group of Solar Cycle 24.

**Columns:** time (peak, UTC), start_time, end_time, flare_class, xrsb_irrad (W/m²), active_region

## global2014_mx_flares.csv
All 368 M/X flares across the whole solar disk in 2014 (all active regions), the
most active year of Solar Cycle 24. Used to test the superposition/clustering
prediction.

**Columns:** time (peak, UTC), flare_class, xrsb_irrad (W/m²), active_region

## Analysis note
Waiting times are computed as differences of consecutive flare peak times (hours),
then normalized to unit mean before computing RMT statistics. See
`code/solar_flare_rmt_pipeline.py`.

## ar11429_mx_flares.csv, ar11515_mx_flares.csv, ar12297_mx_flares.csv
Three additional strong active regions of Solar Cycle 24, selected by *direct*
active-region tag (the 2012/2015 flare reports have dense AR tagging):
- **AR 11429** (Mar 2012): 11 M/X flares (2 X, 9 M)
- **AR 11515** (Jul 2012): 31 M/X flares (0 X, 31 M) — an M-flare factory, most repulsive
- **AR 12297** (Mar 2015): 14 M/X flares (1 X, 13 M)

Pooling all five regions (AR 12673, 12192, 11429, 11515, 12297) gives 112 waiting
times; the pooled ⟨r⟩=0.453 excludes GOE at 95% confidence.

**Columns:** time (peak, UTC), start_time, end_time, flare_class, xrsb_irrad (W/m²), active_region
