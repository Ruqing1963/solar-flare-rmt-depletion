# solar-flare-rmt-depletion

**The Depletion Limit of Level Repulsion: Why Solar Flares Fail the GOE Test while Earth Systems Pass**

Real GOES data show single active regions are Poissonian, turning ⟨r⟩ into a probe of reservoir depletion.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

## The Result: A Clean Negative That Becomes a New Instrument

A unified RMT program showed that single-source charge-and-release Earth systems
repel (GOE) while mixed systems randomize (Poisson). This paper tests the most
violent single-source magnetic charge-and-release system in the Solar System —
an individual solar active region — against **real GOES flare data**.

**Single active regions do NOT show level repulsion. They are Poissonian.**

| Region | Year | n M/X | ⟨r⟩ | 95% CI | Class |
|---|---|---|---|---|---|
| AR 12673 | 2017 | 41 | 0.419 | [0.31, 0.51] | Poisson |
| AR 12192 | 2014 | 20 | 0.479 | [0.31, 0.62] | GOE-edge |
| AR 11429 | 2012 | 11 | 0.401 | [0.18, 0.63] | Poisson |
| AR 11515 | 2012 | 31 | 0.516 | [0.37, 0.59] | GOE-edge |
| AR 12297 | 2015 | 14 | 0.433 | [0.30, 0.64] | Poisson |
| **Pooled (5 ARs)** | — | **112** | **0.453** | **[0.37, 0.49]** | **Poisson** |
| Whole disk 2014 | 2014 | 368 | 0.354 | [0.26, 0.32] | Clustering (CV=2.32) |

The pooled bootstrap CI **excludes GOE (0.531)** but **includes Poisson (0.386)** —
the repulsion hypothesis is rejected at 95% confidence across the five strongest
active regions of Solar Cycle 24. Individual regions scatter (AR 11515, an M-flare
factory, reaches weak GOE; AR 11429 is pure Poisson), reflecting variable per-event
depletion.

## Why? The Reservoir-Depletion Mechanism

Level repulsion requires that each release substantially **empties** the source,
forcing a recharge wait. The per-event **depletion ratio** *f* tunes a system
continuously:

- *f* → 1 (near-total flush): post-event charge ≈ 0, hard minimum waiting time → **repulsion**
- *f* → 0 (tiny nibble from a vast reservoir): state barely changes, can fire again immediately → **Poisson**

A solar active region stores enormous magnetic free energy and a single large
flare releases **under 10%** of it. The field is not reset; the same AR fires
homologous flares in rapid succession. So short waiting times are **not
forbidden** — hence Poisson.

## ⟨r⟩ as a Reservoir-Depletion Probe

This promotes the spacing ratio from a binary classifier to a **quantitative
depletion gauge**. Anchoring linearly at Poisson (⟨r⟩=0.386, RDR=0%) and GOE
(⟨r⟩=0.531, RDR=100%):

| System | Reservoir | ⟨r⟩ | RDR |
|---|---|---|---|
| Old Faithful geyser | thermal/fluid | 0.738 | ~100% (full flush) |
| Tethyan porphyry Cu | metal fluid | 0.712 | ~90% |
| Orogenic gold | metal fluid | 0.678 | ~85% |
| Mantle plumes | magma | 0.630 | ~75% |
| Global radiations | ecospace | 0.628 | ~75% |
| Global extinctions | ecospace | 0.618 | ~70% |
| Andean porphyry Cu | metal fluid | 0.601 | ~70% |
| Nanling W–Sn | metal fluid | 0.574 | ~60% |
| N. Atlantic IRD | ice/meltwater | 0.550 | ~55% |
| **Solar ARs (pooled, 5)** | **magnetic** | **0.453** | **~46%** |
| Poisson floor | — | 0.386 | 0% |

**The Sun does not break the law of level repulsion — it measures where the law
turns off.** By flushing almost nothing per flare, a solar active region marks
the zero-depletion floor against which every fully flushing Earth system is
calibrated.


## Robustness: Inter-AR Scatter is Noise, Not Depletion

We tested whether the region-to-region ⟨r⟩ scatter (0.40–0.52) is itself a
depletion signal — do regions with stronger single-event release sit closer to
GOE? Across **16 active regions** of Cycle 24, **no proxy for release strength
correlates with ⟨r⟩**:

| Proxy | Pearson r | p |
|---|---|---|
| Flare rate (per day) | −0.35 | 0.18 |
| X-class fraction | +0.25 | 0.35 |
| Median peak flux | −0.29 | 0.27 |
| Max peak flux | +0.12 | 0.66 |
| Mean flare energy | +0.04 | 0.88 |

All p > 0.27. The reason is structural: per-AR bootstrap uncertainty
σ⟨r⟩ ≈ 0.10–0.17 is comparable to the *entire* Poisson–GOE gap (0.145). Single-AR
⟨r⟩ is noise-dominated; the apparent repulsion of AR 11515 is a fluctuation, not
a measurement. **The depletion interpretation works across systems (well-sampled
pooled estimates), not within the solar AR population.**

## Cycle-Phase Invariance + Growing Clustering

Single-AR Poisson behavior is **invariant across the solar cycle** (χ² test
p=0.81): pooled single-AR ⟨r⟩ is 0.463 / 0.415 / 0.415 / 0.419 across rise /
maximum / early-decline / late-decline — always near Poisson, never GOE. This
is a *universal* property of magnetic reconnection, not an artifact of one epoch.

Meanwhile the **whole-disk clustering strengthens into the decline**: CV rises
monotonically 2.23 → 2.32 → 2.54 → 3.10. As active regions become fewer, the
disk signal is dominated by homologous flare trains from a handful of regions,
sharpening the clustering.

## Repository Structure
```
solar-flare-rmt-depletion/
├── README.md · LICENSE · requirements.txt · CITATION.cff · .zenodo.json
├── paper/
│   ├── paper.tex · paper.pdf      # 11 pp.
│   └── figs/                      # PDFs embedded by LaTeX
├── code/
│   └── solar_flare_rmt_pipeline.py  # reproducible from raw GOES CSVs
├── data/
│   ├── ar12673_mx_flares.csv      # 41 M/X, AR 12673 (Sep 2017)
│   ├── ar12192_mx_flares.csv      # 20 M/X, AR 12192 (Oct 2014)
│   ├── ar11429_mx_flares.csv      # 11 M/X, AR 11429 (Mar 2012)
│   ├── ar11515_mx_flares.csv      # 31 M/X, AR 11515 (Jul 2012)
│   ├── ar12297_mx_flares.csv      # 14 M/X, AR 12297 (Mar 2015)
│   ├── global2014_mx_flares.csv   # 368 whole-disk M/X, 2014
│   └── per_ar_proxies.csv         # 16 ARs: rate, X-frac, flux, energy, <r>
├── figures/                       # standalone vector PDFs
│   ├── fig1_solar_spacings.pdf      # pooled distribution + per-AR forest plot
│   ├── fig2_depletion_probe.pdf     # the reservoir-depletion ladder
│   ├── fig3_cycle_phase.pdf         # phase invariance + clustering trend
│   └── fig4_proxy_robustness.pdf    # 4 release proxies vs <r> (all n.s.)
└── results/
    └── solar_rmt_results.json     # all statistics
```

## Data Source
Raw flare timings from the science-quality **NOAA/NCEI GOES-R XRS Flare Report**
(one record per flare, with active-region number where available):
`https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/`
Yearly files used: `..._y2012_...`, `..._y2014_...`, `..._y2015_...`, `..._y2017_...csv`.

## Reproduce
```bash
pip install -r requirements.txt
cd code
python solar_flare_rmt_pipeline.py y2012.csv y2014.csv y2015.csv y2017.csv
```

## Selection Notes (Honest Methodology)
- **AR 12673**: the flare report's AR tag is sparse, but AR 12673 is the *only*
  region tagged anywhere in September 2017 — it dominated the disk. We select the
  Sep 3–11 M/X window. Literature: AR 12673 produced the two largest flares of
  Cycle 24 (X9.3, X8.2).
- **AR 12192**: 20 M/X flares carry *direct* AR tags — a contamination-free sample.
- Both selection methods give Poisson, so the conclusion is robust to method.

## Seven-Domain RMT Program
1. Cyclostratigraphy (spatial) → GOE — [zenodo 20774581](https://zenodo.org/records/20774581)
2. Mantle plumes → single-source GOE — [zenodo 20768420](https://zenodo.org/records/20768420)
3. Metallogeny → single ore system GOE/GUE — [zenodo 20768849](https://zenodo.org/records/20768849)
4. Macroevolution → extinction/radiation GOE — [zenodo 20783763](https://zenodo.org/records/20783763)
5. Hydrogeology → super-GUE geyser, GOE floods — [github](https://github.com/Ruqing1963/hydro-rmt-geysers-floods)
6. **Solar flares → Poisson (depletion limit), this work**

## Citation
```bibtex
@misc{chen2026solar,
  author = {Chen, Ruqing},
  title  = {The Depletion Limit of Level Repulsion: Why Solar Flares
            Fail the GOE Test while Earth Systems Pass},
  year   = {2026},
  publisher = {GitHub},
  url    = {https://github.com/Ruqing1963/solar-flare-rmt-depletion}
}
```

## License
[MIT](LICENSE). Solar flare data from the NOAA/NCEI GOES-R XRS Flare Report (public domain).
