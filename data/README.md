# Data

| File | Content | Period |
|---|---|---|
| `bao_co2_pfp_flask.csv` | CO₂ dry-air mole fraction (ppm, WMO CO₂ X2019 scale) of individual flask-air samples collected by programmable flask packages (PFP) at the Boulder Atmospheric Observatory tower | 16 Aug 2007 – 06 Jul 2016 |

**Site:** BAO, Erie, Colorado, USA; 40.05° N, 105.004° W; surface elevation 1578.7 m a.s.l.; sampling intake 305.3 m above ground.

**Columns:**

| Column | Description |
|---|---|
| `datetime_utc` | Sampling time (UTC) |
| `decimal_year` | Sampling time as a decimal year |
| `flask_id` | Flask identifier |
| `co2_ppm` | CO₂ mole fraction (ppm) |
| `intake_height_m` | Intake height above ground (m) |
| `method` | Sampling method code |
| `event_number` | NOAA sample event number |
| `instrument` | Analysis instrument code |
| `qc_flag` | Three-character NOAA quality-control flag |

**Processing:**
- **Valid samples:** a sample is valid when the first character of `qc_flag` is `.`, which means no rejection flag. 3414 of the 3687 samples are valid.
- **Selected columns:** the CSV contains the columns of the NOAA event file needed for this analysis. The values themselves are unchanged.

**Source:** NOAA Global Monitoring Laboratory, Global Greenhouse Gas Reference Network, file `co2_bao_surface-pfp_1_ccgg_event.txt` (version 2026-03-02), https://gml.noaa.gov/aftp/data/trace_gases/co2/pfp/surface/. DOI: [10.15138/gr3w-qm07](https://doi.org/10.15138/gr3w-qm07). Please follow NOAA's citation and fair-use guidelines given in the original file header.
