# Latitude maps of phenological sensitivity

`make_latitude_maps_primary.py` creates the three-panel manuscript figure and
a four-panel supplementary figure from the primary geographic models with a
population random intercept. It reads existing exports and does not fit models.

## Run

Use Python 3.12 and install `requirements_latitude_maps.txt` into an environment:

```powershell
python -m pip install -r requirements_latitude_maps.txt
python make_latitude_maps_primary.py
```

Arial must be installed for the manuscript typography. Figures are exported as
PNG (600 dpi), PDF and editable SVG. The main figure is 180 mm wide; the
four-panel supplementary figure is 240 mm wide.

By default the script reads the run in `model_reviews/geography_primary_models.json`
when available, otherwise the reviewed primary run
`output/phenology_plasticity/geography_mundlak_population_intercept_15km/run_2399ac7eb944`.
An alternative export location and output directory can be supplied explicitly:

```powershell
python make_latitude_maps_primary.py --run "path/to/model_exports" --out "path/to/figures"
```

Model exports and generated figures are intentionally excluded from Git.
The run directory must contain `ALL_geography_within_between_effects.csv` and
`population_geography_components.csv` inside each of `onset_univoltine`,
`offset_univoltine`, `offset_multivoltine` and `onset_multivoltine`.
The script checks coordinate consistency, one record per species-site population,
the equal-population latitude decomposition and the map extent.

## Interpretation

The colour for each population is the within-species latitude interaction
coefficient multiplied by its latitude deviation from the species mean, in units
of 10 degrees. Coefficients and colours are oriented towards advancement for
univoltines and multivoltine onset, and delay for multivoltine offset. Colours
at each site average its species-population contributions with equal weights.
This is only the fitted latitude contribution, not total population plasticity.
Zero denotes no latitude contribution at the species mean latitude. Climate
anomaly units are inherited from source models and are not converted to Celsius.
Panels use separate symmetric colour scales without clipping.

The main figure includes the three within-species latitude interactions with
95% Wald intervals excluding zero. The supplementary figure includes all four
analyses. Insets report the oriented interaction coefficients and 95% intervals
per 10 degrees. Offset models condition on annual onset.

Outputs also include a caption, population contributions, site summaries,
coefficient table and source-file SHA256 hashes. The original automatic numerical
flags are retained in the source exports; practical acceptance of the designated
primary run followed a separate numerical and residual review on 26 October 2026.

