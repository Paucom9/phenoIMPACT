"""Map the within-species latitude contribution from the designated primary models.

No model fitting; site colours average equal-weight species-population contributions.
Run with the bundled Python runtime and figure_python_deps on PYTHONPATH.
"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib import font_manager
from mpl_toolkits.basemap import Basemap

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--run', type=Path, help='Directory containing the four geographic model exports')
parser.add_argument('--out', type=Path, help='Output directory')
args = parser.parse_args()
OUT = args.out.resolve() if args.out else ROOT / 'output/figures/latitude_primary_population_maps'
OUT.mkdir(parents=True, exist_ok=True)
selection_file = ROOT / 'model_reviews/geography_primary_models.json'
selection = json.loads(selection_file.read_text()) if selection_file.exists() else {
    'role': 'primary', 'population_random_intercept': '(1 | SPECIES:SITE_ID)',
    'run': 'output/phenology_plasticity/geography_mundlak_population_intercept_15km/run_2399ac7eb944',
    'review_note': 'Practical numerical and residual review completed on 2026-10-26.'}
RUN = args.run.resolve() if args.run else ROOT / selection['run']
effects_file = RUN / 'ALL_geography_within_between_effects.csv'
effects = pd.read_csv(effects_file)
for font in ['Arial']:
    try:
        font_manager.findfont(font, fallback_to_default=False)
        break
    except ValueError:
        continue
plt.rcParams.update({'font.family': font, 'font.size': 7, 'pdf.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.titleweight': 'normal'})
ANALYSES = ['onset_univoltine', 'offset_univoltine', 'offset_multivoltine', 'onset_multivoltine']
site_tables, coefficients, sources = [], [], [effects_file]
for analysis in ANALYSES:
    file = RUN / analysis / 'population_geography_components.csv'
    sources.append(file)
    p = pd.read_csv(file, dtype={'SITE_ID': str, 'SPECIES': str})
    assert not p.duplicated(['SPECIES', 'SITE_ID']).any()
    assert np.isfinite(p[['latitude_deg', 'longitude_deg', 'latitude_within_species_10deg']]).all().all()
    means = p.groupby('SPECIES').latitude_deg.transform('mean')
    assert np.allclose((p.latitude_deg - means) / 10, p.latitude_within_species_10deg)
    e = effects.loc[(effects.analysis == analysis) & (effects.axis == 'latitude') &
                    (effects.effect_scale == 'within species')]
    assert len(e) == 1
    e = e.iloc[0]
    sign = 1 if analysis == 'offset_multivoltine' else -1
    beta = sign * e['estimate']
    lo, hi = sorted([sign * e['conf.low'], sign * e['conf.high']])
    p['latitude_contribution'] = beta * p.latitude_within_species_10deg
    p['analysis'] = analysis
    p.to_csv(OUT / f'{analysis}_population_contributions.csv', index=False)
    coord_spread = p.groupby('SITE_ID')[['latitude_deg', 'longitude_deg']].agg(lambda x: x.max()-x.min())
    assert (coord_spread.to_numpy() < 1e-8).all()
    site = p.groupby('SITE_ID', as_index=False).agg(
        latitude_deg=('latitude_deg', 'first'), longitude_deg=('longitude_deg', 'first'),
        contribution=('latitude_contribution', 'mean'), n_species=('SPECIES', 'nunique'))
    assert site.latitude_deg.between(38, 66).all() and site.longitude_deg.between(-12, 33).all()
    assert abs(p.latitude_contribution.mean()) < 1e-8
    site['analysis'] = analysis
    site_tables.append(site)
    coefficients.append(dict(analysis=analysis, oriented_beta=beta, low=lo, high=hi,
                             p_Wald=e.p_Wald, n_populations=len(p), n_sites=len(site)))
sites = pd.concat(site_tables, ignore_index=True)
coefs = pd.DataFrame(coefficients).set_index('analysis')
sites.to_csv(OUT / 'site_latitude_contributions.csv', index=False)
coefs.to_csv(OUT / 'latitude_coefficients.csv')

def draw(analyses, name):
    width_mm = 180 if len(analyses) == 3 else 240
    fig, axes = plt.subplots(1, len(analyses), figsize=(width_mm/25.4, 90/25.4))
    fig.subplots_adjust(left=.045, right=.99, bottom=.19, top=.91, wspace=.20)
    for j, (ax, analysis) in enumerate(zip(axes, analyses)):
        site = sites.loc[sites.analysis == analysis]
        e = coefs.loc[analysis]
        event, group = analysis.split('_')
        m = Basemap(projection='merc', llcrnrlon=-12, urcrnrlon=33,
                    llcrnrlat=38, urcrnrlat=66, resolution='l', ax=ax)
        m.drawmapboundary(linewidth=.7, fill_color='white')
        m.drawcoastlines(linewidth=.45, color='#484848')
        m.drawcountries(linewidth=.35, color='#777777')
        m.drawparallels(np.arange(40, 66, 5), labels=[1,0,0,0], linewidth=.25,
                        color='#bbbbbb', dashes=[1,0], fontsize=6)
        m.drawmeridians(np.arange(-10, 31, 10), labels=[0,0,0,1], linewidth=.25,
                        color='#bbbbbb', dashes=[1,0], fontsize=6)
        lim = float(np.abs(site.contribution).max())
        x, y = m(site.longitude_deg.to_numpy(), site.latitude_deg.to_numpy())
        sc = m.scatter(x, y, c=site.contribution, s=2.6, cmap='viridis',
                       norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
                       linewidths=.08, edgecolors='white', alpha=1, zorder=5,
                       rasterized=True)
        ax.set_title(f'{event.capitalize()} — {group}', fontsize=8, pad=7)
        ax.text(-.10, 1.06, 'abcd'[j], transform=ax.transAxes, fontsize=9, fontweight='bold')
        # All estimates use the explicitly oriented biological scale.
        ax.text(.025, .975,
                f'Effect per 10°\n{e.oriented_beta:+.2f} [{e.low:.2f}, {e.high:.2f}]',
                transform=ax.transAxes, va='top', fontsize=6,
                bbox=dict(facecolor='white', edgecolor='#b0b0b0', linewidth=.35, alpha=.97, pad=2.5), zorder=8)
        cb = fig.colorbar(sc, ax=ax, orientation='horizontal', fraction=.034, pad=.095, shrink=.90)
        cb.ax.tick_params(labelsize=6, width=.4, length=2, pad=2)
        cb.outline.set_linewidth(.5)
        cb.set_label('Latitude contribution to ' + ('offset delay' if analysis == 'offset_multivoltine'
                     else f'{event} advancement') + '\n(days per anomaly unit)', fontsize=6, labelpad=3)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / f'{name}.{ext}', dpi=600, facecolor='white')
    plt.close(fig)

draw(ANALYSES[:3], 'latitude_maps_primary_three_panels')
draw(ANALYSES, 'latitude_maps_all_four_supplement')
caption = '''Figure | Latitudinal variation in phenological sensitivity within butterfly species.
Panels show flight onset in univoltine populations (a), flight offset in univoltine populations (b), and flight offset in multivoltine populations (c). Colours represent the modelled contribution of within-species latitude to thermal sensitivity, calculated as the temperature anomaly × within-species latitude coefficient multiplied by each population's latitude deviation from its species mean (in units of 10°). Contributions are oriented so that positive values indicate greater advancement in univoltines and greater offset delay in multivoltines, relative to the response at the species mean latitude. Zero indicates no contribution from latitude, rather than no overall thermal response. Each equally sized point represents a monitoring site; colours average contributions with equal weight across the species populations included at that site. Values show only the latitude component, excluding longitude, between-species components, the baseline response and random slopes. Geographic coefficients were estimated jointly with longitude and between-species components. Inset values give the oriented latitude interaction coefficient and pointwise 95% Wald confidence interval per 10°; anomaly units are inherited from the source models and are not degrees Celsius. Offsets are conditional on annual onset. Separate symmetric colour scales span the full site-level range in each panel and should not be compared as identical magnitudes. Panels display the three within-species latitude interactions whose 95% confidence intervals exclude zero. All four analyses, including multivoltine onset, are provided in the supplementary figure. Models include a population random intercept. Maps visualise fitted linear latitude contributions and do not provide independent evidence of the gradients.
'''
(OUT / 'caption.txt').write_text(caption, encoding='utf-8')
manifest = {'primary_selection': selection, 'map_formula': 'oriented_beta_latitude * latitude_within_species_10deg',
            'site_weighting': 'equal species-population weights within each site',
            'units': 'days per source-anomaly unit; no Celsius conversion',
            'review_note': 'Automatic numerical flags retained; practical acceptance documented in selection review.',
            'source_hashes': {str(f.relative_to(RUN)): hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}}
(OUT / 'provenance.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(coefs.to_string())
print('Saved figures:', OUT)

