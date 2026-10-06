# -*- coding: utf-8 -*-
"""
phenoIMPACT | Manuscript figure: geographic patterns in phenological plasticity
Combined 3-panel layout: each panel contains map (top) + latitude plot (bottom)

Requested design changes:
1) Garamond or closest available fallback.
2) Closed-box scatter panels.
3) Combined column layout: a = onset, b = offset-univoltine, c = offset-multivoltine.
4) Larger text.
5) Remove '(0 = species mean)' from the colourbar label; explain it in the caption.

Run:
    python make_latitude_longitude_figure_combined.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib import font_manager, gridspec
from mpl_toolkits.basemap import Basemap

ROOT = Path(r"E:/phenoIMPACT project/code/phenoIMPACT")
LAT_RUN = (
    ROOT / 'output' / 'phenology_plasticity' /
    'latitude_from_extracted_plasticity_15km' /
    'run_20261005_124834_165c6dd85b83'
)
LONG_RESULTS = ROOT / 'data' / 'longitude_sensitivity_results.csv'
OUT = LAT_RUN / 'manuscript_figure_python_combined'
OUT.mkdir(parents=True, exist_ok=True)

font_candidates = ['Garamond', 'Cormorant Garamond', 'EB Garamond', 'DejaVu Serif', 'serif']
font_family = None
for cand in font_candidates:
    try:
        font_manager.findfont(cand, fallback_to_default=False)
        font_family = cand
        break
    except Exception:
        pass
if font_family is None:
    font_family = 'serif'

plt.rcParams.update({
    'font.family': font_family,
    'font.size': 13,
    'axes.titlesize': 19,
    'axes.labelsize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'axes.linewidth': 0.9,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

coord_files = sorted((LAT_RUN / 'inputs' / 'coordinates').glob('*.csv'))
coords_all = pd.concat([
    pd.read_csv(f, dtype={'SITE_ID': str, 'bms_id': str}) for f in coord_files
], ignore_index=True)
coords = coords_all[['SITE_ID', 'bms_id', 'longitude_deg', 'latitude_deg']].copy()
chk = coords.groupby(['SITE_ID', 'bms_id'])[['longitude_deg', 'latitude_deg']].agg(
    lambda x: x.max() - x.min()
)
if not (chk.to_numpy() < 1e-7).all():
    raise RuntimeError('Coordinate exports disagree.')
coords = coords.drop_duplicates(['SITE_ID', 'bms_id'])

lr = pd.read_csv(LONG_RESULTS)

analyses = ['onset', 'offset_univoltine', 'offset_multivoltine']
titles = {
    'onset': 'Onset',
    'offset_univoltine': 'Offset - univoltine',
    'offset_multivoltine': 'Offset - multivoltine',
}
response_labels = {
    'onset': 'Relative onset advancement',
    'offset_univoltine': 'Relative offset advancement',
    'offset_multivoltine': 'Relative offset delay',
}
orient = {'onset': -1.0, 'offset_univoltine': -1.0, 'offset_multivoltine': 1.0}

fig = plt.figure(figsize=(16.8, 11.2))
outer = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.14,
                          left=0.04, right=0.985, top=0.95, bottom=0.07)

site_outputs = []
partial_outputs = []

for j, a in enumerate(analyses):
    gs_col = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[j],
        height_ratios=[1.08, 0.92], hspace=0.12
    )

    p = pd.read_csv(
        LAT_RUN / a / 'population_plasticity_with_latitude.csv',
        dtype={'SPECIES': str, 'SITE_ID': str, 'bms_id': str}
    )
    cj = coords.rename(columns={'longitude_deg': 'longitude_deg_coord', 'latitude_deg': 'latitude_deg_coord'})
    p = p.merge(cj, on=['SITE_ID', 'bms_id'], how='left', validate='many_to_one')
    if p['longitude_deg_coord'].isna().any():
        raise RuntimeError(f'Missing coordinates for {a}.')
    if not np.allclose(
        p['latitude_deg'].to_numpy(),
        p['latitude_deg_coord'].to_numpy(),
        atol=1e-7, rtol=0
    ):
        raise RuntimeError(f'Latitude mismatch for {a}.')
    p['longitude_deg'] = p['longitude_deg_coord']

    mean_lon = p.groupby('SPECIES')['longitude_deg'].transform('mean')
    lon_w = (p['longitude_deg'] - mean_lon).to_numpy() / 10.0
    lat_w = p['latitude_within_species_deg'].to_numpy() / 10.0
    response_w = orient[a] * p['plasticity_within_species_raw'].to_numpy()

    sub = lr[(lr['analysis'] == a) &
             (lr['weighting'] == 'equal_population') &
             (lr['specification'] == 'latitude_plus_longitude_within_species')]
    latrow = sub[sub['term'] == 'latitude'].iloc[0]
    lonrow = sub[sub['term'] == 'longitude'].iloc[0]

    b_lat = orient[a] * float(latrow['estimate_raw_per_10deg'])
    lat_ci = sorted(orient[a] * np.array([float(latrow['lower_95_raw']), float(latrow['upper_95_raw'])]))
    b_lon = orient[a] * float(lonrow['estimate_raw_per_10deg'])

    p['map_response'] = response_w - b_lon * lon_w
    site = p.groupby(['SITE_ID', 'bms_id', 'longitude_deg', 'latitude_deg'], as_index=False).agg(
        relative_response=('map_response', 'mean'),
        n_species=('SPECIES', 'nunique')
    )
    site['analysis'] = a
    site_outputs.append(site)

    ax_map = fig.add_subplot(gs_col[0])
    m = Basemap(
        projection='merc', llcrnrlon=-12, urcrnrlon=33,
        llcrnrlat=38, urcrnrlat=66, resolution='l', ax=ax_map
    )
    m.drawmapboundary(linewidth=0.85)
    m.drawcoastlines(linewidth=0.48)
    m.drawcountries(linewidth=0.38)
    m.drawparallels(np.arange(40, 66, 5), labels=[1, 0, 0, 0],
                    linewidth=0.22, dashes=[1, 0], fontsize=10)
    m.drawmeridians(np.arange(-10, 31, 10), labels=[0, 0, 0, 1],
                    linewidth=0.22, dashes=[1, 0], fontsize=10)

    x, y = m(site['longitude_deg'].to_numpy(), site['latitude_deg'].to_numpy())
    vals = site['relative_response'].to_numpy()
    lim = float(np.nanquantile(np.abs(vals), 0.98))
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    sizes = 10 + 2.35 * np.sqrt(site['n_species'].to_numpy())

    sc = m.scatter(x, y, c=vals, s=sizes, cmap='viridis', norm=norm,
                   alpha=0.84, linewidths=0, zorder=4)

    ax_map.set_title(titles[a], pad=10, fontweight='semibold')
    ax_map.text(-0.12, 1.07, 'abc'[j], transform=ax_map.transAxes,
                fontsize=16, fontweight='bold', ha='left', va='top', clip_on=False)

    cb = fig.colorbar(sc, ax=ax_map, fraction=0.040, pad=0.018)
    cb.ax.tick_params(labelsize=10, width=0.7, length=3.2)
    cb.outline.set_linewidth(0.7)
    cb.set_label(response_labels[a], fontsize=11.5, labelpad=7)
    cb.ax.text(1.75, 1.01, 'More', transform=cb.ax.transAxes,
               ha='center', va='bottom', fontsize=10)
    cb.ax.text(1.75, -0.03, 'Less', transform=cb.ax.transAxes,
               ha='center', va='top', fontsize=10)

    ax_sc = fig.add_subplot(gs_col[1])
    denom_lon = np.sum(lon_w ** 2)
    lat_res = lat_w - lon_w * (np.sum(lat_w * lon_w) / denom_lon)
    resp_res = response_w - lon_w * (np.sum(response_w * lon_w) / denom_lon)

    partial = pd.DataFrame({
        'analysis': a,
        'SPECIES': p['SPECIES'],
        'SITE_ID': p['SITE_ID'],
        'latitude_within_species_deg': lat_res * 10,
        'relative_response': resp_res,
    })
    partial_outputs.append(partial)

    ax_sc.scatter(partial['latitude_within_species_deg'], partial['relative_response'],
                  s=5.5, alpha=0.12, edgecolors='none', rasterized=True)
    for side in ['top', 'right', 'left', 'bottom']:
        ax_sc.spines[side].set_visible(True)
        ax_sc.spines[side].set_linewidth(0.9)
    ax_sc.axhline(0, linewidth=0.7, linestyle=(0, (2.5, 2.5)), alpha=0.65)
    ax_sc.axvline(0, linewidth=0.7, linestyle=(0, (2.5, 2.5)), alpha=0.65)

    xlo = float(np.nanquantile(partial['latitude_within_species_deg'], 0.005))
    xhi = float(np.nanquantile(partial['latitude_within_species_deg'], 0.995))
    xx = np.linspace(xlo, xhi, 300)
    yy = b_lat * xx / 10.0
    low = lat_ci[0] * xx / 10.0
    high = lat_ci[1] * xx / 10.0
    ax_sc.fill_between(xx, np.minimum(low, high), np.maximum(low, high), alpha=0.18, linewidth=0)
    ax_sc.plot(xx, yy, linewidth=2.05)

    ax_sc.set_xlabel('Latitude within species (°)')
    ax_sc.set_ylabel(response_labels[a])
    ax_sc.tick_params(direction='out', width=0.75, length=3.6)

    slope_text = f"+{b_lat:.2f} [{lat_ci[0]:.2f}, {lat_ci[1]:.2f}] per 10°"
    ax_sc.text(0.03, 0.965, slope_text,
               transform=ax_sc.transAxes, ha='left', va='top', fontsize=11.2,
               bbox=dict(boxstyle='round,pad=0.22', facecolor='white',
                         edgecolor='0.55', linewidth=0.75, alpha=0.95))

pd.concat(site_outputs, ignore_index=True).to_csv(OUT / 'figure_site_values.csv', index=False)
pd.concat(partial_outputs, ignore_index=True).to_csv(OUT / 'figure_partial_values.csv', index=False)

fig.savefig(OUT / 'manuscript_latitude_longitude_combined.png', dpi=420,
            bbox_inches='tight', facecolor='white')
fig.savefig(OUT / 'manuscript_latitude_longitude_combined.pdf',
            bbox_inches='tight', facecolor='white')
print('Saved to:', OUT)
