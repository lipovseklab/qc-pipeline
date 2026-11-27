"""Plot style: UCL colours, palettes with no repeated colours, and publication-ready matplotlib settings."""
from __future__ import annotations

from typing import Dict, Iterable

import matplotlib as mpl
import numpy as np
import scanpy as sc
import seaborn as sb

# palette maths shared with the QC cluster plots and our paper notebooks: bin/ucl_palette.py (bin/ is on sys.path)
from ucl_palette import projection_palette, store_palette, ucl_palette  # noqa: F401

# --- UCL colour scheme ------------------------------------------------------------------------------------------
# Core
DARKPUR = '#361a54'
BRIGHTPUR = '#993bff'
MIDPUR = '#ba82ff'
# Supporting
LIGHTPUR = '#ddbdff'
PALEPUR = '#eedeff'
HERITAGEBLU = '#30d6ff'
# Graph colours, dark
DARKBLU = '#002ea6'
DARKORA = '#781c1c'
DARKFUS = '#9e1a54'
DARKTEA = '#005e5c'
# Graph colours, mid
BLU = '#5487ff'
ORA = '#E36C2A'
FUS = '#ED367D'
TEA = '#57b444'

# Ordered so that neighbours in the list are easy to tell apart
UCL_BASE = [BRIGHTPUR, ORA, TEA, BLU, FUS, HERITAGEBLU, DARKPUR, DARKORA, DARKTEA, DARKBLU, DARKFUS, MIDPUR, LIGHTPUR]


MM = 1 / 25.4                 # inches per millimetre
SINGLE_COLUMN = 85 * MM       # figure widths used by most journals
DOUBLE_COLUMN = 180 * MM


def categorical_palette(categories: Iterable[str]) -> Dict[str, str]:
    """Map each category to a different UCL-derived colour (same input order gives the same colours)."""
    categories = list(categories)
    return dict(zip(categories, ucl_palette(len(categories))))


def age_palette(categories: Iterable[str]) -> Dict[str, str]:
    """The notebook's age colours: seaborn's 'husl' palette with one colour per age (categories in age order)."""
    categories = list(categories)
    return dict(zip(categories, sb.color_palette('husl', len(categories)).as_hex()))


def point_size(n_cells: int, width_in: float) -> float:
    """Scatter point area (pt^2) for a panel `width_in` inches wide. Scanpy's default (120000 / n_cells) is meant for a
    4 inch panel, so it is scaled by panel area; floored so points stay visible for very large datasets."""
    return max(0.3, 120000 / max(n_cells, 1) * (width_in / 4) ** 2)


def usage_cmap():
    """Grey to UCL bright purple, for GEP usage and gene expression."""
    return sb.color_palette(f'blend:gainsboro,{BRIGHTPUR}', as_cmap=True)


def set_publication_style(dpi: int = 300) -> None:
    """Journal-style defaults: small sans-serif text, thin lines, editable text in PDF/SVG, rasterised scatter points."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.titlesize': 8,
        'axes.labelsize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'pdf.fonttype': 42,      # real (editable) text in PDFs
        'svg.fonttype': 'none',  # real text in SVGs
        'savefig.dpi': dpi,
        'savefig.transparent': False,
        'savefig.bbox': 'tight',
        'figure.dpi': 100,
    })
    # scanpy scatter points are drawn as an image inside vector files, so PDFs stay small
    sc.set_figure_params(vector_friendly=True, dpi_save=dpi, frameon=False, fontsize=7)
    # scanpy's set_figure_params changes these again, so set them last (a stray grid line showed in the colour bars)
    mpl.rcParams['font.size'] = 7
    mpl.rcParams['axes.grid'] = False
