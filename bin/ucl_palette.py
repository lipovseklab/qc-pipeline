"""UCL colour palettes for groups without fixed colours, shared by the QC cluster plots (bin/cluster_series.py) and
the cNMF analysis figures (bin/gep_lib/style.py). Needs only numpy, pandas, matplotlib and scikit-learn, so it loads in
every image.

The functions are the same as in our paper's `src/utricle/palettes.py`; keep the two in sync.

    ucl_palette(n)                              n different colours, UCL core colours first, colour-blind aware
    projection_palette(adata, col, basis)       those colours placed so groups touching on a projection differ most
    store_palette(adata, col, palette)          write them to adata.uns['<col>_colors'] (scanpy's convention)
"""
import itertools
import re

import matplotlib as mpl
import numpy as np
import pandas as pd

# UCL colour scheme: bright purple, orange, teal, blue, fuchsia, heritage blue, dark purple, dark orange, dark teal,
# dark blue, dark fuchsia, mid purple, light purple
UCL_BASE = ['#993bff', '#E36C2A', '#57b444', '#5487ff', '#ED367D', '#30d6ff', '#361a54', '#781c1c', '#005e5c',
            '#002ea6', '#9e1a54', '#ba82ff', '#ddbdff']   # bright purple, orange, teal, blue, fuchsia, heritage blue, ...

_AGE = re.compile(r'^([EP])(\d+)$')   # embryonic (E14) before postnatal (P0, P4, P119)


def sorted_categories(values: pd.Series) -> list:
    """Categories in natural order: numbers as numbers (2 before 10), ages E14 < P0 < P4 < P119 chronologically,
    anything else as stored. Categories no cell has are dropped (subsets keep the parent's category list)."""
    cats = list(values.cat.categories) if hasattr(values, 'cat') else sorted(values.unique())
    present = set(values.astype(str))
    cats = [c for c in cats if str(c) in present] or cats
    try:
        return sorted(cats, key=float)
    except (TypeError, ValueError):
        pass
    matches = [_AGE.match(str(c)) for c in cats]
    if all(matches):
        return [c for _, c in sorted(zip([(m.group(1) == 'P', int(m.group(2))) for m in matches], cats))]
    return cats

MIN_LIGHTNESS = 20   # CIELAB L*: leave out near-black ...
MAX_LIGHTNESS = 88   # ... and near-white shades (53 colours left)


def _srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    """(n, 3) sRGB in 0-1 -> CIELAB (D65), for measuring colour differences."""
    rgb = np.asarray(rgb, dtype=float)
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    m = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = lin @ m.T / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.stack([116 * f[:, 1] - 16, 500 * (f[:, 0] - f[:, 1]), 200 * (f[:, 1] - f[:, 2])], axis=1)


def _mix(colour: str, other: str, amount: float) -> str:
    a, b = np.array(mpl.colors.to_rgb(colour)), np.array(mpl.colors.to_rgb(other))
    return mpl.colors.to_hex(a * (1 - amount) + b * amount)


# Colour-vision deficiency simulation (Machado, Oliveira & Fernandes 2009, severity 1.0), applied in linear RGB
_CVD = {
    'deuteranopia': np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]]),
    'protanopia': np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]]),
    'tritanopia': np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]),
}


def _cvd_lab(rgb: np.ndarray) -> list:
    """CIELAB of (n, 3) sRGB colours as seen with normal vision and with each simulated colour-vision deficiency."""
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    views = [rgb]
    for m in _CVD.values():
        sim = np.clip(lin @ m.T, 0, 1)
        views.append(np.where(sim <= 0.0031308, 12.92 * sim, 1.055 * sim ** (1 / 2.4) - 0.055))
    return [_srgb_to_lab(v) for v in views]


def _ucl_pool() -> list:
    """UCL colours plus lighter and darker versions, without the darkest and palest shades."""
    pool = list(UCL_BASE)
    for amount in (0.35, 0.6):
        pool += [_mix(c, '#ffffff', amount) for c in UCL_BASE] + [_mix(c, '#000000', amount) for c in UCL_BASE]
    pool = list(dict.fromkeys(pool))
    lightness = _srgb_to_lab(np.array([mpl.colors.to_rgb(c) for c in pool]))[:, 0]
    return [c for c, l in zip(pool, lightness) if MIN_LIGHTNESS <= l <= MAX_LIGHTNESS]


def ucl_palette(n: int, min_distance: float = 12.0) -> list:
    """`n` different colours: the UCL colours first, then lighter and darker versions of them.

    Up to 13 colours are all core UCL colours; beyond that, tints and shades are added (without near-black, which
    reads as missing data, and near-white, which vanishes on white). Colours are picked one at a time, each the one
    whose smallest difference (CIELAB) from those already chosen is largest, counting a difference as the smallest of
    normal vision and simulated deuteranopia, protanopia and tritanopia, so the order suits colour-blind readers too.
    To place them on a projection so that neighbouring groups get the most different colours, use
    `projection_palette`.
    """
    base = list(UCL_BASE)
    tints = [c for c in _ucl_pool() if c not in base]
    candidates = base + tints
    if n > len(candidates):
        raise ValueError(f'Cannot make {n} distinct colours from the UCL pool ({len(candidates)} available)')
    views = _cvd_lab(np.array([mpl.colors.to_rgb(c) for c in candidates]))
    chosen = [0]
    while len(chosen) < n:
        dist = np.min([np.min(np.linalg.norm(v[:, None, :] - v[chosen][None, :, :], axis=2), axis=1) for v in views],
                      axis=0)
        dist[chosen] = -1
        if len(chosen) < len(base):          # core UCL colours first
            dist[len(base):] = -1
        chosen.append(int(np.argmax(dist)))
    if n > 1:
        lab = views[0]
        closest = min(np.linalg.norm(lab[i] - lab[j]) for i, j in itertools.combinations(chosen, 2))
        if closest < min_distance:
            print(f'WARNING: {n} colours requested; the closest two are only {closest:.1f} apart in CIELAB')
    return [candidates[i] for i in chosen]


def projection_palette(adata, col: str, basis: str, n_neighbors: int = 15) -> dict:
    """UCL colours (ucl_palette) for the groups in `adata.obs[col]`, placed so that groups that touch on the 2-D
    projection `adata.obsm[basis]` get the most different colours (also for colour-blind readers).

    Groups touch when a cell's nearest neighbours on the projection belong to the other group; the more such pairs,
    the more their colours must differ. Groups are coloured most-connected first, each getting the free colour that is
    most different from the colours of the groups it touches (ties: core UCL colours first). Same data, same colours.
    Returns {group: colour}; store it with `store_palette` so every figure uses it.
    """
    from sklearn.neighbors import NearestNeighbors

    key = basis if basis in adata.obsm else f'X_{basis}'
    xy = np.asarray(adata.obsm[key])[:, :2]
    cats = [str(c) for c in sorted_categories(adata.obs[col])]
    index = {c: i for i, c in enumerate(cats)}
    group = np.array([index[str(v)] for v in adata.obs[col]])
    n = len(cats)

    # touching: neighbour pairs across groups, relative to group sizes
    k = min(n_neighbors + 1, len(xy))
    _, nbrs = NearestNeighbors(n_neighbors=k).fit(xy).kneighbors(xy)
    touch = np.zeros((n, n))
    np.add.at(touch, (np.repeat(group, k - 1), group[nbrs[:, 1:]].ravel()), 1)
    np.fill_diagonal(touch, 0)
    size = np.bincount(group, minlength=n).astype(float)
    touch = (touch + touch.T) / np.sqrt(np.outer(size, size))
    weight = touch / touch.max() if touch.max() > 0 else touch

    colours = ucl_palette(n)
    views = _cvd_lab(np.array([mpl.colors.to_rgb(c) for c in colours]))
    diff = np.min([np.linalg.norm(v[:, None, :] - v[None, :, :], axis=2) for v in views], axis=0)

    assigned = {}                                    # group -> colour index
    for g in np.argsort(-weight.sum(axis=1), kind='stable'):
        free = [c for c in range(n) if c not in assigned.values()]
        neighbours = [(h, weight[g, h]) for h in assigned if weight[g, h] > 0]
        if not neighbours:
            assigned[g] = free[0]
            continue
        # most different from touching groups; stronger contact demands a bigger difference
        score = [min(diff[c, assigned[h]] / w for h, w in neighbours) for c in free]
        assigned[g] = free[int(np.argmax(score))]

    # then swap colours between groups while that makes the closest touching pair more different (ties: the contact-
    # weighted mean difference). Pairs with only a few contacts (< 2% of the strongest) don't count as touching.
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n) if weight[a, b] > 0.02]
    if pairs:
        colour_of = np.array([assigned[g] for g in range(n)])

        def quality(col):
            d = np.array([diff[col[a], col[b]] for a, b in pairs])
            return (round(float(d.min()), 6), float(np.sum(d * np.array([weight[a, b] for a, b in pairs]))))

        best = quality(colour_of)
        improved = True
        while improved:
            improved = False
            for a, b in itertools.combinations(range(n), 2):
                trial = colour_of.copy()
                trial[a], trial[b] = trial[b], trial[a]
                q = quality(trial)
                if q > best:
                    colour_of, best, improved = trial, q, True
        assigned = {g: int(colour_of[g]) for g in range(n)}
    return {cats[g]: colours[c] for g, c in assigned.items()}


def store_palette(adata, col: str, palette: dict) -> None:
    """Store {group: colour} as `adata.uns['<col>_colors']` (scanpy's convention), with the column's categories in
    natural order, so scanpy plots and the heatmaps here all use it."""
    cats = [str(c) for c in sorted_categories(adata.obs[col])]
    adata.obs[col] = pd.Categorical(adata.obs[col].astype(str), categories=cats)
    adata.uns[f'{col}_colors'] = [palette[c] for c in cats]
