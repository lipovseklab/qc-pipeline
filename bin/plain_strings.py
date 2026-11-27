"""Store AnnData names and text columns as plain Python strings (anndata's default).

Some AnnData files keep cell/gene names or text columns in pandas' nullable ``string`` type. anndata only writes that
type with ``allow_write_nullable_strings``, and cNMF 1.7.1's ``consensus`` pairs genes wrongly when its gene names are
``string``: the final usage refit divides spectra labelled with those names by gene standard deviations labelled with
plain names, pandas aligns the two by sorting the genes alphabetically, and the usages come out as noise
(docs/cnmf.md). Converting to plain strings when a script reads its input avoids both.
"""
import pandas as pd


def _is_nullable_string(dtype) -> bool:
    return isinstance(dtype, pd.StringDtype)


def _plain_frame(df, prefix: str) -> list:
    """Convert a DataFrame's ``string`` index and columns in place; return what was converted."""
    changed = []
    if _is_nullable_string(df.index.dtype):
        df.index = df.index.astype(object)
        changed.append(f'{prefix}_names')
    for col in df.columns:
        dtype = df[col].dtype
        if _is_nullable_string(dtype):
            # categorical with plain-string categories: anndata writes it (also when every value is missing, which a
            # plain object column of None cannot be) and reads missing values back as missing
            values = df[col].astype(object).where(df[col].notna(), None)
            df[col] = pd.Categorical(values, categories=pd.Index(sorted(values.dropna().unique()), dtype=object))
            changed.append(f'{prefix}.{col}')
        elif isinstance(dtype, pd.CategoricalDtype) and _is_nullable_string(dtype.categories.dtype):
            df[col] = df[col].cat.rename_categories(dtype.categories.astype(object))
            changed.append(f'{prefix}.{col}')
    return changed


def plain_strings(adata) -> list:
    """Convert ``string``-typed obs/var names, columns and categories of ``adata`` (and of ``.raw``) to plain strings,
    in place. Returns the names of what was converted (empty if nothing was), e.g. ``['var_names', 'obs.sample']``.
    """
    changed = []
    if adata.raw is not None:           # edits .raw's own var table in place, without copying its matrix
        changed += [f'raw.{c}' for c in _plain_frame(adata.raw.var, 'var')]
    obs_changed = _plain_frame(adata.obs, 'obs')
    var_changed = _plain_frame(adata.var, 'var')
    # obs_names/var_names setters keep AnnData's internal state consistent
    if 'obs_names' in obs_changed:
        adata.obs_names = adata.obs.index
    if 'var_names' in var_changed:
        adata.var_names = adata.var.index
    return changed + obs_changed + var_changed


def report(path, changed) -> None:
    """Print one line saying what was converted, so the task log shows it."""
    if changed:
        shown = ', '.join(changed[:8]) + (f' and {len(changed) - 8} more' if len(changed) > 8 else '')
        print(f'{path}: stored as plain strings: {shown}')
