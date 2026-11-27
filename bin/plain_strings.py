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


def plain_strings(adata) -> list:
    """Convert ``string``-typed obs/var names, columns and categories of ``adata`` to plain strings, in place.

    Returns the names of what was converted (empty if nothing was), e.g. ``['var_names', 'obs.sample']``.
    """
    changed = []
    if adata.raw is not None:                       # .raw has its own gene names and var table
        raw = adata.raw.to_adata()
        raw_changed = [f'raw.{c}' for c in plain_strings(raw) if not c.startswith('obs')]
        if raw_changed:
            adata.raw = raw
            changed += raw_changed
    for axis in ('obs', 'var'):
        names = getattr(adata, f'{axis}_names')
        if _is_nullable_string(names.dtype):
            setattr(adata, f'{axis}_names', names.astype(object))
            changed.append(f'{axis}_names')
        df = getattr(adata, axis)
        for col in df.columns:
            dtype = df[col].dtype
            if _is_nullable_string(dtype):
                df[col] = df[col].astype(object).where(df[col].notna(), None)
                changed.append(f'{axis}.{col}')
            elif isinstance(dtype, pd.CategoricalDtype) and _is_nullable_string(dtype.categories.dtype):
                df[col] = df[col].cat.rename_categories(dtype.categories.astype(object))
                changed.append(f'{axis}.{col}')
    return changed


def report(path, changed) -> None:
    """Print one line saying what was converted, so the task log shows it."""
    if changed:
        shown = ', '.join(changed[:8]) + (f' and {len(changed) - 8} more' if len(changed) > 8 else '')
        print(f'{path}: stored as plain strings: {shown}')
