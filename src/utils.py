
import yaml
import numpy as np
import pandas as pd


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)


def _is_nullable_dtype(dtype):
    return getattr(dtype, "na_value", None) is pd.NA


def sanitize_nullable_dtypes(adata):
    """Downcast pandas nullable dtypes in obs/var to plain numpy dtypes.

    anndata >=0.11 can write pandas' nullable "string"/"Int64"/"boolean"
    extension-array dtypes (encoding-type 'nullable-string-array' etc.), but
    anndata <0.11 has no registered reader for them at all -- reading such a
    file raises IORegistryError. LINGER pins anndata==0.9.2 (see
    linger_env/pixi.toml), so any h5ad written by this repo's main env (which
    uses a much newer anndata/pandas) and later read by a LINGER script must
    not contain these dtypes anywhere. Three places they can hide:
      - a regular obs/var column with a nullable dtype
      - obs_names/var_names (the DataFrame index) itself
      - a categorical column whose *categories* array has a nullable dtype
        (categoricals are otherwise a normal, always-supported encoding, so
        only the categories array needs downcasting, not the column itself)
    """
    for df in (adata.obs, adata.var):
        for col in df.columns:
            series = df[col]
            dtype = series.dtype
            if isinstance(dtype, pd.CategoricalDtype):
                if _is_nullable_dtype(dtype.categories.dtype):
                    df[col] = series.cat.rename_categories(series.cat.categories.astype(object))
                continue
            if _is_nullable_dtype(dtype):
                if pd.api.types.is_string_dtype(dtype):
                    df[col] = series.astype(object)
                else:
                    df[col] = series.to_numpy(dtype=float, na_value=np.nan)
        if _is_nullable_dtype(df.index.dtype):
            df.index = df.index.astype(object)
    return adata


def split_index(aa):    
    aa.var['source']   = [l[0] for l in aa.var.index.str.split('_', expand=True)]
    aa.var['target']   = [l[1] for l in aa.var.index.str.split('_', expand=True)]
    return aa