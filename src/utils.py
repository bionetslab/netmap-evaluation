
import yaml
import numpy as np
import pandas as pd


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)


def sanitize_nullable_dtypes(adata):
    """Downcast pandas nullable dtype columns in obs/var to plain numpy dtypes.

    anndata >=0.11 can write pandas' nullable "string"/"Int64"/"boolean"
    extension-array columns (encoding-type 'nullable-string-array' etc.), but
    anndata <0.11 has no registered reader for them at all -- reading such a
    file raises IORegistryError. LINGER pins anndata==0.9.2 (see
    linger_env/pixi.toml), so any h5ad written by this repo's main env (which
    uses a much newer anndata) and later read by a LINGER script must not
    contain these dtypes. Categorical columns are left untouched -- those are
    a normal, always-supported anndata encoding.
    """
    for df in (adata.obs, adata.var):
        for col in df.columns:
            dtype = df[col].dtype
            if getattr(dtype, "na_value", None) is pd.NA:
                if pd.api.types.is_string_dtype(dtype):
                    df[col] = df[col].astype(object)
                else:
                    df[col] = df[col].to_numpy(dtype=float, na_value=np.nan)
        # obs_names/var_names themselves can carry the same nullable dtype
        # (e.g. barcodes/gene names re-inferred as pandas "string" dtype by a
        # newer pandas on read), not just regular columns.
        if getattr(df.index.dtype, "na_value", None) is pd.NA:
            df.index = df.index.astype(object)
    return adata


def split_index(aa):    
    aa.var['source']   = [l[0] for l in aa.var.index.str.split('_', expand=True)]
    aa.var['target']   = [l[1] for l in aa.var.index.str.split('_', expand=True)]
    return aa