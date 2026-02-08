from __future__ import annotations

from pathlib import Path

import streamlit as st
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import io
from datetime import datetime
import os, json, hashlib
import requests


import sesame as ssm

PREFERRED_ORDER = [
    "Atmosphere",
    "Lithosphere",
    "Hydrosphere",
    "Biosphere",
    "Cryosphere"
]

@st.cache_data(show_spinner=False)
def list_spheres(atlas_root: str) -> list[str]:
    root = Path(atlas_root).expanduser()
    if not root.exists():
        return []
    spheres = [p.name for p in root.iterdir() if p.is_dir()]
    return sorted(spheres, key=lambda x: PREFERRED_ORDER.index(x) if x in PREFERRED_ORDER else 999)


@st.cache_data(show_spinner=False)
def list_nc_in_sphere(atlas_root: str, sphere: str) -> list[str]:
    root = Path(atlas_root).expanduser() / sphere
    if not root.exists():
        return []
    return sorted([str(p) for p in root.rglob("*.nc")])

@st.cache_data(show_spinner=False)
def open_dataset(path: str) -> xr.Dataset:
    return xr.open_dataset(path)

def pick_variable(ds: xr.Dataset) -> list[str]:
    return sorted(list(ds.data_vars.keys()))

def get_attr(da, key: str, default: str = "") -> str:
    v = da.attrs.get(key, None)
    if v is None:
        return default
    return str(v)

CACHE_DIR = Path("./atlas_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

@st.cache_data(show_spinner=False)
def load_manifest(url: str) -> dict:
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.json()

def cache_path(url: str) -> Path:
    h = hashlib.sha256(url.encode("utf-8")).hexdigest()[:24]
    return CACHE_DIR / f"{h}__{url.split('/')[-1]}"

def download_if_needed(url: str) -> Path:
    local = cache_path(url)
    if local.exists() and local.stat().st_size > 0:
        return local

    tmp = local.with_suffix(local.suffix + ".part")
    with st.spinner(f"Downloading {local.name} ..."):
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
    tmp.replace(local)
    return local


def get_time_options(ds: xr.Dataset) -> tuple[str | None, list[str]]:
    """
    Returns (time_dim_name, time_values_as_strings_for_UI).
    Datetime coordinates are displayed as YYYY only.
    """
    for cand in ["time", "year", "month", "day"]:
        if cand in ds.dims or cand in ds.coords:
            if cand in ds.coords:
                vals = ds[cand].values
                # Datetime handling
                if np.issubdtype(vals.dtype, np.datetime64):
                    years = [str(np.datetime64(v, "Y").astype(int) + 1970) for v in vals]
                    return cand, years
                else:
                    return cand, [str(v) for v in vals]
            else:
                # dim without coord
                return cand, [str(i) for i in range(ds.dims[cand])]
    return None, []


def default_last_coord_value(ds: xr.Dataset, dim: str):
    """Return the last coordinate value (preferred) or last index."""
    if dim in ds.coords:
        vals = ds[dim].values
        if len(vals) == 0:
            return None
        return vals[-1]
    # fallback: integer index (string) if no coord values
    n = ds.sizes.get(dim, 0)
    return (n - 1) if n > 0 else None

def find_dim(ds, candidates):
    for c in candidates:
        if c in ds.dims or c in ds.coords:
            return c
    return None

def fig_to_bytes(fig, fmt: str, dpi: int = 300) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.read()

def parse_levels(s: str):
    s = (s or "").strip()
    if s == "":
        return 10  # default
    # Comma list
    if "," in s:
        parts = [p.strip() for p in s.split(",") if p.strip() != ""]
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            return "INVALID"
        if len(vals) < 2:
            return "INVALID"
        # Must be strictly increasing for clean bins
        if any(vals[i] >= vals[i+1] for i in range(len(vals) - 1)):
            return "INVALID_ORDER"
        return vals
    # Integer
    try:
        n = int(float(s))  # allows "10" or "10.0"
    except ValueError:
        return "INVALID"
    if n < 2:
        return "INVALID"
    return n



st.set_page_config(page_title="The SESAME Human-Earth Atlas Explorer", layout="wide")
st.title("The SESAME Human-Earth Atlas Explorer")


# Sidebar
st.sidebar.header("Data")

MANIFEST_URL = (
    st.secrets.get("ATLAS_MANIFEST_URL", None)
    or os.environ.get("ATLAS_MANIFEST_URL", None)
)

if not MANIFEST_URL:
    st.error("Atlas manifest URL is not configured.")
    st.stop()

manifest = load_manifest(MANIFEST_URL)

spheres = sorted(manifest.keys())
if not spheres:
    st.error("Manifest loaded but contains no spheres.")
    st.stop()

sphere = st.sidebar.selectbox("Sphere", spheres)

urls = manifest.get(sphere, [])
if not urls:
    st.warning(f"No files listed for {sphere}")
    st.stop()

filtered = urls

# If filenames can repeat across subfolders, this mapping can collide.
display_names = [Path(u).name for u in filtered]

file_choice = st.sidebar.selectbox(
    "NetCDF file (.nc)",
    display_names
)

file_url = filtered[display_names.index(file_choice)]

local_nc = download_if_needed(file_url)
file_path = str(local_nc)  # keep compatibility with the rest of your code

ds = open_dataset(file_path)  # reuse your cached open_dataset() wrapper


EXCLUDE_VARS = {"grid_area", "land_frac"}
vars_ = [v for v in pick_variable(ds) if v not in EXCLUDE_VARS]
if not vars_:
    st.error("This NetCDF contains no data variables.")
    st.stop()

var = st.sidebar.selectbox("Variable", vars_)

LAT_DIM = ["lat", "latitude", "y"]
LON_DIM = ["lon", "longitude", "x"]
TIME_DIM = ["time", "year", "month", "day"]
DEPTH_DIM = ["depth", "lev", "level", "z"]

lat_dim = find_dim(ds, LAT_DIM)
lon_dim = find_dim(ds, LON_DIM)
time_dim = find_dim(ds, TIME_DIM)
depth_dim = find_dim(ds, DEPTH_DIM)

da = ds[var]
extra_dims = [d for d in da.dims if d not in {lat_dim, lon_dim} and d is not None]

# Time selection (optional)
time_dim, time_opts = get_time_options(ds)

time_choice = None
time_idx = None

if time_dim and time_opts:
    default_time = str(default_last_coord_value(ds, time_dim))
    time_choice = st.sidebar.selectbox(
        f"Time slice ({time_dim})",
        time_opts,
        index=(time_opts.index(default_time) if default_time in time_opts else len(time_opts) - 1),
        key=f"{Path(file_path).stem}__{var}__time"
    )
    time_idx = time_opts.index(time_choice)
    
# Optional multi-time aggregation
# Optional multi-time aggregation (auto-enabled when user selects times)
time_agg_op = "mean"
time_multi_choices = []

if time_dim and time_opts:
    with st.sidebar.expander("Time aggregation", expanded=False):
        time_multi_choices = st.multiselect(
            f"Select {time_dim} values",
            options=time_opts,
            default=[],
            key=f"{Path(file_path).stem}__{var}__time_multi"
        )

        time_agg_op = st.selectbox(
            "Operation",
            ["mean", "sum", "min", "max", "median", "std"],
            index=0,
            key=f"{Path(file_path).stem}__{var}__time_agg_op"
        )

# Auto-enabled if user chose >=1 entries
time_agg_enabled = len(time_multi_choices) > 0


depth_choice = None
depth_idx = None

if depth_dim and depth_dim in da.dims:
    depth_vals = [str(v) for v in (ds[depth_dim].values if depth_dim in ds.coords else range(ds.sizes[depth_dim]))]
    depth_choice = st.sidebar.selectbox(
        f"Depth slice ({depth_dim})",
        depth_vals,
        index=0,
        key=f"{Path(file_path).stem}__{var}__depth"
    )
    depth_idx = depth_vals.index(depth_choice)  

da0 = ds[var]
title = get_attr(da0, "long_name", default=var)
label = get_attr(da0, "units", default="")

st.sidebar.markdown("---")
st.sidebar.subheader("Data Access")

st.sidebar.markdown(
"""
Faisal, A. A., Kaye, M., Ahmed, M. & Galbraith, E.  
*The SESAME Human-Earth Atlas*.  
figshare (2025).  
https://doi.org/10.6084/m9.figshare.28432499
"""
)

st.sidebar.subheader("Paper")
st.sidebar.markdown(
"""
Faisal, A.A., Kaye, M., Ahmed, M. et al. (2025).
The SESAME Human-Earth Atlas. Scientific Data, 12, 775. 
https://doi.org/10.1038/s41597-025-05087-5
"""
)

st.sidebar.subheader("Code")
st.sidebar.markdown(
"""
SESAME GitHub Repository  
https://github.com/A2Faisal/SESAME
"""
)


# Layout
left, right = st.columns([2, 1], gap="large")

with right:
    st.subheader("Metadata")
    st.write("**Sphere:**", sphere)
    st.write("**File:**", Path(file_path).name)
    st.write("**Variable:**", var)
    st.write("**Dims:**", dict(ds.sizes))
    if var in ds and getattr(ds[var], "attrs", None):
        st.write("**Variable attributes:**")
        st.json({k: str(v) for k, v in ds[var].attrs.items()})

with left:
    st.subheader("Map")

    # ---- Plot controls ----
    with st.expander("Plot options", expanded=False):
        k = lambda name: f"{Path(file_path).stem}__{var}__{name}"  # reset when file/var changes

        cmap = st.text_input("Colormap", value="Blues", key=k("cmap"))
        st.caption(
        "See valid Matplotlib colormaps: "
        "[matplotlib.org colormap reference](https://matplotlib.org/stable/users/explain/colors/colormaps.html)"
        )

        levels_txt = st.text_input(
            "Levels (int or comma-separated list)",
            value="10",
            key=k("levels_txt"),
            help="Examples: 10  |  0, 5, 10, 20, 50"
        )

        remove_ata = st.checkbox("Remove Antarctica", value=False, key=k("remove_ata"))
        out_bound = st.checkbox("Show outer boundary", value=True, key=k("out_bound"))

        st.markdown("#### Value limits (optional)")
        vmin_txt = st.text_input("vmin", value="", key=k("vmin_txt"))
        vmax_txt = st.text_input("vmax", value="", key=k("vmax_txt"))

        # Arrows are only used if limits are set
        extend_max_ui = st.checkbox("Extend above vmax (arrow)", value=True, key=k("extend_max_ui"))
        extend_min_ui = st.checkbox("Extend below vmin (arrow)", value=True, key=k("extend_min_ui"))

    def _parse_float(s: str):
        s = (s or "").strip()
        if s == "":
            return None
        try:
            return float(s)
        except ValueError:
            return "INVALID"

    vmin_arg = _parse_float(vmin_txt)
    vmax_arg = _parse_float(vmax_txt)

    if vmin_arg == "INVALID":
        st.warning("vmin is not a valid number. Leave it blank or enter a numeric value.")
        vmin_arg = None

    if vmax_arg == "INVALID":
        st.warning("vmax is not a valid number. Leave it blank or enter a numeric value.")
        vmax_arg = None

    # Auto-enable arrows only when limits are set
    extend_max_arg = bool(extend_max_ui) if vmax_arg is not None else False
    extend_min_arg = bool(extend_min_ui) if vmin_arg is not None else False

    # Prevent vmin > vmax errors
    if (vmin_arg is not None) and (vmax_arg is not None) and (vmin_arg > vmax_arg):
        st.warning("vmin was greater than vmax. Swapping them to avoid a plotting error.")
        vmin_arg, vmax_arg = vmax_arg, vmin_arg
    
    levels_arg = parse_levels(levels_txt)

    if levels_arg == "INVALID":
        st.warning("Levels is invalid. Use an integer >= 2 or a comma-separated numeric list.")
        levels_arg = 10

    if levels_arg == "INVALID_ORDER":
        st.warning("Levels list must be strictly increasing, like: 0, 5, 10, 20.")
        levels_arg = 10

    fig = None

    try:
        ds_plot = ds

        # ---- Time: aggregate if user selected multiple, else single time ----
        if time_dim and time_dim in ds_plot[var].dims:

            if time_agg_enabled:
                time_indices = [time_opts.index(t) for t in time_multi_choices]
                time_indices = sorted(set(time_indices))  # stable + de-dup
            else:
                time_indices = [time_idx] if time_idx is not None else []

            if time_indices:
                ds_plot = ds_plot.isel({time_dim: time_indices})

                if len(time_indices) > 1:
                    da_tmp = ds_plot[var]

                    if time_agg_op == "mean":
                        da_out = da_tmp.mean(dim=time_dim, skipna=True)
                    elif time_agg_op == "sum":
                        da_out = da_tmp.sum(dim=time_dim, skipna=True)
                    elif time_agg_op == "min":
                        da_out = da_tmp.min(dim=time_dim, skipna=True)
                    elif time_agg_op == "max":
                        da_out = da_tmp.max(dim=time_dim, skipna=True)
                    elif time_agg_op == "median":
                        da_out = da_tmp.median(dim=time_dim, skipna=True)
                    elif time_agg_op == "std":
                        da_out = da_tmp.std(dim=time_dim, skipna=True)
                    else:
                        da_out = da_tmp.mean(dim=time_dim, skipna=True)

                    ds_plot[var] = da_out.squeeze(drop=True)
                else:
                    ds_plot[var] = ds_plot[var].squeeze(drop=True)

        if depth_dim and depth_idx is not None and depth_dim in ds_plot[var].dims:
            if depth_dim in ds_plot.coords:
                selected_depth = ds_plot[depth_dim].values[depth_idx]
                ds_plot = ds_plot.sel({depth_dim: selected_depth})
            else:
                ds_plot = ds_plot.isel({depth_dim: depth_idx})

        plt.close("all")
        ssm.plot_map(
            dataset=ds_plot,
            variable=var,
            title=title,
            label=label,
            color=cmap,
            levels=levels_arg,
            vmin=vmin_arg,
            vmax=vmax_arg,
            extend_min=extend_min_arg,
            extend_max=extend_max_arg,
            remove_ata=remove_ata,
            out_bound=out_bound,
            show=False
        )

        fig = plt.gcf()
        st.pyplot(fig, clear_figure=False)

    except Exception as e:
        st.error(f"plot_map failed: {e}")
        st.info("Fallback: showing xarray plot instead.")

        da_f = ds[var]

        # ---- Time: aggregate if user selected multiple, else single time (fallback) ----
        if time_dim and time_dim in da_f.dims:

            if time_agg_enabled:
                time_indices = [time_opts.index(t) for t in time_multi_choices]
                time_indices = sorted(set(time_indices))
            else:
                time_indices = [time_idx] if time_idx is not None else []

            if time_indices:
                da_f = da_f.isel({time_dim: time_indices})

                if len(time_indices) > 1:
                    if time_agg_op == "mean":
                        da_f = da_f.mean(dim=time_dim, skipna=True)
                    elif time_agg_op == "sum":
                        da_f = da_f.sum(dim=time_dim, skipna=True)
                    elif time_agg_op == "min":
                        da_f = da_f.min(dim=time_dim, skipna=True)
                    elif time_agg_op == "max":
                        da_f = da_f.max(dim=time_dim, skipna=True)
                    elif time_agg_op == "median":
                        da_f = da_f.median(dim=time_dim, skipna=True)
                    elif time_agg_op == "std":
                        da_f = da_f.std(dim=time_dim, skipna=True)
                    else:
                        da_f = da_f.mean(dim=time_dim, skipna=True)

                da_f = da_f.squeeze(drop=True)

            # ---- Download UI ----
            st.markdown("### Save figure")

            fmt = st.selectbox("Format", ["png", "pdf", "svg"], index=0, key="save_fmt")
            dpi = st.number_input("DPI (PNG only)", min_value=72, max_value=600, value=300, step=25, key="save_dpi")

            safe_sphere = sphere.replace(" ", "_")
            safe_file = Path(file_path).stem.replace(" ", "_")
            safe_var = var.replace(" ", "_")
            timestamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
            filename = f"SESAME_{safe_sphere}_{safe_file}_{safe_var}_{timestamp}.{fmt}"

            if fig is not None:
                if fmt in ["pdf", "svg"]:
                    file_bytes = fig_to_bytes(fig, fmt=fmt, dpi=300)
                    mime = "application/pdf" if fmt == "pdf" else "image/svg+xml"
                else:
                    file_bytes = fig_to_bytes(fig, fmt="png", dpi=int(dpi))
                    mime = "image/png"

                st.download_button(
                    label=f"Download {fmt.upper()}",
                    data=file_bytes,
                    file_name=filename,
                    mime=mime,
                    use_container_width=True
                )

st.markdown(
    """
    <hr style="margin-top:2em;margin-bottom:1em;">
    <small>
    © The Surface Earth System Analysis and Modeling Environment (SESAME) Project.<br>
    Developed and maintained by the 
    <a href="https://earthsystemdynamics.org/" target="_blank">
    Earth System Dynamics Lab
    </a>.<br>
    <br>
    The data and visualizations provided here are for research and educational purposes.
    While every effort has been made to ensure accuracy, the authors and affiliated
    institutions assume no responsibility for errors, omissions, or interpretations
    arising from the use of these data.
    </small>
    """,
    unsafe_allow_html=True
)


