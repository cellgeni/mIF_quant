# Introduction

AKOYA pipeline automates low-level spatial proteomics processing from a multichannel image: preprocessing, segmentation (StarDist, or an imported pre-computed mask - see [Imported segmentation mask](#imported-segmentation-mask-optional)), per-cell intensity extraction, and positive/negative marker calls. It is built mostly using [spatialproteomics_cellgeni](https://github.com/cellgeni/spatialproteomics), which was adapted for large images (original package: [spatialproteomics](https://github.com/sagar87/spatialproteomics)). Outputs can include (i) a spatialproteomics object in Zarr format; (ii) CSV tables containing cell positions, per-cell intensity summaries, and marker-presence calls; and (iii) an AnnData file. `adata.X` remains the transformed mean-intensity matrix for backwards compatibility, while the additional per-cell statistics are stored in `adata.obs`. See [Output structure](#output-structure) for a description of every output file and **open_output_files.ipynb** for examples of how to open them.

# Environment

Create Docker container using [Dockerfile](https://github.com/cellgeni/AKOYA_analysis/blob/main/container/Dockerfile). If you run it on farm, please use singularity container: `/nfs/cellgeni/singularity/images/spatialproteomics_nd2.sif`

# Run whole pipeline

## Prepare configuration file

All parameters used in pipeline together with input/output paths should be specified in one configuration file (see as example **conf_AKOYA.yaml**). Below is a description of all parameters of the configuration file.

### Input image and channels

`image_path (str)` - Path to the multiplex image file (e.g. .tif, .qptiff, .ndpi, .nd2) to be processed

`list_of_channels (list[str])` - Channel names in the same order as the channel axis of the loaded image

`channel_for_segmentation (str)` - Which channel to use for StarDist segmentation. Typically the DAPI channel name. Ignored when `segmentation_mask_path` is set (the key must still be present in the config).

`list_of_markers (list[str])` - Marker channels used to compute % positive pixels and generate binary labels / cell-type-like labels via thresholding.

### Cropping

`crop_x (list[int, int])` - X-range to crop as [x_start, x_end].

`crop_y (list[int, int])` - Y-range to crop as [y_start, y_end].

Use `null` for both to process the whole image. If `segmentation_mask_path` is set, the same crop is applied to the mask.

### Segmentation + filtering

`segmentation_label_expansion (int or falsy)` - If truthy, expands segmentation labels by this radius (pixels) using expand_segmentation

`stardist_scale (float, optional)` - Image scale passed to StarDist (default: `3`, preserving the original pipeline behaviour). For cells that are split into multiple objects, try a smaller value such as `1`, `0.75`, or `0.5`; this reduces the apparent cell size seen by the model. This option changes segmentation itself, whereas `segmentation_label_expansion` only expands masks after segmentation. Ignored when `segmentation_mask_path` is set.

`min_area (int or falsy)` - If truthy, filters out segmented objects with area <= min_area.

`max_area (int or falsy)` - If truthy, filters out segmented objects with area >= max_area.

`null`, `None`, an empty value or `False` all disable `min_area` / `max_area` / `segmentation_label_expansion`. The area filters are also applied to an imported mask, so set them to `null` if the mask is already final.

### Imported segmentation mask (optional)

Instead of running StarDist, you can provide your own segmentation (e.g. from Cellpose, QuPath, Mesmer, or a manually curated mask). When `segmentation_mask_path` is set, the StarDist step is skipped and the mask is added to the spatialproteomics object as its `_segmentation` layer (with cell IDs and centroids), so every later step - label expansion, area filtering, intensity statistics, positive/negative calls, plots and all outputs - runs exactly as with StarDist.

`segmentation_mask_path (str or null, optional)` - Path to a label image: `0` = background, `1..N` = one integer ID per cell. Supported formats: `.tif` / `.tiff` / `.ome.tif`, `.npy`, `.png`. Default `null` (use StarDist).

`segmentation_mask_type (str, optional)` - How the imported mask is interpreted. Default `nuclei`.
- `nuclei` - the mask is treated like StarDist output: `segmentation_label_expansion` is applied to grow nuclei into cells, and `split_signal_nuclei_cytoplasm` works as usual.
- `whole_cell` - the mask is used as the final cell outlines: `segmentation_label_expansion` is skipped and `split_signal_nuclei_cytoplasm` is switched off (a warning is printed), because a whole-cell mask carries no nucleus information.

Requirements and behaviour:
- The mask must be 2D and have **the same Y/X size as the original (uncropped) image**; otherwise the pipeline stops with an error. `crop_x` / `crop_y` are then applied to the mask in the same way as to the image.
- Labels must be non-negative integers. They are renumbered to `1..N` internally, so cell IDs in the outputs may differ from the IDs in your mask file.
- A binary mask (only 0/1) is accepted and converted to labels by connected components, but touching cells will be merged into one object - a proper label mask is preferred.
- The log prints `Skipping StarDist; importing <type> mask from <path>` and `Imported segmentation mask: N objects`, which you can use to confirm the mask was used.

Example (Cellpose whole-cell mask, no area filtering):

```
segmentation_mask_path: "/path/to/cellpose_output/image_cp_masks.tif"
segmentation_mask_type: whole_cell
segmentation_label_expansion: 0
min_area: null
max_area: null
```

### Intermediate plots

`save_intermediate_plots (bool)` - Whether to save intermediate QC plots and ROI snapshots.

`number_intermediate_plots (int)` - Number of random ROIs (subimages) to sample for intermediate plotting.

`size_intermediate_plots (int)` - ROI size in pixels (square). Each ROI is size_intermediate_plots × size_intermediate_plots.

`list_of_genes_intermediate_plots (list[str])` - Channel names to display in intermediate ROI plots (alongside DAPI).

Marker-vs-celltype plots saved by `save_individual_marker_presence_plots` include raw marker expression before preprocessing, processed marker expression after thresholding/filtering, and the called positive cells.

### Binary marker presence / label thresholding

`save_individual_marker_presence_plots (bool)` - Whether to plot marker-vs-celltype (binary label) maps.

`fraction_of_positive_pixels (float)` - Backwards-compatible default threshold applied to the fraction of non-zero pixels. It is used for markers/compartments without a matching entry in `positive_cell_rules`.

`manual_marker_thresholds (dict[str, float], optional)` - Per-channel overrides for image preprocessing thresholds, for example `{CD8: 1200, "PD-L1": 0.995}`. Channels not listed here keep the automatic Otsu threshold. Values below 1 are interpreted as relative channel quantiles; values greater than or equal to 1 are interpreted as absolute intensity thresholds.

`split_signal_nuclei_cytoplasm (bool, optional)` - If True, also quantifies nuclei and, when `segmentation_label_expansion` is set, cytoplasm. `whole_cell` means the expanded/combined cell mask (nucleus plus cytoplasm). Whole-cell mean intensity remains in `adata.X`; backwards-compatible transformed compartment means are written to `adata.layers["nuclei_intensity"]` and `adata.layers["cytoplasm_intensity"]`. Not available with an imported mask of `segmentation_mask_type: whole_cell`.

### Per-cell intensity statistics

Statistics are computed from processed pixels after image thresholding and median filtering. For every image channel and each available compartment (`whole_cell`, plus `nuclei` and `cytoplasm` when enabled), the pipeline writes these raw processed-pixel statistics to `adata.obs`:

- `mean`
- `median`
- `std` (population standard deviation, `ddof=0`)
- `variance` (population variance, `ddof=0`)
- `percentage_positive` (fraction of pixels greater than zero)
- every configured quantile

`intensity_quantiles (list[float], optional)` - Quantiles to calculate. Values must be between 0 and 1. The default is `[0.5, 0.75, 0.9, 0.95]`.

Columns use the unambiguous form `<channel>__<compartment>__<metric>`. For example:

```
CD8__whole_cell__mean
CD8__whole_cell__quantile_0.9
CD8__nuclei__median
CD8__cytoplasm__variance
```

The same matrices are exported to Zarr/CSV with names such as `_intensity_median`, `_intensity_median_nuclei`, and `_intensity_quantile_0_9_cytoplasm`.

### Configurable positive/negative calls

`positive_cell_rules (dict, optional)` - Selects the source statistic and threshold independently for every marker. With `{}` (the default), all calls use `fraction_of_positive_pixels`, preserving the previous behaviour. Supported operators are `>=` (default), `>`, `<=`, and `<`.

A rule containing `metric` applies to all available compartments:

```
positive_cell_rules:
  CD8:
    metric: quantile_0.9
    threshold: 25
  CD163:
    metric: median
    operator: ">"
    threshold: 12
```

Rules can differ between whole cell, nucleus, and cytoplasm:

```
positive_cell_rules:
  PD-L1:
    whole_cell: {metric: median, threshold: 12}
    nuclei: {metric: mean, threshold: 8}
    cytoplasm: {metric: percentage_positive, threshold: 0.3}
```

You may also reference an exact `adata.obs` source column. `{marker}`, `{channel}`, and `{compartment}` placeholders are expanded automatically. Existing observation columns such as `area` can also be used:

```
positive_cell_rules:
  CD8:
    column: "CD8__{compartment}__std"
    threshold: 4
  PD-L1:
    column: area
    operator: ">="
    threshold: 150
```

Calls are stored as integer 0/1 columns named `<marker>__<compartment>__positive` in `adata.obs`. They are also available as DataFrames in `adata.obsm["marker_presence_whole_cell"]`, `adata.obsm["marker_presence_nuclei"]`, and `adata.obsm["marker_presence_cytoplasm"]`. The resolved source column, operator, and threshold are recorded in `adata.uns["akoya_positive_cell_rules"]`. When no custom rules are supplied, the historical `*_ct` fields are retained as well.

### Output control

`output_dir (str)` - Path to directory where results and plots are written.

`list_output_formats (list[str])` - Which outputs to save. Supported values in this script: ["zarr", "h5ad", "csv"]

`save_intermediate_zarr (bool)` - If True, saves intermediate sp_object snapshots to `output_dir/sp_object.zarr` after key steps.

`save_omero_segmentation_csv (bool, optional)` - If True, writes OMERO-compatible polygon tables to `output_dir/omero_segmentation_csv/`. One CSV is generated for `whole_cell` and, when their masks are available, `nuclei` and `cytoplasm`. Each file begins with the OMERO fields `object`, `label`, `score`, `confidence_score`, and `polygon`; `polygon` contains a closed WKT string in the form `POLYGON ((x y, ...))`. All `adata.obs` columns are appended, together with one categorical `label-<marker>` column per marker containing `positive` or `negative` for that compartment.

### Optional metadata

`normalise_intensity (bool)` - If True, performs z-score normalization per channel across cells and stores the result alongside the raw intensities.

`pixelsize (float or null, optional)` - Microns per full-resolution pixel. If provided, it is used only to create anndata h5ad object

## Run pipeline

The pipeline depending on image size requires significant amount of memory, so it is recommended for full-tissue crop (with ~(20k x 20k) pixels image and 60 channels) to use 200 Gb of memory or more. Example of submission code can be found in [submit_AKOYA_pipeline.sh](https://github.com/cellgeni/AKOYA_analysis/blob/main/templates/submit_AKOYA_pipeline.sh). Then one can submit a job simply as:

`bsub < submit_AKOYA_pipeline.sh`

## Run pipeline for many samples

In case you want to run the pipeline for many samples, the most annoying part is to prepare all separate configuration files. If you want to keep all parameters the same, you literally need to change in each conf file only "image_path" and "output_dir". You can do it automatically (if you have one example of conf file with all parameters tuned) using notebook [prepare_all_conf](https://github.com/cellgeni/AKOYA_analysis/blob/main/templates/prepare_all_conf.ipynb). And then you can run all of them by submitting number of jobs with [submit_all_AKOYA_jobs.sh](https://github.com/cellgeni/AKOYA_analysis/blob/main/templates/submit_all_AKOYA_jobs.sh)

# Output structure

All results are written to `output_dir` (set in the configuration file). Which files appear depends on the options you enable: `list_output_formats`, `save_intermediate_zarr`, `save_intermediate_plots`, `save_individual_marker_presence_plots`, `save_omero_segmentation_csv`, `split_signal_nuclei_cytoplasm` and `segmentation_label_expansion`. Files that exist only under certain options are marked as *conditional*.

```
output_dir/
├── <sample>.zarr/                    # spatialproteomics object            [list_output_formats: "zarr"]
├── sp_object.zarr/                   # intermediate snapshots of the object [save_intermediate_zarr: True]
├── <sample>.h5ad                     # AnnData object                      [list_output_formats: "h5ad"]
├── <sample>_*.csv                    # CSV tables                          [list_output_formats: "csv"]
├── omero_segmentation_csv/           # OMERO polygon tables                [save_omero_segmentation_csv: True]
│   ├── whole_cell.csv
│   ├── nuclei.csv                    # only if nuclei masks are available
│   └── cytoplasm.csv                 # only if cytoplasm masks are available
└── <plots>                           # QC and marker-presence plots        [save_intermediate_plots / save_individual_marker_presence_plots]
```

<!-- TODO(maintainers): confirm exact file/folder names of the final Zarr, h5ad, CSV and plot outputs against bin/run_AKOYA.py and replace the <placeholders> above. -->

Job logs are not written to `output_dir`: with the provided [submit_AKOYA_pipeline.sh](https://github.com/cellgeni/AKOYA_analysis/blob/main/templates/submit_AKOYA_pipeline.sh) they go to `logs/output<JOBID>.log` (stdout) and `logs/error<JOBID>.log` (stderr) in the directory where you ran `bsub`.

## Description of output files

### Spatialproteomics object (Zarr)

*Enabled by `"zarr"` in `list_output_formats`.* The complete spatialproteomics object, which can be opened with the `spatialproteomics` package (see **open_output_files.ipynb**). It contains the (cropped) multichannel image, the segmentation mask, and the per-cell tables listed below. Intensity matrices are stored with names such as:

| Name pattern | Meaning |
|---|---|
| `_intensity_<metric>` | Per-cell statistic over the **whole cell** (`metric` = `mean`, `median`, `std`, `variance`, `percentage_positive`) |
| `_intensity_<metric>_nuclei` | Same statistic computed only over the **nucleus** (*conditional*: `split_signal_nuclei_cytoplasm: True`) |
| `_intensity_<metric>_cytoplasm` | Same statistic computed only over the **cytoplasm** (*conditional*: `split_signal_nuclei_cytoplasm: True` and `segmentation_label_expansion` set) |
| `_intensity_quantile_<q>[_nuclei/_cytoplasm]` | Quantile `q` of processed pixel intensities per cell, with `.` replaced by `_` (e.g. `_intensity_quantile_0_9_cytoplasm` is the 0.9 quantile). Quantiles come from `intensity_quantiles` (default 0.5, 0.75, 0.9, 0.95) |

Statistics are computed from processed pixels, i.e. after image thresholding and median filtering.

### Intermediate Zarr snapshots (`sp_object.zarr`)

*Conditional: `save_intermediate_zarr: True`.* Snapshots of the spatialproteomics object saved after key pipeline steps (e.g. after preprocessing and segmentation). Useful for debugging or for resuming analysis from an intermediate step without re-running the whole pipeline. They are overwritten at each step and are not the final result.

### AnnData file (h5ad)

*Enabled by `"h5ad"` in `list_output_formats`.* One row per segmented cell, one column per image channel.

| Location | Content |
|---|---|
| `adata.X` | Transformed **mean-intensity** matrix over the whole cell (kept for backwards compatibility) |
| `adata.layers["nuclei_intensity"]` | Transformed mean intensity over nuclei (*conditional*: `split_signal_nuclei_cytoplasm: True`) |
| `adata.layers["cytoplasm_intensity"]` | Transformed mean intensity over cytoplasm (*conditional*: `split_signal_nuclei_cytoplasm: True` and `segmentation_label_expansion` set) |
| `adata.obs["area"]` and other cell properties | Cell geometry and position-related fields (e.g. area, centroid coordinates) |
| `adata.obs["<channel>__<compartment>__<metric>"]` | Per-cell statistics. `compartment` is `whole_cell`, `nuclei` or `cytoplasm`; `metric` is `mean`, `median`, `std`, `variance`, `percentage_positive` or `quantile_<q>` (e.g. `CD8__whole_cell__quantile_0.9`) |
| `adata.obs["<marker>__<compartment>__positive"]` | Integer 0/1 positive/negative call for each marker and compartment, based on `positive_cell_rules` (or `fraction_of_positive_pixels` if no rule is given) |
| `adata.obs["*_ct"]` | Historical binary label fields, retained when no custom `positive_cell_rules` are supplied |
| `adata.obsm["marker_presence_whole_cell"]`, `["marker_presence_nuclei"]`, `["marker_presence_cytoplasm"]` | The same positive/negative calls as DataFrames (cells × markers), one per compartment |
| `adata.uns["akoya_positive_cell_rules"]` | The resolved rule used for every marker: source column, operator and threshold |
| z-score normalised intensities | Per-channel z-score across cells (*conditional*: `normalise_intensity: True`) |
| Pixel size | Microns per pixel is stored only if `pixelsize` is set in the config |

### CSV tables

*Enabled by `"csv"` in `list_output_formats`.* The same per-cell information as plain tables that can be opened without Python:

- **cell positions** - one row per cell with its label/ID and coordinates (and area);
- **per-cell intensity summaries** - one table per statistic and compartment, named like the Zarr matrices (`_intensity_median`, `_intensity_median_nuclei`, `_intensity_quantile_0_9_cytoplasm`, ...), cells × channels;
- **marker-presence calls** - cells × markers tables of 0/1 positive/negative calls.

### OMERO segmentation tables (`omero_segmentation_csv/`)

*Conditional: `save_omero_segmentation_csv: True`.* Segmentation outlines that can be imported into OMERO as ROIs/annotations. One file per compartment: `whole_cell`, and `nuclei` / `cytoplasm` when those masks are available.

| Column | Meaning |
|---|---|
| `object` | Object type in OMERO |
| `label` | Cell ID (label in the segmentation mask) |
| `score`, `confidence_score` | OMERO score fields (placeholders required by the OMERO import format) |
| `polygon` | Closed outline of the cell as WKT: `POLYGON ((x y, ...))` |
| all `adata.obs` columns | Appended after the OMERO fields, so every statistic and call is available per cell |
| `label-<marker>` | One categorical column per marker with value `positive` or `negative` for that compartment |

### Plots

- **Intermediate QC plots and ROI snapshots** - *conditional: `save_intermediate_plots: True`.* `number_intermediate_plots` randomly selected square regions (`size_intermediate_plots` × `size_intermediate_plots` pixels) showing DAPI together with the channels in `list_of_genes_intermediate_plots`, used to check preprocessing and segmentation quality.
- **Marker-vs-celltype plots** - *conditional: `save_individual_marker_presence_plots: True`.* For each marker: the raw marker signal before preprocessing, the processed signal after thresholding/filtering, and the cells called positive. Use them to judge whether the thresholds (`manual_marker_thresholds`, `positive_cell_rules`, `fraction_of_positive_pixels`) are appropriate.

# Run separate steps of the pipeline

To run separately steps from the pipeline (such as image preprocessing, segmentation or intensity extraction) please use as an example notebook [AKOYA_analysis_steps](https://github.com/cellgeni/AKOYA_analysis/blob/main/AKOYA_analysis_steps.ipynb). Please note, that there we use only some of all available from [spatialproteomics](https://github.com/sagar87/spatialproteomics), if you want to find out more about other options of image preprocessing, segmentation, plotting and celltyping please visit [spatialproteomics documentation](https://sagar87.github.io/spatialproteomics/)
