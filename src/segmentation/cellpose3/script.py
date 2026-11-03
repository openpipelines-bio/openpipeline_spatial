import math
import sys
import numpy as np
import spatialdata as sd
from scipy.ndimage import find_objects
from spatialdata.models import Labels2DModel
from spatialdata.transformations import get_transformation
from cellpose import models, transforms
from torch.cuda import is_available as cuda_is_available
from xarray import DataTree


try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False


## VIASH START
par = {
    # Inputs
    "input": "resources_test/xenium_multichannel/xenium_multicellseg_tiny.zarr",
    "input_image": "morphology_focus",
    "cytoplasm_channel": 0,
    "nuclear_channel": 0,
    # Parameters
    "model_type": "nuclei",
    "pretrained_model": None,
    "diameter": None,
    "flow_threshold": 0.4,
    "cellprob_threshold": 0.0,
    "min_size": 15,
    "niter": None,
    "resample": True,
    "augment": False,
    "batch_size": 8,
    "device_type": "cpu",
    "normalize_percentile_low": 1.0,
    "normalize_percentile_high": 99.9,
    "tile_size": 4096,
    "tile_overlap": 128,
    # Outputs
    "output": "foo.zarr",
    "output_labels": "cellpose_labels",
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def global_normalization_bounds(
    image_arr, channels, channel_axis, nchan, percentiles, max_pixels=4_000_000
):
    """Intensity bounds (per model input channel) for the whole image.

    Computed once on an evenly spaced subsample of the full image, so every
    tile is rescaled with the same values. Cellpose's channel selection
    (`convert_image`) is applied to the subsample first, so the bounds are
    computed on exactly the channels the model sees.
    """
    height, width = image_arr.shape[-2:]
    stride = max(1, math.ceil(math.sqrt(height * width / max_pixels)))
    sample = np.asarray(image_arr[..., ::stride, ::stride])
    sample = transforms.convert_image(
        sample, channels, channel_axis=channel_axis, nchan=nchan
    )
    bounds = []
    for c in range(sample.shape[-1]):
        low, high = np.percentile(sample[..., c], percentiles)
        if high <= low:
            # Constant channel (e.g. the all-zero nuclear channel Cellpose
            # adds when none is selected): leave its values unchanged, like
            # Cellpose's own normalization does, instead of dividing by 0.
            low, high = 0.0, 1.0
        bounds.append([low, high])
    return np.array(bounds, dtype=np.float32)


## Read data
logger.info("Reading input data...")
sdata = sd.read_zarr(par["input"])

## Validate inputs
image_key = par["input_image"]
if image_key not in sdata.images:
    raise ValueError(f"Image key '{image_key}' not found in .images.")

if par["normalize_percentile_low"] >= par["normalize_percentile_high"]:
    raise ValueError(
        "'--normalize_percentile_low' "
        f"({par['normalize_percentile_low']}) must be lower than "
        f"'--normalize_percentile_high' ({par['normalize_percentile_high']})."
    )

image_element = sdata.images[image_key]
transformations = get_transformation(image_element, get_all=True)

# Multiscale images (DataTree) expose full-resolution pixel data under the
# "scale0" node; single-scale images expose it directly via .data.
if isinstance(image_element, DataTree):
    image_element = image_element["scale0"]["image"]
# Keep the image lazy (a dask array): a whole-slide image does not fit in
# memory, so only one tile at a time is loaded further down.
image_arr = image_element.data

if image_arr.ndim == 3:
    channels = [par["cytoplasm_channel"], par["nuclear_channel"]]
    channel_axis = 0
    # Cellpose itself either fails with an IndexError that is offset by one
    # from the user-facing 1-based numbering, or (for a single-channel image)
    # silently segments the grayscale image instead of the requested channel.
    n_channels = image_arr.shape[channel_axis]
    for arg_name, channel in zip(["cytoplasm_channel", "nuclear_channel"], channels):
        if channel > n_channels:
            raise ValueError(
                f"'--{arg_name}' ({channel}) is out of range for image "
                f"'{image_key}' with {n_channels} channel(s); valid values are "
                f"0 to {n_channels} (1-based, 0 = none/grayscale)."
            )
else:
    # A 2D image has no channel axis to select from, so a non-zero
    # '--cytoplasm_channel'/'--nuclear_channel' would otherwise be silently
    # dropped here.
    if par["cytoplasm_channel"] != 0 or par["nuclear_channel"] != 0:
        logger.warning(
            "'--cytoplasm_channel' and/or '--nuclear_channel' were set to a "
            f"non-zero value, but the input image '{image_key}' has no "
            f"channel axis (ndim={image_arr.ndim}); ignoring the requested "
            "channel(s)."
        )
    channels = [0, 0]
    channel_axis = None

# Requesting GPU usage (--device_type) on compute resources that don't actually
# have one available would otherwise be silently downgraded deep inside
# Cellpose; resolving and logging it here makes the actual device used
# explicit in this component's own logs.
use_gpu = par["device_type"] == "gpu" and (cuda_is_available() or mps_is_available())
if par["device_type"] == "gpu" and not use_gpu:
    logger.warning(
        "GPU requested via --device_type, but no GPU is available; falling back to CPU."
    )
logger.info(f"GPU enabled? {use_gpu}")

## Load model
if par["pretrained_model"]:
    logger.info(f"Loading custom pretrained model '{par['pretrained_model']}'...")
    model = models.CellposeModel(gpu=use_gpu, pretrained_model=par["pretrained_model"])
else:
    logger.info(f"Loading built-in model '{par['model_type']}'...")
    model = models.CellposeModel(gpu=use_gpu, model_type=par["model_type"])


tile_size = par["tile_size"]
overlap = par["tile_overlap"]
height, width = image_arr.shape[-2:]
tile_origins = [
    (y0, x0) for y0 in range(0, height, tile_size) for x0 in range(0, width, tile_size)
]
percentiles = [par["normalize_percentile_low"], par["normalize_percentile_high"]]

if len(tile_origins) == 1:
    # The image fits in a single tile: let Cellpose compute the normalization
    # percentiles on the image itself (identical to segmenting without tiling).
    normalize = {"percentile": percentiles}
else:
    # Per-tile percentiles would rescale every tile differently (a tile that is
    # mostly background would get its noise stretched into "cells"), so the
    # intensity bounds are computed once for the whole image instead.
    logger.info("Computing normalization bounds on the whole image...")
    lowhigh = global_normalization_bounds(
        image_arr, channels, channel_axis, model.nchan, percentiles
    )
    logger.info(f"Normalization bounds per model input channel: {lowhigh.tolist()}")
    normalize = {"lowhigh": lowhigh}

eval_kwargs = {
    "batch_size": par["batch_size"],
    "channels": channels,
    "channel_axis": channel_axis,
    "diameter": par["diameter"],
    "flow_threshold": par["flow_threshold"],
    "cellprob_threshold": par["cellprob_threshold"],
    "min_size": par["min_size"],
    "niter": par["niter"],
    "augment": par["augment"],
    "resample": par["resample"],
    "normalize": normalize,
}

## Run segmentation
logger.info(
    f"Running Cellpose segmentation on {len(tile_origins)} tile(s) of at most "
    f"{tile_size}x{tile_size} pixels (+{overlap} pixels overlap)..."
)
masks = np.zeros((height, width), dtype=np.uint32)
next_label = 1
for tile_idx, (y0, x0) in enumerate(tile_origins, start=1):
    # Each tile has a "core" (y0:y1, x0:x1); the tiles' cores together cover
    # the image without overlapping. The model is run on the core plus an
    # overlap margin, so that cells crossing the core's border are seen in
    # full. A cell is only kept by the tile whose core contains the centre of
    # the cell's bounding box, so no cell is kept twice.
    y1, x1 = min(y0 + tile_size, height), min(x0 + tile_size, width)
    ty0, tx0 = max(y0 - overlap, 0), max(x0 - overlap, 0)
    ty1, tx1 = min(y1 + overlap, height), min(x1 + overlap, width)

    logger.info(
        f"Segmenting tile {tile_idx}/{len(tile_origins)} "
        f"(y {ty0}-{ty1}, x {tx0}-{tx1})..."
    )
    tile = np.asarray(image_arr[..., ty0:ty1, tx0:tx1])
    tile_masks, _, _ = model.eval(tile, **eval_kwargs)

    for label, bbox in enumerate(find_objects(tile_masks), start=1):
        if bbox is None:
            continue
        cy = ty0 + (bbox[0].start + bbox[0].stop) // 2
        cx = tx0 + (bbox[1].start + bbox[1].stop) // 2
        if not (y0 <= cy < y1 and x0 <= cx < x1):
            continue
        target = masks[
            ty0 + bbox[0].start : ty0 + bbox[0].stop,
            tx0 + bbox[1].start : tx0 + bbox[1].stop,
        ]
        # Pixels already claimed by a cell from a neighbouring tile are kept.
        target[(tile_masks[bbox] == label) & (target == 0)] = next_label
        next_label += 1

# Counting via np.unique would sort the full-size mask array.
n_objects = next_label - 1
logger.info(f"Segmented {n_objects} objects.")

## Store result
out_key = par["output_labels"]
if out_key in sdata.labels:
    logger.warning(
        f"Output labels key '{out_key}' already exists and will be overwritten."
    )

logger.info(f"Storing result in .labels['{out_key}']...")
sdata.labels[out_key] = Labels2DModel.parse(masks, transformations=transformations)

## Write output
logger.info("Saving output data...")
sdata.write(par["output"], overwrite=True)

logger.info("Done!")
