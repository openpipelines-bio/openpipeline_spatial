import sys
import numpy as np
import spatialdata as sd
from spatialdata.models import Labels2DModel
from spatialdata.transformations import get_transformation
from cellpose import models
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
    "input": "resources_test/xenium/xenium_tiny.zarr",
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
    "use_gpu": False,
    "normalize_percentile_low": 1.0,
    "normalize_percentile_high": 99.9,
    # Outputs
    "output": "foo.zarr",
    "output_labels": "cellpose_labels",
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

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
image_arr = np.asarray(image_element.data)

if image_arr.ndim == 3:
    channels = [par["cytoplasm_channel"], par["nuclear_channel"]]
    channel_axis = 0
else:
    channels = [0, 0]
    channel_axis = None

# Requesting GPU usage (--use_gpu) on compute resources that don't actually
# have one available would otherwise be silently downgraded deep inside
# Cellpose; resolving and logging it here makes the actual device used
# explicit in this component's own logs.
use_gpu = par["use_gpu"] and (cuda_is_available() or mps_is_available())
if par["use_gpu"] and not use_gpu:
    logger.warning("GPU requested via --use_gpu, but no GPU is available; falling back to CPU.")
logger.info(f"GPU enabled? {use_gpu}")

## Load model
if par["pretrained_model"]:
    logger.info(f"Loading custom pretrained model '{par['pretrained_model']}'...")
    model = models.CellposeModel(
        gpu=use_gpu, pretrained_model=par["pretrained_model"]
    )
else:
    logger.info(f"Loading built-in model '{par['model_type']}'...")
    model = models.CellposeModel(gpu=use_gpu, model_type=par["model_type"])

## Run segmentation
logger.info("Running Cellpose segmentation...")
masks, flows, styles = model.eval(
    image_arr,
    batch_size=par["batch_size"],
    channels=channels,
    channel_axis=channel_axis,
    diameter=par["diameter"],
    flow_threshold=par["flow_threshold"],
    cellprob_threshold=par["cellprob_threshold"],
    min_size=par["min_size"],
    niter=par["niter"],
    augment=par["augment"],
    resample=par["resample"],
    # Cellpose's own default percentile range (1, 99) can land on a
    # still-zero pixel for mostly-background fluorescence images and
    # degenerate the normalization, so the range is configurable.
    normalize={
        "percentile": [
            par["normalize_percentile_low"],
            par["normalize_percentile_high"],
        ]
    },
)

n_objects = len(np.unique(masks)) - 1
logger.info(f"Segmented {n_objects} objects.")

## Store result
out_key = par["output_labels"]
if out_key in sdata.labels:
    logger.warning(
        f"Output labels key '{out_key}' already exists and will be overwritten."
    )

logger.info(f"Storing result in .labels['{out_key}']...")
sdata.labels[out_key] = Labels2DModel.parse(
    masks.astype(np.uint32), transformations=transformations
)

## Write output
logger.info("Saving output data...")
sdata.write(par["output"], overwrite=True)

logger.info("Done!")
