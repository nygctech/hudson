from scipy.ndimage import shift
from skimage.io import imread,imsave
from skimage.registration import phase_cross_correlation
from skimage.exposure import rescale_intensity
from utils import get_logger,HiSeqImage
import os.path as path
import yaml
import zarr
import os
from os import makedirs,system
import psutil
import subprocess
import numpy as np
import tempfile
from skimage.io import imsave, imread
import importlib.util
from scipy.stats import mode
import sys
import torch
from typing import List,Tuple
sys.path.append(os.path.abspath("/gpfs/commons/home/ecordina/ecordina_innovation/report/standard_comparison/final_figure_src/src/"))
import URNet
#######################################
# Utility Functions
#######################################

def find_translation(ref: np.ndarray, src: np.ndarray) -> Tuple[float, float]:
    """
    Compute the translation needed to align `src` with `ref` using phase cross-correlation.

    Args:
        ref (np.ndarray): Reference image.
        src (np.ndarray): Source image to align.

    Returns:
        Tuple[float, float]: Translation offsets in the y and x directions.
    """
    shift_values, _, _ = phase_cross_correlation(
        ref, src, space='real', disambiguate=True, normalization="phase"
    )
    return shift_values

def compute_crop_region(translations: List[np.ndarray], shape: Tuple) -> Tuple[slice, slice]:
    """Computes the largest valid region after alignment."""
    h, w = shape[-2:]  # Get image dimensions
    x_shifts, y_shifts = zip(*translations)

    min_x, max_x = int(np.floor(min(x_shifts))), int(np.ceil(max(x_shifts)))
    min_y, max_y = int(np.floor(min(y_shifts))), int(np.ceil(max(y_shifts)))

    crop_x_start, crop_x_end = max_x, h + min_x  # Exclude extreme shifts
    crop_y_start, crop_y_end = max_y, w + min_y  # Exclude extreme shifts

    return slice(crop_y_start, crop_y_end), slice(crop_x_start, crop_x_end)

def get_immersion_index(input_string: str) -> float | str:
    """
    Determine the immersion index based on the provided input string.

    Args:
        input_string (str): A string describing the objective or medium, expected to contain 
                            either "Air" or "Oil".

    Returns:
        float: The corresponding immersion index for "Air" (1.0) or "Oil" (1.515).
        str: "Unknown medium" if neither "Air" nor "Oil" is found in the input string.
    """
    immersion_indices = {
        "Air": 1.0,
        "Oil": 1.515
    }
    if "Air" in input_string:
        return immersion_indices["Air"]
    elif "Oil" in input_string:
        return immersion_indices["Oil"]
    else:
        return "Unknown medium"

def generate_psf_file(
    resxy: float,
    resz: float,
    na: float,
    ni: float,
    wavelength: int,
    output_path: str,
    module_command: str = "bash -c 'source /etc/profile && module load deconwolf &&"
) -> None:
    """
    Generates a PSF file using the `dw_bw` command-line tool.

    Args:
        resxy (float): Resolution in XY (nm).
        resz (float): Resolution in Z (nm).
        na (float): Numerical aperture of the lens.
        ni (float): Immersion index.
        wavelength (int): Channel or wavelength used.
        output_path (str): Path where the PSF TIFF file will be saved.
        module_command (str): Shell command to load required module.
    """
    command = (
        f"{module_command} "
        f"dw_bw --resxy {resxy} --resz {resz} "
        f"--NA {na} --ni {ni} --lambda {wavelength} {output_path}'"
    )
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"PSF generation failed: {e}")

def get_psf(
    channel: int,
    temp_folder: str,
    config_path: str = "/commons/instruments/pyseq/machine_settings.yaml"
) -> str:
    """
    Retrieves or generates a PSF (Point Spread Function) file for a given channel.

    Args:
        channel (int): The wavelength or channel value used to parameterize the PSF.
        temp_folder (str): Path to the temporary directory for storing/generated PSF files.
        config_path (str): Path to the YAML config file for microscope settings.

    Returns:
        str: Path to the generated or existing PSF TIFF file.
    """
    psf_path = path.join(temp_folder, f"PSF/PSF_{channel}.tiff")
    if path.isfile(psf_path):
        return psf_path

    try:
        with open(config_path) as stream:
            machines_configs = yaml.safe_load(stream)
    except (yaml.YAMLError, FileNotFoundError) as exc:
        raise RuntimeError(f"Error loading machine config: {exc}")

    objective = machines_configs["Origin"]['Objective']
    pixels = machines_configs["Origin"]['Pixels']

    na = objective['lens_na']
    ni = get_immersion_index(objective['model'])
    resxy = 0.6 * 1000  # nm

    resz = pixels.get('physical_size_z', 0.9) * 1000  # nm, fallback for HiSeq2500

    makedirs(path.join(temp_folder, "PSF"), exist_ok=True)
    fname = path.join(temp_folder, f"PSF/PSF_{resxy}_{resz}_{na}_{ni}_{channel}.tiff")

    if not path.isfile(fname):
        generate_psf_file(resxy, resz, na, ni, channel, fname)

    return fname


def batch_run_deconwolf(image_4d: np.ndarray, psfs: list[str]|str|None = None, n_iter: int = 100, psigma: float = 0.0, gpu: bool = True, tilesize: int = 0,tilepad: int = 0) -> np.ndarray:
    """
    Deconvolves each channel and round of a 4D image using Deconwolf.

    Parameters:
    - image_4d (np.ndarray): 5D array with shape (rounds,channels, y, x)
    - psfs (List[str])|str: list of PSFs, one per channel, if one psf is given, will use the same one for every channels and round, if none, will generate a psf with deconwolf
    - n_iter (int): number of iterations
    - psigma (float): psigma for Deconwolf
    - gpu (bool): whether to use GPU

    Returns:
    - np.ndarray: deconvolved image of same shape as input
    """
    R,C, Y, X = image_4d.shape
    output = np.zeros_like(image_4d)

        
    with tempfile.TemporaryDirectory() as temp_dir:

        if type(psfs) == str:
            psfs = [psfs] * C
        elif psfs is None:
            psfs = [get_psf(channel,temp_dir) for channel in (["587","610","687","740"]*C)[:C]]
        if len(psfs) != C:
            raise ValueError("Number of PSFs must match number of channels.")
            
        input_dir = os.path.join(temp_dir, "input")
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        bash_lines = ["#!/bin/bash", "source /etc/profile", "module load deconwolf"]

        threads_count = psutil.cpu_count() / psutil.cpu_count(logical=False)
        for c in range(C):
            for r in range(R):
                tag = f"c{c}_r{r}"
                image_path = os.path.join(input_dir, f"img_{tag}.tif")
                out_path = os.path.join(output_dir, f"dw_{tag}.tif")
                imsave(image_path, rescale_intensity(image_4d[r,c], out_range=np.float32))
                if tilesize == 0:
                    cmd = f"dw --iter {n_iter} --psigma {psigma} --overwrite --out {out_path} {image_path} {psfs[c]} --threads {threads_count}"
                else:
                    cmd = f"dw --iter {n_iter} --psigma {psigma} --overwrite --out {out_path} {image_path} {psfs[c]} --threads {threads_count} --tilesize {tilesize} --tilepad {tilepad}"
                if gpu and torch.cuda.is_available():
                    cmd += " --gpu"
                bash_lines.append(cmd)

        # Write and run the script
        script_path = os.path.join(temp_dir, "run_dw.sh")
        with open(script_path, "w") as f:
            f.write("\n".join(bash_lines))
        subprocess.run(f"bash {script_path}", shell=True)

        # Read output
        for c in range(C):
            for r in range(R):
                out_path = os.path.join(output_dir, f"dw_c{c}_r{r}.tif")
                if os.path.exists(out_path):
                    output[r,c] = imread(out_path)
                else:
                    print(f"WARNING: Missing output for channel {c}, round {r}")
    
    return output

def pass2Model(
    img: np.ndarray, 
    model: torch.nn.Module, 
    device: torch.device, 
    mask: bool = True, 
    rescale_i: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """
    Pass an image through a deep learning model and return the output.

    Parameters:
    img (np.ndarray): Input image.
    model (torch.nn.Module): PyTorch model for processing the image.
    device (torch.device): Device to run the model on (CPU or GPU).
    mask (bool, optional): Whether to return a segmentation mask. Defaults to True.
    rescale_i (bool, optional): Whether to apply additional intensity rescaling. Defaults to True.

    Returns:
    tuple[np.ndarray, np.ndarray]: Processed image and mask (if applicable).
    """
    # Normalize input image and convert it to float32
    input_img = rescale_intensity(img, out_range=(0, 1)).astype(np.float32)
    
    # Convert image to PyTorch tensor and move to the specified device
    input_torch = torch.from_numpy(input_img[np.newaxis, np.newaxis, ...]).to(device)
    try:
        model.to(device)
        model.module.to(device)
    except:
        model.to(device)
    # Forward pass through the model
    target_out_torch = model(input_torch)

    # Extract output and process accordingly
    if mask:
        # Get the main output and mask from model output
        target_out = target_out_torch[0].detach().cpu().numpy()[0, 0]
        mask_out = target_out_torch[1].detach().cpu().numpy()
    else:
        # Get the main output without mask
        target_out = target_out_torch[0].detach().cpu().numpy()[0, 0]
        mask_out = np.array([[0, 0], [0, 0]])  # Placeholder mask

    # Normalize the output
    target_out = target_out

    if rescale_i:
        # Adjust histogram to remove unwanted intensity peaks
        # target_out = match_histograms(target_out, img)
        target_out = target_out - np.median(target_out)
        target_out[target_out<0] = 0
        target_out = rescale_intensity(target_out,in_range=(target_out.min(),target_out.max()),out_range=(0,1))


    return target_out, mask_out[0][0]

#######################################
# Open Images
#######################################
section_name = snakemake.params.section
# Start logger
smk_logger = get_logger(section_name, filehandler = snakemake.log[0])
smk_logger.info("Input File:"+snakemake.input[0])
hs_image = HiSeqImage(image_path = snakemake.input[0], logger = smk_logger)
zarr_data = np.array(hs_image.im)
smk_logger.info(f"Type of zarr_data {type(zarr_data)}")
if len(zarr_data.shape)==5:
    zarr_data = np.max(zarr_data,axis=2)
smk_logger.info(f"Image Shape {zarr_data.shape}")
#######################################
# Get Parameters
#######################################
method = snakemake.config.get('deconvolution').get('method',"None")
num_iter = snakemake.config.get('deconvolution').get('num_iter',100)
psigma = snakemake.config.get('deconvolution').get('psigma',0)
model_path = snakemake.config.get('deconvolution').get('model_path',None)
machine_name = open(snakemake.config.get("image_path")+"/machine_name.txt").read()
psf_path = snakemake.config.get('deconvolution').get('psf_path',None)
tilepad = snakemake.config.get('deconvolution').get('tilepad',0)
tilesize = snakemake.config.get('deconvolution').get('tilesize',0)
module_location =  snakemake.config.get('deconvolution').get('module_location',"/gpfs/commons/groups/innovation/ecordina/deconvolution/google/Code/")
#######################################
# Pre Processing
#######################################
# Rounds Alignment
# Channel Should Already be aligned in fix_lighting
zarr_data = np.swapaxes(zarr_data,0,1)
out_sum = zarr_data.copy()
out_sum = np.max(out_sum, axis=1)
translations = [find_translation(out_sum[0], out_sum[j]) for j in range(1, zarr_data.shape[0])]
# Apply translations
out_final = np.zeros_like(zarr_data,dtype=np.float32)
out_final[0] = zarr_data[0]
smk_logger.info("Aligning Cycle and Channels")
for cycle, translation in enumerate(translations):
    out_final[cycle + 1] = shift(zarr_data[cycle + 1], [0.0, *translation]) 
# Clip negative values
#out_final = np.maximum(out_final, 0)
# Compute the cropping region
crop_y, crop_x = compute_crop_region(translations, out_final.shape)

# Apply cropping
# out_final = out_final[:,:, crop_x, crop_y]
zarr_data = out_final.copy()
zarr_data = np.swapaxes(zarr_data,0,1)

#######################################
# No Deconvolution
#######################################

if method.lower()=='none':
    hs_image.im.data = zarr_data
    hs_image.im.to_zarr(snakemake.output[0])
    smk_logger.info('None Method specified, Not running deconvolution')

#######################################
# Deconwolf Deconvolution
#######################################

elif method.lower()=='deconwolf' or method.lower()=='dw':
    if num_iter!=0:
        smk_logger.info('Running Deconwolf with '+ str(num_iter)+' iterations')
        gpu_avail = torch.cuda.is_available() 
        smk_logger.info('Deconvoling with GPU' if gpu_avail else "No GPU found, using CPU")
        hs_image.im.data = batch_run_deconwolf(zarr_data, psf_path, n_iter = num_iter, psigma = psigma, gpu = gpu_avail)
        hs_image.im.to_zarr(snakemake.output[0])
        smk_logger.info('Done with Deconvolution')
        #delayed_store = hs_image.write_ome_zarr(snakemake.output[0])\
    else:
        hs_image.im.to_zarr(snakemake.output[0])
        smk_logger.info('Deconwolf specified but no iterations given, Not running deconvolution')

#######################################
# CARE Deconvolution
#######################################

elif method.lower()=='care':

    from csbdeep.models import UpsamplingCARE

    CARE_model = UpsamplingCARE(config=None, name=os.path.basename(model_path), basedir=os.path.dirname(model_path))
    smk_logger.info('Deconvoling with GPU' if torch.cuda.is_available() else "No GPU found, using CPU")
    out_array=np.zeros_like(zarr_data)
    for round_i in range(zarr_data.shape[0]):
        for channel_i in range(zarr_data.shape[1]):
            care_image = CARE_model.predict(zarr_data[round_i,channel_i],"YX", 1,n_tiles=(2,2))
            out_array[round_i,channel_i]=rescale_intensity(care_image[round_i,channel_i],out_range=(0,1)).astype(np.float32)
    smk_logger.info(f"Done with CARE Deconvolution")
    hs_image.im.data=out_array
    hs_image.im.to_zarr(snakemake.output[0])

#######################################
# URNet Deconvolution
#######################################
elif method.lower()=='urnet':

    # Load URNet Model
    device = "cuda" if torch.cuda.is_available() else 'cpu'
    model = torch.load(model_path,map_location=torch.device(device),weights_only=False)
    out_array=np.zeros_like(zarr_data)
    for round_i in range(zarr_data.shape[0]):
        for channel_i in range(zarr_data.shape[1]):
            input_image = np.array(zarr_data[round_i,channel_i])
            deconvolved_image=pass2Model(input_image, model, device, rescale_i=False)[0]
            out_array[round_i,channel_i]=deconvolved_image
    smk_logger.info(f"Done with URNet Deconvolution")
    hs_image.im.data=out_array
    hs_image.im.to_zarr(snakemake.output[0])

else:
    hs_image.im.to_zarr(snakemake.output[0])
    smk_logger.info(f'Unknown Method specified, Got {method} specificied\n Not running Deconvolution')