from zarr.errors import GroupNotFoundError
from utils import get_logger,HiSeqImage
from scipy.ndimage import shift
from skimage.registration import phase_cross_correlation
from skimage import morphology
from skimage.io import imread,imsave
from skimage.exposure import rescale_intensity
import numpy as np
import os.path
import zarr
from os import makedirs,system
import subprocess
import pandas as pd
import glob
import tempfile
from subprocess import run, PIPE
from typing import List,Tuple
#######################################
# Utility Functions
#######################################



def box(point,max_x,max_y,size=8):
    x,y=point
    x2=min(round(x+size),max_x-1)
    y2=min(round(y+size),max_y-1)
    if x2>1 and y2>1:
        x1=min(max(round(x-size),0),x2-1)
        y1=min(max(round(y-size),0),y2-1)
        return x1,x2,y1,y2
    return None,None,None,None

def run_rs_fish(image: np.ndarray,threshold=1e-3,sigma=1e-1,ransac: int = 1,tempdir: str = None,save_output: bool = False,save_name: str =  'output_spots.csv') -> np.ndarray:
    """
    Runs RS-FISH on the given image and returns detected spots.

    Parameters:
    - image (np.ndarray): Input image for spot detection.

    Returns:
    - np.ndarray: Array of detected spot coordinates.
    """

    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, "input_image.tif")
        output_path = os.path.join(temp_dir, "output_spots.csv")
        
        # Save the image
        imsave(input_path, rescale_intensity(image, out_range=np.uint16))
        
        # Run RS-FISH
        _ =run(f"bash -c 'source /etc/profile && module load RS-FISH && rs-fish --image={input_path} --output={output_path} --threshold={threshold} --sigma={sigma} --ransac={ransac}'", stdout=PIPE, stderr=PIPE,shell=True)
        
        # Load and return detected spots
        if os.path.exists(output_path):
            spots = pd.read_csv(output_path)
        else:
            spots = np.array([])
        if save_output:
            _ = run(f"cp {output_path} {tempdir}/{save_name}", stdout=PIPE, stderr=PIPE,shell=True)
    if len(spots)>0:
        return spots[['y','x']].values
    return np.array([])

def filter_close_points(points, min_distance):
    """
    Filters points that are too close to each other.
    
    Parameters:
    points (np.ndarray): An (n,2) array of points.
    min_distance (float): Minimum allowable distance between points.
    
    Returns:
    np.ndarray: Filtered points.
    """
    if len(points) == 0:
        return points
    
    filtered_points = [points[0]]
    
    for point in points[1:]:
        distances = np.linalg.norm(np.array(filtered_points) - point, axis=1)
        if np.all(distances >= min_distance):
            filtered_points.append(point)
    
    return np.array(filtered_points)

#######################################
# Open Images
#######################################
section_name = snakemake.params.section
# Start logger
smk_logger = get_logger(section_name, filehandler = snakemake.log[0])

try:
    # hs_image = HiSeqImage(image_path = snakemake.input[0], logger = smk_logger)
    # zarr_data = hs_image.im
    data = np.array(zarr.open(snakemake.input[0])[section_name])
except GroupNotFoundError:
    # hs_image = HiSeqImage(image_path = snakemake.input[1], logger = smk_logger)
    # zarr_data = hs_image.im
    smk_logger.info("Image not found!")
#######################################
# Get Parameters
#######################################
lowQualityChannels=snakemake.config.get('spotCalling').get('LowQualityChannels',list())
goodQualityChannels=[i for i in range(data.shape[0]) if i not in lowQualityChannels]
footprint=snakemake.config.get('spotCalling').get('footprint',5)
sigma=snakemake.config.get('spotCalling').get('sigma',1)
threshold=snakemake.config.get('spotCalling').get('threshold',1e-3)
ransac=snakemake.config.get('spotCalling').get('ransac',1)
smk_logger.info(f"data is of shape: {data.shape}")
#######################################
# Preprocessing
#######################################

# White Top Hat Filtering on all the individual Channels
footprint = morphology.disk(footprint)
filtered_image = data.copy()
for cycle in range(data.shape[1]):
    for channel in range(data.shape[0]):
        filtered_image[channel,cycle] = morphology.white_tophat(data[channel,cycle,:,:], footprint)


#######################################
# Running RS-FISH
#######################################

# Run RS-FISH on the Channel and Rounds Max projection of the image
max_image=np.max(np.max(filtered_image[goodQualityChannels],axis=0),axis=0)
smk_logger.info(f"OG Data is of shape: {filtered_image.shape}")
smk_logger.info(f"RSFISH Data is of shape: {max_image.shape}")

rs_points = run_rs_fish(max_image,sigma=sigma,threshold=threshold,ransac=ransac)

#######################################
# Spots Post-Processing
#######################################

# DataFrame to Store Results
df_spots = pd.DataFrame(columns=["ID",'x', 'y', 'Cycle', *[f'Channel_{i}' for i in range(filtered_image.shape[0])]])
# Filter spots that are less than 2 pixels apart
filtered=filter_close_points(rs_points, 2)
smk_logger.info(f"x_max,y_max is {np.max(filtered[:,0]),np.max(filtered[:,1])}")
box_size=5
# Get Pixel of Max Intensity in a small region centered around the detected spot in the max projection image
# Correct for any shift left in the image
smk_logger.info(f"Shape is {filtered_image.shape}")
for id,point in enumerate(filtered):
    x1,x2,y1,y2=box(point,filtered_image.shape[-2],filtered_image.shape[-1],box_size)
    if x1 is not None:
        for cycle in range(filtered_image.shape[1]):
            intensities=list()
            for channel in range(filtered_image.shape[0]):
                centered=filtered_image[channel,cycle,x1:x2,y1:y2]
                x_maxes,y_maxes=np.where(centered==centered.max())
                x_max, y_max=x_maxes[0],y_maxes[0]
                intensities.append(centered[x_max, y_max])
            df_spots.loc[len(df_spots)] = [id,x1+x_max, y1+y_max, cycle, *intensities]

smk_logger.info(f"{len(df_spots) / filtered_image.shape[0]} Beads detected")
df_spots = df_spots.sort_values(by=['ID', 'Cycle']).reset_index(drop=True)

#######################################
# Saving
#######################################
df_spots.to_csv(snakemake.output[0],index=False)
if not len(df_spots):
    smk_logger.info("No Spots found in any Channels")
    pd.DataFrame([]).to_csv(snakemake.output[0],index=False)

#######################################
# Function to run RS-FISH on every single channels and rounds
#######################################
# def run_rs_fish_on_folder(
#     data: np.ndarray,
#     threshold: float = 1e-3,
#     sigma: float = 1e-1,
#     ransac: int = 1,
#     logger = None,
# ) -> pd.DataFrame:
#     """
#     Write individual Channels and run RS-FISH on them all .tif images and returns combined results.

#     Parameters:
#     - folder_path (str): Path to folder containing .tif images.
#     - threshold, sigma, ransac: RS-FISH parameters.

#     Returns:
#     - pd.DataFrame: Combined DataFrame of all detected spots.
#     """

#     with tempfile.TemporaryDirectory() as temp_dir:
#         bash_script_path = os.path.join(temp_dir, "run_rsfish.sh")
#         channels_path = os.path.join(temp_dir, "Channels")
#         logger.info("Writing Images in temp dir")
#         # Write Images to Disk and get Threshold
#         for i in range(data.shape[0]):
#             # Rescale Images
#             channel_individual=rescale_intensity(np.max(data[i],axis=0),out_range=np.float32)
#             channel_individual=rescale_intensity(channel_individual,out_range=(0,1))
#             # Get Threshold if not given
#             if threshold is None or threshold == list():
#                 perc=list()
#                 for j in range(101):
#                     perc.append(np.percentile(channel_individual,j))
#                 perc=np.array(perc)
#                 threshold=np.percentile(perc[perc!=0],33)
#             # Write Images to Disk 
#             image_path= os.path.join(channels_path,f"channel_{i}.tiff")
#             imsave(image_path,channel_individual)
#         logger.info("Done Writing Images in temp dir")
#         # Create the bash script
#         with open(bash_script_path, "w") as f:

#             f.write("source /etc/profile")
#             f.write("module load RS-FISH")

#             for i, tif_path in enumerate(tqdm(channels_path, desc="Generating RS-FISH script")):
#                 #Get threshold
#                 if type(threshold)==list:
#                     threshold_i = threshold[i]
#                 else:
#                     threshold_i = threshold
#                 csv_path = os.path.join(channels_path, tif_path.replace(".tif",'.csv'))
#                 cmd = (
#                     f"rs-fish --image="{tif_path}" "
#                     f"--output="{csv_path}" "
#                     f"--threshold={threshold_i} "
#                     f"--sigma={sigma} "
#                     f"--ransac={ransac}n"
#                 )
#                 f.write(cmd)

#                 # Every 10 files, echo progress
#                 if (i + 1) % 10 == 0:
#                     f.write(f'echo "Processed {i + 1} images..."n')

#         # Make script executable and run it
#         run(f"chmod +x {bash_script_path}", shell=True)
#         logger.info("Running RS-FISH on all images...")
#         run(f"bash {bash_script_path}", stdout=PIPE, stderr=PIPE, shell=True)
#         logger.info("Done Running RS-FISH on all images...")
#         # Read all the output CSVs
#         all_spots = []
#         csv_files = sorted(glob.glob(os.path.join(temp_dir, "*_spots.csv")))
#         for csv_file in tqdm(csv_files, desc="Reading RS-FISH outputs"):
#             try:
#                 df = pd.read_csv(csv_file)
#                 df['Channels'] = os.path.splitext(os.path.basename(csv_file))[0].replace('_spots', '')
#                 all_spots.append(df)
#             except Exception as e:
#                 logger.info(f"Error reading {csv_file}: {e}")

#     if all_spots:
#         return pd.concat(all_spots, ignore_index=True)
#     return pd.DataFrame()