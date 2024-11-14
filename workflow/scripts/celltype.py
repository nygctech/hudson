import muon as mu
from scipy.stats import mode
import numpy as np
import yaml
import anndata as ad 
import pandas as pd
import sys
import pickle
from utils import get_cluster, get_logger
from astir.astir import Astir


## Define Normalization Options 

def arcsinh(X, cofactor_dict, var_names):

    
    nrows, ncols = X.shape
    stack = []
    for c in range(ncols):
        cofactor_ = cofactor_dict[var_names[c]]
        stack.append(np.arcsinh(X[:,c]/(cofactor_)))
    return np.stack(stack, axis=1)


def minmax(X, minper = 5, maxper = 99.99):
    
    nrows, ncols = X.shape
    stack = []
    for c in range(ncols):
        minexp, maxexp = np.percentile(X[:,c], [minper, maxper])
        stack.append(((X[:,c]-minexp)/(maxexp-minexp)).clip(min=0))
        
    return np.stack(stack, axis=1)


section_name = snakemake.params.section

# Start logger
smk_logger = get_logger(section_name, filehandler = snakemake.log[0])


## Read in Data
data_path = snakemake.input[0]
smk_logger.info(f'Reading in {data_path}')
mdata = mu.read(data_path)

marker_path = snakemake.config.get('celltype').get('markers', None)
smk_logger.info(f'Reading in {marker_path}')
try: 
    with open(marker_path, 'r') as f: 
        marker_dict = yaml.safe_load(f)
except:
    print("Error: Unable to read in a marker file")
    sys.exit(1)

## Select Protein Data
adata = mdata['protein']

norm = snakemake.config.get('celltype').get('normalization_method', None)

## Perform Normalization 
if norm is None: 
    smk_logger.info('Performing arcsinh nomralization')
    cofactors_dict = {}
    for i,v in enumerate(adata.var.index):
        cofactors_dict[v] =  mode(adata.X[:,i])[0]

    X = arcsinh(adata.X, cofactors_dict, adata.var.index)
else:
    smk_logger.info('Performing minmax normalization')
    minper = snakemake.config.get('celltype').get('minper', 5)
    maxper = snakemake.config.get('celltype').get('maxper', 99.99)
    X = minmax(adata.X, minper, maxper)

## Create Astir Object
smk_logger.info('Initalizing Astir Model')
df_gex = pd.DataFrame(X)
df_gex.columns = adata.var_names
df_gex.index = adata.obs_names
ast = Astir(df_gex, marker_dict)

## Get Training Parameters
N = adata.shape[0]
prop = snakemake.config.get('celltype').get('batch_size_proportion', 100)
batch_size = int(N/prop)
max_epochs = snakemake.config.get('celltype').get('epochs', 1000)
learning_rate = snakemake.config.get('celltype').get('learning_rate', 2e-3)
initial_epochs = snakemake.config.get('celltype').get('initial_epochs', 3)


## Training Model 
smk_logger.info('Training Model')
ast.fit_type(max_epochs = max_epochs,
             batch_size = batch_size,
             learning_rate = learning_rate,
             n_init_epochs = initial_epochs)

cells = ast.get_celltypes()

## Checking to see if output 'Other' Category should be removed
remove = snakemake.config.get('celltype').get('remove_other', False)
if remove:
    probs = ast.get_celltype_probabilities()
    for i, cell in enumerate(cells['cell_type']):
        if cell == 'Other':
            cells['cell_type'][i] = probs.iloc[i].sort_values(ascending=False).index[1]


## Save output 
smk_logger.info('Saving Output')
ast_adata = ad.AnnData(X, obs = cells, var= adata.var)
adata.uns['model'] = ast._type_run_info
mdata.mod['celltypes'] = ast_adata
mdata.write(snakemake.output[0])

smk_logger.info('Writing metrics to Summary')
with open(snakemake.input[1], 'r') as file:
    data = yaml.safe_load(file)

counts = ast.get_celltypes().value_counts()
data['celltyping'] = {}
data['celltyping']['counts'] = counts.to_dict()
data['celltyping']['proportion'] = (counts/N).to_dict()

with open(snakemake.input[1], 'w') as file:
    yaml.dump(data, file)

