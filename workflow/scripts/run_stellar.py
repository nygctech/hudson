from pre.utils import get_config
from utils import get_cluster, get_logger

experiment_config = get_config(snakemake.input[0])
exp_dir = snakemake.config['experiment_directory']
section_name = snakemake.params.section
logger = get_logger(logname = section_name, filehandler = snakemake.log[0])

import argparse
import numpy as np
import torch
import pandas as pd
import anndata
import scanpy as sc
import pickle

from utils import prepare_save_dir
from datasets import GraphDataset, load_data
from STELLAR import STELLAR

labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges, inverse_dict = load_data('/gpfs/commons/groups/nygcfaculty/PySeq/spatial_analysis/tonsil_codex/stellar/data/B004_reg3.csv', '/gpfs/commons/groups/nygcfaculty/PySeq/spatial_analysis/tonsil_codex/stellar/data/B004_reg4.csv', 50, 1)
dataset = GraphDataset(labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges)

stellar = STELLAR(args, dataset)
stellar.train()
temp, results = stellar.pred()

from collections import defaultdict
inverse_dict2 = defaultdict(lambda: "unknown", inverse_dict)
results2 = [inverse_dict2[x] for x in results]

