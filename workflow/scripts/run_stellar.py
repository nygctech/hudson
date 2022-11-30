# from pre.utils import get_config
# from utils import get_cluster, get_logger
# logger = get_logger(logname = section_name, filehandler = snakemake.log[0])
print(snakemake.input[0])
print(list(snakemake.params.ref)[0])
print(list(snakemake.params.typecol)[0])
import numpy as np
import torch
import pandas as pd
import anndata
import scanpy as sc
import pickle
from collections import defaultdict
import argparse
from utils import prepare_save_dir
from datasets import GraphDataset, load_data
from STELLAR import STELLAR
print("settings")
# settings
parser = argparse.ArgumentParser(description='STELLAR')
parser.add_argument('--dataset', default='Hubmap', help='dataset setting')
parser.add_argument('--seed', type=int, default=1, metavar='S', help='random seed (default: 1)')
parser.add_argument('--name', type=str, default='STELLAR')
parser.add_argument('--epochs', type=int, default=20)
parser.add_argument('--lr', type=float, default=1e-3)
parser.add_argument('--wd', type=float, default=5e-2)
parser.add_argument('--input-dim', type=int, default=48)
parser.add_argument('--num-heads', type=int, default=22)
parser.add_argument('--num-seed-class', type=int, default=0)
parser.add_argument('--sample-rate', type=float, default=0.05)
parser.add_argument('-b', '--batch-size', default=1, type=int,
                    metavar='N',
                                    help='mini-batch size')
parser.add_argument('--distance_thres', default=50, type=int)
parser.add_argument('--savedir', type=str, default='./')
args = parser.parse_args(args=[])
args.cuda = torch.cuda.is_available()
args.device = torch.device("cuda" if args.cuda else "cpu")

print("load")
# load
labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges, inverse_dict = load_data(list(snakemake.params.ref)[0], snakemake.input[0], list(snakemake.params.typecol)[0], 50, 1)
dataset = GraphDataset(labeled_X, labeled_y, unlabeled_X, labeled_edges, unlabeled_edges)
# stellar
print("learn")
stellar = STELLAR(args, dataset)
stellar.train()
temp, results = stellar.pred()
# output
print("out")
inverse_dict2 = defaultdict(lambda: "unknown", inverse_dict)
results2 = [inverse_dict2[x] for x in results]
with open(snakemake.output[0], "w") as outfile:
      outfile.write("\n".join(results2))
