#!/usr/bin/env python2.7

import sys

sys.path.append('/opt/hydra/')

import argparse
import bnpy
import logging
import numpy as np
import os
import pandas as pd
import shutil
import tempfile
import uuid
import subprocess

import library.analysis as hydra
from library.utils import mkdir_p
from library.fit import get_assignments

ens_to_hugo = {}
with open('/opt/typhon/data/EnsGeneID_Hugo_Observed_Conversions.txt') as f:
    header = next(f)
    for line in f:
        hugo, ens = line.strip().split('\t')
        ens_to_hugo[ens] = hugo


def fit_enrichment_models(data, diagnosis, output_dir):
    logger = logging.getLogger('root')
    models_pth = os.path.join('/opt/typhon/models/multivariate', diagnosis)

    # Load Enrichment Analysis
    for model in os.listdir(models_pth):
        logger.info("Applying %s model" % model)
        model_pth = os.path.join(models_pth, model, model)
        hmodel = bnpy.ioutil.ModelReader.load_model_at_prefix(model_pth,
                                                              prefix=model)
        # Load original training data
        train_pth = os.path.join(models_pth, model, model, 'training-data.tsv')
        train_data = pd.read_csv(train_pth, sep='\t', index_col=0)

        fit = hydra.PreFitMultivariateModel(hmodel, train_data)
        assignment, subgsea = fit.sub_cluster_gsea(data['TPM'])
        probs = fit.get_probability(data['TPM'])
        for c in range(probs.shape[1]):
            logger.info("Cluster %d Probability: %.3f" % (c, probs[0, c]))

        # Place model in cluster
        logger.info("Place in cluster %d" % assignment)
        if pd.isnull(assignment):
            logger.info("WARNING: Could not classify sample!")
            continue

        output_dir = os.path.join(output_dir, model)
        mkdir_p(output_dir)
        feature_src = os.path.join(models_pth, model, 'features', str(assignment), "cluster-%d-GSEA.tsv" % assignment)
        feature_dest = os.path.join(output_dir, "cluster-%d-GSEA.tsv" % assignment)
        shutil.copyfile(feature_src, feature_dest)

        sub_dest = os.path.join(output_dir, 'subcluster-%d-GSEA.tsv' % assignment)
        subgsea.sort_values("NES", ascending=False).to_csv(sub_dest, sep='\t')


def fit_druggable_genes(data, diagnosis, output):
    logger = logging.getLogger('root')
    models_pth = os.path.join('/opt/typhon/models/univariate/', diagnosis)

    df = pd.DataFrame(columns=["gene", "component", "fraction"])

    drug_genes = pd.read_csv(os.path.join(models_pth, "druggable-genes", "druggable-genes.tsv"),
                             sep='\t')

    drug_genes = drug_genes["search_term"].unique()

    train_pth = os.path.join(models_pth, 'druggable-genes', 'training-data.tsv')
    train_data = pd.read_csv(train_pth, sep='\t', index_col=0)
    for gene in drug_genes:
        model_pth = os.path.join(models_pth, "MultiModalGenes", gene)

        try:
            hmodel = bnpy.ioutil.ModelReader.load_model_at_prefix(model_pth,
                                                                  prefix=gene)

        except IOError:
            continue

        gene_train = train_data.reindex([gene])
        fit = hydra.PreFitMultivariateModel(hmodel, gene_train)
        a = fit.get_assignments(data["TPM"])

        means = []
        for comp in range(len(hmodel.allocModel.get_active_comp_probs())):
            means.append(hmodel.obsModel.get_mean_for_comp(comp)[0])

        means = np.argsort(means)
        if a == means[0]:
            df.loc[len(df), :] = [gene, "LOW", "%.2f" % hmodel.allocModel.get_active_comp_probs()[means[0]]]

        elif a == means[-1]:
            df.loc[len(df), :] = [gene, "HIGH", "%.2f" % hmodel.allocModel.get_active_comp_probs()[means[0]]]

        else:
            df.loc[len(df), :] = [gene, "MIDDLE", "%.2f" % hmodel.allocModel.get_active_comp_probs()[means[0]]]

    output_path = os.path.join(output, "druggable-gene-classification.tsv")
    df.to_csv(output_path, sep='\t', index=False)


def zscore(data, diagnosis, output, constant=0.05):
    logger = logging.getLogger('root')

    rnk_temp = os.path.join(tempfile.gettempdir(), 'RNK' + str(uuid.uuid4()))
    fgsea_temp = os.path.join(tempfile.gettempdir(), 'FGSEA' + str(uuid.uuid4()))

    models_pth = os.path.join('/opt/typhon/genes/', diagnosis)
    train_pth = os.path.join(models_pth, 'druggable-genes', 'training-data.tsv')
    train_data = pd.read_csv(train_pth, sep='\t', index_col=0)

    intersection = data.index.intersection(train_data.index)

    train_data = train_data.reindex(intersection)
    data = data.reindex(intersection)

    center = data.sub( train_data.mean(axis=1), axis=0 ) 
    zscore = center.divide( train_data.std(axis=1) + constant, axis=0 )
    zscore = zscore.sort_values('TPM', ascending=False)
    print(zscore.head())
    zscore.to_csv(rnk_temp,
                  header=False,
                  sep='\t')

    cmd = ['Rscript',
           '/opt/typhon/bin/fgsea.R',
           '/opt/typhon/data/Human_GOBP_AllPathways_no_GO_iea_December_01_2018_symbol.gmt',
           rnk_temp,
           fgsea_temp]

    subprocess.check_call(cmd)

    res = pd.read_csv(fgsea_temp, index_col=0)
    res = res.sort_values("NES", ascending=False)
    pth = os.path.join(output, "cohort-level-GSEA.tsv")
    res.to_csv(pth, sep='\t')


def main():
    """
    Typhon Pipeline
    """
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--diagnosis',
                        help='Patient diagnosis',
                        required=True)

    parser.add_argument('--RSEM',
                        help='Path to N-of-1 RSEM file.',
                        required=True)

    parser.add_argument('--TME',
                        help='Runs xCell tumor microenvironment profiling anlaysis.',
                        action="store_true",
                        required=False)

    parser.add_argument('-o', '--output-dir',
                        help='Output directory',
                        default='typhon-output')

    parser.add_argument('--debug',
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    available_models = ['MYCN-NA-Neuroblastoma', 'Osteosarcoma']

    if args.diagnosis not in available_models:
        msg = "Please select one of the following diagnoses:\n%s" % '\n'.join(available_models)
        raise ValueError(msg)

    # Set up logger
    level = logging.INFO

    # Make the output directory if it doesn't already exist
    mkdir_p(args.output_dir)

    logging.basicConfig(filename=os.path.join(args.output_dir, 'typhon.log'),
                        level=level)
    logging.getLogger().addHandler(logging.StreamHandler())
    logger = logging.getLogger('root')

    # Read in RSEM file
    data = pd.read_csv(args.RSEM, sep='\t')

    # Convert to hugo ids
    data['hugo'] = data['gene_id'].map(ens_to_hugo)
    tpm = data.reindex(['hugo', 'TPM'], axis=1).groupby('hugo').sum()
    exp = np.log2(tpm + 1)

    logger.info("Starting run...")
    fit_enrichment_models(exp, args.diagnosis, args.output_dir)

    fit_druggable_genes(exp, args.diagnosis, args.output_dir)

    zscore(exp, args.diagnosis, args.output_dir)


if __name__ == '__main__':
    main()
