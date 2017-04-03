import collections

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import tables

np.random.seed(0)

GeneBCMatrix = collections.namedtuple('GeneBCMatrix', ['gene_ids', 'gene_names', 'barcodes', 'matrix'])


def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            matrix = sp_sparse.csc_matrix((dsets['data'], dsets['indices'], dsets['indptr']), shape=dsets['shape'])
            return GeneBCMatrix(dsets['genes'], dsets['gene_names'], dsets['barcodes'], matrix)
        except tables.NoSuchNodeError:
            raise Exception("Genome %s does not exist in this file." % genome)
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


def save_matrix_to_h5(gbm, filename, genome):
    flt = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=flt) as f:
        try:
            group = f.create_group(f.root, genome)
            f.create_carray(group, 'genes', obj=gbm.gene_ids)
            f.create_carray(group, 'gene_names', obj=gbm.gene_names)
            f.create_carray(group, 'barcodes', obj=gbm.barcodes)
            f.create_carray(group, 'data', obj=gbm.matrix.data)
            f.create_carray(group, 'indices', obj=gbm.matrix.indices)
            f.create_carray(group, 'indptr', obj=gbm.matrix.indptr)
            f.create_carray(group, 'shape', obj=gbm.matrix.shape)
        except:
            raise Exception("Failed to write H5 file.")


def subsample_matrix(gbm, barcode_indices):
    return GeneBCMatrix(gbm.gene_ids, gbm.gene_names, gbm.barcodes[barcode_indices], gbm.matrix[:, barcode_indices])


def get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()

################## MAIN CODE ###################################
filtered_matrix_h5 = "1M_neurons_filtered_gene_bc_matrices_h5.h5"
genome = "mm10"

tsne = pd.read_csv("analysis/tsne/2_components/projection.csv")
clusters = pd.read_csv("analysis/clustering/graphclust/clusters.csv")

gene_bc_matrix = get_matrix_from_h5(filtered_matrix_h5, genome)

# np.savetxt('gbcm.csv', delimiter=",", X=gene_bc_matrix.matrix.data)
# print ("hello world")
print(gene_bc_matrix.matrix)

# load TSNE and graph clustering
# tsne = pd.read_csv("analysis/tsne/2_components/projection.csv")
# clusters = pd.read_csv("analysis/clustering/graphclust/clusters.csv")
#
# subsample_bcs = 2000
# subset = np.sort(np.random.choice(gene_bc_matrix.barcodes.size, size=subsample_bcs, replace=False))
# subsampled_matrix = subsample_matrix(gene_bc_matrix, subset)
# subsampled_tsne = tsne.loc[subset, :]
# subsampled_clusters = clusters.loc[subset, :]
#
# umis_per_cell = np.asarray(subsampled_matrix.matrix.sum(axis=0)).squeeze()
# genes_per_cell = np.asarray((subsampled_matrix.matrix > 0).sum(axis=0)).squeeze()
#
# print (gene_bc_matrix[:1])