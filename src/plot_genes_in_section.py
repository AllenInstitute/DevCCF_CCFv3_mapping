# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 14:43:18 2023

@author: ashwin.bhandiwad
"""

import os, subprocess
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import requests
import SimpleITK as sitk
import pathlib
import anndata

def plot_section( xx, yy, cc=None, val=None, fig_width = 6, fig_height = 6, cmap=None ) :
    
    plt.style.use('dark_background')
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)
    if cmap is not None :
        sc=plt.scatter(xx,yy,s=0.5,c=val,marker='.',cmap=cmap)
    elif cc is not None :
        sc=plt.scatter(xx,yy,s=0.5,color=cc,marker='.')
        
    cbar = plt.colorbar(sc)
    cbar.set_label('LogCPM')
    ax.set_ylim(11,0)
    ax.set_xlim(0,11)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    return fig, ax

def create_expression_dataframe( ad, gf, section ) :
    gdata = ad[:,gf.index].to_df()
    gdata.columns = gf.gene_symbol
    joined = section.join( gdata )
    return joined

def expression_volume(asubset, df, gene):
    
    gf = asubset.var[asubset.var.gene_symbol == gene]
    gdata = asubset[:,gf.index].to_df()
    gdata.columns = gf.gene_symbol
    joined = df.join(gdata)
    
    return joined.loc[joined[gene]>0]

def plot_gene( asubset, gene, section) :
    
    gf = asubset.var[asubset.var.gene_symbol == gene]
    expression = create_expression_dataframe( asubset, gf, section )
    fig,ax = plot_section(expression['x'], expression['y'],val=expression[gene],cmap='plasma')
    print(max(expression[gene]))
    plt.savefig(f'../figures/{gene}_expression.png')
    
def plot_gene_volume( asubset, cell_joined, gene) :
    
    vol = np.zeros([1140,800,1320],dtype=np.float32)
    joined_filtered = expression_volume(asubset, cell_joined, gene)
    for idx in range(len(joined_filtered)):
        
        x = int(100*(joined_filtered['x_ccf'][idx]))
        y = int(100*(joined_filtered['y_ccf'][idx]))
        z = int(100*(joined_filtered['z_ccf'][idx]))
        
        vol[z,y,x] = joined_filtered[gene][idx]
    
    vol_img = sitk.GetImageFromArray(vol)
    # vol_img = sitk.PermuteAxes(vol_img,[2,1,0])
    sitk.WriteImage(vol_img,f'../figures/{gene}_volume.nrrd',True)
    
def plot_section_series( df, asubset, gene, blist, cmap=None, fig_width = 30, fig_height = 5) :

    plt.style.use('dark_background')
    fig, ax = plt.subplots(1,len(blist))
    fig.set_size_inches(fig_width, fig_height)
    
    for idx,bsl in enumerate(blist) :
        
        section = df[cell_filtered['brain_section_label'] == bsl]
        gf = asubset.var[asubset.var.gene_symbol == gene]
        expression = create_expression_dataframe( asubset, gf, section )
        xx = expression['x']
        yy = expression['y']
        vv = expression[gene]
        
        if cmap is not None :
            sc=ax[idx].scatter(xx,yy,s=1.0,c=vv,marker='.',cmap=cmap)
        else :
            sc=ax[idx].scatter(xx,yy,s=1.0,color=vv,marker=".")
        
        if idx==len(blist)-1:
            cbar = plt.colorbar(sc, ax=ax[len(blist)-1])
            cbar.set_label('LogCPM')
            
        ax[idx].axis('equal')
        ax[idx].set_xlim(0,11)
        ax[idx].set_ylim(11,0)
        ax[idx].set_xticks([])
        ax[idx].set_yticks([])
        
        ax[idx].set_title("%s" % (bsl) )

    plt.subplots_adjust(wspace=0.01, hspace=0.01)
    plt.savefig(f'../figures/{gene}_expression.png')
    return fig, ax
    

def aggregate_by_metadata( df, gnames, value, sort=False ) :
    grouped = df.groupby(value)[gnames].mean()
    if sort :
        grouped = grouped.sort_values(by=gnames[0],ascending=False)
    return grouped

url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/releases/20230630/manifest.json'
manifest = json.loads(requests.get(url).text)

path_dir = '../data/'


metadata = manifest['file_listing']['MERFISH-C57BL6J-638850']['metadata']
view_directory = os.path.join( path_dir, 
                               manifest['directory_listing']['MERFISH-C57BL6J-638850']['directories']['metadata']['relative_path'], 
                              'views')
cache_views = False
if cache_views :
    os.makedirs( view_directory, exist_ok=True )
    
expression_matrices = manifest['file_listing']['MERFISH-C57BL6J-638850']['expression_matrices']
rpath = expression_matrices['C57BL6J-638850']['log2']['files']['h5ad']['relative_path']
file = os.path.join( path_dir, rpath)

adata = anndata.read_h5ad(file,backed='r')
gene = adata.var

#gnames = ['Grik3','Penk','Gad2','Sox2']
# gnames = ['Nr4a3','Crym','Calb1','Gpr139','Otof','Gad2']
gnames = ['Calb1','Calb2','Dcn','Grp','Trhr','Dio3']
pred = [x in gnames for x in gene.gene_symbol]
gene_filtered = gene[pred]

cell = pd.read_csv(view_directory+'/cell_metadata_with_cluster_annotation.csv',dtype={"cell_label":str})

cell_filtered = cell.loc[cell['cell_label']==adata.obs.index]
pred = cell_filtered['low_quality_mapping'] == False
cell_filtered = cell_filtered[pred]
cell_filtered.set_index('cell_label',inplace=True)

rpath = metadata['example_genes_all_cells_expression']['files']['csv']['relative_path']
file = os.path.join( path_dir, rpath )
exp = pd.read_csv(file,dtype={"cell_label":str})
exp.set_index('cell_label',inplace=True)

pred = (cell_filtered['brain_section_label'] == 'C57BL6J-638850.40')
section = cell_filtered[pred]

asubset = adata[cell_filtered.index,gene_filtered.index].to_memory()

adata.file.close()
del adata, cell_filtered, gene_filtered, pred

# # blist = np.unique(cell_filtered['brain_section_label'])
# blist = ['C57BL6J-638850.50','C57BL6J-638850.38','C57BL6J-638850.28']
# fig, ax = plot_section_series( cell_filtered, asubset, 'Penk', blist, cmap='plasma')
# blist = ['C57BL6J-638850.46','C57BL6J-638850.40','C57BL6J-638850.38']
# fig, ax = plot_section_series( cell_filtered, asubset, 'Tac2', blist, cmap='plasma')

# view_directory = os.path.join( path_dir, 
#                                 manifest['directory_listing']['MERFISH-C57BL6J-638850-CCF']['directories']['metadata']['relative_path'], 
#                               'views')

# view_directory = pathlib.Path( view_directory )

# cell_joined = pd.read_csv(os.path.join( view_directory, 'cell_metadata_with_parcellation_annotations.csv'),dtype={'cell_label':str})
# cell_joined.set_index('cell_label',inplace=True)

# for gene_name in gnames:
#     plot_gene_volume(asubset,cell_joined,gene_name)

# plot_gene(asubset,'Penk',section)
# plot_gene(asubset,'Grik3',section)
# plot_gene(asubset,'Gad2',section)
