"""
Quantify cell type subclass distribution by CCF region.
Used to generate Fig 5f"""


import math, os, re
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import SimpleITK as sitk
import requests
import pathlib
import seaborn as sns
from anytree import Node

def ccf_color_generator(order_list=None):
    with open('../data/metadata/687527670_mouse_brain_major_division_set.json', 'r') as file:
        ccf_12levels = json.load(file)
        file.close()

    ccf_color = []
    ccf_name = []
    for item in ccf_12levels:
            ccf_color.append('#'+item.get('color_hex_triplet'))
            ccf_name.append(item.get('acronym'))
    
    if order_list is not None:
        ccf_name = np.array(ccf_name)
        ccf_color_ordered = []
        for item in order_list:
            ccf_color_ordered.append(ccf_color[np.where(ccf_name==item)[0][0]])
        
        ccf_color = ccf_color_ordered
        
    return ccf_color

def dev_color_id(dev_level,dev_metadata):
    
    dev_color = []
    for item in dev_level:
        col = dev_metadata.loc[dev_metadata['Acronym']==item,'Hexidecimal'].values[0]

        dev_color.append('#'+col)
        
    return dev_color

def add_dev_labels(cell_joined):
    
    annotation = sitk.ReadImage('../data/KimLabDevCCFv001.nrrd')
    x = cell_joined['x_ccf'].to_numpy()
    y = cell_joined['y_ccf'].to_numpy()
    z = cell_joined['z_ccf'].to_numpy()
    
    dev_ccf_labels = []
    for idx in range(len(x)):
        dev_ccf_labels.append(annotation[int(100*x[idx]),int(100*y[idx]),int(100*z[idx])])
        
    cell_joined['dev_ccf_id'] = dev_ccf_labels
    
    return cell_joined

def fix_acronym_formatting(acronym):
    
    fixed_acronym = np.copy(acronym)
    for idx,item in enumerate(acronym):
        if (item[0]=="'"):
            fixed_acronym[idx] = item[1:-1]
        elif (item[-1]=="'"):
            fixed_acronym[idx] = item[:-1]
        else:
            fixed_acronym[idx] = item
            
    return fixed_acronym

def csv_to_tree(metadata_filename):
    
    metadata = pd.read_csv(metadata_filename)
    region_id = metadata['ID16'].tolist()
    region_name = metadata['Name'].to_numpy(dtype=str)
    acronym = metadata['Acronym'].to_numpy()

    parent_id = metadata['Parent ID16'].to_numpy(dtype=int)
    
    regions_dict = id_to_acronym(acronym,region_id)
    
    root = initialize_root(acronym,region_id,region_name)
    nodes = build_tree(root,acronym,region_id,regions_dict,parent_id,region_name)
    
    return root,nodes

def id_to_acronym(acronym,region_id):
    regions_dict = {}
    for idx,j in enumerate(region_id):
        regions_dict[j] = acronym[idx]
    return regions_dict

def initialize_root(acronym,region_id,region_name):
    
    return Node(acronym[0],parent=None,region_id=region_id[0],full_name=region_name[0])

def build_tree(root,acronym,region_id,regions_dict,parent_id,region_name):
    
    nodes={}
    nodes[root.name]=root
    for count,name in enumerate(acronym[1:]):
        nodes[name]=Node(
            name,
            parent=nodes[regions_dict[float(parent_id[count+1])]],
            region_id=region_id[count+1],
            full_name=region_name[count+1]
            )
    return nodes

def get_parent_id(cell_joined,dev_nodes,dev_level):
    
    dev_id = cell_joined['dev_ccf_id'].to_numpy()
    parent_acronym = np.empty_like(dev_id,dtype='U8')
    for node_name in dev_level:
        
        descendant_nodes = dev_nodes[node_name].descendants
        descendant_list = [x.region_id for x in descendant_nodes]
        
        matching_descendants = np.isin(dev_id,descendant_list)
        
        parent_acronym[matching_descendants] = dev_nodes[node_name].name
    
    return parent_acronym
        

def celltype_heatmap(heatmap_table,x,y,col,xticks,order_list,save_filename):
    
    plt.style.use('dark_background')
    ax = sns.heatmap(heatmap_table, annot=False,vmin=0,vmax=1, cmap='plasma',rasterized=True)

    plt.yticks(range(len(order_list)),order_list,fontsize=24)
    plt.xticks(xticks,[],fontsize=24)
    plt.ylabel('CCF region',fontsize=34)
    plt.xlabel('Subclass',fontsize=34)
    
    ax = plt.bar(x,y,color=col)

    plt.gcf().set_size_inches(45,15)
    plt.savefig(save_filename, dpi=100)
    plt.clf()
    
def remove_groups(cell_count,column_name,groups_list):
    
    for group in groups_list:
        unassigned_idx = cell_count.loc[cell_count[column_name]==group].index
        cell_count.drop(list(unassigned_idx),axis=0,inplace=True)
    
    return cell_count

def class_column_order(cell_joined):
    
    subclass = np.unique(cell_joined['subclass'])
    class_id = np.empty_like(subclass)
    for idx,item in enumerate(subclass):
        class_level = cell_joined.loc[cell_joined['subclass']==item,'class'].values[0]
        class_id[idx] = class_level

    class_idx = np.argsort(class_id,kind='stable')
    class_order = subclass[class_idx]

    return class_order,class_id
        
    
def grouped_xticks(df,low_level,high_level):
    
    ndivs = np.unique(df[low_level])
    class_ticks = []
    for div in ndivs:
        n = df.loc[df[low_level]==div,high_level].to_list()[0]
        class_ticks.append(n)
        
    tick_label,xticks = np.unique(class_ticks, return_counts=True)
    tick = np.cumsum(xticks)
    return tick,tick_label

def bar_data(df,bar_list,level,color_level):
    
    color = []
    class_name = []
    for bar in bar_list:
        n = df.loc[df[level]==bar,color_level].to_list()[0]
        color.append(n)
        m = df.loc[df[level]==bar,'class'].to_list()[0]
        class_name.append(m)
    x_val = list(range(len(bar_list)))
    y_val = np.ones(len(x_val))
    return x_val,y_val,color,class_name
    
def reorder_dataframe(df,ccf_level,order_list):
    
    df[ccf_level] = pd.Categorical(df[ccf_level], categories=order_list, ordered=True)
    
    return df.reset_index(drop=True)

def generate_order_list(region_list):
    
    ont = pd.read_csv(path_dir+'MouseAtlas_ontologies_notree.csv')
    order = np.zeros(len(region_list),dtype=np.int32)
    for count,j in enumerate(region_list):
        
        order[count] = ont.loc[(ont['Acronym']==j),'InDel'].values[0]

    return region_list[order.argsort()]

def plot_dev(cell_joined):
    
    dev_metadata_filename = os.path.join(path_dir,'metadata/AllenDevMouseOntologyStructure_DevCCFv004.csv')
    dev_tree,dev_nodes = csv_to_tree(dev_metadata_filename)
    dev_level = ['NeoCx','MesoCx','AlloCx','SPall','THy','PHy','p3','p2','p1','m1','m2','PPH','PH','PMH','MH']
    dev_column = 'dev_parent_name'

    parent_acronym = get_parent_id(cell_joined,dev_nodes,dev_level)
    cell_joined[dev_column] = parent_acronym
    
    dev_metadata = pd.read_csv(dev_metadata_filename)
    dev_color = dev_color_id(dev_level,dev_metadata)
    plot_legend(dev_color,dev_level,'../figures/dev_colors.svg',annotate=False)

    dev_cells = remove_groups(cell_joined,dev_column,[''])
    dev_cells.sort_values(by=['class','subclass'])

    cells_by_subclass = dev_cells.groupby('subclass').mean()
    cells_by_subclass = cells_by_subclass.iloc[:,0].to_list()
    cell_count = dev_cells.groupby([dev_column,'subclass']).count()
    cell_count.reset_index(inplace=True)

    region_list= np.unique(cell_count[dev_column])

    dev_order = np.array(dev_level)

    cell_county = reorder_dataframe(cell_count,dev_column,dev_order)
    cell_county.sort_values(by='dev_parent_name',inplace=True)

    dff = cell_county.pivot(index=dev_column,columns='subclass',values='cluster')
    total_by_subclass = dff.sum(axis=0).to_numpy()
    dff_norm = dff.div(total_by_subclass)
    class_order,class_id = class_column_order(cell_joined)
    dff_norm=dff_norm.reindex(class_order,axis=1)
    xticks,xtick_labels = grouped_xticks(dev_cells,'subclass','class')
    x,y,col,class_name = bar_data(dev_cells,list(dff_norm.columns),'subclass','class_color')

    celltype_heatmap(dff_norm,x,y,col,xticks,dev_level,'../figures/devccf_quantification.svg')
    
    dff_norm.to_csv('../data/Fig5f_DevCCF_subclass_normalized.csv')

def plot_ccf(cell_joined):
    
    ccf_level = 'parcellation_division'

    remove_group_list = ['AQ','V3','V4','VL','brain-unassigned','cbf','cm','eps','fiber tracts-unassigned','lfbs','mfbs','scwm','unassigned']
    cell_joined = remove_groups(cell_joined,ccf_level,remove_group_list)

    cell_count = cell_joined.groupby([ccf_level,'subclass']).count()
    cell_count.reset_index(inplace=True)

    region_list= np.unique(cell_count[ccf_level])
    # order_list = generate_order_list(region_list)
    order_list = ['Isocortex','CTXsp','HPF','OLF','STR','PAL','HY','TH','MB','CB','P','MY']
    
    ccf_color = ccf_color_generator(order_list)
    plot_legend(ccf_color,order_list,'../figures/ccf_colors.svg',annotate=False)

    cell_county = reorder_dataframe(cell_count,ccf_level,order_list)

    dff = cell_county.pivot(index=ccf_level,columns='subclass',values='cluster')
    total_by_subclass = dff.sum(axis=0).to_numpy()
    dff_norm = dff.div(total_by_subclass)
    class_order,class_id = class_column_order(cell_joined)
    dff_norm=dff_norm.reindex(class_order,axis=1)

    xticks,xtick_labels = grouped_xticks(cell_joined,'subclass','class')
    x,y,col,class_name = bar_data(cell_joined,list(dff_norm.columns),'subclass','class_color')
    
    names,col_idx = np.unique(class_name,return_index=True)
    colors = np.array(col)[col_idx]
    
    plot_legend(colors,names,'../figures/class_legend.svg')

    celltype_heatmap(dff_norm,x,y,col,xticks,order_list,'../figures/ccf_quantification.svg')
    dff_norm.to_csv('../data/Fig5f_CCF_subclass_normalized.csv')
    
def plot_legend(colors, names, save_filename, ncols=1,annotate=True):

    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 3

    n = len(names)
    nrows = math.ceil(n / ncols)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 100

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7
        
        if annotate==True:
            ax.text(text_pos_x, y, name, fontsize=14,
                    horizontalalignment='left',
                    verticalalignment='center')

        ax.add_patch(
            mpatch.Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                      height=18, facecolor=colors[i], edgecolor='0.7')
        )
        
        plt.savefig(save_filename)

    return fig

path_dir = '../data/'
url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/releases/20230630/manifest.json'
manifest = json.loads(requests.get(url).text)

view_directory = os.path.join( path_dir, 
                                manifest['directory_listing']['MERFISH-C57BL6J-638850-CCF']['directories']['metadata']['relative_path'], 
                              'views')

view_directory = pathlib.Path( view_directory )

cell_joined = pd.read_csv(os.path.join( view_directory, 'cell_metadata_with_dev_CCF_parcellation_annotations.csv'),dtype={'cell_label':str})

cells_df = add_dev_labels(cell_joined)

dev_ccf_data = plot_dev(cells_df)
ccf_data = plot_ccf(cells_df)

