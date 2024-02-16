"""
Set of functions for generating DevCCF to CCFv3 mappings as a Sankey flow diagram
Used to generate Fig 5c

"""

import nrrd, json
import SimpleITK as sitk
import numpy as np
import pandas as pd
from anytree import Node
import plotly.graph_objects as go


def read_metadata_excel(metadata_filename):
    
    dtype_types = {'ID16':np.int32,
                   'Acronym':str,
                   'r':np.int32,
                   'g':np.int32,
                   'b':np.int32,
                   'Graph Order':float,
                   'Parent ID16':np.int32
                   }
    
    return pd.read_excel(metadata_filename,dtype=dtype_types)
    
    
def filter_metadata(volume,metadata_filename,dev=False):
    
    metadata = read_metadata_excel(metadata_filename)
    
    volume_ids = np.unique(volume)
    if dev==True:
        filtered_metadata = metadata[metadata['ID16'].isin(volume_ids)]
    else:
        filtered_metadata = metadata[metadata['ID'].isin(volume_ids)]
    
    return filtered_metadata,metadata

def json_to_tree(json_filename):
    
    with open(json_filename, 'r') as file:
        dataset = json.load(file)
        file.close()
        
    # initialize root node
    root_data = dataset['msg'][0]
    root = Node(name='997', 
                full_name=root_data.get('name'),
                acronym=root_data.get('acronym'),
                parent=None, 
                color=root_data.get('color_hex_triplet'))
    nodes = {'997' : root}
    tree_builder(root_data.get('children'),nodes)
    
    return root,nodes

def tree_builder(d,nodes):

    if isinstance(d,list):
        for item in d:
            if isinstance(item, dict):
                node_uid = str(item.get('id'))
                node = nodes[node_uid] = Node(name=node_uid,
                                       full_name=item.get('name'),
                                       acronym=item.get('acronym'),
                                       parent=nodes[str(item.get('parent_structure_id'))], 
                                       color=item.get('color_hex_triplet'))
                tree_builder(item.get('children'),nodes)
            
    return nodes

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
    for count,region_acronym in enumerate(acronym[1::]):
        name = acronym[count+1]
        nodes[name]=Node(
            name,
            parent=nodes[regions_dict[float(parent_id[count+1])]],
            region_id=region_id[count+1],
            full_name=region_name[count+1]
            )
    return nodes

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

def csv_to_tree(csv_filename):
    
    dev_ccf_info = pd.read_csv(csv_filename,sep=',')
    region_id = dev_ccf_info['ID16'].tolist()
    region_name = dev_ccf_info['Name'].tolist()
    acronym = dev_ccf_info['Acronym'].tolist()
    acronym = fix_acronym_formatting(acronym)
    parent_id = dev_ccf_info['Parent ID16'].tolist()
    
    regions_dict = id_to_acronym(acronym,region_id)
    
    root = initialize_root(acronym,region_id,region_name)
    nodes = build_tree(root,acronym,region_id,regions_dict,parent_id,region_name)
    
    return root,nodes

def combine_by_level(ccf,nodes,level,dev=True):
    
    volume = np.zeros(np.shape(ccf),dtype=np.int16)
    for count,area in enumerate(level):
        all_descendants=nodes[area].descendants

        for des in all_descendants:
            if dev==False:
                volume[np.where(ccf==float(des.name))]=count+1
            else:
                volume[np.where(ccf==float(des.region_id))]=count+1
            
    return volume

def generate_mapping(ccf,dev_ccf,cutoff=250):
    
    brain_mask = np.nonzero(ccf)
    mapping = []

    for n in range(np.shape(brain_mask)[1]):
        point = tuple([brain_mask[0][n],brain_mask[1][n],brain_mask[2][n]])
        
        if ((ccf[point]>0) & (dev_ccf[point]>0)):
            mapping.append(ccf[point] + 1j*dev_ccf[point]) #Hack: combine ccf and devCCF together as a complex number to compute unique pairs

    unique_pairs, counts = np.unique(mapping, return_counts=True)
    
    mapping = np.array([unique_pairs.real,unique_pairs.imag,counts]).transpose()
    
    mapping = mapping[mapping[:,2]>cutoff] #Filter mappings with less than 250 voxels

    return mapping


def ccf_level_generator(ccf_json_filename):
    with open(ccf_json_filename, 'r') as file:
        ccf_12levels = json.load(file)
        file.close()
        
    ccf_level = []
    ccf_color = []
    for item in ccf_12levels:
            ccf_level.append(str(item.get('id')))
            ccf_color.append('#'+item.get('color_hex_triplet'))
    return ccf_level,ccf_color

def dev_color_id(dev_level,dev_full):
    
    dev_color = []
    for item in dev_level:
        rgb = dev_full.loc[dev_full['Acronym']==item]
        r = rgb['R'].values[0]
        g = rgb['G'].values[0]
        b = rgb['B'].values[0]
        
        dev_color.append(f'rgba({r},{g},{b},1)')
        
    return dev_color

def get_terminals(ccf_nodes):
    
    terminal_node = np.array(['ID','Acronym','Full name'])
    greys = ccf_nodes['8'].descendants
    for node in greys:
        if node.is_leaf == True:
            print(node.acronym)
            terminal_node = np.append(terminal_node,np.array([node.name,node.acronym,node.full_name]),axis=0)

    return np.reshape(terminal_node,[-1,3])

ccf = sitk.ReadImage('../data/P56_CCFv3_Annotation_50um.nrrd',sitk.sitkFloat64)
ccf = sitk.GetArrayFromImage(sitk.PermuteAxes(ccf,[2,1,0]))
dev_ccf = sitk.ReadImage('../data/P56_KimLabDevCCFv002_Annotations_20um_50um.nii.gz',sitk.sitkFloat64)
dev_ccf = sitk.GetArrayFromImage(sitk.PermuteAxes(dev_ccf,[2,1,0]))

ccf_metadata,ccf_full = filter_metadata(ccf,'../data/CCF_metadata.xlsx')
dev_metadata,dev_full = filter_metadata(dev_ccf,'../data/AllenDevMouseOntologyStructure_DevCCFv004.xlsx',dev=True)

ccf_tree,ccf_nodes = json_to_tree('../data/adult_mouse_ontology.json')
dev_tree,dev_nodes = csv_to_tree('../data/AllenDevMouseOntologyStructure_DevCCFv004.csv')

dev_level = ['NeoCx','MesoCx','AlloCx','SPall','THy','PHy','p3','p2','p1','m1','m2','PPH','PH','PMH','MH']
dev_combined = combine_by_level(dev_ccf,dev_nodes,dev_level)
dev_color = dev_color_id(dev_level,dev_full)

ccf_level,ccf_color = ccf_level_generator('../data/687527670_mouse_brain_major_division_set.json')
ccf_label = ['Isocortex','OLF','HPF','CTXsp','STR','PAL','TH','HY','MB','P','MY','CB']
        
ccf_combined = combine_by_level(ccf,ccf_nodes,ccf_level,False)

mapping = generate_mapping(ccf_combined,dev_combined)

nodes2=mapping[:,0]-1
nodes1=[int(max(nodes2)+point) for point in mapping[:,1]]
nodes1_order = list(np.argsort(nodes1))
weights=np.log(mapping[:,2])

dev_label = [dev_level[int(x)] for x in mapping[:,1]-1]
ccf_label_df = [ccf_label[int(x)] for x in nodes2]

lineweights = [dev_color[int(x-12)] for x in nodes1]
alpha = 0.8
lineweights_alpha = [f'{x[:-2]}{alpha})' for x in lineweights]

fig = go.Figure(data=[go.Sankey(
    arrangement="perpendicular",
    node = dict(
        pad = 10,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label = ccf_label+dev_level,
        color = ccf_color+dev_color
    ),
    link = dict(
      source = nodes1, 
      target = nodes2,
      value = weights,
      color = lineweights_alpha
  ))])

fig.update_layout(width=600,
                  height=600,
                  font=dict(size = 12, color = 'white'),
                  plot_bgcolor='black',
                  paper_bgcolor='black')
fig.show(renderer='png')
fig.write_image('../figures/dev_ccf_sankey_250voxelthreshold.svg')

df = pd.DataFrame({'DevCCF':dev_label,
                    'DevCCF_id':nodes1,
                    'CCFv3':ccf_label_df,
                    'CCF_id':nodes2,
                    'pixel_overlap':mapping[:,2],
                    'Log2_overlap':weights})
df.to_csv('../data/Fig_3c_devCCF16_CCFv3_voxel_mapping.csv',index=False)

