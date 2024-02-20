"""
Set of functions to access Allen Brain cell (ABC) atlas, download MERFISH data, and annotate with DevCCF labels.
"""
import os, json, requests, pathlib, subprocess, time
import pandas as pd
import matplotlib.pyplot as plt

def download_file( file_dict ) :
    
    print(file_dict['relative_path'],file_dict['size'])
    local_path = os.path.join( path_dir, file_dict['relative_path'] )
    local_path = pathlib.Path( local_path )
    remote_path = manifest['resource_uri'] + file_dict['relative_path']

    command = "aws s3 cp --no-sign-request %s %s" % (remote_path, local_path)
    start = time.process_time()
    result = subprocess.run(command.split(),stdout=subprocess.PIPE)
    print("time taken: ", time.process_time() - start)



url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/releases/20230630/manifest.json'
manifest = json.loads(requests.get(url).text)

path_dir = '../data/'

for r in manifest['directory_listing'] :
    
    r_dict =  manifest['directory_listing'][r]
    
    for d in r_dict['directories'] :
        
        if d != 'metadata' :
            continue
        
        d_dict = r_dict['directories'][d]
        local_path = os.path.join( path_dir, d_dict['relative_path'])
        local_path = pathlib.Path( local_path )
        remote_path = manifest['resource_uri'] + d_dict['relative_path']
        
        command = "aws s3 sync --no-sign-request %s %s" % (remote_path, local_path)
        result = subprocess.run(command.split(),stdout=subprocess.PIPE)


view_directory = os.path.join( path_dir, 
                                manifest['directory_listing']['MERFISH-C57BL6J-638850-CCF']['directories']['metadata']['relative_path'], 
                              'views')

view_directory = pathlib.Path( view_directory )
cache_views = True
if cache_views :
    os.makedirs( view_directory, exist_ok=True )
    

metadata = manifest['file_listing']['MERFISH-C57BL6J-638850']['metadata']


rpath = metadata['cell_metadata_with_cluster_annotation']['files']['csv']['relative_path']
file = os.path.join( path_dir, rpath)
cell = pd.read_csv(file,dtype={"cell_label":str})
cell.rename(columns={'x':'x_section','y':'y_section','z':'z_section'},inplace=True)
cell.set_index('cell_label',inplace=True)

metadata = manifest['file_listing']['MERFISH-C57BL6J-638850-CCF']['metadata']
rpath = metadata['reconstructed_coordinates']['files']['csv']['relative_path']
file = os.path.join( path_dir, rpath)
reconstructed_coords = pd.read_csv(file,dtype={"cell_label":str})
reconstructed_coords.rename(columns={'x':'x_reconstructed','y':'y_reconstructed','z':'z_reconstructed'},inplace=True)
reconstructed_coords.set_index('cell_label',inplace=True)
cell_joined = cell.join(reconstructed_coords,how='inner')

metadata = manifest['file_listing']['MERFISH-C57BL6J-638850-CCF']['metadata']
rpath = metadata['ccf_coordinates']['files']['csv']['relative_path']
file = os.path.join( path_dir, rpath)
ccf_coords = pd.read_csv(file,dtype={"cell_label":str})
ccf_coords.rename(columns={'x':'x_ccf','y':'y_ccf','z':'z_ccf'},inplace=True)
ccf_coords.drop(['parcellation_index'],axis=1,inplace=True)
ccf_coords.set_index('cell_label',inplace=True)
cell_joined = cell_joined.join(ccf_coords,how='inner')

metadata = manifest['file_listing']['Allen-CCF-2020']['metadata']
rpath = metadata['parcellation_to_parcellation_term_membership_acronym']['files']['csv']['relative_path']
file = os.path.join( path_dir, rpath)
parcellation_annotation = pd.read_csv(file)
parcellation_annotation.set_index('parcellation_index',inplace=True)
parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]

rpath = metadata['parcellation_to_parcellation_term_membership_color']['files']['csv']['relative_path']
file = os.path.join( path_dir, rpath)
parcellation_color = pd.read_csv(file)
parcellation_color.set_index('parcellation_index',inplace=True)
parcellation_color.columns = ['parcellation_%s'% x for x in  parcellation_color.columns]

cell_joined = cell_joined.join(parcellation_annotation,on='parcellation_index')
cell_joined = cell_joined.join(parcellation_color,on='parcellation_index')

if cache_views :
    file = os.path.join( view_directory, 'cell_metadata_with_dev_CCF_parcellation_annotations.csv')
    cell_joined.to_csv( file )
 