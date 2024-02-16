"""
Set of functions used to generate MERFISH section images with DevCCF/CCFv3 labels.
Used in Fig 5d-e.
"""

import os, subprocess
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import requests
import SimpleITK as sitk
import pathlib


def plot_section( xx=None, yy=None, cc=None, val=None, pcmap=None, 
                 overlay=None, extent=None, bcmap=plt.cm.Greys_r, alpha=1.0,
                 fig_width = 6, fig_height = 6 ) :
    plt.style.use('dark_background')
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)

    if xx is not None and yy is not None and pcmap is not None :
        plt.scatter(xx,yy,s=0.5,c=val,marker='.',cmap=pcmap)
    elif xx is not None and yy is not None and cc is not None :
        plt.scatter(xx,yy,s=0.5,color=cc,marker='.',zorder=1)   
        
    if overlay is not None and extent is not None and bcmap is not None :
        plt.imshow(overlay, cmap=bcmap, extent=extent,alpha=alpha,zorder=2)
        
    ax.set_ylim(11,0)
    ax.set_xlim(0,11)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    
    return fig, ax

def image_info( img ) :
    print('size: ' + str(img.GetSize()) + ' voxels')
    print('spacing: ' + str(img.GetSpacing()) + ' mm' )
    print('direction: ' + str(img.GetDirection()) )
    print('origin: ' + str(img.GetOrigin()))

path_dir = '../data/'
url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/releases/20230630/manifest.json'
manifest = json.loads(requests.get(url).text)

view_directory = os.path.join( path_dir, 
                                manifest['directory_listing']['MERFISH-C57BL6J-638850-CCF']['directories']['metadata']['relative_path'], 
                              'views')

view_directory = pathlib.Path( view_directory )

cell_joined = pd.read_csv(os.path.join( view_directory, 'cell_metadata_with_parcellation_annotations.csv'),dtype={'cell_label':str})
cell_joined.set_index('cell_label',inplace=True)
z = 100*cell_joined['x_ccf'].to_numpy()
y = 100*cell_joined['y_ccf'].to_numpy()
x = 100*cell_joined['z_ccf'].to_numpy()

dev_ccf = sitk.ReadImage(path_dir+'KimLabDevCCFv001.nrrd',8)
dev_ccf_np = sitk.GetArrayViewFromImage(dev_ccf)

parcellation = []
for idx in range(len(x)):
    if x[idx]>=1017:
        point = dev_ccf_np[1140-int(x[idx]),int(y[idx]),int(z[idx])]
    else:
        point = dev_ccf_np[int(x[idx]),int(y[idx]),int(z[idx])]
    parcellation.append(point)

parcellation = np.array(parcellation,dtype=np.int32)
parcellation[parcellation>20000] = 126651558

cell_joined['dev_ccf_id'] = parcellation

dev_ccf_metadata = pd.read_csv(os.path.join( path_dir,'metadata', 'AllenDevMouseOntologyStructure_DevCCFv003.csv'))

dev_ccf_color = []
dev_ccf_name = []
for idx in range(len(parcellation)):
    
    color = dev_ccf_metadata.loc[dev_ccf_metadata['ID']==parcellation[idx]]
    if len(color)==0:
        color = '000000'
    else:
        color = color['Hexidecimal'].to_numpy()[0]
    dev_ccf_color.append('#'+color)
    
    if parcellation[idx]>0:
        name = dev_ccf_metadata.loc[dev_ccf_metadata['ID']==parcellation[idx],'Acronym'].to_list()[0]
        dev_ccf_name.append(name[2:-2])
    else:
        dev_ccf_name.append([])
    
cell_joined['dev_color'] = dev_ccf_color
cell_joined['dev_acronym'] = dev_ccf_name


brain_section = 'C57BL6J-638850.40'
pred = (cell_joined['brain_section_label'] == brain_section )
section = cell_joined[pred]

fig, ax = plot_section(xx=section['z_ccf'], yy=section['y_ccf'], 
                                cc=section['parcellation_structure_color'])
res = ax.set_title("Neuortransmitter - Reconstructed Coordinates")

volumes = manifest['file_listing']['MERFISH-C57BL6J-638850-CCF']['image_volumes']
print("reading resampled_average_template")
rpath = volumes['resampled_average_template']['files']['nii.gz']['relative_path']
file = os.path.join( path_dir, rpath)
average_template_image = sitk.ReadImage( file )
average_template_array = sitk.GetArrayViewFromImage( average_template_image )

print("reading resampled_annotation")
rpath = volumes['resampled_annotation']['files']['nii.gz']['relative_path']
file = os.path.join( path_dir, rpath)
annotation_image = sitk.ReadImage( file )
annotation_array = sitk.GetArrayViewFromImage( annotation_image )

print("reading resampled_annotation_boundary")
rpath = volumes['resampled_annotation_boundary']['files']['nii.gz']['relative_path']
file = os.path.join( path_dir, rpath)
annotation_boundary_image = sitk.ReadImage( file )
annotation_boundary_array = sitk.GetArrayViewFromImage( annotation_boundary_image )

size = average_template_image.GetSize()
spacing = average_template_image.GetSpacing()
extent = (-0.5 * spacing[0], (size[0]-0.5) * spacing[0], (size[1]-0.5) * spacing[1], -0.5 * spacing[1] )

zindex = int(section.iloc[0]['z_reconstructed'] / 0.2)
boundary_slice = annotation_boundary_array[zindex,:,:]

file = os.path.join( view_directory, 'cell_metadata_with_dev_CCF_parcellation_annotations.csv')
cell_joined.to_csv( file )

cell_joined = pd.read_csv(file,dtype={"cell_label":str})

brain_section = 'C57BL6J-638850.40'
pred = (cell_joined['brain_section_label'] == brain_section )
section = cell_joined[pred]

fig, ax = plot_section(section['x_reconstructed'], section['y_reconstructed'], 
                        cc=section['parcellation_structure_color'],
                        # overlay=boundary_slice, extent=extent, 
                        # bcmap=plt.cm.Greys, alpha = 1.0*(boundary_slice>0),
                        fig_width = 10, fig_height = 10 )
plt.savefig('../figures/CCF_section.png')

fig, ax = plot_section(section['x_reconstructed'], section['y_reconstructed'], 
                        cc=section['parcellation_structure_color'],
                        # overlay=boundary_slice, extent=extent, 
                        # bcmap=plt.cm.Greys, alpha = 1.0*(boundary_slice>0),
                        fig_width = 10, fig_height = 10 )
plt.savefig('../figures/CCF_section.png')

fig, ax = plot_section(section['x_reconstructed'], section['y_reconstructed'], 
                        cc=section['neurotransmitter_color'],
                        # overlay=boundary_slice, extent=extent, 
                        # bcmap=plt.cm.Greys, alpha = 1.0*(boundary_slice>0),
                        fig_width = 6, fig_height = 6 )
plt.savefig('../figures/neurotransmitter_section.png')

fig, ax = plot_section(section['x_reconstructed'], section['y_reconstructed'], 
                        cc=section['class_color'],
                        # overlay=boundary_slice, extent=extent, 
                        # bcmap=plt.cm.Greys, alpha = 1.0*(boundary_slice>0),
                        fig_width = 6, fig_height = 6 )
plt.savefig('../figures/class_section.png')

fig, ax = plot_section(section['x_reconstructed'], section['y_reconstructed'], 
                        cc=section['subclass_color'],
                        # overlay=boundary_slice, extent=extent, 
                        # bcmap=plt.cm.Greys, alpha = 1.0*(boundary_slice>0),
                        fig_width = 6, fig_height = 6 )
plt.savefig('../figures/subclass_section.png')