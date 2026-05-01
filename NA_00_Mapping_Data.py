# %%
'''
Network Analysis - 00 - Mapping Data to Correspondence Particles

This script maps data to the correspondence particles for downstream Network Analysis scripts.

Needed Files:
    project directory
    group folders and subsequent nested subject-specific folders
    mean surface (.vtk or .stl) and related correspondence model (nx3 particles file, .txt or .particles)
    feature map files (.txt, .csv, or .xlsx)
        1st column with correspondence particle id# and other columns with mapped data (each column represents a subsequent timepoint)
        

Multiple analyses and different data can be ran from the same project file, hence the necessary architecture
The script is looking for (from bottom to top) feature maps in their respective surface folders,
which are in their subject-specific folder, in their associated group folder, in one main directory named 'Data_Needs_Mapped'.
The architecture makes it quickly and readily available to import the data for each subject. So long as they all 
have surface and feature map folders with the same named schema. The individual files themselves (.stl, .csv, etc)
can be named whatever you want. It is their location within the folders that is what the code will look for.

Any folder with a name that needs to be specific is denoted with an asterisk (*)
Folders are denoted with a forward slash (/) and files will have a file extension (.stl, .csv, etc)

How the project directory should look...
Project Directory                                   (can be located anywhere, named anything, and this is what you will select when you first run the script)
    └─/Data_Needs_Mapped                            (./Data_Needs_Mapped contains all necessary data, structure is required for importing properly)    
        ├─/Subject_01                               (subject-specific folders)
        ├─/Subject_02
        ...
        └─Norm_12                                   [example]
            ├─/Surface_01
            ├─/Surface_02
            ...
            └─/Pelvis                               [example]
                ├─/Particles_File
                ├─/Surface_File
                ...
                ├─/norm12_pelvis_local.particles    [example]
                └─/norm12_pelvis.stl                [example]
            └─Norm_02_FE_Mean.csv                   [example] (sits at same level as surface folders, data to be mapped)

Created by:     Rich J. Lisonbee, MS
Organization:   University of Utah, Orthopaedic Research Laboratory
PI:             Andrew E. Anderson, PhD
Date:           4/1/2026

Modified by:
Date:

'''

import numpy as np
import pandas as pd
import trimesh
import scipy
import scipy.spatial
import scipy
import tkinter as tk
import FreeSimpleGUI as sg
import os
import warnings

os.system('cls||clear')
np.random.seed(0)
print('Network Analysis - 00')


# %%
# Project Directory
sg.popup('Please select the Project Directory with the Data_Needs_Mapped folder and nested subject-specific folders containing the particles file (.txt or .particles), feature map (.csv), and surface (.stl) files', keep_on_top=True)
proj_dir = tk.filedialog.askdirectory(title="Project Directory (don't select Data_Needs_Mapped)", initialdir=os.getcwd())
# nmap_dir = os.path.join(proj_dir,'Data_Needs_Mapped')
sg.popup('Please select the Data_Needs_Mapped folder', keep_on_top=True)
nmap_dir = tk.filedialog.askdirectory(title="Folder of data that needs mapped", initialdir=proj_dir)

if not os.path.exists(nmap_dir):
    warnings.warn("Required folder '\Data_Needs_Mapped' does not exist!")
    
out_dir = os.path.join(proj_dir,'Data')


# %%
subject_folders = os.listdir(nmap_dir)
surf_option     = [surf for surf in os.listdir(os.path.join(nmap_dir,subject_folders[0])) if os.path.isdir(os.path.join(nmap_dir,subject_folders[0],surf))]
part_option     = [part for part in os.listdir(os.path.join(nmap_dir,subject_folders[0])) if not os.path.isdir(os.path.join(nmap_dir,subject_folders[0],part))]

if len(part_option) > 1:
    warnings.warn("Found more than one feature map file")

# ***User Inputs***
layout  = [
    [sg.Text('How does the data need to be mapped?')],
    [sg.Radio('SurfaceA to ParticlesA', group_id=1, default=True),
     sg.Radio('SurfaceB to ParticlesA (front facing)', group_id=1),
     sg.Radio('SurfaceB to ParticlesA (back facing)', group_id=1)],
    [sg.Text('example:  SurfaceA = Pelvis, ParticlesA = Pelvis Particles, SurfaceB = Cartilage')],
    [sg.Text('          SurfaceB to ParticlesA (front facing) will map with surfaces facing one another (i.e., osteochondral boundary)')],
    [sg.Text('          SurfaceB to ParticlesA (back facing) will map with surfaces facing opposite one another (i.e., articular cartilage surface back to bone)')],
    [sg.Text('Enter the name of the Feature Map (for outputs)')],
    [sg.Input('FE_Data',
        key='-FEAT-NAME-',
        size=(40,1))],
    [sg.Text('Process specific frame? (leave blank if you want to process full activity, otherwise enter frame number starting at 1 for the first frame)')],
    [sg.Input('',
        key='-FRAME-',
        size=(40,1))],
    # [sg.Text('Move surface file to the output folder?')],
    # [sg.Checkbox('[check for yes to move, leave blank for no to make a copy instead]',
    #     key='-MOVE_SURF-')],
    [sg.Ok(), sg.Cancel()]
]
window = sg.Window("Mapping Inputs", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()

# move_surf   = values['-MOVE_SURF-']
feat_name   = values['-FEAT-NAME-']
if values['-FRAME-'] != '':
    frame_sel   = abs(int(values['-FRAME-']))
else:
    frame_sel = -1
# ***

map_type = [k for k, v in values.items() if v is True][0]

if map_type == 0:
    # ***Surface Selection UI***
    layout  = [
        [sg.Text('Pick SurfaceA (to be mapped to):')],
        [sg.Listbox(surf_option,
            key='-SURFA-',
            select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
            size=(40, 3))],
        [sg.Text('Distance limit from face centers to include particle in analysis (set to 0 to grab all particles, ignoring distance):')],
        [sg.Input('0',
            key='-DISTA-',
            size=(40, 3))],  
        [sg.Ok(), sg.Cancel()]
    ]
    
    window = sg.Window("Surface Selection", layout, modal=True, keep_on_top=True)
    event, values = window.read()
    window.close()
    
    surfA_id = values['-SURFA-'][0]
    dist_roi = float(values['-DISTA-'][0])
    #  ***
else:
    # ***Surface Selection UI***
    layout  = [
        [sg.Text('Pick SurfaceA (to be mapped to):')],
        [sg.Listbox(surf_option,
            key='-SURFA-',
            select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
            size=(40, 3))],  
        [sg.Text('Pick SurfaceB (to be mapped from):')],
        [sg.Listbox(surf_option,
            key='-SURFB-',
            select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
            size=(40, 3))],          
        [sg.Ok(), sg.Cancel()]
    ]
    
    window = sg.Window("Surface Selection", layout, modal=True, keep_on_top=True)
    event, values = window.read()
    window.close()
    
    surfA_id = values['-SURFA-'][0]
    surfB_id = values['-SURFB-'][0]
    #  ***


# %%
cp_list = []
loc     = []
meshA   = []
meshB   = []
partA   = []
cp_list = []
id_list = []
data    = []
subject_folders = os.listdir(nmap_dir)
n       = 0
# todo: save mappings so users can load them rather than needing to process multiple times for different features...
# it is a lot of bookkeeping and file management
for subject in subject_folders:
    print(subject)
    
    data_name = [name for name in os.listdir(os.path.join(nmap_dir,subject)) if name.endswith(('.csv','.xlsx'))]
    if len(data_name) > 1:
        warnings.warn(f"{subject} has more than one feature data in its folder")
    if data_name[0].endswith('.csv'):
        data.append(np.genfromtxt(os.path.join(nmap_dir,subject,data_name[0]),delimiter=','))
    elif data_name[0].endswith('.xlsx'):
        data.append(pd.read_excel(os.path.join(nmap_dir,subject,data_name[0])).to_numpy())
    
    if map_type == 0:
        # SurfaceA to ParticlesA
        # Load Surfaces
        meshA_name = [name for name in os.listdir(os.path.join(nmap_dir,subject,surfA_id)) if name.endswith(('.stl','.obj','.ply'))]
        if len(meshA_name) > 1:
            warnings.warn(f"{subject} has more than one surface in the {surfA_id} folder")
            
        # Load Particles
        partA_name = [name for name in os.listdir(os.path.join(nmap_dir,subject,surfA_id)) if name.endswith(('.particles','.txt'))]
        if len(partA_name) > 1:
            warnings.warn(f"{subject} has more than one particles file in the {surfA_id} folder")
        partA.append(np.loadtxt(os.path.join(nmap_dir, subject, surfA_id, partA_name[0])))
        
        # switches to mesh B for later logic consistency...
        meshB.append(trimesh.load(os.path.join(nmap_dir, subject, surfA_id, meshA_name[0])))
        
        if np.abs(dist_roi) > 0:
            tree    = scipy.spatial.cKDTree(partA[n])
            ind1    = tree.query_ball_point(meshB[n].triangles_center,r=2)
    
            cp_list.append(np.unique(np.array(np.concatenate(ind1),dtype=int)))
        else:
            cp_list.append(list(range(0,len(partA[n]),1)))
        
    else:
        # Load Surfaces
        meshA_name = [name for name in os.listdir(os.path.join(nmap_dir,subject,surfA_id)) if name.endswith(('.stl','.obj','.ply'))]
        if len(meshA_name) > 1:
            warnings.warn(f"{subject} has more than one surface in the {surfA_id} folder")
            
        meshB_name = [name for name in os.listdir(os.path.join(nmap_dir,subject,surfB_id)) if name.endswith(('.stl','.obj','.ply'))]
        if len(meshA_name) > 1:
            warnings.warn(f"{subject} has more than one surface in the {surfB_id} folder")
        
        meshA.append(trimesh.load(os.path.join(nmap_dir, subject, surfA_id, meshA_name[0])))
        meshB.append(trimesh.load(os.path.join(nmap_dir, subject, surfB_id, meshB_name[0])))
        
        # Load Particles
        partA_name = [name for name in os.listdir(os.path.join(nmap_dir,subject,surfA_id)) if name.endswith(('.particles','.txt'))]
        if len(partA_name) > 1:
            warnings.warn(f"{subject} has more than one particles file in the {surfA_id} folder")
        partA.append(np.loadtxt(os.path.join(nmap_dir, subject, surfA_id, partA_name[0])))
        
        # ray tracing intersection
        face_centers = meshA[n].triangles_center
        face_normals = meshA[n].face_normals
        
        loc1, i_ray, i_tri = meshB[n].ray.intersects_location(ray_origins=face_centers, ray_directions=face_normals)
    
        tree    = scipy.spatial.cKDTree(partA[n])
        ind1    = tree.query_ball_point(loc1,r=1)
        
        if map_type == 1:
            # SurfaceB to ParticlesA (front facing)
            loc2, i_ray, i_tri = meshB[n].ray.intersects_location(ray_origins=face_centers, ray_directions=face_normals)
        elif map_type == 2:
            # SurfaceB to ParticlesA (back facing)
            loc2, i_ray, i_tri = meshB[n].ray.intersects_location(ray_origins=face_centers, ray_directions=-face_normals)
        tree    = scipy.spatial.cKDTree(partA[n])
        ind2    = tree.query_ball_point(loc2,r=1)
        
        indices = np.concatenate((ind1, ind2))
        
        locations = np.vstack((loc1,loc2))

        cp_list.append(np.unique(np.array(np.concatenate(indices),dtype=int)))
        loc.append(locations)
    n += 1

sets = [set(arr) for arr in cp_list]
temp = set.intersection(*sets)
i_cp = np.array(list(temp))

print(len(i_cp))


# %% Data
max_data    = []
mean_data   = []
near_data   = []
for n in range(len(subject_folders)):
    # remove rows with nans in case there was header information 
    data[n] = data[n][~np.isnan(data[n]).all(axis=1)]
    face_centers        = meshB[n].triangles_center
    face_ids            = np.array(data[n][:,0],dtype=int)
    tree                = scipy.spatial.cKDTree(partA[n][i_cp])
    distances, indices  = tree.query(face_centers[face_ids])
    
    if frame_sel > 0:
        local_max   = np.empty((0, 1))
        local_mean  = np.empty((0, 1))
    else:
        local_max   = np.empty((0, len(data[n][0,1:])))
        local_mean  = np.empty((0, len(data[n][0,1:])))
    i = np.unique(indices)
    for idx in range(len(i)):
        i_rows = np.where(indices == i[idx])
        x_max = []
        x_mean = []
        if frame_sel > 0:
            col = frame_sel
            x_mean.append(np.mean(data[n][i_rows,col]))
            x_max.append(np.max(data[n][i_rows,col]))
        elif frame_sel == -1:
            for col in range(1,len(data[n][0,:])):
                x_mean.append(np.mean(data[n][i_rows,col]))
                x_max.append(np.max(data[n][i_rows,col]))
        x_max   = np.array(x_max)
        x_mean  = np.array(x_mean)
      
        local_max   = np.vstack((local_max,   x_max))
        local_mean  = np.vstack((local_mean, x_mean))
    max_data.append(local_max)
    mean_data.append(local_mean)
    
    face_centers        = meshB[n].triangles_center
    face_ids            = np.array(data[n][:,0],dtype=int)
    tree                = scipy.spatial.cKDTree(face_centers[face_ids])
    distances, indices  = tree.query(partA[n][i_cp])
    
    t = data[n][indices,1:]
    t = np.array(t)   
    
    near_data.append(t)


# %% Save Data
for n in range(len(subject_folders)):
    subject = subject_folders[n]
    
    subj_out = os.path.join(out_dir,subject,surfA_id,feat_name + '_Mean')
    os.makedirs(subj_out,exist_ok=True)
    df = pd.DataFrame(np.hstack((i_cp.reshape(-1,1),mean_data[n])))
    df.to_excel(os.path.join(subj_out, feat_name + '_Mean_Results.xlsx'),index = False, header=False)
    
    subj_out = os.path.join(out_dir,subject,surfA_id,feat_name + '_Max')
    os.makedirs(subj_out,exist_ok=True)
    df = pd.DataFrame(np.hstack((i_cp.reshape(-1,1),max_data[n])))
    df.to_excel(os.path.join(subj_out, feat_name + '_Max_Results.xlsx'),index = False, header=False) 
    
    subj_out = os.path.join(out_dir,subject,surfA_id,feat_name + '_Near')
    os.makedirs(subj_out,exist_ok=True)
    df = pd.DataFrame(np.hstack((i_cp.reshape(-1,1),near_data[n])))
    df.to_excel(os.path.join(subj_out, feat_name + '_Near_Results.xlsx'),index = False, header=False)  

print(feat_name)
print(f'Results Output: {out_dir}')
print('Complete!')
print('')
print('Consolidate the output subject folders into their respective group folders and you are ready for NA_01')


