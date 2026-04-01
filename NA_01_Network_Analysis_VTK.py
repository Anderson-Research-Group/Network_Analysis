# %%
'''
Network Analysis - 01 - Statistical Analysis of Spatiotemporal Data

This script performs Network Analysis on user feature maps and correspondence models.
Important to note that this workflow will not perform pairing or mapping of the
feature maps data for you, and requires that you accomplish this before utilizing 
these packages. This is due to the highly variable and complex nature of mapping 
these results. 
Once data has been mapped to the particles (via nearest neighbor, focal operations,
local maxima or minima, etc.) in a way that you would like to investigate, this
script will allow for spatiotemporal comparison of the mapped data.

Types of tests available:
    paired t-tests      (either between two groups, or between two feature maps)
    two-sample t-tests  (between groups)
    *if you choose to combine groups (select multiple group folders) it assumes
    that you are performing a paired test between two feature maps
    *all tests assume two-tailed
    *it is highly encouraged that when using circular shift it is appropriate for your application

Needed Files:
    project directory
    group folders and subsequent nested subject-specific folders
    mean surface (.vtk or .stl) and related correspondence model (nx3 particles file, .txt or .particles)
    feature map files (.txt, .csv, or .xlsx)
        1st column with correspondence particle id# and other columns with mapped data (each column represents a subsequent timepoint)
        

Multiple analyses and different data can be ran from the same project file, hence the necessary architecture
The script is looking for (from bottom to top) feature maps in their respective surface folders,
which are in their subject-specific folder, in their associated group folder, in one main directory named 'Data'.
The architecture makes it quickly and readily available to import the data for each subject. So long as they all 
have surface and feature map folders with the same named schema. The individual files themselves (.stl, .csv, etc)
can be named whatever you want. It is their location within the folders that is what the code will look for.

Any folder with a name that needs to be specific is denoted with an asterisk (*)
Folders are denoted with a forward slash (/) and files will have a file extension (.stl, .csv, etc)

How the project directory should look...
Project Directory                               (can be located anywhere, named anything, and this is what you will select when you first run the script)
    └─/Data*                                    (./Data contains all necessary data, structure is required for importing properly)    
        ├─/Group_A                              (group folders)
        ├─/Group_B
        ...
        └─/Controls                             [example]
            ├─/Subject_01                       (subject-specific folders)
            ├─/Subject_02
            ...
            └─Norm_12                           [example]
                ├─/Surface_01
                ├─/Surface_02
                ...
                └─/Pelvis                       [example]
                    ├─/Mapped_Data_01
                    ├─/Mapped_Data_02
                    ...
                    └─/FE_Data                  [example] (only has one file in it)
                        └─Norm_02_FE_Mean.xlsx  [example]

Original Work:  Penny R. Atkins, PhD
                doi: 10.1007/s10439-023-03270-6
Adapted by:     Rich J. Lisonbee, MS and Bergen Braun, MS
Organization:   University of Utah, Orthopaedic Research Laboratory
PI:             Andrew E. Anderson, PhD
Date:           1/15/2025

Modified by:
Date:

'''

import numpy as np
import pandas as pd
import trimesh
import pyvista as pv
import itertools
import open3d as o3d
import scipy
import spm1d
import scipy.spatial
import scipy
import tkinter as tk
import FreeSimpleGUI as sg
import time
from pathlib import Path
import os, sys

os.system('cls||clear')
np.random.seed(0)
print('Network Analysis - 01')


# %%
# ***User Inputs***
proj_dir_choice = False
layout  = [
    [sg.Text('Run script from project file?')],
    [sg.Radio('Yes', group_id=1),
     sg.Radio('No', group_id=1, default=True)],
    [sg.Ok(), sg.Cancel()]
]
window = sg.Window("Project File/Directory", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()
# ***

if event == 'Cancel' or event == sg.WIN_CLOSED:
    sys.exit()
    
if values[0]:
    proj_dir_choice = True
#  ***

if proj_dir_choice:
    root = tk.Tk()
    # root.window()
    root.attributes('-topmost',True)
    proj_file = tk.filedialog.askopenfile(mode='r',title='Select project file', initialdir=os.getcwd(),parent=root)
    root.destroy()
    
    proj_nam = os.path.basename(proj_file.name)
    proj_dir = proj_file.name.replace(proj_nam,'')
    proj_file.close()
else:
    # Project Directory
    sg.popup('Please select the Project Directory', keep_on_top=True)
    proj_dir = tk.filedialog.askdirectory(title="Select a folder", initialdir=os.getcwd())


# %% Group Selection
if proj_dir_choice:
    # find available group names (project file)
    if proj_nam.split('.')[-1] == 'xlsx':
        project_file = pd.read_excel(os.path.join(proj_dir,proj_nam))
    elif proj_nam.split('.')[-1] == 'csv':
        project_file = pd.read_csv(os.path.join(proj_dir,proj_nam))
    
    # groups to include
    group_ids = project_file.columns[project_file.columns.str.contains('group_')].tolist()
    
    # ***Group IDs Selection UI***
    layout  = [
        [sg.Text('Pick Group ID from project:')],
        [sg.Listbox(group_ids,
            key='-GROUP-ID-',
            select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
            size=(40, 6))],   
        [sg.Ok(), sg.Cancel()]
    ]

    window = sg.Window("Group ID Selection", layout, modal=True, keep_on_top=True)
    event, values = window.read()
    window.close()

    group_id = values['-GROUP-ID-']
    #  ***
    
    group_names_avail = list(set(project_file[group_id[0]].to_list()))
    
else:
    # find available group names (directories)
    group_names_avail = os.listdir(os.path.join(proj_dir,'Data'))
    

# ***Group Selection UI***
group_ui_text = 'Pick Group(s) to include'
if proj_dir_choice:
    group_ui_text += f' from ID {group_id[0]}'
group_ui_text += ':'

layout  = [
    [sg.Text(group_ui_text)],
    [sg.Listbox(group_names_avail,
        key='-GROUP-',
        select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,    # multi-select
        size=(40, 6))],   
    [sg.Ok(), sg.Cancel()]
]

window = sg.Window("Group Selection", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()

group_names = values['-GROUP-']
#  ***


# %% Surfaces Available
if proj_dir_choice:
    # find available surfaces names (project file)
    surface_names_avail = project_file.columns[project_file.columns.str.contains('shape_')].tolist()

else:
    # find available surfaces names (directories)
    surface_names = []
    for group_count in range(len(group_names)):
        print(group_names[group_count])
        subj_names = []
        subj_names = [name for name in os.listdir(os.path.join(proj_dir,'Data',group_names[group_count])) if os.path.isdir(os.path.join(proj_dir,'Data',group_names[group_count],name))]
        for subj_count in range(len(subj_names)):
            surf_names = [name for name in os.listdir(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count])) if os.path.isdir(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],name))]
            surface_names.append(surf_names)
            
    surface_names_avail = list(set(itertools.chain.from_iterable(surface_names)))


# %% User Inputs
tests_avail = ['two-tailed t-test', 'paired t-test']
# ***User Inputs***
layout  = [
    [sg.Text('Alpha Value:')],
    [sg.Input('0.05',
        key='-ALPHA-',
        size=(40,1))],
    [sg.Text('Number of Iterations (null distribution):')],
    [sg.Input('10000',
        key='-NUM-ITER-',
        size=(40,1))],
    [sg.Text('Would you like to add regularization?')],
    [sg.Checkbox('Regularization:',
        key='-REGU-OPT-')],
    [sg.Text('Regularization (if yes above)[10^n]')],
    [sg.Input('-12',
        key='-REGU-',
        size=(40,1))],    
    [sg.Text('Test Type:')],
    [sg.OptionMenu(tests_avail,
        key='-TEST-TYPE-',
        default_value=tests_avail[0],               # single-select
        size=(40, 1))],
    [sg.Text('Surface Selection:')],
    [sg.OptionMenu(surface_names_avail,
        key='-SURFACE-',
        default_value=surface_names_avail[0],       # single-select
        size=(40, 1))],
    [sg.Text('Is it a volumetric correspondence model?')],
    [sg.Checkbox('Volumetric:',
        key='-VOLUMETRIC-')],
    [sg.Text('Circular Shift permutations? (temporal data only, verify if appropriate)')],
    [sg.Checkbox('Circular Shift:',
        key='-SHIFT-OPT-')],    
    [sg.Ok(), sg.Cancel()]
]

window = sg.Window("Select Options", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()
if event == 'Cancel' or event == sg.WIN_CLOSED:
    sys.exit()

alpha_val       = float(values['-ALPHA-'])
num_iterations  = int(values['-NUM-ITER-'])
test_type       = tests_avail.index(values['-TEST-TYPE-'])
surf_name       = values['-SURFACE-']
shift_opt       = values['-SHIFT-OPT-']
regu_opt        = values['-REGU-OPT-']
if regu_opt:
    regu_val    = float(values['-REGU-'])
vol_opt         = values['-VOLUMETRIC-']
# ***

# check if the groups selected are paired, or wish to combine groups (would then need two feature maps for paired comparison)
combine_groups = False
if test_type == 1 and len(group_names) > 1:
    # ***Combine Groups Check***
    layout  = [
        [sg.Text('Combine Groups?')],
        [sg.Radio('Yes', group_id=1),
         sg.Radio('No', group_id=1, default=True)],
        [sg.Ok(), sg.Cancel()]
    ]
    window = sg.Window("Combination", layout, modal=True)
    event, values = window.read()
    window.close()
    # ***
    
    if event == 'Cancel' or event == sg.WIN_CLOSED:
        sys.exit()
        
    if values[0]:
        combine_groups = True

print(tests_avail[test_type])

    
# %% Mean Surface and Points Selection
sg.popup('Please select the Mean Surface file', keep_on_top=True)
mean_model_name = tk.filedialog.askopenfilename(title="Select Mean Surface", initialdir=proj_dir,filetypes=[("stl",'*.stl'),("vtk",'*.vtk')])

sg.popup('Please select the Mean Particles file for:\n' f'{mean_model_name}', keep_on_top=True)
mean_point_name = tk.filedialog.askopenfilename(title="Select Mean Particles", initialdir=mean_model_name.split('/')[:-1],filetypes=[("particles",'*.particles'),("txt",'*.txt')])


# %% Feature Map Selection
if proj_dir_choice:
    # find available feature names (project file)
    # project_file.columns[project_file.columns.str.contains('shape_')].tolist()
    # make the assumption the only scalar field
    p = project_file[surf_name].tolist()[project_file[group_id[0]].tolist().index(group_names[0])]
    p = Path(str(p).strip()).expanduser()
    if p.is_absolute() is False:
        p = (proj_dir / p).resolve()
    
    mesh = pv.read(p)
    values['feature1'] = mesh.active_scalars_name

else:
    # find available feature names (directories)
    feature_names = []
    for group_count in range(len(group_names)):
        subj_names = []
        subj_names = [name for name in os.listdir(os.path.join(proj_dir,'Data',group_names[group_count])) if os.path.isdir(os.path.join(proj_dir,'Data',group_names[group_count],name))]
        for subj_count in range(len(subj_names)):
            feat_names = [name for name in os.listdir(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name)) if os.path.isdir(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name,name))]
            feature_names.append(feat_names)
    feature_names_avail = list(set(itertools.chain.from_iterable(feature_names)))

    # ***Feature Selection UI***
    layout  = [
        [sg.Text('Pick Feature Map:')],
        [sg.OptionMenu(feature_names_avail,
            key='feature1',
            default_value=feature_names_avail[0],      # single-select
            size=(40, 1))],
    ]
    
    # if needing to grab two feature maps
    if combine_groups or len(group_names) == 1 and len(feature_names_avail) > 1:
        layout.append(
            [[sg.Text('Pick 2nd Feature Map (paired data):')],
             sg.OptionMenu(feature_names_avail,
                key='feature2',
                default_value=feature_names_avail[1],      # single-select
                size=(40, 1))]
            )
    layout.append([sg.Ok(), sg.Cancel()])
    
    window = sg.Window("Group Selection", layout, modal=True, keep_on_top=True)
    event, values = window.read()
    window.close()

feature_selection = [values['feature1']]
print(f'Feature Selected: {feature_selection[0]}')
if combine_groups or len(group_names) == 1:
    feature_selection.append(values['feature2'])
    print(f'Additional Feature Selected: {feature_selection[1]}')
#  ***


# %% Compile Data
compiled_data   = []
group_ids       = []
if proj_dir_choice:
    subj_names = project_file[group_id[0]].tolist()
    for subj_count in range(len(subj_names)):
        for group_count in range(len(group_names)):
            if subj_names[subj_count] == group_names[group_count]:
                p = project_file[surf_name].tolist()[subj_count]
                p = Path(str(p).strip()).expanduser()
                if p.is_absolute() is False:
                    p = (proj_dir / p).resolve()
                mesh = pv.read(p)
                
                temp_data = np.array(mesh[mesh.active_scalars_name]).reshape(-1,1)
                
                if regu_opt:
                    temp_data += np.power(10,regu_val) * np.random.rand(*temp_data.shape)
                    
                compiled_data.append(temp_data)
                group_ids.append(group_count)
    compiled_data = np.array(compiled_data)
    group_ids   = np.array(group_ids)
                
else:
    for group_count in range(len(group_names)):
        subj_names = []
        subj_names = [name for name in os.listdir(os.path.join(proj_dir,'Data',group_names[group_count])) if os.path.isdir(os.path.join(proj_dir,'Data',group_names[group_count],name))]
        for subj_count in range(len(subj_names)):
            feat_name = os.listdir(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name,feature_selection[0]))
            if feat_name[0].split('.')[-1] == 'xlsx':
                temp_data = pd.read_excel(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name,feature_selection[0],feat_name[0]), header=None).to_numpy()
            elif feat_name[0].split('.')[-1] == 'csv':
                temp_data = pd.read_csv(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name,feature_selection[0],feat_name[0]), header=None).to_numpy()
                
            if regu_opt:
                temp_data[:,1:] += np.power(10,regu_val) * np.random.rand(*temp_data[:,1:].shape)
            
            compiled_data.append(temp_data)
            if combine_groups or len(group_names) == 1 and len(feature_names_avail) > 1:
                group_ids.append(0)
                feat_name = os.listdir(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name,feature_selection[1]))
                temp_data = pd.read_excel(os.path.join(proj_dir,'Data',group_names[group_count],subj_names[subj_count],surf_name,feature_selection[1],feat_name[0]), header=None).to_numpy()
                # ADD TEXT AND CSV FILETYPES
                if regu_opt:
                    temp_data[:,1:] += np.power(10,regu_val) * np.random.rand(*temp_data[:,1:].shape)
                
                compiled_data.append(np.array(temp_data))
                group_ids.append(1)
            else:
                group_ids.append(group_count)
    
    group_ids   = np.array(group_ids)

def common_cp(arrays):
    common = np.unique(arrays[0][:, 0])
    
    for a in arrays[1:]:
        common = np.intersect1d(common, np.unique(a[:, 0]))
        
    first_flags = arrays[0][:, 0]
    pos = {flag: np.where(first_flags == flag)[0][0] for flag in common}
    cp_ids = np.array(sorted(common, key=lambda f: pos[f]),dtype=int)

    filtered_array = []
    for a in arrays:
        flags = a[:, 0]
        mask = np.isin(flags, cp_ids)
        a_filtered = a[mask]

        idx_map = {}
        filt_indices = np.nonzero(mask)[0]
        for i, orig_idx in enumerate(filt_indices):
            f = flags[orig_idx]
            idx_map.setdefault(f, []).append(i)

        parts = []
        for f in cp_ids:
            if f in idx_map:
                parts.append(a_filtered[idx_map[f]])
        if parts:
            a_reordered = np.vstack(parts)
        else:
            a_reordered = a_filtered

        filtered_array.append(a_reordered)
    return np.array(cp_ids), np.array(filtered_array)

if proj_dir_choice:
    cp_ids = np.arange(0,len(compiled_data))
    filtered_array = compiled_data
    
    all_data = []
    for s in range(filtered_array.shape[0]):
        all_data.append(filtered_array[s,:,:])
    
else:
    cp_ids, filtered_array = common_cp(compiled_data)

    all_data = []
    for s in range(filtered_array.shape[0]):
        all_data.append(filtered_array[s,:,1:])
    
all_data = np.stack(all_data, axis=2)
print(all_data.shape)
num_pts     = all_data.shape[0]
particles   = np.ones((num_pts,0),dtype=bool)
timepoints  = all_data.shape[1]


#%% Traditional SPM 
A = np.copy(np.transpose(all_data[:, :, np.where(group_ids == 0)[0]].reshape(num_pts*timepoints, len(np.where(group_ids == 0)[0]))))
B = np.copy(np.transpose(all_data[:, :, np.where(group_ids == 1)[0]].reshape(num_pts*timepoints, len(np.where(group_ids == 1)[0])))) 
if test_type == 0:
    # nonparametric two-tailed t-test
    snpm = spm1d.stats.nonparam.ttest2(A,B)
    # degree of freedom
    df = len(np.where(group_ids == 0)[0]) + len(np.where(group_ids == 1)[0]) - 2
    snpmi = snpm.inference(alpha_val, two_tailed=True, iterations=10000, force_iterations=True)

if test_type == 1:
    # paired two-tailed t-test 
    snpm = spm1d.stats.ttest_paired(A, B)
    # degree of freedom
    df = len(np.where(group_ids == 0)[0]) - 1
    snpmi = snpm.inference(alpha_val, two_tailed=True)
    
Z           = snpmi.z       # flattened test statistic  (i.e., t value) over only non-zero-variance nodes
tradzstar   = snpmi.zstar   # critical test statistic  (i.e., critical t value)
fvalues     = snpm.z.reshape(len(particles), timepoints)

tradzsig    = Z.copy()
tradzsig[np.abs(Z) < tradzstar] = 0
tradzsig = tradzsig.reshape(len(particles),timepoints)


#%% Connectivity
mean_mesh       = trimesh.load(mean_model_name)

mesh_normals    = mean_mesh.face_normals
mesh_points     = mean_mesh.vertices

mean_point      = np.loadtxt(mean_point_name)
mean_shape      = mean_point.reshape(-1, 3)

if not vol_opt:
    kdtree          = scipy.spatial.KDTree(mesh_points)
    dist, pts_index = kdtree.query(mean_shape)
    pts_index.mean()
    
    pcd             = o3d.geometry.PointCloud()
    pcd.points      = o3d.utility.Vector3dVector(mean_shape)
    pcd.normals     = o3d.utility.Vector3dVector(mesh_normals[pts_index, :])
    
    distances       = pcd.compute_nearest_neighbor_distance()
    avg_dist        = np.mean(distances)
    radius          = avg_dist * 1.5
    
    mesh            = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))
    tris            = np.asarray(mesh.triangles)
elif vol_opt:
 # todo: update when we have example dataset
     def subdivide_cube_to_tetrahedra(corners):
         return [
             [corners[0], corners[1], corners[3], corners[4]],
             [corners[1], corners[2], corners[3], corners[6]],
             [corners[1], corners[4], corners[5], corners[6]],
             [corners[3], corners[4], corners[6], corners[7]],
             [corners[1], corners[3], corners[4], corners[6]],
             [corners[3], corners[6], corners[7], corners[4]],
         ]    
     
     hexahedrons         = mean_mesh.faces
     volume_connectivity = mean_mesh.faces # todo: what are both of these, and where do we get them? 
     tetrahedra = []
     for cube in hexahedrons:  # each cube is a list of 8 point indices
         tetrahedra.extend(subdivide_cube_to_tetrahedra(cube))

     target_set     = set(cp_ids)
     
     tetrahedrons   = [tetra for tetra in volume_connectivity if all(pt in target_set for pt in tetra)]
     tris           = np.array(tetrahedrons)   


# %% Volumetric Function
# n = 50 # this defines the density of the grid and should be carried over from when particles are warped -- don't know where the best location for this is or how to label it
# reader  = vtk.vtkPolyDataReader()
# reader.SetFileName(mean_model_name)
# reader.Update()
# mesh    = reader.GetOutput()

# x = np.linspace(mesh.GetBounds()[0], mesh.GetBounds()[1], n)
# y = np.linspace(mesh.GetBounds()[2], mesh.GetBounds()[3], n)
# z = np.linspace(mesh.GetBounds()[4], mesh.GetBounds()[5], n)
# grid = np.array(np.meshgrid(x, y, z)).reshape(3, -1).T

# implicit_function = vtk.vtkImplicitPolyDataDistance()
# implicit_function.SetInput(mesh)

# points_inside = []
# points_index = []
# for k,pt in enumerate(grid):
#     if implicit_function.FunctionValue(pt) < 0:  # Negative means inside
#         points_inside.append(pt)
#         points_index.append(k)

# point_cloud = pv.PolyData(points_inside)

# nx, ny, nz = n,n,n  # same as your linspace resolution
# grid_3d = grid.reshape((nx, ny, nz, 3))  # shape: (nx, ny, nz, 3)
# points_flat = grid_3d.reshape(-1, 3)
# def corner_indices(i, j, k, nx, ny):
#     base = lambda a, b, c: a + b * nx + c * nx * ny
#     return [
#         base(i, j, k),
#         base(i+1, j, k),
#         base(i+1, j+1, k),
#         base(i, j+1, k),
#         base(i, j, k+1),
#         base(i+1, j, k+1),
#         base(i+1, j+1, k+1),
#         base(i, j+1, k+1)
#     ]

# hexahedrons = []
# for i in range(nx - 1):
#     for j in range(ny - 1):
#         for k in range(nz - 1):
#             corners = corner_indices(i, j, k, nx, ny)
#             hexahedrons.append(corners)

# hexa_np = np.array(hexahedrons, dtype=np.int64)
# cell_sizes = np.full((hexa_np.shape[0], 1), 8)
# cell_array = np.hstack((cell_sizes, hexa_np)).flatten()
# cell_types = np.full(hexa_np.shape[0], pv.CellType.HEXAHEDRON, dtype=np.uint8)

# def subdivide_cube_to_tetrahedra(corners):
#     return [
#         [corners[0], corners[1], corners[3], corners[4]],
#         [corners[1], corners[2], corners[3], corners[6]],
#         [corners[1], corners[4], corners[5], corners[6]],
#         [corners[3], corners[4], corners[6], corners[7]],
#         [corners[1], corners[3], corners[4], corners[6]],
#         [corners[3], corners[6], corners[7], corners[4]],
#     ]


#%% Permutations
# Establish initial significance before clustering based on connectivity
print('Permutations...')
start = time.time()
cluster_samples, cluster_time_points, cluster_subjects = all_data.shape
perms           = np.zeros((num_iterations, cluster_subjects), dtype=int)
F_perms         = np.zeros((len(particles), timepoints, num_iterations))

# degree of freedom
if test_type == 0:
    df = len(np.where(group_ids == 0)[0]) + len(np.where(group_ids == 1)[0]) - 2
elif test_type == 1:
    df = len(np.where(group_ids == 0)[0]) - 1

it_count = 0


for it in range(num_iterations):
    if it == 0:
        perms[it, :] = range(cluster_subjects)
    else:
        perms[it, :] = np.random.permutation(cluster_subjects)
        # https://link.springer.com/article/10.1007/s11222-012-9370-4
        # Either form of permuting for paired tests are asymptotically valid
        # for paired testing. 
        # 1) swapping condition labels, not the subjects labels (maintaining dependency structures)
        # 2) permutations of all observations (neglecting dependency structure)
        
    # Copy original data
    permuted = all_data.copy()
    
    ### Circular Shift ###
    # random time shift 
    # PMCID: PMC6625832
    # https://stats.stackexchange.com/questions/599147/can-i-use-a-permutation-test-on-timeseries-data
    # maintaining temporal dependency structure
    if shift_opt == True:
        if it > 0:
            for s in range(cluster_subjects):
                shift = np.random.randint(cluster_time_points)
                permuted[:,:,s] = np.roll(permuted[np.random.permutation(cluster_samples),:,s],shift,axis=1)
                
    A = np.copy(np.transpose(permuted[:, :, perms[it, np.where(group_ids == 0)[0]]].reshape(num_pts*timepoints, len(np.where(group_ids == 0)[0]))))
    B = np.copy(np.transpose(permuted[:, :, perms[it, np.where(group_ids == 1)[0]]].reshape(num_pts*timepoints, len(np.where(group_ids == 1)[0]))))
    
    if test_type == 0:
        snpm_perms = spm1d.stats.nonparam.ttest2(A, B)
        
    if test_type == 1: 
        snpm_perms = spm1d.stats.ttest_paired(A, B) 
    
    F_perms[:, :, it]       = snpm_perms.z.reshape(len(particles), timepoints)
        
end = time.time()
time_mins = (end-start)/60
print(f"Elapsed time: {time_mins:.4f} minutes")    
print('Permutations: done')


#%% Clustering
print('Clustering...')
start = time.time()
# initialize variables
fclust_size = np.zeros((1, 1))
lc_it       = np.zeros((num_iterations))
# (assume one-sided)
if test_type == 0:
    ns = len(np.where(group_ids == 0)[0]) + len(np.where(group_ids == 1)[0])
    df = len(np.where(group_ids == 0)[0]) + len(np.where(group_ids == 1)[0]) - 2
elif test_type == 1:
    ns = len(np.where(group_ids == 0)[0])
    df = len(np.where(group_ids == 0)[0]) - 1
    
# vvv Cluster F_critical vvv    
# Fcrit       = scipy.stats.t.ppf(1 - p_value, df) # one-tailed
Fcrit       = scipy.stats.t.ppf(1 - alpha_val/2, df) # two-tailed

# iterate
for it in range(num_iterations):
    cnum = 1
    all_nodes = np.zeros([len(particles), timepoints])
    comp_size = []

    components = [None] * timepoints
    # for each timepoint
    for t in range(timepoints):
        # find particles with signficant test statistic
        nodes = np.where(abs(F_perms[:, t, it]) > Fcrit)[0]  
        nn = list()
        # for each particle find respective triangle
        for nd in nodes:
            nodeinx = np.where(tris == nd)[0] # Volumetric Connectivity tris --> tetrahedrons
            # if nodeinx[0].size != 0:
            #     nodeinx = np.array(nodeinx)
            #     nn.append(nodeinx.tolist())
            nn.extend(nodeinx)  # find tris containing significant nodes
        nodetris = np.unique(nn)
        # create graph (each edge of mesh)
        if nodetris.size > 0:
            edges = []
            for tri in tris[nodetris, :]: # Volumetric Connectivity tris --> tetrahedrons
                edges.extend(list(itertools.combinations(tri, 2)))

            # compute connected components and add to components list
            components[t] = list(trimesh.graph.connected_components(edges, min_len=1,nodes=nodes))  # less than three will not have a triangle
            for component in components[t]:
                all_nodes[component, t] = cnum
                comp_size.append(len(component))
                cnum = cnum + 1

    # find connected component across time by particle
    nzind       = np.nonzero(all_nodes)
    checkcomp   = []
    duplicates  = list(set([i for i in list(nzind[0]) if list(nzind[0]).count(i) > 1]))  # particles that are duplicated in non-zero
    ccdupind    = []
    for dp in duplicates:
        nxind       = np.nonzero(all_nodes[dp, :].ravel())[0]
        sequences   = np.split(nxind, np.array(np.where(np.diff(nxind) > 1)[0]) + 1)
        for seq in sequences:
            if len(seq) > 1:  # find components that span more than one timepoint and add to list
                dataadded = 0
                if len(checkcomp) == 0:
                    checkcomp.append(np.array(all_nodes[dp, seq]))
                else:
                    for check in range(len(checkcomp)):
                        if len(set.intersection(set(checkcomp[check]), set(all_nodes[dp, seq]))) > 0:
                            checkcomp[check] = np.array(np.unique(np.append(checkcomp[check], all_nodes[dp, seq]).ravel()))
                            if dataadded == 1:
                                # remove duplicate entry
                                ccdupind.append(check)
                            dataadded = 1
                    if dataadded == 0:
                        checkcomp.append(np.array(all_nodes[dp, seq]))

    # remove duplicate components by index
    for ccdi in reversed(np.unique(np.array(ccdupind))):
        del checkcomp[ccdi]

    # final check for components that overlap
    removecc = []
    for ind0 in range(len(checkcomp)):
        for ind1 in range(ind0 + 1, len(checkcomp)):
            if len(set.intersection(set(checkcomp[ind0]), set(checkcomp[ind1]))) > 0:
                checkcomp[ind0] = np.unique(np.append(checkcomp[ind0], checkcomp[ind1]))
                removecc        = np.append(removecc, ind1)
    for rcc in reversed(np.unique(removecc)):
        del checkcomp[rcc.astype(int)]

    # for all larger components - add the size and sum of all components
    comp_size       = np.array(comp_size)
    connectedcomps  = {}
    evaluatedcomps  = []
    sizecount       = 0
    for cc in checkcomp:
        sizecount       = sum(comp_size[cc.astype(int) - 1]) + sizecount
        evaluatedcomps  = np.append(evaluatedcomps, cc.astype(int) - 1)
        lc_it[it]       = max(lc_it[it], sum(comp_size[cc.astype(int) - 1]))
    for x in range(0, len(comp_size)):
        if x not in evaluatedcomps:
            sizecount = comp_size[x] + sizecount
            lc_it[it] = max(lc_it[it], comp_size[x])
    if sizecount != sum(comp_size):
        print('  Total summed size', sizecount, 'compared to', sum(comp_size), len(np.nonzero(all_nodes)[0]), '!!!!!!!!!')
        print(evaluatedcomps)
        print(checkcomp, duplicates)

# fclust_size     = int(np.ceil(np.quantile(lc_it, 1 - p_value, method='linear')))
fclust_size     = int(np.ceil(np.percentile(lc_it,(1-alpha_val)*100, method='linear')))

end = time.time()
time_mins = (end-start)/60
print(f"Elapsed time: {time_mins:.4f} minutes")
print('Clustering: done')
print(f"{(1-alpha_val)*100}th percentile cluster size: {fclust_size}")


#%% Network Analysis
print('Network Analysis...') 
# start network-based analysis
fvalues_size    = np.ones([len(particles), timepoints])
components      = [None] * timepoints
p_values_network    = 2 * (1 - scipy.stats.t.cdf(abs(fvalues), df)) # two-tailed p-value

cnum = 1
all_nodes = np.zeros([len(particles), timepoints])
comp_size = []

# for each timepoint
for t in range(timepoints):
    nodes   = np.where(abs(p_values_network[:, t]) < alpha_val)[0]  # find particles with signficant test statistic
    # nodes   = np.where(abs(fvalues[:, t]) > Fcrit)[0]
    nn      = list()
    # for each particle find respective triangle
    for nd in nodes:
        nodeinx = np.where(tris == nd)[0] # Volumetric Connectivity tris --> tetrahedrons
        nn.extend(nodeinx)  # tris containing significant particles
    nodetris = np.unique(nn)
    # create graph (each edge of mesh)
    if nodetris.size > 0:
        edges = []
        for tri in tris[nodetris, :]: # Volumetric Connectivity tris --> tetrahedrons
            edges.extend(list(itertools.combinations(tri, 2)))

        # compute connected components and add to components list
        components[t] = list(trimesh.graph.connected_components(edges, min_len=1, nodes=nodes))  # less than three will not have a triangle
        for component in components[t]:
            all_nodes[component, t] = cnum
            comp_size.append(len(component))
            cnum = cnum + 1

# find connected component across time by particle
nzind       = np.nonzero(all_nodes)
uniquepairs = dict()
checkcomp   = []
ccdupind    = []
duplicates  = list(set([i for i in list(nzind[0]) if list(nzind[0]).count(i) > 1]))  # particles that are duplicated in non-zero

for dp in duplicates:
    nxind = np.nonzero(all_nodes[dp, :].ravel())[0]
    sequences = np.split(nxind, np.array(np.where(np.diff(nxind) > 1)[0]) + 1)
    for seq in sequences:
        if len(seq) > 1:
            dataadded = 0
            # add array of components to a list for reviewing (more than two adjacent times)
            if len(checkcomp) == 0:
                checkcomp.append(np.array(all_nodes[dp, seq]))
            else:
                for check in range(len(checkcomp)):
                    if len(set.intersection(set(checkcomp[check]), set(all_nodes[dp, seq]))) > 0:
                        checkcomp[check] = np.array(np.unique(np.append(checkcomp[check], all_nodes[dp, seq]).ravel()))
                        if dataadded == 1:
                            # remove duplicate entry
                            ccdupind.append(check)
                        dataadded = 1
                if dataadded == 0:
                    checkcomp.append(np.array(all_nodes[dp, seq]))

# remove duplicate components by index
for ccdi in reversed(np.unique(np.array(ccdupind))):
    del checkcomp[ccdi]

# final check for components that overlap
removecc = []
for ind0 in range(len(checkcomp)):
    for ind1 in range(ind0 + 1, len(checkcomp)):
        if len(set.intersection(set(checkcomp[ind0]), set(checkcomp[ind1]))) > 0:
            checkcomp[ind0] = np.unique(np.append(checkcomp[ind0], checkcomp[ind1]))
            removecc        = np.append(removecc, ind1)
for rcc in reversed(np.unique(removecc)):
    del checkcomp[rcc.astype(int)]

# for all larger components - add the size and sum of all components
comp_size       = np.array(comp_size)
sig_nodes       = np.zeros([len(particles), timepoints])
connectedcomps  = {}
evaluatedcomps  = []
cccount         = 0
sizecount       = 0
for cc in checkcomp:
    connectedcomps[cccount]                 = {}
    connectedcomps[cccount]['components']   = cc

    ccpts       = []
    starttime   = timepoints
    endtime     = 0
    for component in connectedcomps[cccount]['components']:
        [pts, ti]   = np.where(all_nodes == component)
        starttime   = min(starttime, min(ti))
        endtime     = max(endtime, max(ti))
        ccpts       = np.unique(np.append(pts, ccpts))
    connectedcomps[cccount]['particles']    = ccpts
    connectedcomps[cccount]['size']         = sum(comp_size[cc.astype(int) - 1])
    connectedcomps[cccount]['starttime']    = starttime
    connectedcomps[cccount]['endtime']      = endtime

    # accumulators for checking
    evaluatedcomps  = np.append(evaluatedcomps, connectedcomps[cccount]['components'])
    sizecount       = sizecount + connectedcomps[cccount]['size']
    cccount         = cccount + 1

for x in range(1, len(comp_size) + 1):
    if x not in evaluatedcomps:
        sizecount                               = comp_size[x - 1] + sizecount
        connectedcomps[cccount]                 = {}
        connectedcomps[cccount]['components']   = x
        [pts, ti]                               = np.where(all_nodes == x)
        connectedcomps[cccount]['particles']    = pts
        connectedcomps[cccount]['size']         = comp_size[x - 1]
        connectedcomps[cccount]['starttime']    = min(ti)
        connectedcomps[cccount]['endtime']      = max(ti)
        cccount = cccount + 1

# if sizecount != sum(comp_size):
    # print(evaluatedcomps, cccount, duplicates)
    # print('  Total summed size', sizecount, 'compared to', sum(comp_size),
    #       len(np.nonzero(all_nodes)[0]), '!!!!!!!!!')

for component in connectedcomps.keys():
    complen = connectedcomps[component]['size']
    if complen >= int(fclust_size):
        # print('    Cluster of', complen, 'nodes from time', connectedcomps[component]['starttime'],
        #       'to', connectedcomps[component]['endtime'], ', pval =',
        #       1 - scipy.stats.percentileofscore(lc_it, complen) / 100)
        try:
            for comp in connectedcomps[component]['components']:
                ccloc = np.where(all_nodes == comp)
                if int(complen) >= int(fclust_size):
                    fvalues_size[ccloc[0], ccloc[1]] = p_values_network[ccloc[0], ccloc[1]] #fvalues[ccloc[0], ccloc[1]]#
        except:
            comp = connectedcomps[component]['components']
            ccloc = np.where(all_nodes == comp)
            # print('  ', comp, len(ccloc[0]))
            if int(complen) >= int(fclust_size):
                fvalues_size[ccloc[0], ccloc[1]] = p_values_network[ccloc[0], ccloc[1]] #fvalues[ccloc[0], ccloc[1]]#
    else:
        print('      Insignificant cluster of', complen, 'nodes from time',
              connectedcomps[component]['starttime'], 'to', connectedcomps[component]['endtime'],', pval =', 1 - scipy.stats.percentileofscore(lc_it, complen) / 100)
        
print('Network Analysis: done')

p_values_cluster    = fvalues_size
temp                = p_values_cluster[p_values_cluster < 1]
print(f"Number of significant nodes: {len(temp)}")


# %% Grab Timestamp
os.makedirs(os.path.join(proj_dir,'Network_Analysis_Results'), exist_ok=True)
if test_type == 0:
    test_add = 'tTest'
elif test_type == 1:
    test_add = 'paired'
    
timestamp       = time.strftime("%Y%m%d_%H%M%S")
out_dir_name    = f"Analysis_{test_add}_{timestamp}"

os.makedirs(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name), exist_ok=True)
out_dir = os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name)
os.makedirs(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Stats'), exist_ok=True)
os.makedirs(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Maps'), exist_ok=True)

# %% Save Results
out_arr = np.concatenate((cp_ids.reshape(-1,1),p_values_cluster), axis=1)

out_feat_name = surf_name
for n in group_names:
    if combine_groups:
        out_feat_name += '_' + n
        out_group_name = out_feat_name
    else:
        if n == group_names[0]:
            out_feat_name += '_' + n
        else:
            out_feat_name += '_vs_' + n
out_feat_name += '-' + feature_selection[0]
if len(feature_selection) > 1:
    out_feat_name = out_feat_name + '_vs_' + feature_selection[1]
    
# save to PolyData
out_points = []
out_points = pv.PolyData(mean_point[cp_ids,:])
for t in range(timepoints):
    out_points[f'Network_Analysis_{t}']  = p_values_cluster[:,t]
out_points.save(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Stats','Stats_' + out_feat_name + '.vtk'), binary=False)

# Output Stats
out_arr = np.concatenate((cp_ids.reshape(-1,1),p_values_cluster), axis=1)
df = pd.DataFrame(out_arr)
df.to_excel(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Stats','Stats_' + out_feat_name + '.xlsx'),index = False, header=False)

# Output Maps Results as .vtk
for g in range(2):
    if combine_groups:
        for f in range(len(feature_selection)):
            out_arr = np.mean(all_data[:, :, np.where(group_ids == f)[0]], axis=2)
            out_points = []
            out_points = pv.PolyData(mean_point[cp_ids,:])
            for t in range(timepoints):
                out_points[f'{feature_selection[f]}_{t}']  = out_arr[:,t]
            out_points.save(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Maps','Maps_' + out_group_name + '-' + feature_selection[f] + '.vtk'), binary=False)             
            
            out_arr = np.concatenate((cp_ids.reshape(-1,1),out_arr), axis=1)
            df = pd.DataFrame(out_arr)        
            df.to_excel(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Maps','Maps_' + out_group_name + '-' + feature_selection[f] + '.xlsx'),index = False, header=False)
    else:  
        out_arr = np.mean(all_data[:, :, np.where(group_ids == g)[0]], axis=2)
        out_points = []
        out_points = pv.PolyData(mean_point[cp_ids,:])
        for t in range(timepoints):
            out_points[f'{feature_selection[0]}_{t}']  = out_arr[:,t]
        out_points.save(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Maps','Maps_' + surf_name + '_' + group_names[g] + '-' + feature_selection[0] + '.vtk'), binary=False)        
        
        out_arr = np.concatenate((cp_ids.reshape(-1,1),out_arr), axis=1)
        df = pd.DataFrame(out_arr)         
        df.to_excel(os.path.join(proj_dir,'Network_Analysis_Results',out_dir_name,'Maps','Maps_' + surf_name + '_' + group_names[g] + '-' + feature_selection[0] + '.xlsx'),index = False, header=False)

print(f'Results Output: {out_dir}')
print('Complete!')

