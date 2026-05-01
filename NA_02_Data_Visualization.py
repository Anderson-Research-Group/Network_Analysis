# %%
'''
Network Analysis - 02 - Visualization of Results

This script creates figures of the Network Analysis results. Allows for some
user customization (colors, colormaps, glyph settings, etc.) but by no means
is a comprehensive package for every scenario. 

Needed Files:
    project directory
    surface (*.vtk or *.stl) and related correspondence model (nx3 particles file, *.txt or *.particles)
    results from NA_01 script (*.xlsx, must be in created architecture)

Created by:     Rich J. Lisonbee, MS
Organization:   University of Utah, Orthopaedic Research Laboratory
PI:             Andrew E. Anderson, PhD
Date:           1/15/2026

Modified by:
Date:

'''

import numpy as np
import pandas as pd
import matplotlib
import pyvista as pv
import tkinter as tk
import FreeSimpleGUI as sg
import re
import imageio.v2 as imageio
import io
import os, sys

os.system('cls||clear')
np.random.seed(0)

sg.set_options(use_ttk_buttons=False)


# %% Project Directory
sg.popup('Please select the Project Directory', keep_on_top=True)
proj_dir = tk.filedialog.askdirectory(title="Select a folder", initialdir=os.getcwd())


# %% Results Selection
result_dir_names = [name for name in os.listdir(os.path.join(proj_dir,'Network_Analysis_Results'))]

# remove figure settings from choices if it is there
try:
    del result_dir_names[result_dir_names.index('Figure_Settings')]
except ValueError:
    pass

# ***Result Selection UI***
layout  = [
    [sg.Text('Result Selection:')],
    [sg.Listbox(result_dir_names,
        key='result',
        select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
        size=(40, 10))],   
    [sg.Ok(), sg.Cancel()]
]

window = sg.Window("Result Selection", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()

result_dir = values['result'][0]
#  ***


# %% Select the Stats
result_names = []
result_names = [name for name in os.listdir(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Stats')) if name.split('.')[-1] == 'xlsx']


# ***Stats Result Selection UI***
layout  = [
    [sg.Text('Result Selection:')],
    [sg.Listbox(result_names,
        key='result',
        select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
        size=(80, 10))],   
    [sg.Ok(), sg.Cancel()]
]

window = sg.Window("Result Selection", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()

result_selection = values['result']
#  ***

data = pd.read_excel(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Stats',result_selection[0]), header=None).to_numpy()
cp_idx = data[:,0]


# %% Select the Maps
sg.popup('Select what feature map (data) you would like to visualize, if two are selected it will take the difference.', keep_on_top=True)

result_names = []
result_names = [name for name in os.listdir(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Maps')) if name.split('.')[-1] == 'xlsx']


# ***Stats Result Selection UI***
layout  = [
    [sg.Text('Maps Selection:')],
    [sg.Listbox(result_names,
        key='result',
        select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,    # multi-select
        size=(80, 10))],   
    [sg.Ok(), sg.Cancel()]
]

window = sg.Window("Feature Maps Selection", layout, modal=True, keep_on_top=True)
event, values = window.read()
window.close()

maps_selection = values['result']
#  ***

diff_maps = False
pick0 = 0
if len(maps_selection) > 1:
    diff_maps = True
    # ***Stats Result Selection UI***
    layout  = [
        [sg.Text('First Map (will be subtracted by the other):')],
        [sg.Listbox(maps_selection,
            key='result',
            select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # single-select
            size=(80, 10))],   
        [sg.Ok(), sg.Cancel()]
    ]

    window = sg.Window("Difference Feature Map Selection", layout, modal=True, keep_on_top=True)
    event, values = window.read()
    window.close()

    diff_selection = values['result']
    #  ***
    
    pick0 = maps_selection.index(diff_selection[0])
    pick1 = int(maps_selection.index(diff_selection[0]) - 1)
    map0_data = pd.read_excel(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Maps',maps_selection[pick0]), header=None).to_numpy()[:,1:]
    map1_data = pd.read_excel(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Maps',maps_selection[pick1]), header=None).to_numpy()[:,1:]
    
    map_data = map0_data - map1_data
    
    map_name = maps_selection[pick0].split('-')[-1].split('.')[0] + ' - ' + maps_selection[pick1].split('-')[-1].split('.')[0]
    print(f"Feature Maps Difference: {map_name}")
    
else:
    map_data = pd.read_excel(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Maps',maps_selection[pick0]), header=None).to_numpy()[:,1:]
    map_name = maps_selection[pick0].split('-')[-1].split('.')[0]
    print(f"Feature Map: {map_name}")


# %% Function Definitions
def unit(v):
    v = np.asarray(v, dtype=float)
    nrm = np.linalg.norm(v)
    return v / nrm if nrm != 0 else v

def align_axis_to_target(mesh, center, src_axis, dst_axis):
    src = unit(np.asarray(src_axis, dtype=float))
    dst = unit(np.asarray(dst_axis, dtype=float))

    v = np.cross(src, dst)
    s = np.linalg.norm(v)
    c = np.dot(src, dst)

    if s < 1e-12:
        # Already aligned (c>0) or opposite (c<0); handle 0° or 180°.
        angle_deg = 0.0 if c > 0 else 180.0
        # Any axis works for 180°—choose z for determinism
        axis = np.array([0, 0, 1])
    else:
        axis = v / s
        angle_deg = np.degrees(np.arctan2(s, c))

    mesh.rotate_vector(axis, angle_deg, point=center, inplace=True)
    return mesh

def Particle_PointCloud_Bead(model,point,radius=1):
    bead_merge = []
    for p in range(point.shape[0]):
        i = model.find_closest_point(point[p,:])
        n = model.point_data['Normals'][i]
        bead = pv.Sphere(center=point[p,:], direction=n, radius=radius)
        bead = align_axis_to_target(bead, center=point[p,:], src_axis=[0, 0, 1], dst_axis=n)
        bead.cell_data['cp_id'] = np.full(bead.n_cells, p, dtype=int)
        bead_merge.append(bead)
            
        bead_merged = pv.merge(bead_merge)
    return bead_merged

def Particle_PointCloud_Disc(model,point,point_sig,radius=1,rad_scalez=0.5):
    disc_merge  = []
    for p in range(point.shape[0]):        
        if np.isin(p,point_sig):
            i = model.find_closest_point(point[p,:])
            n = model.point_data['Normals'][i]
            disc = pv.ParametricSuperToroid(center=point[p,:], xradius=radius, yradius=radius, zradius=radius*rad_scalez)
            disc = align_axis_to_target(disc, center=point[p,:], src_axis=[0, 0, 1], dst_axis=n)
            disc.cell_data['cp_id'] = np.full(disc.n_cells, p, dtype=int)
            disc_merge.append(disc)
            
    disc_merged = pv.PolyData()
    if len(disc_merge) > 0:
        disc_merged = pv.merge(disc_merge)
    return disc_merged

def hex_to_rgb01(hex_color: str):
    # Ensure the '#' is removed
    if hex_color.startswith('#'):
        hex_color = hex_color[1:]
        
    # The input to fromhex must be a 6-digit hex string
    if len(hex_color) == 3:
        hex_color = "".join([c*2 for c in hex_color])

    return tuple(bytes.fromhex(hex_color))

def ColorMap_Bin(map_data,merged_mesh,cmap_choice='jet',clim=None):
    cmap = matplotlib.colormaps[cmap_choice]
    if clim == None:
        feat_min = np.array(np.floor(np.min(map_data)))
        feat_max = np.array(np.ceil(np.max(map_data)))
    else:
        feat_min = np.array(clim[0])
        feat_max = np.array(clim[1])
    
    map_data_norm = map_data / (feat_max - feat_min)
    
    cmap_list = []
    for time_point in range(map_data_norm.shape[1]):
        cmap_time = []
        for c in range(map_data_norm.shape[0]):
            # np.full(cmap(map_data_norm[c,time_point])[:3],merged_mesh.n_cells)
            cmap_time.append([float(x) for x in cmap(map_data_norm[c,time_point])[:3]])
        cmap_list.append(cmap_time)
    return cmap_list, [feat_min, feat_max]
    
def ColorMap_Assign(cp_idx,merged_mesh,color_map_list):
    updated_mesh = []
    for cp in range(len(cp_idx)):
        ext_mesh = merged_mesh.extract_values(cp)
        ext_mesh.cell_data['colors'] = np.ones((int(ext_mesh.n_cells),3)) * color_map_list[cp]
        
        updated_mesh.append(ext_mesh)
        
    mesh_colored = pv.merge(updated_mesh)
    return mesh_colored

def ColorPicker(surf_color='#c0c0c0',surf_trans=1,glyph_scale=1,disc_color='#000000',cmap_choice='jet',cbar_inc=True,cbar_dir=True,clim=(0,10),window_size=(900,900)):
    def make_cmap(cmap_choice):
        # makes the swatch for the colormap
        cmap    = matplotlib.colormaps[cmap_choice]
        width   = 120
        height  = 30
        x = np.linspace(0, 1, width)
        img = np.tile(x, (height, 1))
        fig, ax = matplotlib.pyplot.subplots(figsize=(width/100, height/100), dpi=100)
        ax.imshow(img, cmap=cmap, aspect='auto')
        ax.axis('off')
        buf = io.BytesIO()
        matplotlib.pyplot.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        matplotlib.pyplot.close(fig)
        buf.seek(0)
        return buf.getvalue()
    
    cmap_choices = matplotlib.colormaps()#['jet','RdBu','PuOr','managua']
    # ***Initial Figure UI***
    # colorbar inclusion logic
    cbar_inc_def = [True, False]
    if not cbar_inc:
        cbar_inc_def = [False, True]
    # colorbar direction logic
    cbar_dir_def = [True, False]
    if not cbar_dir:
        cbar_dir_def = [False, True]    
    
    layout = [
        [sg.Text('Color Choices', justification='center')],
        # SURFACE COLOR
        [sg.Text("", key="-SWATCH-SURF-COLOR-", size=(11, 1), background_color=surf_color),
         sg.Button("Pick Surface Color",    key="-PICK-SURF-COLOR-", use_ttk_buttons=False)],
        # SURFACE TRANSPARENCY
        [sg.Text('Surface Transparency'),
         sg.Input(default_text=surf_trans,    key='-SURF-TRANS-', size=(10, 1), justification='center')],     
        # RING COLOR - SIGNIFICANCE
        [sg.Text("", key="-SWATCH-RING-COLOR-", size=(11, 1), background_color=disc_color), 
         sg.Button("Pick Ring Color (if statistically significant)", key="-PICK-RING-COLOR-", use_ttk_buttons=False)],
        # GLYPH SCALE
        [sg.Text('Particle Scale'),
         sg.Input(default_text=glyph_scale, key='-GLYPH-SCALE-', size=(10, 1), justification='center')],        
        # COLORMAP
        [sg.Button("", key="-SWATCH-COLORMAP-", image_data=make_cmap(cmap_choice)),
         sg.Button("Pick Colormap",         key="-PICK-COLORMAP-", use_ttk_buttons=False)],
        # SHOW COLORBAR
        [sg.Radio('Include Colorbar',       key='-CBAR0-', group_id=1, default=cbar_inc_def[0]),
         sg.Radio('Hide Colorbar',          key='-CBAR1-', group_id=1, default=cbar_inc_def[1])],
        # COLORBAR DIRECTION
        [sg.Radio('Vertical Colorbar',      key='-CBAR-DIR0-', group_id=2, default=cbar_dir_def[0]),
         sg.Radio('Horizontal Colorbar',    key='-CBAR-DIR1-', group_id=2, default=cbar_dir_def[1])],    
        # COLORMAP LIMITS
        [sg.Text('Colormap Limits'),
         sg.Input(default_text=clim[0],     key='-CLIM1-', size=(10, 1), justification='center'),
         sg.Input(default_text=clim[1],     key='-CLIM2-', size=(10, 1), justification='center')],
        # FIGURE SETTINGS SECTION
        [sg.Text('Figure Settings', justification='center')],    
        # WINDOW
        [sg.Text('Window Size    '),
         sg.Input(default_text=window_size[0],key='-WIN1-', size=(10, 1), justification='center'),
         sg.Input(default_text=window_size[1],key='-WIN2-', size=(10, 1), justification='center')],
        # PROCEED/CANCEL
        [sg.Button("OK")],
    ]
    win = sg.Window("Figure Settings", layout, finalize=True, keep_on_top=True)
    
    k = 0
    while k == 0:
        event, values = win.read()
        if event == "-PICK-SURF-COLOR-":
            # returns ((r,g,b), '#RRGGBB') or (None, None) if cancelled
            rgb, hex_color = tk.colorchooser.askcolor(parent=win.TKroot)
            if hex_color:
                win["-SWATCH-SURF-COLOR-"].update(background_color=hex_color)
                win.refresh()
                surf_color = hex_color
        if event == "-PICK-RING-COLOR-":
            # returns ((r,g,b), '#RRGGBB') or (None, None) if cancelled
            rgb, hex1_color = tk.colorchooser.askcolor(parent=win.TKroot)
            if hex1_color:
                win["-SWATCH-RING-COLOR-"].update(background_color=hex1_color)
                win.refresh()
                disc_color = hex1_color
        if event == "-PICK-COLORMAP-":
            print('Colormap Choices:')
            # ***
            layout1 = [
                [sg.Text('Color Map Selection', justification='center')],
                [sg.Input(default_text='',key='-MAN-INPUT-')],
                [sg.Listbox(cmap_choices,
                    key='-CMAP-CHOICE-',
                    select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,    # multi-select
                    size=(80, 10))],                
                [sg.Button("OK")]
                ]
            win1 = sg.Window("Color Choices", layout1, finalize=True, keep_on_top=True)
            event1, values1 = win1.read()
            if event1 == "OK":
                cmap_choice1 = values1['-CMAP-CHOICE-']
                cmap_choice2 = values1['-MAN-INPUT-']
                
                if values1['-CMAP-CHOICE-']:
                    cmap_choice = cmap_choice1[0]
                if values1['-MAN-INPUT-']:
                    if cmap_choice2 in matplotlib.colormaps():
                        cmap_choice = cmap_choice2
                print(cmap_choice)
                win["-SWATCH-COLORMAP-"].update(image_data=make_cmap(cmap_choice))
                win1.close()
            # ***
        if event == "OK":
            surf_trans  = float(values['-SURF-TRANS-'])
            glyph_scale  = float(values['-GLYPH-SCALE-'])
            window_size = (int(values['-WIN1-']), int(values['-WIN2-']))
            clim        = (float(values['-CLIM1-']), float(values['-CLIM2-']))
            cbar_inc = False
            if values['-CBAR0-']:
                cbar_inc = True
            cbar_dir = False
            if values['-CBAR-DIR0-']:
                cbar_dir = True            
            k = 1
            break
        if event == sg.WIN_CLOSED:
            print('exiting application...')
            sys.exit()
    
    win.close()
        # ***
        # all of these outputs could have been a dictionary... oh well
    return(surf_color, surf_trans, glyph_scale, disc_color, cmap_choice, cbar_inc, cbar_dir, clim, window_size)


# %% Mean Surface and Points Selection
sg.popup('Please select the Surface file (results will be visualized on)', keep_on_top=True)
mean_model_name = tk.filedialog.askopenfilename(title="Select Surface", initialdir=proj_dir,filetypes=[("stl",'*.stl'),("vtk",'*.vtk')])

sg.popup('Please select the Particles file for:\n' f'{mean_model_name}', keep_on_top=True)
mean_point_name = tk.filedialog.askopenfilename(title="Select Particles", initialdir=mean_model_name.split('/')[:-1],filetypes=[("particles",'*.particles'),("txt",'*.txt')])


# %% Load Previous Camera Orientation
# window init
pos     = []
focal   = []
vup     = []

window_size = (900, 900)

load_cam = False

if os.path.isdir(os.path.join(proj_dir,'Network_Analysis_Results','Figure_Settings')):
    setting_names = []
    setting_names = [name for name in os.listdir(os.path.join(proj_dir,'Network_Analysis_Results','Figure_Settings')) if name.split('.')[-1] == 'txt']
    
    # ***Break Loop UI***
    layout  = [
        [sg.Text('Load previous camera settings?')],
        [sg.Checkbox('Check to Load:',
                     key='-CAM-')],        
        [sg.Ok()]
    ]
    window = sg.Window("Figure Settings", layout, modal=True)
    event, values = window.read()
    window.close()

    if event == sg.WIN_CLOSED:
        sys.exit()
    # ***
    
    load_cam = values['-CAM-']
    
    if load_cam: 
        # ***Stats Result Selection UI***
        layout  = [
            [sg.Text('Camera Settings:')],
            [sg.Listbox(setting_names,
                key='-CAM-CHOICE-',
                select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,
                size=(30, 5))],   
            [sg.Ok(), sg.Cancel()]
        ]
    
        window = sg.Window("Feature Maps Selection", layout, modal=True, keep_on_top=True)
        event, values = window.read()
        window.close()
    
        cam_selection = values['-CAM-CHOICE-'][0]
        #  ***
        
        temp = np.loadtxt(os.path.join(proj_dir,'Network_Analysis_Results','Figure_Settings',cam_selection),delimiter=',')
        pos     = temp[0,:]
        focal   = temp[1,:]
        vup     = temp[2,:]
        
        window_size = (int(temp[3,0]), int(temp[3,1]))


# %% Load Beads and Discs
model = pv.read(mean_model_name)
model.compute_normals(inplace=True)
point = np.loadtxt(mean_point_name)

# randomly generate significance
ri              = np.random.default_rng().choice(np.arange(0,data[:,0].shape[0]), size=(1, np.ceil(data[:,0].shape[0] * 0.10).astype(int)), replace=False)
point_sig_rand  = data[ri,0]

point_model     = point[data[:,0].astype(int),:]
bead_merged     = Particle_PointCloud_Bead(model,point_model)
disc_merged     = Particle_PointCloud_Disc(model,point,point_sig_rand)
# point_sig = data[np.where(data[:,1] < 1)[0],0]
# disc_merged     = Particle_PointCloud_Disc(model,point,point_sig)
# print(disc_merged)


# colorbar init
cmap_array, clim    = ColorMap_Bin(map_data,bead_merged,clim=None)
color_map_array     = cmap_array[0]
bead_colored        = ColorMap_Assign(cp_idx,bead_merged,color_map_array)


# %% Plot - Initial
sg.set_options(use_ttk_buttons=False)

surf_color, surf_trans, glyph_scale, disc_color, cmap_choice, cbar_inc, cbar_dir, clim, window_size = ColorPicker(clim=clim,window_size=window_size)
cmap_array, clim    = ColorMap_Bin(map_data,bead_merged,clim=clim,cmap_choice=cmap_choice)
color_map_array     = cmap_array[0]
bead_colored        = ColorMap_Assign(cp_idx,bead_merged,color_map_array)
prev_scale          = glyph_scale

sg.popup('Reminder that the particles shown are randomly generated for visualization purposes!', keep_on_top=True)

k = 0
adjust_settings = 1
append_text = ''
while adjust_settings == 1: 
    # grab previous values...
    prev_pos, prev_focal, prev_vup = pos, focal, vup
    
    if prev_scale - glyph_scale !=0:
        print('restitching mesh together...')
        disc_merged     = Particle_PointCloud_Disc(model,point,point_sig_rand,radius=glyph_scale)
        bead_merged     = Particle_PointCloud_Bead(model,point_model,radius=glyph_scale)
        bead_colored    = ColorMap_Assign(cp_idx,bead_merged,color_map_array)
        prev_scale = glyph_scale
        
    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(model, color=surf_color, show_edges=False, opacity=surf_trans)
    plotter.add_mesh(bead_colored, scalars='colors', rgb=True, show_edges=False, opacity=1)
    plotter.mapper.lookup_table = pv.LookupTable(cmap=cmap_choice)
    if cbar_inc:
        cbar = plotter.add_scalar_bar(vertical=cbar_dir, label_font_size=30, interactive=True,position_x=0.8, position_y=0.1)    #,fmt='%10.1f'          
    plotter.update_scalar_bar_range(clim=clim)
    if disc_merged is not None:
        plotter.add_mesh(disc_merged, color=disc_color, show_edges=False, opacity=1)
    headlight = pv.Light(light_type='headlight')
    headlight.intensity = 0.25
    plotter.add_light(headlight)
    
    if k == 0:
        k = 1
        if load_cam:
            plotter.camera_position = (pos, focal, vup)  
    elif k == 1:
        plotter.camera_position = (pos, focal, vup)   
    plotter.render()
    plotter.show() #auto_close=False
    
    # Grab the camera position, focal point, and view-up
    pos, focal, vup = plotter.camera_position
        
    # ***Break Loop UI***
    layout  = [
        [sg.Text('Change figure settings?')],
        [sg.Radio('Yes (Camera/Orientation)', group_id=1, default=True),
         sg.Radio('Yes (Figure Settings)', group_id=1),
         sg.Radio('No (Proceed)', group_id=1)],
        [sg.Text('Text to append to file name:')],
        [sg.Input(append_text,
            key='-APPEND-',
            size=(40,1))],
        [sg.Checkbox('Save camera orientation:',
                     key='-CAM-')],        
        [sg.Ok()]
    ]
    window = sg.Window("Continue changing setttings...", layout, modal=True)
    event, values = window.read()
    window.close()

    if event == sg.WIN_CLOSED:
        sys.exit()
        
    if values[1] == True:       
        surf_color, surf_trans, glyph_scale, disc_color, cmap_choice, cbar_inc, cbar_dir, clim, window_size = ColorPicker(surf_color=surf_color,surf_trans=surf_trans,glyph_scale=glyph_scale,disc_color=disc_color,cmap_choice=cmap_choice,cbar_inc=cbar_inc,cbar_dir=cbar_dir,clim=clim,window_size=window_size)
    
    if values[2]:
        adjust_settings = 0
    
    append_text = values['-APPEND-']
    save_set    = values['-CAM-']
    # ***
    
    cmap_array, clim    = ColorMap_Bin(map_data,bead_merged,clim=clim,cmap_choice=cmap_choice)
    color_map_array     = cmap_array[0]
    bead_colored        = ColorMap_Assign(cp_idx,bead_merged,color_map_array)
    
    
# %% Create Figures / Images
out_dir     = os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Figures',map_name)
out_cmap    = os.path.join(out_dir.replace(map_name,map_name+'_ColorMap'))
os.makedirs(out_dir, exist_ok=True)
os.makedirs(out_cmap, exist_ok=True)
    
for time_point in range(len(data[0,1:])):
    disc_merged = None
    point_sig   = data[np.where(data[:,time_point+1] < 1)[0],0]
    if len(point_sig) > 0:
        disc_merged = Particle_PointCloud_Disc(model,point,point_sig,radius=glyph_scale)
    
    color_map_array     = cmap_array[time_point]
    bead_colored        = ColorMap_Assign(cp_idx,bead_merged,color_map_array)
    
    plotter = pv.Plotter(off_screen=True,window_size=window_size)
    plotter.add_mesh(model, color=surf_color, show_edges=False, opacity=surf_trans,smooth_shading=False)
    plotter.add_mesh(bead_colored, scalars='colors', rgb=True, show_edges=False, opacity=1,smooth_shading=False)
    plotter.mapper.lookup_table = pv.LookupTable(cmap=cmap_choice)
    if cbar_inc:
        cbar = plotter.add_scalar_bar(vertical=cbar_dir, label_font_size=30, interactive=False,position_x=0.90, position_y=0.05)
        plotter.update_scalar_bar_range(clim=clim)
    if disc_merged is not None:
        plotter.add_mesh(disc_merged, color=disc_color, show_edges=False, opacity=1,smooth_shading=False)
    headlight = pv.Light(light_type='headlight')
    headlight.intensity = 0.25
    plotter.add_light(headlight)
    
    plotter.camera_position = (pos, focal, vup)

    plotter.disable()
    plotter.show()
    print(time_point)
    if append_text == '':
        plotter.screenshot(os.path.join(out_dir,f'fig_{time_point}.tiff'))
    else:
        plotter.screenshot(os.path.join(out_dir,f'fig_{time_point}_{append_text}.tiff'))
    plotter.close()

# save colorbar just in case
plotter = pv.Plotter(off_screen=True,window_size=window_size)
plotter.add_mesh(model, color=surf_color, show_edges=False, opacity=0,smooth_shading=False)
plotter.add_mesh(bead_colored, scalars='colors', rgb=True, show_edges=False, opacity=0,smooth_shading=False)
plotter.mapper.lookup_table = pv.LookupTable(cmap=cmap_choice)
plotter.add_scalar_bar(vertical=cbar_dir, label_font_size=30,fmt='%10.2f', interactive=False, position_x=0.5, position_y=0.25)
plotter.update_scalar_bar_range(clim=clim)
plotter.disable()
plotter.show()
plotter.screenshot(os.path.join(out_cmap,f'Color_Map_{cmap_choice}.tiff'))
plotter.close()


# %% Save Figure Settings
if save_set:
    fig_set_dir = os.path.join(proj_dir,'Network_Analysis_Results','Figure_Settings')
    os.makedirs(fig_set_dir, exist_ok=True)
    
    # ***Figure Settings UI***
    layout  = [
        [sg.Text('Name of Figure Settings File')],
        [sg.Input('cam_position_0',
            key='-FILE-NAME-',
            size=(40,1))],     
        [sg.Ok()]
    ]
    window = sg.Window("Figure Settings", layout, modal=True)
    event, values = window.read()
    window.close()

    if event == sg.WIN_CLOSED:
        sys.exit()
    fig_set_name = values['-FILE-NAME-']
    # ***
    
    np.savetxt(os.path.join(fig_set_dir,fig_set_name + '.txt'), np.array((pos, focal, vup, (window_size[0], window_size[1] ,0))), delimiter=',')    


# %% Video Creation (stitching together the .tiff files from last section)
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

# Create Video (only if more than one timepoint)
if len(data[0,1:]) > 1:
    if append_text == '':
        vid_output = os.path.join(out_dir,f'Video_{map_name}.mp4')
    else:
        vid_output = os.path.join(out_dir,f'Video_{map_name}_{append_text}.mp4')
    frames = [name for name in os.listdir(out_dir) if name.split('.')[-1] == 'tiff']
    
    frames = sorted_alphanumeric(frames)
    
    with imageio.get_writer(
        vid_output,
        fps=20,
        codec="libx264",
        quality=8,            # 0–10, higher = better quality (subject to codec)
        macro_block_size=None # allow odd sizes
    ) as writer:
        for fp in frames:
            img = imageio.imread(os.path.join(out_dir,fp))
            # Ensure 3-channel (some TIFFs may be grayscale or RGBA)
            if img.ndim == 2:             # gray → RGB
                img = imageio.core.util.Array.stack([img]*3, axis=-1)
            elif img.shape[-1] == 4:      # RGBA → RGB
                img = img[..., :3]
            writer.append_data(img)

print('Output Results Located:')
print(os.path.join(proj_dir,'Network_Analysis_Results',result_dir,'Figures',map_name))
print('')
print('complete!')

