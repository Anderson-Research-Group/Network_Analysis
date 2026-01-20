# Network_Analysis
Code for computing and visualizing Network Analysis results for spatiotemporal data. 
 
##### NA_01_Network_Analysis.py  
Performs spatiotemporal statistical analysis on mapped data at correspondence locations (particles) between different groups or between two different mapped data (sometimes referred to as feature maps).  
  
The results of this analysis will be saved in the project directory under /*proj_dir*/Network_Analysis_Results/Analysis_*date*_*time*/Stats and the combined map data (group mean at particle) under ../Maps
  
Output data includes .vtk surface files and .xlsx spreadsheets that can be used independently or with Paraview for visualization. If found to be not statistically significant at a point, it will be labeled as 1 in the stats files, otherwise it will have the adjusted p-value stored if it is significant. This makes querying the data easier.  
##### NA_02_Data_Visualization.py  
Creates figures (and video if temporal data) of the Network Analysis results and mapped feature data. Allows for some user customization of the figures.
   
## Resources    
[Network Analysis](https://link.springer.com/article/10.1007/s10439-023-03270-6)  
[ShapeWorks](https://github.com/SCIInstitute/ShapeWorks)  
[Python](https://www.python.org)  
[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  

## Installation
1. Clone the repository  
   >cd my/GitHubRepo/path  
   >git clone https://github.com/Anderson-Research-Group/Network_Analysis.git
     
2. Run installation bash script (Windows)
   Open Anaconda Prompt/Teriminal:
   >cd my/GitHubRepo/path/Network_Analysis  
   >install_network_analysis.bat
   
3. Ready to run the scripts in your preferred IDE or command line  
     *for Sypder users you will need to install spyder-kernels (example below)*
     >pip install spyder-kernels=3.0  
  
## File Organization    
Types of tests available:  
    paired t-tests      (either between two groups, or between two feature maps)  
    two-sample t-tests  (between groups)  
    *if you choose to combine groups (select multiple group folders) it assumes  
    that you are performing a paired test between two feature maps  
    *all tests assume two-tailed  
    *it is highly encouraged that when using circular shift it is appropriate for your application  

#### Needed Files  
    project directory  
       'Data' folder
           group folders and subsequent nested subject-specific folders  
           feature map files (.txt, .csv, or .xlsx)  
           1st column with correspondence particle id# and other columns with mapped data (each column represents a subsequent timepoint)  
    mean surface (.vtk or .stl) and related correspondence model (nx3 particles file, .txt or .particles  
  
Multiple analyses and different data can be ran from the same project file, hence the necessary architecture. The script is looking for (from bottom to top) feature maps in their respective surface folders, which are in their subject-specific folder, in their associated group folder, in one main directory named 'Data'. The architecture makes it quickly and readily available to import the data for each subject. So long as they all   
have surface and feature map folders with the same named schema. The individual files themselves (.stl, .csv, etc) can be named whatever you want. It is their location within the folders that is what the code will look for.  
  
Any folder with a name that needs to be specific is denoted with an asterisk (*)  
Folders are denoted with a forward slash (/) and files will have a file extension (.stl, .csv, etc)  

### How the project directory should look...  
Project Directory (can be located anywhere, named anything, and this is what you will select when you first run the script)  
&emsp;└─/Data* (./Data contains all necessary data, structure is required for importing properly)     
&emsp;&emsp;├─/Group_A &emsp;(group folders)  
&emsp;&emsp;├─/Group_B  
&emsp;&emsp;...  
&emsp;&emsp;└─/Controls [example]  
&emsp;&emsp;&emsp;├─/Subject_01    (subject-specific folders)  
&emsp;&emsp;&emsp;├─/Subject_02  
&emsp;&emsp;&emsp;...  
&emsp;&emsp;&emsp;└─Norm_12 [example]  
&emsp;&emsp;&emsp;&emsp;├─/Surface_01  
&emsp;&emsp;&emsp;&emsp;├─/Surface_02  
&emsp;&emsp;&emsp;&emsp;...  
&emsp;&emsp;&emsp;&emsp;└─/Pelvis [example]  
&emsp;&emsp;&emsp;&emsp;&emsp;├─/Mapped_Data_01  
&emsp;&emsp;&emsp;&emsp;&emsp;├─/Mapped_Data_02  
&emsp;&emsp;&emsp;&emsp;&emsp;...  
&emsp;&emsp;&emsp;&emsp;&emsp;└─/FE_Data [example]        (only has one file in it)  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;└─Norm_02_FE_Mean.xlsx [example]  
  
### Python Dependencies (Linux/macOS/Windows 
   numpy  
   pandas  
   trimesh  
   pyvista  
   open3d  
   spm1d  
   scipy  
   FreeSimpleGUI  
   matplotlib  
   imageio 
