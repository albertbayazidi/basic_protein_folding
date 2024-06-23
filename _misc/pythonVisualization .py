#%% 
from utils import *
from plot import *  
label = True
# %%  INITAL PROTEIN
path = 'data/monte_carlo_sim/task_2.5/inital_and_last/inital.txt'
size, kords, monomers = convert_array_n_3(path)

path = 'data/monte_carlo_sim/task_2.5/inital_and_last/NN_inital.txt'
_, _, nearestMonomers = convertDataMatrix(path)

obj2 = findMonomerparKoords(nearestMonomers) #findes the nearestneighbour for each monomer  

if kords.shape[0] > 100:
    label = False

ax, fig = plotProtein(kords,obj2,labels = label)

#%% LAST MONTE CARLO ITERAION

path = 'data/monte_carlo_sim/task_2.5/inital_and_last/last.txt'
size, kords, monomers = convert_array_n_3(path)

path = 'data/monte_carlo_sim/task_2.5/inital_and_last/NN_last.txt'
_, _, nearestMonomers = convertDataMatrix(path)

obj2 = findMonomerparKoords(nearestMonomers) #findes the nearestneighbour for each monomer  

ax, fig = plotProtein(kords,obj2,labels = label)

#%% ZOOM in
if kords.shape[0] > 100:
    ax, fig = plotProtein(kords,obj2, sizeX = [155,160], sizeY = [148,150])
    ax, fig = plotProtein(kords,obj2, sizeX = [240,240], sizeY = [148,150])
else:
    print('To few monomers to zoom in')
# %% parse logger
path = 'data/monte_carlo_sim/task_2.5/main_sim/Protein2D.txt'
nr_mc_steps, nr_mc_sweeps, Koords_tensor, AminoAcid_matrix, NN_tensor, EndToEndDistanceArray, Energy_array = ReadLogger2D(path)

# %% make images
if nr_mc_steps*nr_mc_sweeps > 5000:
    print('To many images')
else:
    png_path = "data/png/image"
    make_images_Protein2d(nr_mc_steps, nr_mc_sweeps, Koords_tensor, NN_tensor, png_path,skip = size, sizeX = [50,65], sizeY = [75,93])   

#%% Animate

if nr_mc_steps*nr_mc_sweeps > 5000:
    print('To many images to create a gif')
else:
    gif_path = "data/monte_carlo_sim/task_2.5/gif/main_sim.gif"
    createGif(nr_mc_steps, nr_mc_sweeps, png_path, gif_path, gif_length = 500, skip = size)   

#%% End to end distance
plotEndToEndDistance(EndToEndDistanceArray)
## NEED TO IMPLEMENT SAVE IN FUNCTION

#%% Energy
plotEnergy(Energy_array, window_size=300)
## NEED TO IMPLEMENT SAVE IN FUNCTION
# %%
