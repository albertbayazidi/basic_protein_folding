import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from utils import *

font_step = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 25,
        }

font_sweep = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 18,
        }


def plotProtein(obj, obj2, sizeX = [-1,-1], sizeY = [-1,-1], labels = True):
    if sizeX[0] == -1:
        xMin = np.min(obj[:,0])
        xMax = np.max(obj[:,0])
        yMin = np.min(obj[:,1])
        yMax = np.max(obj[:,1])
    else:
        xMin, yMin = sizeX[0],sizeY[0]
        xMax, yMax = sizeX[1],sizeY[1]

    # plot settings
    fig = plt.figure(figsize=(8,8))
    ax  = fig.add_subplot(111)
    chartBox = ax.get_position() 
    ax.set_position([chartBox.x0, chartBox.y0 + chartBox.y0*1.2, chartBox.width, chartBox.height *0.8])

    # protein
    ax.plot(obj[0,0], obj[0,1], 'o', color = 'red' ,label = 'Start')
    ax.plot(obj[:2,0], obj[:2,1], color = 'black') 
    ax.plot(obj[1:,0], obj[1:,1], 'o-', color = 'black', label = 'Covalent Bonds') 
    ax.plot(obj[-1,0], obj[-1,1], 'o', color = 'green', label = 'End')

    # non-covalent bonds
    if len(obj2.shape) == 2:  
        for i in range(obj2.shape[1]):
            X = obj[obj2[:,i],0]
            Y = obj[obj2[:,i],1]
            ax.plot(X,Y, '--', color = 'red', label = f'non-Covalent nr{i+1}')

    # misc plot settings
    if labels: 
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=3)
    else:
        pass
    
    ax.set_xlim(-5 + xMin ,5 + xMax)
    ax.set_ylim(-5 + yMin ,5 + yMax)
    ax.grid()
    return ax, fig

# works for small proteins and frew X
def make_images_Protein2d(nr_mc_steps, nr_mc_sweeps, Koords_tensor, NN_tensor,
                     png_path, skip = 1, sizeX = [-1,-1], sizeY = [-1,-1]):

    tot = nr_mc_steps * nr_mc_sweeps
    
    curr_mc_step = 0
    curr_mc_sweep = 0

    for i in range(0,tot,skip):
        obj2 = findMonomerparKoords(NN_tensor[:,:,i])
        obj1 = Koords_tensor[:,:,i]
        ax, fig = plotProtein(obj1, obj2, sizeX, sizeY)
  
        ax.text(25, 32, f'MC step nr = {curr_mc_step}', fontdict=font_step)
        ax.text(25, 31.5, f'MC sweep nr = {curr_mc_sweep}', fontdict=font_sweep)
        
        fig.savefig(f'{png_path + str(i)}.png',pad_inches=0.1)
        plt.close(fig)
        
        curr_mc_sweep += skip
        if curr_mc_sweep == nr_mc_sweeps:
            curr_mc_sweep = 0
            curr_mc_step += 1
       
    print('Images created')

def createGif(nr_mc_steps, nr_mc_sweeps, png_path, gif_path, gif_length = 100, skip = 1,):
    tot = nr_mc_steps * nr_mc_sweeps
    # Create a gif
    images = [Image.open(f'{png_path + str(n)}.png') for n in range(0,tot,skip)]
    images[0].save(f'{gif_path}', save_all=True,
            append_images=images[1:], duration=gif_length, loop=0)

    print('Gif created')

def plotEndToEndDistance(EndToEndDistanceArray):
    plt.plot(EndToEndDistanceArray, label='End to End Distance')
    plt.xlabel('MC step')
    plt.ylabel('End to End Distance (Eucledian distance)')
    plt.title('End to End Distance of the protein over time')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True)
    plt.grid()
    plt.show()

def plotRog(Rog_array):
    plt.plot(Rog_array)
    plt.xlabel('MC step')
    plt.ylabel('Rog')
    plt.title('Rog of the protein over time')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True)
    plt.grid()
    plt.show()

def plotEnergy(Energy_array, window_size=10):
    # Calculate the running average using a moving window
    running_avg = np.convolve(Energy_array, np.ones(window_size)/window_size, mode='valid')

    # Plot the running average of the energy
    plt.plot(range(window_size, len(Energy_array) + 1), running_avg, label=f'Running Average (window size={window_size})')

    plt.xlabel('MC step')
    plt.ylabel('Energy')
    plt.title('Energy of the protein over time')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True)
    plt.grid()
    plt.show()
