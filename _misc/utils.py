import numpy as np
import matplotlib.pyplot as plt

def convert_array_n_3(file_name):
    """
    Convert a file with the following format:
    n
    x1 y1 m1
    x2 y2 m2
    ...
    xn yn mn
    into two numpy array with the following format:
    [[x1, y1], [[m1],
     [x2, y2],  [m2],
     ...
     [xn, yn]], [mn]]
    return the size of the array and the array itself
    """
    
    data = np.loadtxt(file_name ,delimiter=" ", dtype=str)
    mod_data = np.char.split(data[1:], sep=",")
    size = int(int(data[0])/2)

    matrix = np.zeros((size, 3), dtype=float)
    for i in range(int(size)):
        matrix[i,:] = np.array(mod_data[i],dtype=float)
        
    kords = matrix[:,:2]
    monomers = matrix[:,2]
    return size, np.array(kords,dtype=int), monomers


def convertDataMatrix(file_name):
    with open(file_name, 'r') as file:
        data = file.read()
        lines = data.strip().split('\n')
        size1 = int(lines[0])
        size2 = int(lines[1])

    matrix = np.zeros((size1, size2))
    
    for i in range(2, len(lines)):
        mod_lines = np.array(lines[i].split(',')[:-1])
        matrix[i-2] = mod_lines

    return size1, size2, matrix


def findMonomerparKoords(nearestMonomers):


    nearestMonomers = nearestMonomers.astype(int)
    firstNeigbour = nearestMonomers > 0
    indexStart = np.argwhere(firstNeigbour)[:,0] # might want 0 to be cahnged
    indexEnd = nearestMonomers[firstNeigbour]
    otherBoundsIndex = np.vstack((indexStart, indexEnd))
    
    return otherBoundsIndex


def ReadLogger2D(path):
    # parsing the data
    main_data = np.loadtxt(path ,delimiter=" ", dtype=str)
    mod_data = np.char.split(main_data[3:], sep=",")
    
    # size of the arrays
    nr_mc_steps = int(main_data[0])
    size_of_koords_array = int(int(main_data[2])/2)  # +1 since we know keep track of mcsweps     
    start_of_NN_matrix = size_of_koords_array + 2
    size1_of_NN_matrix =  int(mod_data[size_of_koords_array][0])
    size2_of_NN_matrix =  int(mod_data[size_of_koords_array+1][0])
    data_len = size_of_koords_array + size1_of_NN_matrix + 6

    # 6 is based on nr of monomers, nr of monommer neighbours, size of these arrays, the end to end distance. and the current mc_sweep will ned +1 for ROG 

    # makeing storage arrays
    Koords_matrix = np.zeros((size_of_koords_array, 2,nr_mc_steps*size_of_koords_array)) 
    NN_matrix = np.zeros((size1_of_NN_matrix, size2_of_NN_matrix,nr_mc_steps*size_of_koords_array))
    EndToEndDistance = np.zeros(nr_mc_steps*size_of_koords_array)
    Energy_array = np.zeros(nr_mc_steps*size_of_koords_array)
    AminoAcid_array = np.zeros((size_of_koords_array))

    for j in range(nr_mc_steps*size_of_koords_array): 
        # temporary variables
        time_shift = data_len * j 
        koords_line = mod_data[time_shift:size_of_koords_array+time_shift]
        NN_lines = mod_data[start_of_NN_matrix + time_shift :start_of_NN_matrix+size1_of_NN_matrix + time_shift]
        
        EndToEndDistance[j] = mod_data[data_len*(j+1)-4][0] # will be modified when including ROG
        Energy_array[j] = mod_data[data_len*(j+1)-3][0]# will be modified when including ROG

        for i in range(size_of_koords_array):
            mod_lines = np.array(NN_lines[i][:-1])
            Koords_matrix[i,:,j] = np.array(koords_line[i][:-1])
            AminoAcid_array[i] = koords_line[i][2]
            NN_matrix[i,:,j] = mod_lines

    
    return nr_mc_steps, size_of_koords_array, Koords_matrix, AminoAcid_array, NN_matrix, EndToEndDistance, Energy_array