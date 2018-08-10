import math, random, tarfile
import h5py
from tables import *
import numpy as np



PathIN = "/gpfs/cfel/cmi/labs/cryo-2/simulation_BGC/2018-08-03/"
PathOUT = "/gpfs/cfel/cmi/labs/cryo-2/simulation_BGC/2018-08-03/"
particle_size = 500
flowfields = [5,10,20,25,30,35,40,45,50,60,70]

for a in range(len(flowfields)):
    print('Sorting through '+str(flowfields[a])+' SCCM data.')
    files = np.genfromtxt(PathIN+'/'+str(particle_size)+'nm_'+str(flowfields[a]).zfill(2)+'sccm_maxwell/FileList.txt',dtype=str)
    output_file = str(flowfields[a]).zfill(2)+'sccm_'+str(particle_size)+'nm.h5'
    output = []

    for i in range(len(files)):
        output_temp = []
        print("Sorting file ",i," of ",len(files))
        FileIN = h5py.File(PathIN+'/'+str(particle_size)+'nm_'+str(flowfields[a]).zfill(2)+'sccm_maxwell/'+files[i],'r')
        initial = np.transpose(np.array(FileIN['/source/source']))
        final = np.transpose(np.array(FileIN['/detectors/-0.052']))
        initial_ids = np.ndarray.tolist(initial[:,-1])
        final_ids = np.ndarray.tolist(final[:,-1])
        for j in range(len(initial_ids)):
            try:
                matching_id=final_ids.index(initial_ids[j])
            except ValueError:
                pass
            else:
                output.append(np.concatenate((initial[j,0:6], final[matching_id,0:7]), axis=0))
        #output.append(output_temp)
        FileIN.close()
    #print(output)
    output = np.reshape(output,(-1,13))
    #print(output)
    #print(final_ids)
    #print(final[:,-1])

    FileOUT = h5py.File(PathOUT+'/'+str(particle_size)+'nm_'+str(flowfields[a]).zfill(2)+'sccm_maxwell/'+output_file, 'w')
    FileOUT.create_dataset('/particles', data=output)
    FileOUT.close()
