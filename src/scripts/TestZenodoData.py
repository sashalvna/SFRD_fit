# Trying to get Zenodo data to work
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

save_loc    =  '/Users/lieke/surfdrive/Documents'+'/SFRD_fit/src/tex/figures/' #/n/home04/lvanson/

data_dir = 'src/data/vanSon21/'




# import tarfile

# tar = tarfile.open(data_dir + 'COMPAS_Output_wRates.tar.gz', "r:gz")

# for tarinfo in tar:
#     print(tarinfo.name, "is", tarinfo.size, "bytes in size and is")
#     if tarinfo.isreg():
#         print("a regular file.")
#     elif tarinfo.isdir():
#         print("a directory.")
#     else:
#         print("something else.")
# tar.close()


data = h5.File(data_dir + 'COMPAS_Output_wRates.tar.gz', 'r')#'COMPAS_Output_wWeights.h5', 'r')

print(data.keys())

print('!!! Je bent nu in deze file!!!!')

data.close

fix,ax = plt.subplots()
ax.plot(np.linspace(0,1,10),np.linspace(0,1,10))
plt.savefig(save_loc + 'testdata.jpg')

