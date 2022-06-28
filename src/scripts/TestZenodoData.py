# Trying to get Zenodo data to work
import h5py as h5
import numpy as np

save_loc    =  '/Users/lieke/surfdrive/Documents'+'/SFRD_fit/src/tex/figures/' #/n/home04/lvanson/

data_dir = 'src/data/vanSon21/COMPAS_Output_wRates/'
data = h5.File(data_dir + 'COMPAS_Output_wWeights.h5', 'r')

print(data.keys())

data.close

fix,ax = plt.subplots()
ax.plot(np.linspace(0,1,10),np.linspace(0,1,10))
plt.savefig(save_loc + 'testdata.jpg')
