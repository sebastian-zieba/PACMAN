import numpy as np
import os
from astroquery.mast import Observations

t_min_obs =  56364.52973075 - 0.0001
t_max_obs =  56378.980103359994 + 0.0001

proposal_obs = Observations.query_criteria(proposal_id=13021,  instrument_name='WFC3/IR', project='HST')
print("Number of observations:",len(proposal_obs))

select = (t_min_obs <= proposal_obs['t_min'].value.data) & (proposal_obs['t_min'].value.data <= t_max_obs)

proposal_obs_select = proposal_obs[select]

data_products = Observations.get_product_list(proposal_obs_select)
print("Number of results:",len(data_products))

data_products_ima = data_products[data_products['productSubGroupDescription'] == 'IMA']
print("Number of results with ima extension:",len(data_products_ima))

Observations.download_products(data_products_ima,mrp_only=False)


file_dir = os.path.dirname(os.path.realpath(__file__))

root_dir = file_dir + '/mastDownload/HST' # Specify root directory to be searched for .sav files.
move_dir = file_dir
filelist = []
for tree,fol,fils in os.walk(root_dir):
    filelist.extend([os.path.join(tree,fil) for fil in fils if fil.endswith('.fits')])
for fil in filelist:
    name = fil.split('/')[-1]
    os.rename(fil,move_dir + '/' + name)

# delete the mastDownload directory
os.system("rm -r {0}".format(file_dir + '/mastDownload'))
