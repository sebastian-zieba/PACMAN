from astropy.time import Time
import os
from astroquery.mast import Observations

#start and end time of interest
t_obs_utc = ['2013-03-13T12:34:53', '2013-03-15T08:34:55']
t_obs = Time(t_obs_utc, format='isot', scale='utc')

t_min_obs =  t_obs.mjd[0] - 0.0001
t_max_obs =  t_obs.mjd[1] + 0.0001

#seach from observations meeting our criteria
proposal_obs = Observations.query_criteria(proposal_id=13021, instrument_name='WFC3/IR', project='HST')
print("Number of observations:",len(proposal_obs))
print(proposal_obs)

#only choose the ones in our timewindow of interest
select = (t_min_obs <= proposal_obs['t_min'].value.data) & (proposal_obs['t_min'].value.data <= t_max_obs)
proposal_obs_select = proposal_obs[select]

#lists all files
data_products = Observations.get_product_list(proposal_obs_select)
print("Number of results:",len(data_products))

#only keep the IMA files
data_products_ima = data_products[data_products['productSubGroupDescription'] == 'IMA']
print("Number of results with ima extension:",len(data_products_ima))

#download the data
Observations.download_products(data_products_ima,mrp_only=False)

#the files were all saved in separate directories. lets move all IMA files into a common directory.
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

