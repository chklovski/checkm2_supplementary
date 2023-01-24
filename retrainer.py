
import numpy as np
import pandas as pd

import zenodo_backpack


import pickle
from sklearn.preprocessing import MinMaxScaler
import lightgbm as lgb
import scipy
from scipy import sparse
import os
import warnings

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow import keras
from tensorflow.keras import layers

from numpy import unique
from numpy import argmax
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.layers import Lambda
from tensorflow.keras.layers import Flatten
from tensorflow.keras import Model
from tensorflow.keras.layers import Input
from tensorflow.keras import regularizers
from tensorflow.keras.layers import Concatenate
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import Reshape
from tensorflow.keras import optimizers
from tensorflow.keras import utils


from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix, vstack, hstack
from scipy import sparse
import scipy
import subprocess
import sys
import os
from os import listdir
import pickle
import logging
import shutil

loglevel = logging.INFO
logging.basicConfig(level=loglevel, format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


nthreads = sys.argv[1]
new_GTDB_version = sys.argv[2]
DOI = sys.argv[3]

if not new_GTDB_version.startswith('r') or not len(new_GTDB_version) == 4:
    logging.info('New GTDB version must be in the format rXXX where XXX are numbers')
    sys.exit(1)
    
current_dir = os.getcwd()
good_genome_folder = f'{current_dir}/new_genomes/'



#let us acquire external data

backpack_downloader = zenodo_backpack.ZenodoBackpackDownloader()
backpack = backpack_downloader.download_and_extract(f'{current_dir}', f'{DOI}')
shutil.move(f'{current_dir}/old_vectors.zb/payload_directory', f'{current_dir}/old_vectors/')

#tensorflow can't handle sparse matrices without using batch generators
def batch_generator(X, y, batch_size, samples_per_epoch, CONVOLUTIONAL=False):
    number_of_batches = samples_per_epoch/batch_size
    counter=0
    shuffle_index = np.arange(np.shape(y)[0])
    np.random.shuffle(shuffle_index)
    X =  X[shuffle_index, :]
    y =  y[shuffle_index]
    #z = weights[shuffle_index]
    while 1:
        index_batch = shuffle_index[batch_size*counter:batch_size*(counter+1)]
        X_batch = X[index_batch,:].todense()
        y_batch = y[index_batch]
        #z_batch = z[index_batch]
        counter += 1
        if not CONVOLUTIONAL:
            yield(np.array(X_batch),y_batch)
        else:
            X_batch = np.array(X_batch)
            #yield(X_batch, y_batch)
            yield(X_batch.reshape(X_batch.shape[0], X_batch.shape[1], 1),y_batch)
        if (counter < number_of_batches):
            np.random.shuffle(shuffle_index)
            counter=0


def returnPredictors(dataset, mode):
    dataset_refined = dataset.copy()
    del dataset_refined['Accession']

    del dataset_refined['Domain']
    del dataset_refined['Phylum']
    del dataset_refined['Class']
    del dataset_refined['Order']
    del dataset_refined['Family']
    del dataset_refined['Genus']
    del dataset_refined['Species']
    contamination = dataset_refined['Contamination'].values
    del dataset_refined['Contamination']
    del dataset_refined['Name']

    completeness = np.array(dataset_refined['Completeness']) 

    predictors = np.array(dataset_refined.drop('Completeness', axis=1))
    if mode == 'Comp':
        return predictors, completeness
    elif mode == 'Cont':
        return predictors, contamination
    else:
        print('Please specify mode: Comp or Cont')
        return None




# IMPORTANT: your new good genomes must be protein files (.faa) ideally called by prodigal
logging.info('Making synthetic genomes.')
try:
    #this will create a series of fragged_* folders in good_genome_folder
    cmd = f"cd {current_dir}/{good_genome_folder}; mkdir output; readlink -f *.faa | while read -r genome; do bash ../faa_fragger.sh $genome; done"

    logging.debug(cmd)
    subprocess.call(cmd, shell=True)
except Exception as e:
    logging.error('An error occured while running {}: {}'.format(cmd, e))
    sys.exit(1)
    
    #now run checkm2 to dump input vectors for each 
logging.info('Running CheckM2 and dumping input vectors.')

try:    
    cmd = f"cd {current_dir}/{good_genome_folder}; for i in fragged_*; do checkm2 predict -i $i -o output/$i -t {nthreads} -x faa --genes --dbg_vectors --remove_intermediates; done"
    
    logging.debug(cmd)
    subprocess.call(cmd, shell=True)
except Exception as e:
    logging.error('An error occured while running {}: {}'.format(cmd, e))
    sys.exit(1)

    
#remove intermediate synthetic genomes
logging.info('Removing generated synthetic genomes')

try:    
    cmd = f"cd {current_dir}/{good_genome_folder}; for i in fragged_*; do rm -rf $i; done"
    
    logging.debug(cmd)
    subprocess.call(cmd, shell=True)
except Exception as e:
    logging.error('An error occured while running {}: {}'.format(cmd, e))
    sys.exit(1)


#now read in pickled vectors and concatenate them

#contamination model starts giving weird results if genome is extremely small due to very few annotated proteins, 
#so exclude AALength < 100,000 

cont_threshold = 100000

#lists to hold vectors
comp_scaled_vectors, comp_scaled_labels = [],[]
comp_raw_vectors, comp_raw_labels = [],[]
cont_raw_vectors, cont_raw_labels = [],[]

scalerfile = f'{current_dir}/scaler.sav'
scaler = pickle.load(open(scalerfile, 'rb'))


logging.info(f'Going through CheckM2 output and concatenating vectors. Excluding from contamination model genomes with AA length < {cont_threshold}.')
    
cm2_folders = listdir(f'{good_genome_folder}/output/')
comp_files = sorted([f for f in cm2_folders if f.startswith('fragged_')])
pickled_vectors = []
#we read in each vector; transform it for scaled use; 
for folder in cm2_folders:
    exclude_from_cont = False
    pickled_location = f'{good_genome_folder}/output/{folder}'
    pickles = listdir(pickled_location)
    pickles = sorted([f for f in pickles if f.endswith('.pkl')])
    interim_pickled_list = []
    for pickle in pickles:
        pickled_vector = pd.read_pickle(f'{good_genome_folder}/output/{folder}/{pickle}')
        pickled_vector['Completeness'] = pickled_vector['Name'].apply(lambda x: x.split('_pc_complete')[0].split('_')[-1]).astype(float)
        pickled_vector['Contamination'] = pickled_vector['Name'].apply(lambda x: x.split('pc_contaminated')[0].split('_pc_completeXXX')[-1])
        interim_pickled_list.append(pickled_vector)
    
    full_pickles = pd.concat(interim_pickled_list)
#    full_pickles['ActualName'] = full_pickles['Name'].apply(lambda x: x.split('pc_completeXXX')[0].split('_')[-1])
    if full_pickles[(full_pickles['Completeness'] == 100) & (full_pickles['Completeness'] == 100)]['AALength'].values[0] < cont_threshold:
        exclude_from_cont = True
    #make raw comp vectors
    
    del full_pickles['Name']
    comp_raw_vectors.append(csr_matrix(full_pickles.iloc[:, :-2].values))
    comp_raw_labels.append(np.array(full_pickles['Completeness'].values))

    #make raw cont vectors
    if not exclude_from_cont:
        cont_raw_vectors.append(csr_matrix(full_pickles.iloc[:, :-2].values))
        cont_raw_labels.append(np.array(full_pickles['Contamination'].values))

    
    #make scaled vectors
    comp_scaled_vectors.append(csr_matrix(scaler.transform(full_pickles.iloc[:, :-2].values))[:, :20021])
    comp_scaled_labels.append(np.array(full_pickles['Completeness'].values))

#now concatenate all vectors with previously existing vectors
logging.info(f'Concatenating new feature vectors.')

comp_scaled_vectors = vstack(comp_scaled_vectors)
comp_scaled_labels = np.concatenate(comp_scaled_labels)
comp_raw_vectors = vstack(comp_raw_vectors)
comp_raw_labels = np.concatenate(comp_raw_labels)
cont_raw_vectors = vstack(cont_raw_vectors)
cont_raw_labels = np.concatenate(cont_raw_labels)



#now concatenate all vectors with previously existing vectors


logging.info(f'Concatenating new feature vectors with old existing vectors in {current_dir}/old_vectors')


old_release_loc = f'{current_dir}/old_vectors/'

previous_release_version = listdir(old_release_loc)
previous_release_version = [f for f in previous_release_version if f.endswith('.npy')][0].split('_')[1]

logging.info(f'Previous release version was: {previous_release_version}')




old_comp_scaled_vectors = scipy.sparse.load_npz(f'{old_release_loc}/all_{previous_release_version}_scaled.npz')[:, :20021]
old_comp_scaled_labels = np.load(f'{old_release_loc}/all_{previous_release_version}_scaled_complabels.npy')
old_comp_raw_vectors = scipy.sparse.load_npz(f'{old_release_loc}/all_{previous_release_version}_raw_comps.npz')
old_comp_raw_labels = np.load(f'{old_release_loc}/all_{previous_release_version}_raw_compabels.npy')
old_cont_raw_vectors = scipy.sparse.load_npz(f'{old_release_loc}/all_{previous_release_version}_raw_conts.npz')
old_cont_raw_labels = np.load(f'{old_release_loc}/all_{previous_release_version}_raw_contlabels.npy')




comp_scaled_vectors = vstack([comp_scaled_vectors, old_comp_scaled_vectors])
comp_scaled_labels = np.concatenate([comp_scaled_labels, old_comp_scaled_labels])
comp_raw_vectors = vstack([comp_raw_vectors, old_comp_raw_vectors])
comp_raw_labels = np.concatenate([comp_raw_labels, old_comp_raw_labels])
cont_raw_vectors = vstack([cont_raw_vectors, old_cont_raw_vectors])
cont_raw_labels = np.concatenate([cont_raw_labels, old_cont_raw_labels])



compparams = {
    'boosting_type': 'gbdt',
    'objective': 'regression',
    'metric': {'mae'},
    'num_leaves': 11,
    'min_data_in_leaf': 150, 
    'learning_rate': 0.2,
    'feature_fraction': 0.5,
    'bagging_fraction': 0.5,
    'bagging_freq': 3,
    'verbose': 1, 
    'reg_sqrt': True,
    'num_threads': nthreads, 
    'deterministic': True,
    'min_child_weight': 180,
}


contparams = {
    'boosting_type': 'gbdt',
    'objective': 'regression',
    'metric': {'mae'},
    'num_leaves':211,
    'min_data_in_leaf': 150, 
     'learning_rate': 0.2,
     'feature_fraction': 0.9,
     'bagging_fraction': 0.8,
     'bagging_freq': 5,
     'verbose': 1, 
     'reg_sqrt': True,
     'num_threads': nthreads, 
 }


gb_comp_train = lgb.Dataset(comp_raw_vectors, label=comp_raw_labels)
gb_cont_train = lgb.Dataset(cont_raw_vectors, label=cont_raw_labels)

logging.info('Training General model for Completeness using lightgbm.')
new_gb_model_comp = lgb.train(compparams,
               gb_comp_train,
               num_boost_round=450, valid_sets=gb_comp_train)
               

new_gb_model_comp.save_model(f'{current_dir}/new_models/general_model_COMP.gbm')

logging.info('Training General model for Contamination using lightgbm.')

new_gb_model_cont = lgb.train(contparams,
               gb_cont_train,
               num_boost_round=450, valid_sets=gb_cont_train)

new_gb_model_cont.save_model(f'{current_dir}/new_models/model_CONT.gbm ')


model_checkpoint_folder = f'{current_dir}/model_checkpoint/'
logging.info(f'Training Specific model for Completeness using tensorflow. Checkpointing models for later assessment at {model_checkpoint_folder}.')


#set up NN model:

model = Sequential()
model.add(Conv1D(180, kernel_size=10, strides=10,  input_shape=(20021,1)))
model.add(BatchNormalization())
model.add(Conv1D(180, kernel_size=10, strides=10, activation='relu'))
model.add(BatchNormalization())
model.add(Conv1D(180, kernel_size=10, strides=10,  activation='relu'))
model.add(BatchNormalization())
model.add(Conv1D(100, kernel_size=10, strides=10, activation='relu'))
model.add(BatchNormalization())
model.add(Flatten())
model.add(Dense(100, activation="relu"))
model.add(Dense(1, activation='sigmoid'))
model.compile(optimizer='adam', loss='mse', weighted_metrics=['mean_absolute_error'], loss_weights=comp_scaled_labels * 500) 

#last layer is sigmoid so need labels to be between 0 and 1
train_labels = comp_scaled_labels / 100
epochs_num = 10
batch_size = 786 #786

my_callbacks = [
    keras.callbacks.EarlyStopping(patience=10),
    keras.callbacks.ModelCheckpoint(filepath=os.path.join(f'{model_checkpoint_folder}', 'model.{epoch:02d}.h5')),
    keras.callbacks.TensorBoard(log_dir='{}/'.format(model_checkpoint_folder)),
]

#reshape feature vectors for CONV1D layer:

#train_vectors = comp_scaled_vectors.shape(comp_scaled_vectors.shape[0], comp_scaled_vectors.shape[1], 1)

train_samples = comp_scaled_vectors.shape[0]

print(comp_scaled_vectors.shape)
print(train_labels.shape)

model.fit(batch_generator(comp_scaled_vectors, train_labels, batch_size,train_samples, True),
                    steps_per_epoch=(int(train_samples/batch_size)), epochs=epochs_num, callbacks=my_callbacks, verbose=1)



#because we used BatchNormalization, results can fluctuate a lot between model training runs. We run a test on randomly sampled RefSeq r202 synthetic genomes to find the best models across all MIMAG cut-offs to ensure no bias 

#first we need to load in test data:

#change to csv in the future
test_data = f'{current_dir}/old_vectors/testdata.csv'
test_data = pd.read_csv(test_data, sep='\t')

nnmodels = listdir(f'{model_checkpoint_folder}/')
nnmodels = [n for n in nnmodels if n.endswith('.h5')]

logging.info(f'Selecting from {len(nnmodels)} generated nn models')
dfs = []




bestH, bestM, bestL, bestName = 0,0,0,None
for mod in nnmodels:    
    
    logging.info(f'Now processing: {mod}')
    
    picked_model = keras.models.load_model('{}/{}'.format(model_checkpoint_folder, mod))


    preds, comps = returnPredictors(test_data[test_data['Completeness'] > 89], 'Comp')
    #comps = comps/100


    AR_train = scaler.transform(preds)
    AR_train = AR_train[:, :20021]
    newAR_train = AR_train.reshape(AR_train.shape[0], AR_train.shape[1], 1)

    predictions = picked_model.predict(newAR_train)
    predictions[predictions > 100] = 100
    predictions[predictions < 0] = 0

    predictions = predictions.flatten()
    predictions = predictions * 100
    
    error = predictions - comps

    x1 = pd.DataFrame({'Error': error})
    x1['Model'] = mod
    x1['MIMAG'] = 'HQ'
    


    
    
    mquees = test_data[test_data['Completeness'] > 50]
    mquees = test_data[test_data['Completeness'] < 90]
    
    preds, comps = returnPredictors(mquees, 'Comp')
    #comps = comps/100


    AR_train = scaler.transform(preds)
    AR_train = AR_train[:, :20021]
    newAR_train = AR_train.reshape(AR_train.shape[0], AR_train.shape[1], 1)

    predictions = picked_model.predict(newAR_train)
    predictions[predictions > 100] = 100
    predictions[predictions < 0] = 0

    predictions = predictions.flatten()
    predictions = predictions * 100
    
    error = predictions - comps
    #MQs.append(error)
    
    x3 = pd.DataFrame({'Error': error})
    x3['Model'] = mod
    x3['MIMAG'] = 'MQ'

    
    
    preds, comps = returnPredictors(test_data[test_data['Completeness'] < 50], 'Comp')
    #comps = comps/100

    AR_train = scaler.transform(preds)
    AR_train = AR_train[:, :20021]
    newAR_train = AR_train.reshape(AR_train.shape[0], AR_train.shape[1], 1)

    predictions = picked_model.predict(newAR_train) #verbose=1
    predictions[predictions > 100] = 100
    predictions[predictions < 0] = 0

    predictions = predictions.flatten()
    predictions = predictions * 100
    
    error = predictions - comps
    x2 = pd.DataFrame({'Error': error})
    x2['Model'] = mod
    x2['MIMAG'] = 'LQ'

    if bestH + bestL + bestM == 0:
        bestH = abs(x1['Error'].mean())
        bestM = abs(x3['Error'].mean())
        bestL = abs(x2['Error'].mean())
        bestName = mod
        
        
    else:
        if (abs(x1['Error'].mean()) + abs(x3['Error'].mean()) + abs(x2['Error'].mean())) < (bestH + bestM + bestL):
            bestH = abs(x1['Error'].mean())
            bestM = abs(x3['Error'].mean())
            bestL = abs(x2['Error'].mean())
            bestName = mod
            
    logging.info('Current model is {} with MAE for HQ: {} MQ: {} LQ: {}'.format(mod, abs(x1['Error'].mean()), abs(x2['Error'].mean()), abs(x3['Error'].mean())))

            
logging.info(f'Best model was {bestName} with MAE for HQ: {bestH} MQ: {bestM} LQ: {bestL}')

#copy best model
os.replace(f'{model_checkpoint_folder}/{bestName}', f'{current_dir}/new_models/specific_model_COMP.hd5')

#now save all new feature vectors

vectors_out_path = f'{current_dir}/new_vectors/'

scipy.sparse.save_npz(f'{vectors_out_path}/all_{new_GTDB_version}_scaled.npz', comp_scaled_vectors)
np.save(f'{vectors_out_path}/all_{new_GTDB_version}_scaled_complabels.npy', comp_scaled_labels)
scipy.sparse.save_npz(f'{vectors_out_path}/all_{new_GTDB_version}_raw_comps.npz', comp_raw_vectors)
np.save(f'{vectors_out_path}/all_{new_GTDB_version}_raw_compabels.npy', comp_raw_labels)
scipy.sparse.save_npz(f'{vectors_out_path}/all_{new_GTDB_version}_raw_conts.npz', cont_raw_vectors)
np.save(f'{vectors_out_path}/all_{new_GTDB_version}_raw_contlabels.npy', cont_raw_labels)




#add new complete genomes to the reference data for specific/general model selection


old_ref_data = scipy.sparse.load_npz(f'{old_release_loc}/min_ref_rsdata_{previous_release_version}.npz')

refdata = full_pickles[full_pickles['Completeness'] == 100]
refdata = refdata[refdata['Contamination'] == 0]
refdata = csr_matrix(scaler.transform(refdata.iloc[:, :-2].values))
new_ref_vectors =  vstack([old_ref_data, refdata])

scipy.sparse.save_npz(f'{vectors_out_path}/min_ref_rsdata_{new_GTDB_version}.npz', new_ref_vectors)

#need a copy for the database as well as a copy for the CheckM2 update
shutil.copy(f'{vectors_out_path}/min_ref_rsdata_{new_GTDB_version}.npz', f'{new_models}/min_ref_rsdata_{new_GTDB_version}.npz')


#move test data from old_vectors to new_vectors
os.replace(f'{current_dir}/old_vectors/testdata.csv', f'{vectors_out_path}/testdata.csv')



#finally, use zenodobackpack to package up new reference data

creator = zenodo_backpack.ZenodoBackpackCreator()#
creator.create(f"{vectors_out_path}", f"{current_dir}/{new_GTDB_version}_database.zb.tar.gz", f"{new_GTDB_version}")