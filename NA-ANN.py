import os,sys,gc
import pandas as pd
import numpy as np
import timeit
from datetime import datetime
import keras
import tensorflow as tf
from keras import backend as K
from keras import layers
from keras.models import Sequential
from keras.layers import Dense, Dropout, GaussianNoise, AlphaDropout, Reshape
from keras.layers import Flatten, LocallyConnected1D, LocallyConnected2D
from keras.optimizers import Adam, Adagrad, Adadelta
from keras.backend.tensorflow_backend import set_session 
from keras.utils import to_categorical
##set default values

learning_rate = 0.01
JobID = 1
ps = 25
optim = "adam"
X_file = "KE.geno.csv"
Y_file = "KE_pheno.csv"
CV_file = "KE_cv_pw.csv"
label = "DtSILK"
start_time = timeit.default_timer()
act="relu"
drop_rate = str('0.5,0.5,0.5')
arc = str('63,63')
DG = 'D,D,D,D,D,G'
LC = True
training_epochs = 25
hyp = False

### parse command line arguments

for i in range (1,len(sys.argv),2):
    if sys.argv[i] == "-x":
        X_file = sys.argv[i+1]
    elif sys.argv[i] == "-y":
        Y_file = sys.argv[i+1]
    elif sys.argv[i] == "-cv":
        CV_file = sys.argv[i+1]
    elif sys.argv[i] == "-JobID":
        JobID = int(sys.argv[i+1])
    elif sys.argv[i] == "-label":
        label = sys.argv[i+1]
    elif sys.argv[i] == "-act":
        act = str(sys.argv[i+1])
    elif sys.argv[i] == "-epochs":
        training_epochs = int(sys.argv[i+1])
    elif sys.argv[i] == "-lr":
        learning_rate = float(sys.argv[i+1])
    elif sys.argv[i] == "-arc":
        arc = sys.argv[i+1]
    elif sys.argv[i] == "-ps":
        ps = int(sys.argv[i+1])
    elif sys.argv[i] == "-dr":
        drop_rate=str(sys.argv[i+1])
    elif sys.argv[i] == "-LC":
         LC = bool(sys.argv[i+1])
    elif sys.argv[i] == "hyp":
        hyp = bool(sys.argv[i+1])
    else:
        print('unknown option ' + str(sys.argv[i]))
        quit()

        
        
## change dir to data location


#os.chdir('/home/jaf81qa/jan_storage/tens')
x = pd.read_csv(X_file, index_col = 0, sep="\t") #modification for reading csv-file, seperation is tab
#os.chdir("/storage/full-share/genoPred/maze")
y = pd.read_csv(Y_file, index_col = 0)
cv_folds = pd.read_csv(CV_file,index_col=0)

## select column of phenotype file via columnname

y = y[[label]]
## activity_regularizer=regularizers.l1(0.01)))

def build_network(arc,drop_rate,LC,DG):
    def add_drops(model,drop_out,k):
        if DG[k].upper() == 'D':
            model.add(Dropout(drop_out[0]))
        elif DG[k].upper() == 'G':
            model.add(GaussianNoise(drop_out[k]))
        elif DG[k].upper() == "A":
            model.add(AlphaDropout(drop_out[k]))
        else:
            pass
        return model    
    DG = DG.strip().split(",")
    arc = arc.strip().split(",")
    archit = []
    for layer in  arc:
        archit.append(int(layer))
    layer_number = len(archit)        
    drop_rate = drop_rate.strip().split(",")
    drop_out = []
    for drops in drop_rate:
        drop_out.append(float(drops)) 
    model = Sequential()
    if LC == True:
        model.add(Reshape(input_shape=(x_train.shape[1],x_train.shape[2]),target_shape=(x_train.shape[1],x_train.shape[2]))) #x.shape[2] the different shape of the encoding data 
        model.add(LocallyConnected1D(1,10,strides=7,input_shape=(x_train.shape[1],x_train.shape[2]))) #same as line 107
        model.add(Flatten())
        start = 0
        model = add_drops(model,drop_out,start)
    elif LC == False:
        model.add(Dense(archit[0], kernel_initializer='truncated_normal', activation=act, input_shape=(x_train.shape[1],x_train.shape[2])))  #same as line 107 
        model = add_drops(model,drop_out,start)
    start = 1
    for k in range(start,len(archit)):
        model.add(Dense(archit[k], kernel_initializer='truncated_normal', activation=act))
        model = add_drops(model,drop_out,k)
    model.add(Dense(1, kernel_initializer='truncated_normal'))
    return(model)

     
config = tf.ConfigProto()
#config.gpu_options.per_process_gpu_memory_fraction = 0.1
config.gpu_options.allow_growth = True
set_session(tf.Session(config=config))
#Reshape Data replace NA and encode
y = np.asarray(y)
x = np.asarray(x)
for i in range(x.shape[1]):
    x[np.argwhere(np.isnan(x[:,i])).flatten(),i] = 2
x = to_categorical(x)


if not os.path.isfile("RESULTScv50.txt"):
    out2 = open("RESULTScv50.txt",'w')
    out2.write('DateTime\tCompTime\tDF\tGenos\tPhenos\tCV_fold\tArchit\tConv\tActFun\tEpochs\tdrop_rate\tAccuracy\n' )

    
for k in range(1,51):
    print("Training on cv fold "+ str(k))
    cv = cv_folds['cv_' + str(k)]
    num_cvs = np.ptp(cv) + 1
    
    i = 1
    x_train = x[cv != i] 
    x_test = x[cv == i] 
    y_train = y[cv != i]
    y_test = y[cv == i]

    yhat = np.zeros(shape = y_test.shape)

    model = build_network(arc,drop_rate,LC,DG)
    model.compile(loss='mse', optimizer=Adam(lr=0.01,decay = 0.001),metrics=['accuracy'])
    model.fit(x_train,y_train, epochs=training_epochs , verbose=0) 
#    score = model.evaluate(x_test, y_test, verbose=0)
    bla = model.predict(x_test)
    y_sub= y[np.asarray(cv == i)]
    
    print(model.summary())
    print('\n')
    print(label)        

    comp_time = int(round(timeit.default_timer() - start_time,0))

    DateTime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    acc = np.corrcoef(bla[:,0],np.asarray(y_sub)[:,0])[0,1]

    out2 = open("RESULTScv50.txt", 'a')
    out2.write('%s\t%i\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%i\t%s\t%0.5f\n' % (
        DateTime, comp_time, label, X_file, Y_file, int(k), arc, LC,  act,int(training_epochs), drop_rate, round(acc,4)))

    del model,bla, x_train, x_test, y_train, y_test 
    K.clear_session() 
    gc.collect()
    
    config = tf.ConfigProto()
    #config.gpu_options.per_process_gpu_memory_fraction = 0.1
    config.gpu_options.allow_growth = True
    set_session(tf.Session(config=config))
