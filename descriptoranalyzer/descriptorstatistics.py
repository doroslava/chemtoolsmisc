import pandas
import numpy
import matplotlib.pyplot as plt
import seaborn
import os

def read_descriptor_file(file_prefix, file_type, n_rep = 2):
    """
    Function to read  descriptor files from MD simulations. It can
    also read multiple replicas. For the input file it assumes that they
    contain two columns - first for the time step and second for the measured 
    descriptor. It stores the replicas in the dictionary.
    
    Arguments:
        file_prefix (char): file prefix for the files to be read
        file_type (char) : file sufix (type) for the file to be read
        n_rep (int) : number of the replicas to be read
        
    Returns:
        descriptor_dict (dictionary): dictionary containing descriptor information
        for every timestep of MD simulation. It can also contain information
        for multiple replicas
        
    """
    descriptor_list = []
    for i in range(n_rep):
        file_name = file_prefix + str(i) + "." + file_type
        descriptor_dict = pandas.read_csv(os.path.join("data/", file_name), header=None, sep=" ")
        descriptor_dict[0] = descriptor_dict[0].astype(int)
        descriptor_dict = descriptor_dict.rename(columns={1:file_prefix})
        
        descriptor_list.append(descriptor_dict)
        
    return descriptor_list

def merge_descriptors(descriptor_list):
    
    final_list = []
    
    num_descriptors = len(descriptor_list)
    num_replicas = len(descriptor_list[0])
    
    for i in range(num_replicas):
        
        df = descriptor_list[0][i][[0]]
        
        
        for j in range(num_descriptors):
                
            df = df.merge(descriptor_list[j][i], how = "left")
                 
        final_list.append(df)
        
    return (final_list)        
    
def define_conformation_states (descriptor_list, treshold = [75, 0.75]):
    
    tmp_descriptor_list = merge_descriptors(descriptor_list)    
    
    final_conformation_states = []
    
    length = len(tmp_descriptor_list)
    
    for i in tmp_descriptor_list:
        
            define_state = lambda x: x[1:].abs() > treshold
            
            conformation_states = i.apply(define_state, axis = 1)
             
            conformation_states = pandas.to_numeric(conformation_states.apply(sum, axis = 1))
            
            final_conformation_states.append(conformation_states)
    
    return (final_conformation_states)

def define_conformation_labels(conformation_states_list, old_labels_list, new_labels_list):
    conformation_labels_list = []
    for i in conformation_states_list:
        conformation_labels_list.append(i.replace(old_labels_list, value=new_labels_list))
    return conformation_labels_list

def barplot_conformation_states(conformation_states_list):
    f, axes = plt.subplots(1, len(conformation_states_list), figsize=(5, 10), sharey=True)
    
    for i in range(len(conformation_states_list)):
        seaborn.countplot(conformation_states_list[i], ax=axes[i])
        
def timeseries_conformation_states(conformation_states_list):  
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(numpy.array(conformation_states_list), aspect='auto', cmap=plt.cm.gray,interpolation='nearest')
        





