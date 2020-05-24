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
        descriptor. It stores the replicas for the specified descriptor in the list.

        Arguments:
            file_prefix (char): file prefix for the files to be read
            file_type (char) : file sufix (type) for the file to be read
            n_rep (int) : number of the replicas to be read

        Returns:
            descriptor_list (list): list containing descriptor information
            for every timestep of MD simulation. It can also contain information
            for multiple replicas
        
    """
    descriptor_list = []
    
    for i in range(n_rep):
        file_name = file_prefix + str(i) + "." + file_type
        
        try:
            descriptor = pandas.read_csv(os.path.join("data/", file_name), header=None, sep=" ")
        except:
            print ("File not found")
        
        descriptor[0] = descriptor[0].astype(int)
        descriptor = descriptor.rename(columns={1:file_prefix})
        
        descriptor_list.append(descriptor)
        
    return descriptor_list

def merge_descriptors(descriptor_list):
    """
        Function to merge different descriptors based on the molecular dynamics simulation time frames. As input
        it takes list of the specified descriptors, and merges them in one data frame for each replica. 
        
        Arguments:
            descriptor_list (list): List of pandas dataframes for different descriptors. Each entity in 
            desciptor_list is an additional list with a data frame for each replica. 
        
        Returns:
            final_list (list): list containing pandas dataframe with merged descriptors for replicas
    """

    
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
    """
        Function to define different conformational states through MD based on treshold values for the
        specified descriptors. 
        
        Arguments:
            descriptor_list (list): List of pandas dataframes with merged descriptors for each replica.
            treshold (float) : List of treshold values for the definition of conformation states
        
        Returns:
            final_conformation_states (list): list containing pandas Series with defined conformational states
            through MD for each replica
    """
    
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
     
    """
        
        Function to assign user-specified labels for different conformational states through MD.
        
        Arguments:
            conformation_states_list (list): List of pandas Series with the conformation states
            old_labels_list (list) : List of labels for conformation states in integer format
            new_labels_list (list) : List of labels of conformation states in user - defined format
        
        Returns:
            conformation_labels_list (list): list containing pandas Series with user-defined defined labels for
            conformational states
     """
        
    conformation_labels_list = []
    for i in conformation_states_list:
        conformation_labels_list.append(i.replace(old_labels_list, value=new_labels_list))
    return conformation_labels_list

def barplot_conformation_states(conformation_states_list):
    
    """        
        Function to plot barplot of occurences of different states through MD simulations
        
        Arguments:
            conformation_labels_list (list): list containing pandas Series with defined labels for
            conformational states
        
        Returns:
            None
               
     """
    f, axes = plt.subplots(1, len(conformation_states_list), figsize=(5, 10), sharey=True)
    
    for i in range(len(conformation_states_list)):
        seaborn.countplot(conformation_states_list[i], ax=axes[i])
        
def timeseries_conformation_states(conformation_states_list):  
    """        
        Function to plot time occurences of different states through MD simulations
        
        Arguments:
            conformation_labels_list (list): list containing pandas Series with defined labels for
            conformational states
        
        Returns:
            None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(numpy.array(conformation_states_list), aspect='auto', cmap=plt.cm.gray,interpolation='nearest')
    ax.yaxis.set_ticks([])
    

    
        





