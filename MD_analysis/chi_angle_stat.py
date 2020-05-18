import pandas
import numpy
import matplotlib.pyplot as plt
import seaborn
import os
from scipy.stats import ttest_ind, levene
from __future__ import division

from scipy.special import stdtr
from numpy import inf

#reading Dimer2 files

chiDimer2A= {}
distanceDimer2A = {}

path = "/mdspace/dora/DIO3/MD/DIO3_dimer/DIO3_1QQ2/2/Analysis"

for i in range(5):
    name1 = 'ChiPHEA_' + str(i) + ".xvg"
    name2 = 'distA_' + str(i) + ".xvg"
    chiDimer2A[i] = pandas.read_csv(os.path.join(path,name1),header=None,sep=" ")
    distanceDimer2A[i] = pandas.read_csv(os.path.join(path,name2),header=None,sep=" ")

chiDimer2B= {}
distanceDimer2B= {}

for i in range(5):
    name1 = 'ChiPHEB_' + str(i) + ".xvg"
    name2 = 'distB_' + str(i) + ".xvg"
    chiDimer2B[i] = pandas.read_csv(os.path.join(path,name1),header=None,sep=" ")
    distanceDimer2B[i] = pandas.read_csv(os.path.join(path,name2),header=None,sep=" ")

# reading Monomer files
chiMonomer = {}
distanceMonomer = {}



path = "/mdspace/dora/DIO3/MD/DIO3_dimer/DIO3_mono/Analysis"

for i in range(5):
    name1 = 'ChiPHEA_' + str(i) + ".xvg"
    name2= 'dist_' + str(i) + ".xvg"
    chiMonomer[i] = pandas.read_csv(os.path.join(path,name1),header=None,sep=" ")
    distanceMonomer[i] = pandas.read_csv(os.path.join(path,name2),header=None,sep=" ")

#Defining the functions

tresh1=75
tresh2=0.75

def toSeabornCountplotData (dataFrame1,dataFrame2,treshold1=tresh1,treshold2=tresh2):
    array1 = pandas.to_numeric(dataFrame1.iloc[:, 1]).abs()
    array2 = pandas.to_numeric(dataFrame2.iloc[:, 1]).abs()
    array = (array1 < treshold1) & (array2 < treshold2)
    array = array.replace(True, value="Open")
    array = array.replace(False, value="Closed")
    return (array)

def toAverageBarPlotData (dataFrame1,dataFrame2,treshold1=tresh1,treshold2=tresh2):
    array1 = pandas.to_numeric(dataFrame1.iloc[:, 1]).abs()
    array2 = pandas.to_numeric(dataFrame2.iloc[:, 1]).abs()
    array = (array1 < treshold1) & (array2 < treshold2)
    array = array.replace(True, value="Open")
    array = array.replace(False, value="Closed")
    A = array.value_counts(sort=False)
    if ((A[0] == 752) | (A[0] == 721)):
        if (A.index[0]=="Open" ):
            A = A.append(pandas.Series([0], index=["Closed"]))
        else:
            A = A.append(pandas.Series([0], index=["Open"]))
    return (A)

def toTimeSeriesData(dataFrame1,dataFrame2,treshold1=tresh1,treshold2=tresh2):
    array1 = pandas.to_numeric(dataFrame1.iloc[:, 1]).abs()
    array2 = pandas.to_numeric(dataFrame2.iloc[:, 1]).abs()
    array = (array1 < treshold1) & (array2 < treshold2)
    array = array [0:721]
    matrix = array.as_matrix()
    return(matrix)

###Seperate barplots

f, axes = plt.subplots(1, 15, figsize=(1, 10), sharey=True)

for i in range(5):
    seaborn.countplot(toSeabornCountplotData(chiMonomer[i], distanceMonomer[i]), ax=axes[i], order=["Closed", "Open"])
    seaborn.countplot(toSeabornCountplotData(chiDimer2A[i], distanceDimer2A[i]), ax=axes[i+5], order=["Closed", "Open"])
    seaborn.countplot(toSeabornCountplotData(chiDimer2B[i], distanceDimer2B[i]), ax=axes[i+10], order=["Closed", "Open"])

#Average barplots

mono = [toAverageBarPlotData(chiMonomer[i],distanceMonomer[i])["Open"]/toAverageBarPlotData(chiMonomer[i],distanceMonomer[i])["Closed"] for i in range(5)]

dimer2A = [toAverageBarPlotData(chiDimer2A[i],distanceDimer2A[i])["Open"]/toAverageBarPlotData(chiDimer2A[i],distanceDimer2A[i])["Closed"] for i in range(5)]

dimer2B = [toAverageBarPlotData(chiDimer2B[i],distanceDimer2B[i])["Open"]/toAverageBarPlotData(chiDimer2B[i],distanceDimer2B[i])["Closed"] for i in range(5)]

dimer = []

for i in range(5):
    if dimer2A[i]> dimer2B[i]:
        dimer.append(dimer2A[i])
    else:
        dimer.append(dimer2B[i])
from numpy import inf

mono= numpy.asarray(mono)
dimer = numpy.asarray(dimer)
dimer = numpy.nan_to_num(dimer)



f, axes = plt.subplots(1, 1, figsize=(1,5), sharey=True)

axes.set_ylabel('Open/Closed')


seaborn.stripplot(mono,orient="v",color="blue")

seaborn.stripplot(dimer,orient="v",color="red")
axes.legend(["monomer","dimer"])






#Time series data


#plotarray = numpy.append([toTimeSeriesData(dataMonomer_0,distanceMonomer_0)],[toTimeSeriesData(dataMonomer_1,distanceMonomer_1)],axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataMonomer_2,distanceMonomer_2)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataMonomer_3,distanceMonomer_3)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataMonomer_4,distanceMonomer_4)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer1_0,distanceDimer1_0)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer1_1,distanceDimer1_1)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer1_2,distanceDimer1_2)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer1_3,distanceDimer1_3)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer1_4,distanceDimer1_4)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer2_0,distanceDimer2_0)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer2_1,distanceDimer2_1)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer2_2,distanceDimer2_2)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer2_3,distanceDimer2_3)]),axis=0)
#plotarray = numpy.concatenate((plotarray,[toTimeSeriesData(dataDimer2_4,distanceDimer2_4)]),axis=0)






####ploting through time
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(plotarray, aspect='auto', cmap=plt.cm.gray,interpolation='nearest')
labels = [ "monomer0","monomer1","monomer2","monomer3","monomer4","dimer1_0","dimer1_1","dimer1_2","dimer1_3","dimer1_4","dimer2_0",
           "dimer2_1","dimer2_2","dimer2_3","dimer2_4"]

ax.set_yticks(numpy.arange(len(labels)))
ax.set_yticklabels(labels)




