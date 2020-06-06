import sys
from abaqus import *
from abaqusConstants import *
from odbAccess import *
from abaqusConstants import *
from types import IntType
import numpy as np
import math
import os
######################################################################################
#
#               SDV extraction script
#                       Created by Shlomo Spizter
#                           17.04.2020
#
#######################################################################################
#
#   1)  Make sure your ABAQUS working dir has the .odb requested. the script will look there.
#   2)  Update lines 23-34 and run. 
#   
#######################################################################################

odb_name = "6by4_SL_L1_10042020"#'6by4_SL_L1_10042020' # to be filled <-----------------------------------------------------------------------------------
#### MAKE SURE YOUR WORKING DIR HAS THE .ODB#####

# desired name for the results file
FileName = odb_name + '_SDVs.csv'                 # to be filled <------------------------------------------------------------

# Custom variables names

TheStep = 'LVI'                #CAPS!                       # to be filled <-------------------------------------------
ThePart = 'SL_PLATE_6BY4-1' #remember, assembly adds "-1"   # to be filled <-------------------------------------------
TheElementSet = "PLATE_SECTION" #(as in the part)           # to be filled <-------------------------------------------
TheSDVindexList = np.array([7,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54])#([7,8,11,12,13,14,15,16,22,23,24,25,26,27,33,34,35,36,37,38,44,45,46,47,48,49]) # to be filled <-------------------------------------------
              # ^-----------------------------    # to be filled <-------------------------------------------

#######################################################################################

# Opening the odb
odb = openOdb(odb_name + '.odb', readOnly=True)

############## Select the last frame for damage parameters
Lastframe = odb.steps[TheStep].frames[-1]

plate_elements_subset = odb.rootAssembly.instances[ThePart].elementSets[TheElementSet]

############## Create matrix to store values

tempValueList = Lastframe.fieldOutputs["SDV1"].getSubset(region=plate_elements_subset).values
numElements = len(tempValueList)
NumValues = len(TheSDVindexList)+1
SDVmatrix = np.zeros([numElements, NumValues])

indexNum = 0 # index 0 left for label
print("SDV extraction initiated")
for i in TheSDVindexList:

    index = "SDV{}".format(str(i))

    #SDVs = Lastframe.fieldOutputs[index]
    #SDVs_name = Lastframe.fieldOutputs[index].name
    #SDVs_value = Lastframe.fieldOutputs[index].values
    Plate_SDV_value = Lastframe.fieldOutputs[index].getSubset(region=plate_elements_subset).values
    if indexNum == 0:
        for ele in Plate_SDV_value:
            SDVmatrix[ele.elementLabel-1][0] = ele.elementLabel   #here the is a possible problem about element labels and numbering, consider different label and #ing     
    indexNum = indexNum +1
    for element in Plate_SDV_value:
        SDVmatrix[element.elementLabel-1][indexNum] = element.data#
        
print("Finished creating SDV Matrix, Creating .csv file......")    
try:
    np.savetxt(FileName, SDVmatrix, delimiter=",")
    print("CSV File was created under the name of "+FileName)
except:
    print("Original name is taken.")
    try:
        np.savetxt(odb_name + '_SDVs_temp.csv', SDVmatrix, delimiter=",")
        print("CSV File was created under the temporary name of "+ odb_name + "_SDVs_temp.csv \n")
        print("because the original name was taken, please rename the created file name.")
        print("Have a good day!")
    except:
        print("The file failed to save as the name is taken, please solve file naming and versions to avoid data loss and try again.")
        print("Have a great day!")


