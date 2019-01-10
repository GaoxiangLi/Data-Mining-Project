#KDD dataset:
#   Determining protein class (3rd element) by clusters of native sequences (block ID)
#   Protein class is 1 if homologous to the native sequence, 0 if non-homologous(decoys)

import sys
import random
import time

#Each Data_Point corresponds to each line in the dataset
class Data_Point:
  def __init__(self, block_id, example_id, protein_class, features):
    self.block = block_id
    self.id = example_id
    self.protein_class = protein_class
    self.features = features #Array of feature data, not data point metadata

def load_kdd(filename):
  with open(filename, 'r') as kdd_set:
    print("Reading in dataset...")
    for line in kdd_set:
      values = line.split()
      
      #string values are converted into proper types and inserted
      block_id = int(values[0])        #First value is Block ID, multiple lines per Block ID
      example_id = int(values[1])      #Second value is Example ID, unique per line
      protein_class = int(values[2])   #Third value is Protein Class (the desired answer)
      feature_values = values[3:]        #Remaining 74 values are feature values/data (all floats)
      for i in range(len(feature_values)):
        feature_values[i] = float(feature_values[i])
      
      #Create a new data_point class and add it to the list
      datapoints.append(Data_Point(block_id,example_id,protein_class,feature_values))
  print("Reading complete.\n")    

#Creates a lightweight coreset of size m according to the paper's algorithm
def create_coreset(m):
  print("Creating lightweight coreset...")
  print("---Finding the mean of the data...")
  #-------------------------------------
  #1. Find mu = mean of data points X
  #-------------------------------------
  mean = [] #mean value for each of the values associated with a datapoint
  total_number = 0
  
  for datapoint in datapoints:
    #if mean is empty, create initial values
    if total_number == 0:
      for i in range(len(datapoint.features)):
        mean.append(float(0))
        
    #add up feature values to mean array
    for i in range(len(datapoint.features)):
      mean[i] += datapoint.features[i]
      
    total_number += 1
    
  #Once mean array is populated, calculate means
  for i in range(len(mean)):
    mean[i] /= total_number

  #-------------------------------------
  #2. Create summation of distances d(x,mean)^2 
  #-------------------------------------
  print("---Finding differences squared sum between the mean and data values...")
  distances_sum = [] #each value represents the summation of differences between the mean and each value in the dataset squared
  distances = []     #Separate array keeping d(x,mu)^2 for each of the datapoints
  initialized = 0
  
  #Sum up differences for each example, and add to array
  for datapoint in datapoints:  
    #Initialize distances_sum array
    if initialized == 0:
      for i in range(len(datapoint.features)):
        distances_sum.append(float(0))
      initialized = 1
    
    temp_distance = 0.0 #Reset distance value for each of the datapoints
    
    #Add new distance to distances_sum array's existing value and to temp_distance
    for i in range(len(datapoint.features)):
      distances_sum[i] += abs(mean[i]-datapoint.features[i])
      temp_distance += abs(mean[i]-datapoint.features[i])
    
    distances.append(temp_distance**2) #Add the current datapoint's d(x,mu)^2, useful in second part of probability equation creation
  
  #Calculate total distance from mean for second part of probability equation and then square it
  total_distance = 0.0
  for value in distances_sum:
    total_distance += value
  total_distance = total_distance**2
  
  #-------------------------------------
  #3. Create q(x) probability for each example in dataset
  #-------------------------------------
  print("---Creating q(x) probability array...")
  q = []
  uniform_distribution = 0.5*(1/float(total_number)) #Used in first part of probability equation
  
  #Algorithm step 3: create probability distribution
  for i in range(len(datapoints)):
    q.append(uniform_distribution+0.5*(distances[i]/total_distance))
  
  #-------------------------------------
  #4. Sample m weighted points for coreset using probability distribution
  #-------------------------------------
  print("---Sampling",m,"points to be used in lightweight coreset")
  
  #Generate weights 1/(m*q(x))
  for i in range(len(q)):
    weight = 1/(float(m)*q[i])
    q[i] = weight
  #Error checking so m is not more than total dataset
  if int(m) >= total_number:
    for i in range(len(q)):
      coreset.append(i)
  else:
    #Sample m datapoints based on weights in q
    
    id_list = list(range(len(datapoints))) #Create a list of ids for all datapoints
    for i in range(int(m)):
      #From a list of all ids of datapoints, use weights q to sample a single datapoint id and add it to the coreset id list
      selected_id = random.choices(id_list, q)[0] #returns a list with only 1 value, so we take the first value
      coreset.append(selected_id)
      
      #Make the chosen id's weight zero so it doesn't get chosen again
      q[selected_id] = 0.0
  print("Coreset creation complete.\n")
  
  
#Takes the coreset and writes it to a new file called export.dat  
def export_coreset():
  new_line = '' #prevents a newline being printed at the beginning
  export_file = open("export.dat", "w+")
  print("Exporting lightweight coreset...")
  for id_value in coreset:
    datapoint = datapoints[id_value]
    export_file.write(new_line+str(datapoint.block)+'\t'+str(datapoint.id)+'\t'+str(datapoint.protein_class))
    for value in datapoint.features:
      export_file.write('\t'+str(value))
    new_line = '\n' #every line after the first has a newline printed
  export_file.close()
  print("Export complete.")  
  
#-------------------------------------
#-----------Main Code Block-----------
#-------------------------------------
start_time = time.time()
datapoints = [] #Master list of all data points
coreset = []    #ids of selected datapoints for coreset

#parameter checking and assignment
if len(sys.argv) != 3:
  print("Incorrect number of parameters.")
  print("usage: python coresets.py dataset_filename m")
  exit()
else:
  m = sys.argv[2]
  filename = sys.argv[1]
	
load_kdd(filename)  #load dataset
create_coreset(m)   #create corresponding lightweight coreset with size m
export_coreset()    #export finished lightweight coreset
print('\n'+"Elapsed Program Time:",(time.time()-start_time),"seconds.")