import matplotlib.pyplot as plt
import numpy as np



class ReadFile():

    def __init__(self,filename : str): 

        self.x = [] 
        self.y = []

        with open(filename,'r') as f:
            lines = f.readlines()
            for line in lines:
                val = [float(s) for s in line.split() ] 
                self.x.append(val[0])
                self.y.append(val[1])

        self.x = np.array(self.x,dtype=np.double)
        self.y = np.array(self.y,dtype=np.double)

    def getData(self, get_experimental_data: bool = True):
        return self.x, self.y

#awk '!/    nan/' dmf_theta_vs_T_0.95sat.txt >> temp && mv temp file.text
#awk '!/   inf/' file.text >> temp && mv temp dmf_theta_vs_T_0.95sat.txt

#file1 = ReadFile("dediff_theta_vs_T_0.1sat.txt")
file2 = ReadFile("data_0.31.txt")
