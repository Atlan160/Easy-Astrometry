import numpy as np
import os.path
from os import path

def save_text(text, filename="output", ending=".txt", header=""):

    #if !path.exists(filename+ending):
        

    file1= open(filename+ending,"w")
    file1.write(text)
    file1.close()

def save_to_file_2D(x, y, filename="output", acuracy=3, header=""):
    
    output=np.array([x,y])
    output=output.transpose()
    np.savetxt(filename+".dat", output, fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")

def save_to_file_3D(x, y,z, filename="output", acuracy=3, header=""):
    
    output=np.array([x,y,z])
    output=output.transpose()
    np.savetxt(filename+".dat", output, fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")

def save_to_file_4D(x, y,z,t, filename="output", acuracy=3, header=""):
    
    output=np.array([x,y,z,t])
    output=output.transpose()
    np.savetxt(filename+".dat", output, fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")

def save_to_file_5D(a,b,c,d,e, filename="output", acuracy=3, header=""):
    
    output=np.array([a,b,c,d,e])
    output=output.transpose()
    np.savetxt(filename+".dat", output, fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")

def save_to_file_8D(a,b,c,d,e,f,g,h, filename="output", acuracy=3, header=""):
    
    output=np.array([a,b,c,d,e,f,g,h])
    output=output.transpose()
    np.savetxt(filename+".dat", output, fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")

def save_to_file_1D(x, filename="output", acuracy=3, header=""):
    
    output=np.array([x])
    output=output.transpose()
    np.savetxt(filename+".dat", output, fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")

def save_to_file_ND(x, filename="output", acuracy=3, header=""):
    np.savetxt(filename+".dat", x.transpose(), fmt='%.'+str(acuracy)+'e', header=header, delimiter='   ')
    print("file saved")