import fileinput
#import NumPy as np
import sys
o = open("RR22.f","w").close    # reset the file
o = open("RR22.f","a")        # open for append
org = open("RR21.f","r")       # the file with the permission of read
#org = fileinput.FileInput("test.f",inplace=1)
#skipline = np.arange(122,139,1)   # to generate an array
for line_number,line in enumerate( org ):    
       if line.strip()=="E(1:(Estep/9*2)) = 0.0" :
            o.write("E(1:(Estep/9*2)) = -1.4 \n")
       else:
            o.write(line)     
o.close()
fileinput.close()



















