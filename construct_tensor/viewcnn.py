import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

dimension=64
input_tensor=open("image_mo.txt",'r')

yy=np.empty([dimension,dimension],dtype=float)

count=-1
for line in input_tensor:
    count+=1
    line=line.strip().split()
    yy[count][:]=np.asarray(line,dtype=float)

input_tensor_s1=open("image_set.txt",'r')
input_tensor_s2=open("image_seb.txt",'r')

se_t=np.empty([dimension,dimension],dtype=float)
se_b=np.empty([dimension,dimension],dtype=float)

count=-1
for line in input_tensor_s1:
    count+=1
    line=line.strip().split()
    se_t[count][:]=np.asarray(line,dtype=float)

count=-1
for line in input_tensor_s2:
    count+=1
    line=line.strip().split()
    se_b[count][:]=np.asarray(line,dtype=float)

print("Mo max val:",np.max(yy))
print("Mo max val:",np.max(se_t))
print("Mo max val:",np.max(se_b))
#max_val=np.max(se_t)
#print(max_val)
#se_t=se_t/max_val
#print(np.max(se_t))
#img = Image.fromarray(yy)
plt.imshow(yy,cmap=plt.cm.cool)
plt.show()
plt.imshow(se_t,cmap=plt.cm.binary)
plt.show()
plt.imshow(yy+se_t/2.0+se_b/2.0,cmap=plt.cm.binary)
plt.show()

