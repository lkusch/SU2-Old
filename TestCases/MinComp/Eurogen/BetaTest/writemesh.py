import sys
import numpy
import string
import os

#path=sys.argv[1]

a=10
b=50
file=open("mesh.txt","w") 
file.write("NDIME= 2"+"\n")
file.write("NELEM= "+str(a*b)+"\n") 
counter=0        
ps1=0
ps2=4
ps3=4+2*(a-1)+2*(b-1)
ps4=4+2*(a-1)+2*(b-1)
for i in range(0,a):
    p1=ps2+(i-1)
    p2=p1+1
    p3=ps3+i*(b-1)
    p4=ps3+(i-1)*(b-1)
    if(i==0):
        p1=ps1
        p2=ps2
        p4=p3-1
    if(i==a-1):
        p2=ps1+1
        p3=ps2+(a-1)
    file.write("9        "+str(p1)+"         "+str(p2)+"         "+str(p3)+"         "+str(p4)+"         "+str(counter)+"\n")
    counter=counter+1
    for j in range(1,b-1):
        p1=p4
        p2=p3
        p3=p2+1
        p4=p1-1
        if(i!=0): 
            p4=p1+1
        file.write("9        "+str(p1)+"         "+str(p2)+"         "+str(p3)+"         "+str(p4)+"         "+str(counter)+"\n")
        counter=counter+1
    p1=p4
    p2=p3
    p3=4+2*(a-1)+(b-1)-i-1
    p4=p3+1
    if(i==0): 
        p4=3
    if(i==a-1):
        p3=2
    file.write("9        "+str(p1)+"         "+str(p2)+"         "+str(p3)+"         "+str(p4)+"         "+str(counter)+"\n")
    counter=counter+1
file.write("NPOIN= "+str(a*b+a+b+1)+"\n")
counter=0
p1=0
file.write(str(0)+" "+str(a)+" "+str(counter)+"\n")
counter=counter+1
file.write(str(0)+" "+str(0)+" "+str(counter)+"\n")
counter=counter+1
file.write(str(b)+" "+str(0)+" "+str(counter)+"\n")
counter=counter+1
file.write(str(b)+" "+str(a)+" "+str(counter)+"\n")
counter=counter+1
for i in range(1,a):
    file.write(str(0)+" "+str(a-i)+" "+str(counter)+"\n")
    counter=counter+1
for i in range(1,b):
    file.write(str(i)+" "+str(0)+" "+str(counter)+"\n")
    counter=counter+1
for i in range(1,a):
    file.write(str(b)+" "+str(i)+" "+str(counter)+"\n")
    counter=counter+1
for i in range(1,b):
    file.write(str(b-i)+" "+str(a)+" "+str(counter)+"\n")
    counter=counter+1
for i in range(1,a):
    for j in range(1,b):
        file.write(str(j)+" "+str(a-i)+" "+str(counter)+"\n")
        counter=counter+1
file.write("NMARK=3"+"\n")
file.write("MARKER_TAG= left\n")
file.write("MARKER_ELEMS= "+str(a)+"\n")
p1=0
p2=4
for i in range(1,a+1):
    file.write("3 "+str(p1)+" "+str(p2)+"\n")
    p1=p2
    p2=p1+1
    if(i==a-1):
        p2=1
file.write("MARKER_TAG= upperleft\n")
file.write("MARKER_ELEMS= 1\n")
file.write("3 "+str(0)+" "+str(4+2*(a-1)+2*(b-1)-1)+"\n")
file.write("MARKER_TAG= lowerright\n")
file.write("MARKER_ELEMS= 1\n")
file.write("3 "+str(2)+" "+str(4+(a-1)+(b-1)-1)+"\n")




        
