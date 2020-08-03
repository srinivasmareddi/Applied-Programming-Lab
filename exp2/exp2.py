def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
import numpy as np
import sys
import cmath as cm
circuit = sys.argv[1]
line=[]
listy=[]
J=complex(1j)
try:
    with open(circuit) as myfile:
        line=myfile.readlines();
except IOError:
    print("File doesnot exist")
    sys.exit()

line_list=[]
for i in line:
    line_list.append(i.split(' '))

b = False
for i in line:
    if i == '.circuit\n':
        b=True
    if i in ('.end','.end\n'):
        b=False
    if b :
        listy.append((i.split('#')[0]).split());
for i in line_list:
    if i[0]=='.ac':
        f = 2*(np.pi)*float(i[2])
listy.pop(0)
nodes={}
comp={}
srcnds={}
rad=0
two_terminal={'R','L','C','V','I'}
for i in listy:
    if i[0][0] in two_terminal:
        if i[0][0] in ('V','I'):
            if i[3] == 'dc':
                comp.update({i[0]:float(i[4])})
            else:
                rad=np.deg2rad(float(i[5]))
                comp.update({i[0]:((cm.rect(float(i[4]),rad))/2)})
            srcnds.update({i[0]:[i[1],i[2]]})
        elif i[0][0] == 'R':
            comp.update({i[0]:float(i[3])})
        elif i[0][0] == 'L':
            comp.update({i[0]:(float(i[3])*J*f)})
        elif i[0][0] == 'C':
            comp.update({i[0]:(1/(float(i[3])*J*f))})
    if i[1] not in nodes:
        nodes.update({i[1]:[i[0]]})
    else:
        nodes[i[1]].append(i[0])
    if i[2] not in nodes:
        nodes.update({i[2]:[i[0]]})
    else:
        nodes[i[2]].append(i[0])
n=0
for i in listy:
    if i[0][0]=='V':
        n+=1
nodes.pop('GND')
print('List of nodes:',nodes,'\n')
print('Values of components:',comp,'\n')
print('List of source nodes:',srcnds,'\n')
c=len(nodes)+n
x = np.array([np.zeros(c) for i in range(c)],dtype=complex)
I = np.array(np.zeros(c),dtype=complex)
for i in range(len(nodes)):
    for m in nodes[list(nodes)[i]]:
        if m[0] == 'I':
            if srcnds[m][0]== list(nodes)[i] :
                I[i]=-comp(m);
            else:
                I[i]=comp(m);
    for j in range(len(nodes)):
            if i == j:
                for k in nodes[list(nodes)[i]]:
                    if k[0] in ('R','L','C'):
                        x[i][j]+=1/comp[k]
            else:
                for l in intersection(nodes[list(nodes)[i]],nodes[list(nodes)[j]]):
                    if l[0] in ('R','L','C'):
                        x[i][j] =x[i][j]-1/comp[l]
for i in range(len(nodes)):
    for j in range(len(nodes),c):
        for k in nodes[list(nodes)[i]]:
            if k[0] == 'V':
                if srcnds[k][0] == list(nodes)[j-len(nodes)]:
                    x[i][j]=1;
                elif srcnds[k][1]==list(nodes)[j-len(nodes)]:
                    x[i][j]=-1;
for i in range(len(nodes),c):
    I[i]=comp[list(srcnds)[i-len(nodes)]]
    for j in range(len(nodes)):
        if list(srcnds)[i-len(nodes)][0]=='V':
            if srcnds[list(srcnds)[i-len(nodes)]][0]==list(nodes)[j]:
                x[i][j]=1
            elif srcnds[list(srcnds)[i-len(nodes)]][1]==list(nodes)[j]:
                x[i][j]=-1
print('Current Matrix:',I,'\n')
print('Impedance Matrix:\n',x,'\n')
y=np.linalg.solve(x,I)
print('Solution Matrix:',y,'\n')
for i in range(len(nodes)):
    print('V',list(nodes)[i],'=',y[i])
for i in range(len(nodes),c):
    print('I through',list(srcnds)[i-len(nodes)],'=',y[i])
#adjust for v peak to peak
#adjust for frequency
