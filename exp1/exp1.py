import sys
circuit = sys.argv[1]
line=[]
list=[]
try:
    with open(circuit) as myfile:
        line=myfile.readlines();
except IOError:
    print("File doesnot exist")
    sys.exit()
b=False

for i in line:
    if i == '.circuit\n':
        b=True
    if i =='.end\n':
        b=False
    if b :
        list.append((i.split('#')[0]).split());
list.reverse()
list.pop()

for k in list:
    k.reverse()
    for x in range(len(k)):
        print(k[x],end=" ")
    print()
