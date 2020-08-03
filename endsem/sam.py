from pylab import *
import sys

Lx=0.10
Ly=0.20
Vs=1
M=25
N=50
k=6
delta=0.0000001
N0=1500

x=linspace(int(0),int(M),M)
y=linspace(int(0),int(N),N)
tank=array([zeros(M)]*N)

karr=linspace(1*(N/10),9*(N/10),9)
error=array([zeros(N0)]*9)

def tankvol(k):
    prev=tank.copy()
    tank[1:-(k+1),1:-1]=0.25*(tank[1:-(k+1),0:-2]+tank[1:-(k+1),2:]+tank[0:-(k+2),1:-1]+tank[2:-k,1:-1])
    tank[-(k+1),1:-1]=(1/3)*(2*tank[-k,1:-1]+tank[-(k+2),1:-1])
    tank[-k:-1,1:-1]=0.25*(tank[-k:-1,0:-2]+tank[-k:-1,2:]+tank[-(k+1):-2,1:-1]+tank[-(k-1):,1:-1])
    tank[1:-1,0]=tank[1:-1,1]
    tank[1:-1,-1]=tank[1:-1,-2]
    tank[-1,1:-1]=tank[-2,1:-1]
    tank[1:,0],tank[1:,-1],tank[-1,:],tank[0,:]=0,0,0,Vs
    return [tank,prev]

def func(k):
    kind=int(k*10/N)-1
    for i in range(N0):
        if i==0:
            tankvol(k)
            error[kind][i]=abs(tankvol(k)[0]-tankvol(k)[1]).max()
        else:
            if error[kind][i-1]>=delta:
                tankvol(k)
                error[kind][i]=abs(tankvol(k)[0]-tankvol(k)[1]).max()
            else:
                maxi=i
                break
    return [error,maxi]
for i in range(9):
    maxi=func(int(karr[i]))[1]
    print(maxi)
