from pylab import *

#Values given in question
Lx=0.10
Ly=0.20
Vs=1

color=['blue','red','green','hotpink','orange','black','purple','cyan','olive']

#Option of chosing default or custom values
#Accepts both capital and small letters
opt=input('Enter Y/N for default values: ')
if opt=='Y' or opt=='y':
    M=25
    N=50
    k=6
    delta=10**-13
    N0=1500
elif opt=='N' or opt=='n':
    M=int(input("Number of nodes along x(M): "))
    N=int(input("Number of nodes along y(N): "))
    k=int(input("Index height of liquid: "))
    delta=float(input("Least value of error: "))
    N0=int(input('Max no of iterations: '))
else:
    print('Invalid Input')
    quit()
#Creating arrays for use in loop
x=linspace(int(0),int(M),M)
y=linspace(int(0),int(N),N)
Ex=array([zeros(M)]*(N))
Ey=array([zeros(M)]*(N))
tank=array([zeros(M)]*N)
#Creating Arrays for storing values for multiple h/Ly values
karr=linspace(1*(N/10),9*(N/10),9)
Qtoparr=array(zeros(9))
Qfluidarr=array(zeros(9))
Exarr=array([[zeros(M)]*(N)]*9)
Eyarr=array([[zeros(M)]*(N)]*9)
tankarr=array([[zeros(M)]*(N)]*9)
#error array was created with ones so that when it is called in the for loop
#it continues evaluating for error[-1] for the first time
error=array([ones(N0)]*9)
maxiarr=array(zeros(9))
#Function to be called for evaluating potential for one iteration
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

#Algorithm for calculating final values of E,tank,error,Q for different values of h/Ly
for k in karr:
    k=int(k)
    kind=int(k*10/N)-1
    #Algorithm for calculating upto certain accuracy value
    for i in range(N0):
        if error[kind][i-1]>=delta:
            tankvol(k)
            error[kind][i]=abs(tankvol(k)[0]-tankvol(k)[1]).max()
        else:
            maxiarr[kind]=i
            break
    #calculating e field
    # Ex[:,1:-1]=-0.5*(tank[:,2:]-tank[0:,0:-2])*(M/Lx)
    # Ey[1:-1,:]=-0.5*(tank[2:,:]-tank[0:-2,:])*(N/Ly)
    Ex[:,0:-1]=-1*(tank[:,1:]-tank[0:,0:-1])*(M/Lx)
    Ey[0:-1,:]=-1*(tank[1:,:]-tank[0:-1,:])*(N/Ly)
    #Calculating charge on top and in fluid walls
    Qtop,Qfluid=0,0
    for i in range(M):
        Qtop+=(Lx/M)*Ey[1,i]*8.854*10**(-12)
        Qfluid-=(Lx/M)*Ey[-2,i]*8.854*10**(-12)*2
    for i in range(N-k-1,N):
        Qfluid+=(Ly/N)*(Ex[i,1]-Ex[i,-2])*8.854*10**(-12)*2
    #assigning induvidual array elements to a larger array
    Qtoparr[kind]=Qtop
    Qfluidarr[kind]=Qfluid
    Exarr[kind]=Ex
    Eyarr[kind]=Ey
    tankarr[kind]=tank

#array of no of iterations done for each value of h/Ly
print('The Number of iterations done for h/Ly=',karr)
print(maxiarr,'\n')

#verification of continuity at the surface of the liquid
k1=int(karr[1])
E1x=Exarr[1][-k1-2,int(M/2-1)]
E1y=Eyarr[1][-k1-2,int(M/2-1)]
print('Ex in air: ',E1x)
print('Ey in air: ',E1y)
E2x=Exarr[1][-k1-1,int(M/2-1)]
E2y=Eyarr[1][-k1-1,int(M/2-1)]
print('Ex in fluid: ',E2x)
print('Ey in fluid: ',E2y)
print('ratio of normal components of electric fields =',E1y/E2y,'\n')

#verification of snells law
E1=complex(Exarr[4][-k1-1,13],Eyarr[4][-k1-1,13])
E2=complex(Exarr[4][-k1,13],Eyarr[4][-k1-2,13])
snell=cos(angle(E1))/cos(angle(E2))
print('The ratio of sine of incident,refracted angle is:',snell)

#Polyfit of Qtoparr
coeff=polyfit(karr,Qtoparr,4)
pfit=poly1d(coeff)
t=linspace(3,50,5000)

#------Plots------
#quiver and contour plot
figure(0)
title('Voltage Contour Plot + E field')
cont=contour(x,y,tankarr[4])
colorbar(cont)
clabel(cont,inline=True,fontsize=8)
xlabel('x \u2192',size=10)
ylabel('y \u2192',size=10)
q1=quiver(x,y,Exarr[2],Eyarr[2],label='E Vectors')
quiverkey(q1,0.8,0.9,0.5,label='E Vectors',labelcolor='grey',color='grey')
#Charge plots
#Qtop plot
figure(1)
title('Qtop vs k')
xlabel('k \u2192')
ylabel('Qtop \u2192')
grid(True)
plot(karr,Qtoparr,'ro',label='Qtop per m')
# plot(karr,Qtoparr,'r:')
#Polyfit plot
figure(1)
plot(t,pfit(t),'b-',label='Polyfit of Qtop')
legend(loc='upper left')
#Qfluid plot
figure(2)
title('Qfluid vs k')
xlabel('k \u2192')
ylabel('Qfluid \u2192')
grid(True)
plot(karr,Qfluidarr,'go',label='Qfluid per m')
plot(karr,Qfluidarr,'g:')
legend(loc='upper right')
#error vs iteration plots
figure(3)
title('Error in Potential vs Iterations')
xlabel('k \u2192')
ylabel('Error \u2192')
grid(True)
plot(linspace(0,int(maxiarr[2]),int(maxiarr[2])),error[2][:int(maxiarr[2])],color=color[2])
figure(4)
title('Error in Potential vs Iterations(multiple cases)')
xlabel('k \u2192')
ylabel('Error \u2192')
grid(True)
for i in range(1,9):
    plot(linspace(0,int(maxiarr[i]),int(maxiarr[i])),error[i][:int(maxiarr[i])],color=color[i])
show()
