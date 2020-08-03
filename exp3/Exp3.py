from pylab import *
import scipy.special as sp
A=1.05
B=-0.105
with open('fitting1.dat') as datfile:
    dat = datfile.readlines();
N=len(dat)
k=len(dat[0].split())
vect=array([zeros(k)]*N)
for i in range(0,N):
    dat[i]=dat[i].split()
    for j in range(0,k):
        vect[i][j]=float(dat[i][j])
def g(t,A,B):
    return sp.jn(2,t)*A+B*t

t=vect[:,0]
y=array([zeros(N)]*(k-1))
stdev=array(zeros(k-1))
for i in range(0,9):
    y[i]=vect[:,i+1]
    stdev[i]=std(g(t,A,B)-y[i])
print(stdev)
yx=vect[:,1:]
figure(0)
plot(t,g(t,A,B),color='black')
plot(t,yx)
title('Data to be fitted to theory')
xlabel('t \u2192',size=10)
ylabel('f(t)+noise \u2192',size=10)
grid(True)

figure(1)
plot(t,g(t,A,B),color='black',label='f(t)')
errorbar(t[::5],y[0][::5],yerr=stdev[0],fmt='ro',label='errorbar')
title('Data points for \u03C3 = 0.101 with exact function')
xlabel('t \u2192',size=10)
legend(loc='upper right')
grid(True)

func=c_[sp.jn(2,t),t]
a=linspace(0,2,21)
b=linspace(-0.2,0,21)
error=array([zeros(21)]*21)
for i in range(0,21):
    for j in range(0,21):
        error[i][j]=(square((y[0]-g(t,a[i],b[j])))).mean()
x_vals=linspace(0,20,21)
y_vals=linspace(0,20,21)
figure(2)
cont=contour(a,b,error)
colorbar(cont)
clabel(cont,[0.06,0.12,0.18,0.24],fontsize=8)
scatter(1.05,-0.105,color='red')
annotate('Exact Value',(1.05,-0.105),color='blue')
title('Contour plot of errors')
xlabel('a \u2192',size=10)
ylabel('b \u2192',size=10)

Aerr=array(zeros(k-1))
Berr=array(zeros(k-1))
Aerrabs=array(zeros(k-1))
Berrabs=array(zeros(k-1))
for i in range(0,9):
    p,resid,rank,sig=lstsq(func,y[i])
    Aerrabs[i]=abs(A-p[0])
    Berrabs[i]=abs(B-p[1])
    Aerr[i]=A-p[0]
    Berr[i]=B-p[1]
figure(3)
plot(stdev,Aerrabs,'ro',label='Aerr',linestyle='--',linewidth=1)
plot(stdev,Berrabs,'go',label='Berr',linestyle='--',linewidth=1)
legend(loc='upper left')
xlabel('Noise Standard Deviation \u2192',size=10)
ylabel('$MS Error \u2192$',size=10)
title('Variation of error with noise')
grid(True)

figure(4)
errorbar(stdev,Aerrabs,yerr=std(Aerr),fmt='ro',label='Aerr')
errorbar(stdev,Berrabs,yerr=std(Berr),fmt='go',label='Berr')
title('Variation of error with noise')
xlabel('$\u03C3_n \u2192$',size=10)
ylabel('$MS Error \u2192$',size=10)
legend(loc='upper left')
yscale('log')
xscale('log')
grid(True)
show()
