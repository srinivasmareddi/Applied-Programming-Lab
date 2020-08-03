from pylab import *
import sys
import mpl_toolkits.mplot3d.axes3d as p3

if len(sys.argv)==1:
    Nx=25
    Ny=25
    r=8
    Niter=1500
elif len(sys.argv)==5:
    Nx=int(sys.argv[1])
    Ny=int(sys.argv[2])
    r=float(sys.argv[3])
    Niter=int(sys.argv[4])
else:
    print("Invalid Arguments")

x=linspace(int(-Nx/2),int(Nx/2),Nx)
y=linspace(int(-Ny/2),int(Ny/2),Ny)
xe=arange(0,Niter,1)
phi=array([zeros(Nx)]*Ny)
jx=array([zeros(Nx)]*(Ny))
jy=array([zeros(Nx)]*(Ny))
error=array(zeros(Niter))
iter=arange(1,Niter+1,1)

Y,X=meshgrid(y,x)
ii=where(X*X+Y*Y<=r*r)
phi[ii]=1.0

for k in range(Niter):
    oldphi=phi.copy()
    phi[1:-1,1:-1]=0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])
    phi[:,0]=phi[:,1]
    phi[0,:]=phi[1,:]
    phi[:,-1]=phi[:,-2]
    phi[ii]=1.0
    error[k]=(abs(phi-oldphi)).max()

jx[:,1:-1]=-0.5*(phi[:,2:]-phi[0:,0:-2])
jy[1:-1,:]=-0.5*(phi[2:,:]-phi[0:-2,:])

figure(1)
title('Contour Plot Of Voltage')
xlabel('x \u2192',size=10)
ylabel('y \u2192',size=10)
cont=contour(x,y,phi)
colorbar(cont)
clabel(cont,inline=True,fontsize=8)
grid(True)

figure(2)
ax=p3.Axes3D(figure(2))
title('The 3-D surface plot of the potential')
surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=cm.jet)
fit1=polyfit(iter,log(error),1)
fit2=polyfit(iter[499:],log(error[499:]),1)

figure(3)
title('Error in Voltage - Semilog Plot')
xlabel('Iteration Count \u2192',size=10)
ylabel('Error in log \u2192',size=10)
semilogy(xe,error,'b',label='error')
semilogy(iter,abs(exp(fit1[1])*exp(fit1[0]*iter)),'r',label='fit for all')
semilogy(iter,abs(exp(fit2[1])*exp(fit2[0]*iter)),'g',label='fit after 500')
grid(True)
legend(loc='upper right')

figure(4)
title('Error in Voltage - Loglog Plot')
xlabel('Iteration Count \u2192',size=10)
ylabel('Error in log \u2192',size=10)
loglog(xe,error,'b',label='error')
loglog(iter,abs(exp(fit1[1])*exp(fit1[0]*iter)),'r',label='fit for all')
loglog(iter,abs(exp(fit2[1])*exp(fit2[0]*iter)),'g',label='fit after 500')
grid(True)
legend(loc='upper right')

figure(5)
title('1V Region On Plate')
xlabel('x \u2192',size=10)
ylabel('y \u2192',size=10)
scatter(ii[0]-12,ii[1]-12,color='r',marker='o')

figure(6)
title('Current Density Profile On Plate')
xlabel('x \u2192',size=10)
ylabel('y \u2192',size=10)
quiver(x,y,jx,jy,scale=1.2,scale_units='inches',label='J Vectors')#decrease size of arrows
scatter(ii[0]-12,ii[1]-12,color='r',marker='o')

show()
