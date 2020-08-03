from pylab import *
rcParams['axes.grid']=True
# Problem 1
n=2**7
x=linspace(0,2*pi,n+1);x=x[:-1]
y=sin(x)**3
Y=(fftshift(fft(y))/n)
w=linspace(-64,63,n)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin^3(t)$")
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii])*180/pi,'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
y=cos(x)**3
Y=fftshift(fft(y))/n/2
w=linspace(-64,63,n)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos^3(t)$")
subplot(2,1,2)
plot(w,angle(Y)*180/pi,'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii])*180/pi,'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
#Problem 2
n=2**10
x=linspace(-2*pi,2*pi,n+1);x=x[:-1]
y=cos(20*x+5*cos(x))
Y=fftshift(fft(y)/n)*2
w=linspace(-n/2,n/2,n+1);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-70,70])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos(20x+5\cos(x))$")
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii])*180/pi,color='orange',lw=2)
plot(w[ii],angle(Y[ii])*180/pi,'ro',lw=2)
xlim([-70,70])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
#Problem 3
#Accuracy to 6 Digits?
n=2**20
x=linspace(-8*pi,8*pi,n+1);x=x[:-1]
y=exp(-x**2/2)
Y=(fftshift(fft(y))/n)*8
w=linspace(-n/2,n/2,n+1);w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-30,30])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $e(-x^2/2)$")
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii])*180/pi,color='hotpink',lw=2)
xlim([-30,30])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
show()
