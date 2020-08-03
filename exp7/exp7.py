from __future__ import division
from sympy import *
import pylab as p
import scipy.signal as sp
p.rcParams['axes.grid'] = True
def sympy_to_lti(xpr, s=Symbol('s')):
	""" Convert Sympy transfer function polynomial to Scipy LTI """
	num, den = simplify(xpr).as_numer_denom()  # expressions
	p_num_den = poly(num, s), poly(den, s)  # polynomials
	c_num_den = [expand(p).all_coeffs() for p in p_num_den]  # coefficients
	l_num, l_den = [lambdify((), c)() for c in c_num_den]  # convert to floats
	return sp.lti(l_num, l_den)
s=symbols("s")
def lowpass(R1,R2,C1,C2,G,Vi):
	A=Matrix([ [0,0,1,-1/G], [-1/(1+s*R2*C2),1,0,0] , [0,-1*G,G,1], [(-1/R1)-(1/R2)-(s*C1),1/R2,0,s*C1] ])
	b=Matrix([0,0,0,Vi/R1])
	V=A.inv()*b
	return (A,b,V)
def highpass(R1,R3,C1,C2,G,Vi):
	A=Matrix([ [1,-G,0,0] ,[-1,-G,G,0] ,[0,0,s*C2+(1/R3),-C2*s] ,[-1/R1,0,-C2*s,s*(C1+C2)+1/R1] ])
	b=Matrix([0,0,0,s*C1*Vi])
	V=A.inv()*b
	return (A,b,V)

A,b,V=lowpass(10000,10000,1e-11,1e-11,1.586,1)
Vo=V[3]
ww=p.logspace(0,8,801)
ss=1j*ww
hf=lambdify(s,Vo,'numpy')
v=hf(ss)
# Transfer Function Of Lowpass Filter
p.figure(0)
p.title("Transfer Function Of Lowpass Filter")
p.loglog(ww,abs(v),lw=2)
# Step Response of Lowpass Filter
p.figure(1)
p.title("Step Response of Lowpass Filter")
Y=sympy_to_lti(Vo)
t1=p.linspace(0,0.05,1000)
t,x=sp.step(Y,None,t1)
p.plot(x,t)
# Output of Lowpass Filter
p.figure(2)
w_1=2e3*p.pi
w_2=2e6*p.pi
Vin=(w_1/(s**2+w_1**2))+(s/(s**2+w_2**2))
A0,b0,V0=lowpass(10000,10000,1e-11,1e-11,1.586,Vin)
Vout=V0[3]
t2=p.linspace(0,10,1000)
# t2,x=sp.impulse(sympy_to_lti(Vout),None,t2)
u1=p.sin(w_1*t)+p.cos(w_2*t)
t,y,svec=sp.lsim(Y,u1,t2)
p.title('Output of Lowpass Filter')
p.plot(t,y)
#Transfer Function Of Highpass Filter
a1,b1,V1=highpass(10000,10000,1e-9,1e-9,1.586,1)
hf1=lambdify(s,V1[3],'numpy')
v1=hf1(ss)
p.figure(3)
p.title('Transfer Function Of Highpass Filter')
p.loglog(ww,abs(v1),lw=2)
#Output To A Damped Sinusoid
Vda=V1[3]
Vd=sympy_to_lti(Vda)
a=0.5
# t6=p.linspace(1,1.1,1000)
# u2=p.sin(w_2*p.pi*t6)*p.exp(-a*t6)
u=p.sin(w_2*p.pi*t)*p.exp(-a*t)
t=p.linspace(1.1,1.2,1000)
td,y,svec=sp.lsim(Vd,u,t)
p.figure(4)
p.title("Response To A Damped Sinusoid")
p.plot(td,y,label='Output')
# p.plot(t6,u2,label='Input',color='orange',marker='o')
p.legend(loc='upper right')
#Step Response Of Highpass Filter
a2,b2,V2=highpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vst=V2[3]
p.figure(5)
p.title("Step Response Of Highpass Filter")
t4=p.linspace(0,30.07,100000)
Y1=sympy_to_lti(Vst)
t,y=sp.impulse(Y1,None,t4)
p.plot(t,y)
p.show()
