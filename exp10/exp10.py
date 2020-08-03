from pylab import *
from scipy import signal
import csv
rcParams['axes.grid']=True
#Importing Data From h.csv File.
#Problem 1
from csv import reader
def file(filename):
    l=[]
    with open(filename) as inp:
        for row in inp:
            l.append(float(row))
    return l
l=file('h.csv');l=array(l)
#Problem 2
w,h=signal.freqz(l)
figure()
subplot(2,1,1)
title('Magnitude and Phase Response of Filter')
ylabel('|Y| \u2192')
plot(w,20*log10(abs(h)),color='green')
subplot(2,1,2)
ylabel('Phase of Y')
plot(w,angle(h),'ro',markersize=1)
xlabel('w \u2192')
#Problem 3
N=2**10
n=arange(1,N,1)
x=cos(0.2*pi*n)+cos(0.85*pi*n)
figure()
title('Cosine Sequence : $cos(0.2\pi t) + cos(0.85\pi t)$')
xlabel('n \u2192')
ylabel('x \u2192')
plot(n,x,color='red')
xlim([1,100])
#Linear Convolution
#Problem 4
y=signal.convolve(x,l)
figure()
title('Linearly Convolved Cosine Sequence')
xlabel('n \u2192')
ylabel('y \u2192')
plot(range(len(n)+len(l)-1),y,color='darkviolet')
xlim([1,100])
#Using the DFTs
#Problem 5
x_=concatenate((x,zeros(len(l)-1)))
l_=concatenate((l,zeros(len(x)-1)))
y1=ifft(fft(x_)*fft(l_))
figure()
title('Convolvution of Cosine Sequence using DFT')
xlabel('n \u2192')
ylabel('y1 \u2192')
plot(range(len(y1)),real(y1),color='hotpink')
xlim([1,100])
#Circular Convolution
#Problem 6
p=len(l)
n_=int(ceil(log2(p)))
l_=concatenate((l,zeros(2**n_-p)))
p=len(l_)
n1=int(ceil(len(x)/2**n_))
x_=concatenate((x,zeros(n1*int(2**n_)-len(x))))
y2=zeros(len(x_)+len(l_)-1)
for i in range(n1):
    temp=concatenate((x_[i*p:(i+1)*p],zeros(p-1)))
    y2[i*p:(i+1)*p+p-1]+=real(ifft(fft(temp)*fft(concatenate((l_,zeros(len(temp)-len(l_)))))))
figure()
title('Circular Convolvution of Cosine Sequence')
xlabel('n \u2192')
ylabel('y2 \u2192')
plot(range(len(y2)),y2,color='royalblue')
xlim([0,100])
#Zadoff-Chu Sequence
#Problem 7
rows=[]
with open('x1.csv') as zadf:
    csvreader = csv.reader(zadf)
    for row in csvreader:
        rows.append(row)
rows2=[]
for row in rows:
    row = list(row[0])
    try:
        row[row.index('i')]='j'
        rows2.append(row)
    except ValueError:
        rows2.append(row)
        continue
x1=[complex(''.join(line)) for line in rows2]
x2=roll(x1,5)
zer=ifftshift(correlate(x2,x1,'full'))
figure()
title('Correlation of Zadoff-Chu Sequence')
xlabel('n \u2192')
ylabel('Correlation \u2192')
stem(arange(0,len(zer),1),abs(zer),'g',use_line_collection=True)
xlim([1,20])
show()
