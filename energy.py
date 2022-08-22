from cProfile import label
import string
import numpy as np
import math
from numpy.linalg import inv
from numpy.linalg.linalg import norm
from sys import argv
from matplotlib import pyplot as plt
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial
import matplotlib as mpl
from matplotlib import cm
mpl.rcParams['font.family'] = 'serif'
#plt.rcParams['font.size'] = 18
#plt.rcParams['axes.linewidth'] = 2
plt.rcParams['text.usetex'] = True

from cycler import cycler
line_cycler   = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                 cycler(marker=["4", "2", "3", "1", "+", "x", "."]))
plt.rc("axes", prop_cycle=line_cycler)

def cheb(n,m):
    wyn=[]
    for k in range(0,n+1):
        wyn.append(math.cos((2*k-1)/(4*n)*math.pi)*m)
    wyn.reverse()
    return wyn

class Cords:
    def __init__(self,angles,r):    
        self.angles=angles  
        self.rtab=r
        self.coefs=None #interpolating polynomial coefs
        self.x=[]   #x coordinates of points
        self.y=[]   #y coordinates of points

    def convert(self):  
        temp=con(self.rtab[0],self.angles[0]) #convert points from polar coordinates to cartesian ones
        self.x.append(-temp[0])
        self.y.append(temp[1])
        for i in range(1,len(self.rtab)):
            temp=con(self.rtab[i],self.angles[i])
            self.x.append(temp[0])
            self.x.append(-temp[0])
            self.y.append(temp[1])
            self.y.append(temp[1])
        te=[]
        n=2*len(self.angles)-1
        for i in range(0,n): #create van der mond matrix used to coefs calculation
            row=[]
            for j in range(0,n):
                row.append(self.x[i]**j)
            te.append(row)
        te2=[]
        for i in range(0,n):
            te2.append([self.y[i]])
        #print(np.array(te))
        temp=np.dot(inv(np.array(te)),np.array(te2)) #temp is vector containing coefs of polynomial 
        self.coefs=[]
        for i in temp:
            self.coefs.append(i[0])

    def convert2(self): #same as convert but uses lagrange interpolation from scipy
        temp=con(self.rtab[0],self.angles[0]) #convert points from polar coordinates to cartesian ones
        self.x.append(-temp[0])
        self.y.append(temp[1])
        for i in range(1,len(self.rtab)):
            temp=con(self.rtab[i],self.angles[i])
            self.x.append(temp[0])
            self.x.append(-temp[0])
            self.y.append(temp[1])
            self.y.append(temp[1])
        
        self.coefs=Polynomial(lagrange(self.x,self.y)).coef

    def normal(self,pos):   #calculating normal vector
        temp=0
        for a in range(1,len(self.coefs)):
            temp+=pos[0]**(a-1)*self.coefs[a]*a
        #print(temp)
        return np.array([-temp,1])

    def fun(self,x):    #calculating position y coordinate of surface based on x coordinate
        wyn=0
        for a in range(0,len(self.coefs)):
            wyn+=(x**a)*self.coefs[a]
        return np.array([x,wyn])
    
    def fun2(self,x): #same as fun but output is on bottom side of strar
        wyn=0
        for a in range(0,len(self.coefs)):
            wyn+=x**a*self.coefs[a]
        return np.array([x,-wyn])

    def normal2(self,pos):  #normal vector on the second side of star
        temp=0
        for a in range(1,len(self.coefs)):
            temp+=pos[0]**(a-1)*self.coefs[a]*a
        return np.array([-temp,-1])

    def len(self,pos):  #prepared but not needed 
        temp=self.normal(pos)
        return math.sqrt(pos[0]**2+pos[1]**2)

    def integrate(self,n,mag):  #integrate using evenly spaced grid
        X=np.linspace(0,max(self.x),n)
        dx=max(self.x)/n
        wyn=0
        for i in range(0,n):    #integrating upper part of star
            x=X[i]
            pos=self.fun(x)
            wyn+=pot(pos,mag)*np.vdot(self.normal(pos),mag_field(pos,mag))*x*dx
        for i in range(0,n):    #integrating bottom part of strar
            x=X[i]
            pos=self.fun2(x)    
            wyn+=pot(pos,mag)*np.vdot(self.normal2(pos),mag_field(pos,mag))*x*dx
        return  math.pi*wyn
    
    def integrate_cheb(self,n,mag): #integrating using chebishev nodes
        X=cheb(n,max(self.x))
        wyn=0
        for i in range(0,n):    #upper side
            x=X[i]
            if i==0:
                dx=(X[1]+X[2])/2
            elif i==n-1:
                dx=max(self.x)-(X[-1]+X[-2])/2
            else:
                dx=(X[i+1]-X[i-1])/2
            pos=self.fun(x)
            wyn+=pot(pos,mag)*np.vdot(self.normal(pos),mag_field(pos,mag))*x*dx
        for i in range(0,n):    #bottom side
            x=X[i]
            if i==0:
                dx=(X[1]+X[2])/2
            elif i==n-1:
                dx=max(self.x)-(X[-1]+X[-2])/2
            else:
                dx=(X[i+1]-X[i-1])/2
            pos=self.fun2(x)
            wyn+=pot(pos,mag)*np.vdot(self.normal2(pos),mag_field(pos,mag))*x*dx
        return  math.pi*wyn


    def show(self,color): #show star and plot normal vector
        x=[*self.x[::-2],*self.x[1::2],*((self.x[1::2])[::-1]),*self.x[::2]]
        y=[*self.y[::-2],*self.y[1::2],*([-i for i in self.y[1::2]])[::-1],*[-i for i in self.y[::2]]]
        x=[*x,x[0]]
        y=[*y,y[0]]
        plt.plot(x,y,c=color)

def pot(pos,mag):   #magnetic potential
    r=math.sqrt(pos[0]**2+pos[1]**2)
    return mag*pos[1]/(4*math.pi*r**3)

def mag_field(pos,mag): #magnetic induction vector (in 2D cartesian)
    r=math.sqrt(pos[0]**2+pos[1]**2)
    return 1/(math.pi*4*r**5)*np.array([3*pos[1]*mag*pos[0],3*pos[1]*mag*pos[1]-mag*r**2])

def con(r,phi): #convert polar cords to cartesian
    return np.array([r*math.sin(phi),r*math.cos(phi)])

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
def floating(X):
    wyn=[]
    for x in X:
        wyn.append(float(x))
    return wyn

def get(name):
    wyn=[]
    with open(name) as f:
        for t in f:
            line=t.split()
            wyn.append(floating(line))
    return wyn

def app(file1,file2,tab):
    ina=open(file1,"r")
    out=open(file2,"w")
    temp=ina.readline()
    out.write(temp)
    temp=ina.readline()
    out.write(temp[:-1]+" External energy\n")
    for i in tab:
        temp=ina.readline()
        out.write(temp[:-1]+" {}\n".format(i))
    ina.close()
    out.close()
def printujemy(x):
    wyn=0
    for a in x:
        wyn+=a*a
    return wyn

if __name__=="__main__":
    w=get("psurface.txt")
    wyn=[]
    k=1
    color=[]
    p_arr=[]
    for c in w:
        mm=c[0]
        color.append(c[1])
        a=c[2:]
        print(printujemy(a),min(a),max(a))
        X=np.linspace(0,math.pi,len(a))
        #print(type(X))
        p=Cords(X[:len(a)//2:k],a[:len(a)//2:k])
        #print(a,X)
        p.convert()
        wyn.append(1.25663706*10**(-6)*(10**32*mm)**2*p.integrate_cheb(1000,1)/1000**3)
        p_arr.append(p)
    print(len(wyn))
    print(wyn)
    app("results.txt",argv[1],wyn)