from pylab import * 
from numpy import * 
import time 

w = 1 #plot size in diagonal 
subplot(111, aspect = 'equal') 
ion(); 

def euler( X, t, h, df, pars ): 
	return X + h*df( X, t, pars) 

def impEuler( X, t, h, df, pars ): 
	X1 = X + h*df( X, t, pars) 
	X2 = X1 + h*df( X1, t, pars) 
	return (X1+X2)/2. 

def rk4step(X, t, h, ratesofchange, pars):  # "Runge-Kutta 4th order" step 
	halfh = h/2. 
	slope1 = ratesofchange( X               , t      , pars ) 
	slope2 = ratesofchange( X + halfh*slope1, t+halfh, pars ) 
	slope3 = ratesofchange( X + halfh*slope2, t+halfh, pars ) 
	slope4 = ratesofchange( X +     h*slope3, t+h    , pars ) 
	slope  = ( slope1 + 2.*slope2 + 2.*slope3 + slope4 )/6. 
	return X + h*slope 
	 
def df( X, t, pars): #rates of change 
	x,y,u,v = X 
	m1, m2, r1, r2, G, a = pars #a is angular frequency of orbit around 2 suns 
	x1 = r1*cos(a*t);	y1 = r1*sin(a*t) 
	x2 = -r2*cos(a*t);	y2 = -r2*sin(a*t)		 
	k1 = (-G*m1)/((((x-x1)**2)+((y-y1)**2))**1.5) 
	k2 = (-G*m2)/((((x-x1)**2)+((y-y1)**2))**1.5) 
	du = k1*(x-x1)+k2*(x-x2);	dv = k1*(y-y1)+k2*(y-y2) 
	rates = array([ u, v, du, dv ]) 
	return rates	 

T = 2*pi # time we want to get to 
h = 0.01 # number of Euler steps to take/Time-step 
n = int(1010*ceil(T/h)) # number of steps to take 

m1 = 0.75 
m2 = 1.-m1 
r1 = m2 
r2 = m1 
G = 1. 
a = 1. # rotational velocity 
#pars = array([ m1, m2, r1, r2, G, a ]) # pack up our pars in an incredibly functional array suitcase 
p = array([ m1, m2, r1, r2, G, a ]) 
d = 7.125 
X = array([ d, 0., 0., -.25 ]) 

#s1, = plot( [r1], [0], 'ro', markersize = 10) 
#s2, = plot( [-r2], [0], 'go', markersize = 10) 
pl, = plot([X[0]], [X[1]], 'bo', markersize = 5) 

dp = d*2. 
plot( [-dp,dp], [-dp,dp], 'w') 


for i in range(n): 
	t = i*h 
	#print X, t, n 
#	s1.set_xdata(r1*cos(a*t));	s1.set_ydata(r1*sin(a*t)) 
#	s2.set_xdata(-r2*cos(a*t));	s2.set_ydata(-r2*sin(a*t))	 
	X = rk4step( X, t, h, df, p) 
	#X = impEuler(X, t, h, df, p)		 
	x, y = X[:2] 
	F = df( X, t, p) 
	xr = F[2];	yr = F[3]; 
	if i%1000 == 0: 
		print X, t, n 
		plot([x,x+xr], [y, y+yr], 'g-') 
		pl.set_xdata(X[0])	 
		pl.set_ydata(X[1])	 
		draw() 
ioff();		raw_input()
