import numpy as np
import matplotlib.pyplot as plt

# model /parameters
steps = 1000.0; bmin = 0.0; bmax = 1.40; alpha=0.5; delta=0.45; tau=10.0; T=10.0
b=np.linspace(bmin,bmax,steps+1); Eeq=(T-delta/alpha)/(1+tau*delta)
a1=alpha*Eeq+delta+1
a2=alpha*Eeq+delta+alpha*Eeq*delta-alpha*T+alpha*Eeq+alpha*Eeq*tau*delta
a3=-alpha*Eeq*delta; a4=a3; a5=-a4

# functions 
def calcFunc(iv,par1,par2,par3,par4,par5,par6): # function to be evaluated
   dv=((par1*par4-par3)*iv**3+(par2*par3-par4*par5)*iv)/(par4*(iv**4
-par2*iv**2)+par3*(par1*iv**2-par5))-np.tan(par6*iv)
   soln=np.column_stack((iv,dv))
   return soln

def checkFunc(par1,par2,par3,par4,par5):
   b=np.sqrt(0.5*(-(par1**2-2.0*par2)+np.sqrt((par1**2-2.0*par2)**2
-4*(par2**2-2*par1*par5-par4**2))))
   return b

fig1 = plt.figure() # plot population and carrying capacity
# create data to plot
answer=calcFunc(b,a1,a2,a3,a4,a5,tau)
# specify axis labels, size and location
ax = plt.gca()
label = ax.set_xlabel('b', fontsize = 20)
ax.xaxis.set_label_coords(0.5, -0.06)
label = ax.set_ylabel(r'''$G(b,\tau)$''', fontsize = 20)
ax.yaxis.set_label_coords(-0.05, 0.5)
# plot data
plt.plot(answer[:,0],answer[:,1],'-g', linewidth=1)
# define plot boundaries & specify legend size/location
plt.xlim(bmin, bmax); plt.ylim(-10, 10)
#plt.xticks(range(0)," ")
#plt.yticks(range(0)," ")
# save as pdf and display plot
plt.savefig("bifurcationPlot.pdf")
plt.show()
