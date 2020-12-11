import matplotlib.pyplot as plt
import numpy as np
import math as m
import time
from scipy.stats import linregress

#Integrateurs :

def euler(W0,pot,n,dt,F):
	"""W0 is the initial value vector of structure [x,y,vx,vy]
		n is the number of iterations of the time step
		dt is the time step"""
	W=W0
	iter=n
	for i in range(iter):
			np.savetxt(F,W)
			Wp=pot(W)
			W[0]= W[0]+Wp[0]*dt
			W[1]= W[1]+Wp[1]*dt
			W[2]= W[2]+Wp[2]*dt
			W[3]= W[3]+Wp[3]*dt
	return()

def RK2(W0,pot,n,dt,F):
	"""W0 is the initial value vector of structure [x,y,vx,vy]
		n is the number of iterations of the time step
		dt is the time step"""
	W=W0
	W12=np.array([0.,0.,0.,0.])
	iter=n
	for i in range(iter):
		np.savetxt(F,W)
		Wp=pot(W)
		W12=W+0.5*dt*Wp
		Wp12=pot(W12)
		W=W+dt*Wp12
	return()

def RK4(W0,pot,n,dt,F):
	"""W0 is the initial value vector of structure [x,y,vx,vy]
		n is the number of iterations of the time step
		dt is the time step"""
	W=W0
	iter=n
	for i in range(iter):
		np.savetxt(F,W)
		k1=pot(W)
		Wk2=W+0.5*dt*k1
		k2=pot(Wk2)
		Wk3=W+0.5*dt*k2
		k3=pot(Wk3)
		Wk4=W+dt*k3
		k4=pot(Wk4)
		W=W+dt*(k1+2*k2+2*k3+k4)/6
	return()

#Potentiels :

def Kepler(W):
	'''W is the vector [x,y,vx,vy]
		this function returns the time derivative of W
	'''
	Wpoint=np.array([0.,0.,0.,0.])
	tmp0=W[0]
	tmp1=W[1]
	Wpoint[0]=W[2]
	Wpoint[1]=W[3]
	Wpoint[2]=-tmp0/(m.pow(tmp0**2+tmp1**2,1.5))
	Wpoint[3]=-tmp1/(m.pow(tmp0**2+tmp1**2,1.5))
	return(Wpoint)

def HenonHeiles(W):
	'''W is the vector [x,y,vx,vy]
		this function returns [vx,vy,ax,ay]
		with a=F=-grad(V)
	'''
	Wpoint=np.array([0.,0.,0.,0.])
	tmp0=W[0]
	tmp1=W[1]
	Wpoint[0]=W[2]
	Wpoint[1]=W[3]
	Wpoint[2]=-tmp0*(1+2*tmp1)
	Wpoint[3]=-tmp0*tmp0-tmp1*(1-tmp1)
	return(Wpoint)

# Main, affichage :

def integre(func,pot,W0,n,dt,nb_orbit):
    '''func is the integrator (euler, RK2 or RK4)
    pot is the potential (Kepler, HenonHeiles)
    W0 is the [x,y,vx,vy] initial vector
    n*dt is the integration time with dt the timestep
    '''
    F=open(".\CoordVector.txt", "w")
    for j in range(nb_orbit):
        func(W0[j],pot,n,dt,F)
    F.close()
    return()


def affichage(filename,orbit_nb):
    '''Display the trajectory Y=f(X) based on the [x,y,vx,vy] vectors saved in the file "filename"
        Possibility of displaying total energy as a function of time'''
    f=open(filename,"r")
    lu=f.readline
    W=np.loadtxt(filename)
    print("length of W =",len(W))
    X,Y,Vx,Vy,E_kinetic,E_kepler,V_henon_heiles,T=np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb)),np.zeros((orbit_nb,len(W)//orbit_nb))
    i=0
    for j in range(orbit_nb):
        while(i*4+3<=len(W)/nb_orbit):
            X[j][i]=W[(j+i)*4]
            Y[j][i]=W[(j+i)*4+1]
            Vx[j][i]=W[(j+i)*4+2]
            Vy[j][i]=W[(j+i)*4+3]
            #E_kinetic[j][i]=(Vx[i]**2+Vy[i]**2)/2
            #E_kepler[j][i]=-1/m.sqrt(X[i]**2+Y[i]**2
            #V_henon_heiles[j][i]=0.5*(X[i]**2+Y[i]**2+2*X[i]**2*Y[i]-2*Y[i]**3/2
            #T[j][i]=i*len(W)/1000
            #np.append(X[j],W[i*4])
            #np.append(Y[j],[i*4+1])
            #np.append(Vx[j],W[i*4+2])
            #np.append(Vy[j],W[i*4+3])
            #np.append(E_kinetic[j],(Vx[i]**2+Vy[i]**2)/2)
            #np.append(E_kepler[j],-1/m.sqrt(X[i]**2+Y[i]**2)
            #np.append(V_henon_heiles[j],0.5*(X[i]**2+Y[i]**2+2*X[i]**2*Y[i]-2*Y[i]**3/2)
            #np.append(T[j],i*len(W)/1000)
            i+=1
        #plt.subplot(2,2,1)
        plt.plot(X[j],Y[j],'r')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Trajectoire y=f(x)')
        #plt.subplot(2,2,3)
        #plt.plot(T,E_kinetic+E_kepler) #or V_henon_heiles
        #plt.xlabel('Temps')
        #plt.ylabel('Energie')
    t_affichage = time.time()
    plt.show()
    return(t_affichage)

# def E_dt(func,pot,W0,n,DT):
# 	E=[]
# 	for j in range(len(DT)):
# 		integre(func,pot,W0,n,DT[j],1)
# 		f=open("D:\\Documents\\TPS\\3A\\Simulations numeriques\\CoordVector.txt","r")
# 		lu=f.readline
# 		W=np.loadtxt("D:\\Documents\\TPS\\3A\\Simulations numeriques\\CoordVector.txt")
# 		X,Y,Vx,Vy=[],[],[],[]
# 		i=0
# 		while(i*4+3<=len(W)):
# 			X.append(W[i*4])
# 			Y.append(W[i*4+1])
# 			Vx.append(W[i*4+2])
# 			Vy.append(W[i*4+3])
# 			i+=1
# 		E.append((Vx[-1]**2+Vy[-1]**2)/2-1/m.sqrt(X[-1]**2+Y[-1]**2))
# 	return(E)
#
# def affichageE_dt(pot,W0,n,DT):
# 	E_euler=E_dt(euler,pot,W0,n,DT)
# 	E_RK2=E_dt(RK2,pot,W0,n,DT)
# 	E_RK4=E_dt(RK4,pot,W0,n,DT)
# 	E_euler=abs(E_euler)
# 	E_euler_log=np.log10(E_euler)
# 	E_RK2=abs(E_RK2)
# 	E_RK2_log=np.log10(E_RK2)
# 	E_RK4=abs(E_RK2)
# 	E_RK4_log=np.log10(Sigma_E_RK4[i])
# 	E=[E_euler_log,E_RK2_log,E_RK4_log]
# 	plt.plot(np.log10(DT),E[0],'b',np.log10(DT),E[10],'r',np.log10(DT),E[2],'g')
# 	plt.xlabel('log10(DT)')
# 	plt.ylabel('log10(E)')
# 	plt.show()
# 	return()


def PoincareSection(filename):
    '''Display the Poincare section Vy=f(y)
        based on the [x,y,vx,vy] vectors saved in the file "filename"
        in the x=0 plane, vx>0
    '''
    f=open(filename,"r")
    lu=f.readline
    W=np.loadtxt(filename)
    print("length of W =",len(W))
    i=0
    Y,Vy=[],[]
    while((i+1)*4+3<=len(W)):
        if(W[i*4]*W[(i+1)*4]<=0): # vx > 0 and we cross x=0 (sign changes) #W[i*4+2]>0 and
			#Linear interpolation
            y=W[i*4+1]-W[i*4]*(W[(i+1)*4+1]-W[i*4+1])/(W[(i+1)*4]-W[i*4])
            # y(x=0) = y_a - x_a * (y_b - y_a) / (x_b - x_a)
            vy=W[i*4+3]-W[i*4]*(W[(i+1)*4+3]-W[i*4+3])/(W[(i+1)*4]-W[i*4])
			# vy(x=0) = vy_a - x_a * (vy_b - vy_a) / (x_b - x_a)
            Y.append(y)
            Vy.append(vy)
            i+=1
        else:
            i+=1
    plt.plot(Y,Vy,'r.')
    plt.xlabel('Y')
    plt.ylabel('Vy')
    plt.title('Poincare section Vy=f(y)')
    t_affichage = time.time()
    plt.show()
    return([Y,Vy])


def go(nb_orbit):
    start = time.time()
    W0=np.array([1.,0.,0.,1.])
    integre(RK2,Kepler,W0,4000,0.01,nb_orbit)
    end_integre = time.time()
    end_affichage=affichage(".\CoordVector.txt",nb_orbit)
    # DT=np.array([0.0001,0.001,0.01,0.1])
    # affichageE_dt(Kepler,W0,100000,DT)
    print("Integration time in s = ", end_integre - start)
    print("Displaying time in s = ", end_affichage - end_integre)
    return(0)

def goHenonHeiles(E,nb_orbit):
    start = time.time()
    xlim=m.sqrt(2*E) #xlim is x / E = V(x,0)
    x0=[np.random.uniform(0,xlim) for _ in range(nb_orbit)] #we choose a random initial x in the authorized range
    Ec=[E-0.5*x0[i]**2 for i in range(len(x0))]
    vy=[m.sqrt(2*Ec[i]) for i in range(len(Ec))]
    W0=[np.array([x0[i],0.,0.,vy[i]]) for i in range(len(x0))]
    integre(RK4,HenonHeiles,W0,50000,0.1,nb_orbit)
    end_integre = time.time()
	#end_affichage=affichage("D:\\Documents\\TPS\\3A\\Simulations numeriques\\CoordVector.txt")
    end_affichage=PoincareSection(".\CoordVector.txt")
    print("Integration time in s = ", end_integre - start)
	#print("Displaying time in s = ", end_affichage - end_integre) #change return value of PoincareSection by t_affichage


def orbite_reguliere_ou_pas(E):
    '''E total energy
    This function simulates one orbit and returns 1 if it is regular
    or 0 if it is chaotic.
    '''
    #start = time.time()
    #mu=0
    Delta_W=[]
    xlim=m.sqrt(2*E) #xlim is x / E = V(x,0)
    x0_1=np.random.uniform(0,xlim) #we choose a random initial x in the authorized range
    x0=[x0_1,x0_1-1e-7]
    Ec=[E-0.5*x0[0]**2,E-0.5*x0[1]**2]
    vy=[m.sqrt(2*Ec[0]),m.sqrt(2*Ec[1])]
    w01=[x0[0],0.,0.,vy[0]]
    w02=[x0[1],0.,0.,vy[1]]
    W0=np.array([w01,w02])
	#print("Conditions initiales =",W0[0])
    integre(RK4,HenonHeiles,W0,100000,0.01,2) #adjust parameters ?
    W=np.loadtxt("./CoordVector.txt")
    i=0
    X,Y,Vx,Vy=[],[],[],[]
    W1,W2=[],[]
    while(i<len(W)//2):
        W1.append(W[i])
        W2.append(W[len(W)//2+i])
        i+=1
    i=0
    while(i*4+3<=len(W1)):
        Delta_W.append((W1[i*4]-W2[i*4])**2+(W1[i*4+1]-W2[i*4+1])**2+(W1[i*4+2]-W2[i*4+2])**2+(W1[i*4+3]-W2[i*4+3])**2) #distance dans l'espace des phases
		#mu+=Delta_W[i]
        i+=1
    lnDelta_W=[np.log(Delta_W[i]) for i in range(len(Delta_W))]
    x=[i for i in range(len(lnDelta_W))]
    a,b,r,p_value, std_err=linregress(x,lnDelta_W)
    if(r**2>0.95): #a verifier, suffisant ou pas ? lineaire si Delta_W exponentiel donc si chaotique
        res=0
    else:
        res=1
    end_integre = time.time()
    #print("Integration time in s = ", end_integre - start)
    #print("mu = ", mu)
    return(res)


if __name__=="__main__":
    start = time.time()
    nb_orbit=1
    E=1/100
    #goHenonHeiles(E,nb_orbit)
    E=[1/100,1/12,1/10,1/8,1/6]
    aire=[]
    nb_orbits=2
    avancement=0
    for k in range(len(E)):
        A=0
        for i in range(nb_orbits):
            A+=orbite_reguliere_ou_pas(E[k])
            current_orbit = k*nb_orbits+i
            avancement = current_orbit*100/(nb_orbits*len(E))
            print("Avancement à {} %".format(avancement))
        aire.append(A)
    fig7=[aire[i]/aire[0] for i in range(len(aire))]
    end_integre = time.time()
    print("Integration time in s = ", end_integre - start)
    plt.scatter(E,fig7)
    plt.xlabel("E")
    plt.ylabel("Ratio of chaotic orbits over regular orbits (normalized w.r.t E=1/100)")
    plt.show()

