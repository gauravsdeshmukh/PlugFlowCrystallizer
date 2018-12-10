# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 23:20:56 2018

@author: Gaurav
"""

import scipy as sci
import numba as nb
import pandas as pd
import os
############################CLASS SPACE########################################
class Boundary:
    def __init__(self,boundary_type,boundary_value,constant=0):
        self.DefineBoundary(boundary_type,boundary_value,constant)
        
    def DefineBoundary(self,boundary_type,boundary_value,constant):
        self.type=boundary_type
        self.value=boundary_value
        self.constant=constant

class Space:
    def __init__(self):
        pass
    
    def CreateMesh(self,rowpts,colpts):
        self.rowpts=rowpts
        self.colpts=colpts
        self.u=sci.zeros((self.rowpts+2,self.colpts+2))
        self.v=sci.zeros((self.rowpts+2,self.colpts+2))
        self.u_star=sci.zeros((self.rowpts+2,self.colpts+2))
        self.v_star=sci.zeros((self.rowpts+2,self.colpts+2))
        self.p=sci.zeros((self.rowpts+2,self.colpts+2))
        self.T=sci.zeros((self.rowpts+2,self.colpts+2))
        self.C=sci.zeros((self.rowpts+2,self.colpts+2))
        self.SetEmptyMesh()
        
    def SetDeltas(self,breadth,length):
        self.dx=length/self.colpts
        self.dy=breadth/self.rowpts

    def SetEmptyMesh(self):
        self.u_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.v_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.T_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.C_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.p_c=sci.zeros((self.rowpts,self.colpts))
        self.u_c=sci.zeros((self.rowpts,self.colpts))
        self.v_c=sci.zeros((self.rowpts,self.colpts))
        self.T_c=sci.zeros((self.rowpts,self.colpts))
        self.C_c=sci.zeros((self.rowpts,self.colpts))

class Fluid:
    def __init__(self,rho,mu,k,cp,D):
        self.SetFluidProperties(rho,mu)
        self.SetThermalProperties(k,cp)
        self.SetTransportProperties(D)
    
    def SetFluidProperties(self,rho,mu):
        self.rho=rho
        self.mu=mu
        
    def SetThermalProperties(self,k,cp):
        self.k=k
        self.cp=cp
        self.alpha=self.k/(self.rho*self.cp)
        
    def SetTransportProperties(self,D):
        self.D=D
        
        
class Crystal:
    def __init__(self,name,space):
        self.name=name
        self.space=space
        self.rowpts=self.space.rowpts
        self.colpts=self.space.colpts
        self.GetSolubility()
        self.SetMeshes()
        self.SetGeometricParameters()
        
    def GetSolubility(self):
        cwdir=os.getcwd()
        if(os.path.basename(os.getcwd())=="Result"):
            os.chdir("..")
        filename="SolubilityData.xlsx"
        path=cwdir+"\\"+filename
        sol=pd.read_excel(path)
        
        row=sol[sol["Name"]==self.name]
        self.a0=row["a0"].values[0]
        self.a1=row["a1"].values[0]
        self.a2=row["a2"].values[0]
        
    def SetConstants(self,MW,kg,g,rho,dH):
        self.MW=MW
        self.kg=kg
        self.g=g
        self.rho=rho
        self.dH=dH
        
    def SetGeometricParameters(self,K_a=sci.pi,K_v=sci.pi/6):
        self.K_a=K_a
        self.K_v=K_v
        
    def SetMeshes(self):
        self.N=sci.zeros((self.rowpts+2,self.colpts+2))
        self.L=sci.zeros((self.rowpts+2,self.colpts+2))
        self.A=sci.zeros((self.rowpts+2,self.colpts+2))
        self.V=sci.zeros((self.rowpts+2,self.colpts+2))
        self.G=sci.zeros((self.rowpts+2,self.colpts+2))
        self.Cstar=sci.zeros((self.rowpts+2,self.colpts+2))
        self.N_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.L_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.A_next=sci.zeros((self.rowpts+2,self.colpts+2))
        self.V_next=sci.zeros((self.rowpts+2,self.colpts+2))
    
##########################BOUNDARY SPACE#######################################

def SetUBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.u[:,0]=left.value
    elif(left.type=="N"):
        space.u[:,0]=-left.value*space.dx+space.u[:,1]
    
    if(right.type=="D"):
        space.u[:,-1]=right.value
    elif(right.type=="N"):
        space.u[:,-1]=right.value*space.dx+space.u[:,-2]
        
    if(top.type=="D"):
        space.u[-1,:]=2*top.value-space.u[-2,:]
    elif(top.type=="N"):
        space.u[-1,:]=top.value*space.dy+space.u[-2,:]
     
    if(bottom.type=="D"):
        space.u[0,:]=2*bottom.value-space.u[1,:]
    elif(bottom.type=="N"):
        space.u[0,:]=-bottom.value*space.dy+space.u[1,:]
        

def SetVBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.v[:,0]=2*left.value-space.v[:,1]
    elif(left.type=="N"):
        space.v[:,0]=-left.value*space.dx+space.v[:,1]
    
    if(right.type=="D"):
        space.v[:,-1]=2*right.value-space.v[:,-2]
    elif(right.type=="N"):
        space.v[:,-1]=right.value*space.dx+space.v[:,-2]
        
    if(top.type=="D"):
        space.v[-1,:]=top.value
    elif(top.type=="N"):
        space.v[-1,:]=top.value*space.dy+space.v[-2,:]
     
    if(bottom.type=="D"):
        space.v[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.v[0,:]=-bottom.value*space.dy+space.v[1,:]
    
def SetPBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.p[:,0]=left.value
    elif(left.type=="N"):
        space.p[:,0]=-left.value*space.dx+space.p[:,1]
    
    if(right.type=="D"):
        space.p[1,-1]=right.value
    elif(right.type=="N"):
        space.p[:,-1]=right.value*space.dx+space.p[:,-2]
        
    if(top.type=="D"):
        space.p[-1,:]=top.value
    elif(top.type=="N"):
        space.p[-1,:]=top.value*space.dy+space.p[-2,:]
     
    if(bottom.type=="D"):
        space.p[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.p[0,:]=-bottom.value*space.dy+space.p[1,:]
    
 
def SetTBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.T[:,0]=left.value
    elif(left.type=="N"):
        space.T[:,0]=-left.value*space.dx+space.T[:,1]
    elif(left.type=="C"):
        space.T[:,0]=left.value*space.T[:,1]+left.constant
    
    if(right.type=="D"):
        space.T[1,-1]=right.value
    elif(right.type=="N"):
        space.T[:,-1]=right.value*space.dx+space.T[:,-2]
    elif(right.type=="C"):
        space.T[:,-1]=right.value*space.T[:,-2]-right.constant
        
    if(top.type=="D"):
        space.T[-1,:]=top.value
    elif(top.type=="N"):
        space.T[-1,:]=top.value*space.dy+space.T[-2,:]
    elif(top.type=="C"):
        space.T[-1,:]=top.value*space.T[-2,:]-top.constant
     
    if(bottom.type=="D"):
        space.T[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.T[0,:]=-bottom.value*space.dy+space.T[1,:]
    elif(bottom.type=="C"):
        space.T[0,:]=bottom.value*space.T[1,:]+bottom.constant
        
def SetCBoundary(space,left,right,top,bottom):
    if(left.type=="D"):
        space.C[:,0]=left.value
    elif(left.type=="N"):
        space.C[:,0]=-left.value*space.dx+space.C[:,1]
    elif(left.type=="C"):
        space.C[:,0]=left.value*space.C[:,1]+left.constant
    
    if(right.type=="D"):
        space.C[1,-1]=right.value
    elif(right.type=="N"):
        space.C[:,-1]=right.value*space.dx+space.C[:,-2]
    elif(right.type=="C"):
        space.C[:,-1]=right.value*space.C[:,-2]-right.constant
        
    if(top.type=="D"):
        space.C[-1,:]=top.value
    elif(top.type=="N"):
        space.C[-1,:]=top.value*space.dy+space.C[-2,:]
    elif(top.type=="C"):
        space.C[-1,:]=top.value*space.C[-2,:]-top.constant
     
    if(bottom.type=="D"):
        space.C[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.C[0,:]=-bottom.value*space.dy+space.C[1,:]
    elif(bottom.type=="C"):
        space.C[0,:]=bottom.value*space.C[1,:]+bottom.constant 
        
########################FUNCTION SPACE#########################################
def SetTimeStep(CFL,space,fluid):
    if(space.u.any()!=0):
        dt_hyper=CFL/max(sci.amax(space.u)/space.dx,sci.amax(space.v)/space.dy)
    else:
        dt_hyper=CFL*space.dx 
        
    dt_para=min(space.dx**2/(2*fluid.mu),space.dy**2/(2*fluid.mu))
    dt_temp=min(space.dx**2/(2*fluid.alpha),space.dy**2/(2*fluid.alpha))
    dt_conc=min(space.dx**2/(2*fluid.D),space.dy**2/(2*fluid.D))
    dt_min=min(dt_hyper,dt_para,dt_temp,dt_conc)
    space.dt=dt_min
 
def GetStarredVelocities(space,fluid):
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u=space.u.astype(float)
    v=space.v.astype(float)
    u_star=space.u_star.astype(float)
    v_star=space.v_star.astype(float)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    mu=float(fluid.mu)
    
    u_star=u.copy()
    v_star=v.copy()
    
    u1_y=(u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy)
    u1_x=(u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx)
    u2_y=(u[2:,1:cols+1]-2*u[1:rows+1,1:cols+1]+u[0:rows,1:cols+1])/(dy**2)
    u2_x=(u[1:rows+1,2:]-2*u[1:rows+1,1:cols+1]+u[1:rows+1,0:cols])/(dx**2)
    v_face=(v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols]+v[2:,1:cols+1]+v[2:,0:cols])/4
    u_star[1:rows+1,1:cols+1]=u[1:rows+1,1:cols+1]-dt*(u[1:rows+1,1:cols+1]*u1_x+v_face*u1_y)+(dt*(mu/rho)*(u2_x+u2_y))   

    v1_y=(v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy)
    v1_x=(v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx)
    v2_y=(v[2:,1:cols+1]-2*v[1:rows+1,1:cols+1]+v[0:rows,1:cols+1])/(dy**2)
    v2_x=(v[1:rows+1,2:]-2*v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols])/(dx**2)
    u_face=(u[1:rows+1,1:cols+1]+u[1:rows+1,2:]+u[0:rows,1:cols+1]+u[0:rows,2:])/4
    v_star[1:rows+1,1:cols+1]=v[1:rows+1,1:cols+1]-dt*(u_face*v1_x+v[1:rows+1,1:cols+1]*v1_y)+(dt*(mu/rho)*(v2_x+v2_y))
    
    space.u_star=u_star.copy()
    space.v_star=v_star.copy()   
        
@nb.jit    
def SolvePressurePoisson(space,fluid,left,right,top,bottom):
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u_star=space.u_star.astype(float)
    v_star=space.v_star.astype(float)
    p=space.p.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    factor=1/(2/dx**2+2/dy**2)
    
    error=1
    tol=1e-3

    i=0
    while(error>tol):
        i+=1
        p_old=p.astype(float,copy=True)
        
#        term_1=(1/dt)*(((v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy))+((u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx)))
#        term_2=((u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx))**2
#        term_3=((v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy))**2
#        term_4=2*((u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy))*((v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx))
        p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2+(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2
        ustar1_x=(u_star[1:rows+1,2:]-u_star[1:rows+1,0:cols])/(2*dx)
        vstar1_y=(v_star[2:,1:cols+1]-v_star[0:rows,1:cols+1])/(2*dy)
        p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y)
        error=sci.amax(abs(p-p_old))
        #Apply Boundary Conditions
        SetPBoundary(space,left,right,top,bottom)
        
        if(i>500):
            tol*=10
            
    
@nb.jit
def SolveMomentumEquation(space,fluid):
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u_star=space.u_star.astype(float)
    v_star=space.v_star.astype(float)
    p=space.p.astype(float)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    u_next=space.u_next.astype(float,copy=False)
    v_next=space.v_next.astype(float,copy=False)

    p1_x=(p[1:rows+1,2:]-p[1:rows+1,0:cols])/(2*dx)
    u_next[1:rows+1,1:cols+1]=u_star[1:rows+1,1:cols+1]-(dt/rho)*p1_x

    p1_y=(p[2:,1:cols+1]-p[0:rows,1:cols+1])/(2*dy)
    v_next[1:rows+1,1:cols+1]=v_star[1:rows+1,1:cols+1]-(dt/rho)*p1_y
            
def AdjustUV(space):
    space.u[1:-1,1:-1]=space.u_next[1:-1,1:-1].copy()
    space.v[1:-1,1:-1]=space.v_next[1:-1,1:-1].copy()
    
def SetCentrePUVTCNLAV(space,crystal):
    space.p_c=space.p[1:-1,1:-1]
    space.u_c=space.u[1:-1,1:-1]
    space.v_c=space.v[1:-1,1:-1]
    space.T_c=space.T[1:-1,1:-1]
    space.C_c=space.C[1:-1,1:-1]
    crystal.N_c=crystal.N[1:-1,1:-1]
    crystal.L_c=crystal.L[1:-1,1:-1]
    crystal.A_c=crystal.A[1:-1,1:-1]
    crystal.V_c=crystal.V[1:-1,1:-1]
   
def SetInitialT(space,T_int):
    space.T[:,:]=T_int    
        
@nb.jit
def SolveEnergyEquation(space,crystal,fluid):    
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u=space.u.astype(float,copy=False)
    v=space.v.astype(float,copy=False)
    T=space.T.astype(float,copy=False)
    T_next=space.T_next.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    alpha=float(fluid.alpha)
    cp=float(fluid.cp)
    rhoc=float(crystal.rho)
    dH=float(crystal.dH)
    K_v=float(crystal.K_v)
    G=crystal.G.astype(float)
    A=crystal.A.astype(float)
    
    T1_y=(T[2:,1:cols+1]-T[0:rows,1:cols+1])/(2*dy)
    T1_x=(T[1:rows+1,2:]-T[1:rows+1,0:cols])/(2*dx)
    T2_y=(T[2:,1:cols+1]-2*T[1:rows+1,1:cols+1]+T[0:rows,1:cols+1])/(dy**2)
    T2_x=(T[1:rows+1,2:]-2*T[1:rows+1,1:cols+1]+T[1:rows+1,0:cols])/(dx**2)
    source=-3*(rhoc)*(dH/cp)*K_v*G[1:rows+1,1:cols+1]*A[1:rows+1,1:cols+1]
    T_next[1:rows+1,1:cols+1]=T[1:rows+1,1:cols+1]-dt*(u[1:rows+1,1:cols+1]*T1_x+v[1:rows+1,1:cols+1]*T1_y)+(dt*alpha*(T2_x+T2_y))+(dt*source)
    
def AdjustT(space):
    space.T[1:-1,1:-1]=space.T_next[1:-1,1:-1].copy()
    
def SetInitialC(space,C_int):
    space.C[:,:]=C_int    
 
@nb.jit
def SolveSpeciesEquation(space,crystal,fluid):
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u=space.u.astype(float,copy=False)
    v=space.v.astype(float,copy=False)
    C=space.C.astype(float,copy=False)
    C_next=space.C_next.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    D=float(fluid.D)
    rhoc=float(crystal.rho)
    MW=float(crystal.MW)
    K_v=float(crystal.K_v)
    G=crystal.G.astype(float)
    A=crystal.A.astype(float)
    
    C1_y=(C[2:,1:cols+1]-C[0:rows,1:cols+1])/(2*dy)
    C1_x=(C[1:rows+1,2:]-C[1:rows+1,0:cols])/(2*dx)
    C2_y=(C[2:,1:cols+1]-2*C[1:rows+1,1:cols+1]+C[0:rows,1:cols+1])/(dy**2)
    C2_x=(C[1:rows+1,2:]-2*C[1:rows+1,1:cols+1]+C[1:rows+1,0:cols])/(dx**2)
    source=-3*(rhoc/MW)*K_v*G[1:rows+1,1:cols+1]*A[1:rows+1,1:cols+1]
    C_next[1:rows+1,1:cols+1]=C[1:rows+1,1:cols+1]-dt*(u[1:rows+1,1:cols+1]*C1_x+v[1:rows+1,1:cols+1]*C1_y)+(dt*D*(C2_x+C2_y))+(dt*source)
    
def AdjustC(space):
    space.C[1:-1,1:-1]=space.C_next[1:-1,1:-1].copy()  
    
def GetCstar(space,crystal):
    T=space.T.astype(float,copy=True)
    T-=273.16
    Cstar_unmodded=crystal.a0+crystal.a1*T+crystal.a2*T**2
    Cstar=(Cstar_unmodded/crystal.MW)*10*1000/1000 #kmol/m3
    crystal.Cstar=Cstar.copy()

def GetG(space,crystal):
    kg=float(crystal.kg)
    g=float(crystal.g)
    C=space.C.astype(float)
    Cstar=crystal.Cstar.astype(float)
    G=kg*sci.sign(C-Cstar)*(abs(C-Cstar)/Cstar)**g
    crystal.G=G.copy()

def SetInitialNLAV(crystal,N_init,L_init):
    K_a=crystal.K_a
    K_v=crystal.K_v
    crystal.N[:,0]=N_init
    crystal.L[:,0]=L_init
    crystal.A[:,0]=K_a*L_init**2
    crystal.V[:,0]=K_v*L_init**3
    crystal.N_next[:,0]=N_init
    crystal.L_next[:,0]=L_init
    crystal.A_next[:,0]=K_a*L_init**2
    crystal.V_next[:,0]=K_v*L_init**3

@nb.jit
def SolvePopulationBalanceEquation(space,crystal):
    rows=int(space.rowpts)
    cols=int(space.colpts)
    u=space.u.astype(float,copy=False)
    dx=float(space.dx)
    dt=float(space.dt)
    G=crystal.G.astype(float)
    N=crystal.N.astype(float)
    L=crystal.L.astype(float)
    A=crystal.A.astype(float)
    V=crystal.V.astype(float)
    N_next=crystal.N_next.astype(float,copy=False)
    L_next=crystal.L_next.astype(float,copy=False)
    A_next=crystal.A_next.astype(float,copy=False)
    V_next=crystal.V_next.astype(float,copy=False)
    
    N1_x=(N[1:rows+1,1:cols+1]-N[1:rows+1,0:cols])/dx
    N_next[1:rows+1,1:cols+1]=N[1:rows+1,1:cols+1]-(dt*u[1:rows+1,1:cols+1]*N1_x)
    
    L1_x=(L[1:rows+1,1:cols+1]-L[1:rows+1,0:cols])/dx
    L_next[1:rows+1,1:cols+1]=L[1:rows+1,1:cols+1]-(dt*u[1:rows+1,1:cols+1]*L1_x)+(dt*G[1:rows+1,1:cols+1]*N[1:rows+1,1:cols+1])

    A1_x=(A[1:rows+1,1:cols+1]-A[1:rows+1,0:cols])/dx
    A_next[1:rows+1,1:cols+1]=A[1:rows+1,1:cols+1]-(dt*u[1:rows+1,1:cols+1]*A1_x)+2*(dt*G[1:rows+1,1:cols+1]*L[1:rows+1,1:cols+1])
    
    V1_x=(V[1:rows+1,1:cols+1]-V[1:rows+1,0:cols])/dx
    V_next[1:rows+1,1:cols+1]=V[1:rows+1,1:cols+1]-(dt*u[1:rows+1,1:cols+1]*V1_x)+2*(dt*G[1:rows+1,1:cols+1]*A[1:rows+1,1:cols+1])
    
def AdjustNLAV(crystal):
    crystal.N=crystal.N_next.copy()
    crystal.L=crystal.L_next.copy()
    crystal.A=crystal.A_next.copy()
    crystal.V=crystal.V_next.copy()
    
def WriteToFile(space,crystal,iteration,interval):
    if(iteration%interval==0):
        cwdir=os.getcwd()
        if(iteration==0):
            if(os.path.isdir("Result")==False):
                os.mkdir("Result")
                os.chdir("Result")
                cwdir=os.getcwd()
            elif(os.path.isdir("Result")==True):
                os.chdir("Result")
                cwdir=os.getcwd()
                filelist=os.listdir()
                for file in filelist:
                    os.remove(file)
        if(os.path.basename(os.getcwd())!="Result"):
            os.chdir("Result")
            cwdir=os.getcwd()
        filename=f"PUV{iteration}.txt"
        path=os.path.join(cwdir,filename)
        with open(path,"w") as f:
            for i in range(space.rowpts):
                for j in range(space.colpts):
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(space.p_c[i,j],space.u_c[i,j],space.v_c[i,j],space.T_c[i,j],space.C_c[i,j],crystal.N_c[i,j],crystal.L_c[i,j],crystal.A_c[i,j],crystal.V_c[i,j]))
    
#######################END OF FILE#############################################        
