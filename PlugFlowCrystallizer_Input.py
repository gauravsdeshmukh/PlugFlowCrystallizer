# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 17:33:26 2018

@author: Gaurav
"""

import scipy as sci
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from PlugFlowCrystallizer import *

###############################################################################
###########################USER INPUT BEGINS###################################
############################DEFINE SPATIAL AND TEMPORAL PARAMETERS#############
length=10
breadth=0.4
colpts=150
rowpts=40
time=100
crystal_timedelay=20 #The time after which crystals are introduced in the system
###############################MISC############################################
CFL_number=0.1 #Do not touch this unless solution diverges
file_flag=1 #Keep 1 to print results to file
file_interval=50 #Regular time interval after which files are written
plot_flag=1 #Keep 1 to plot results at the end
###########################DEFINE PHYSICAL PARAMETERS##########################
rho=1
mu=0.01
k=0.6
cp=41.8
D=1.5e-02
###########################DEFINE CRYSTAL PARAMETERS###########################
MW=101.1032
kg=1.1612*10**(-4)
g=1.32
rhoc=1400
dH=-50
##########################DEFINE INITIAL MOMENTUM PARAMETERS###################
u_in=2
v_wall=0
p_out=0
##########################DEFINE INITIAL TEMPERATURE PARAMETERS###############
h=4000
T_jacket=303
T_in=363
T_init=293
##########################DEFINE INITIAL CONCENTRATION CONDITIONS##############
C_init=0
C_in=150
#########################DEFINE INITIAL CRYSTAL CONDITIONS#####################
N_init=1000
L_cryst=1e-3
L_init=N_init*L_cryst

###############################################################################
########################CREATE SPACE OBJECT####################################
pfc=Space()
pfc.CreateMesh(rowpts,colpts)
pfc.SetDeltas(breadth,length)
water=Fluid(rho,mu,k,cp,D)

###############################################################################
#########################BOUNDARY DEFINITIONS##################################
########################HEAT CONVECTION BOUNDARY###############################
jacketvaluetop=1/(1-(h*pfc.dy/water.k))
jacketconstanttop=(h*pfc.dy/water.k)/(1-(h*pfc.dy/water.k))*T_jacket
jacketvaluebottom=1/(1+(h*pfc.dy/water.k))
jacketconstantbottom=(h*pfc.dy/water.k)/(1+(h*pfc.dy/water.k))*T_jacket
########################CREATE BOUNDARY OBJECTS################################
###########################VELOCITY############################################
flow=Boundary("D",u_in)
noslip=Boundary("D",v_wall)
zeroflux=Boundary("N",0)
############################PRESSURE###########################################
pressureatm=Boundary("D",p_out)
############################TEMPERATURE########################################
jackettop=Boundary("C",jacketvaluetop,jacketconstanttop)
jacketbottom=Boundary("C",jacketvaluebottom,jacketconstantbottom)
intemp=Boundary("D",T_in)
############################CONCENTRATION######################################
inconc=Boundary("D",C_in)

###############################################################################
#########################CREATE CRYSTAL OBJECT#################################
kno3=Crystal("Potassium nitrate",pfc)
kno3.SetConstants(MW,kg,g,rhoc,dH)

###############################################################################
#########################SET INITIAL CONDITIONS IN DOMAIN######################
SetInitialT(pfc,T_init)
SetInitialC(pfc,C_init)
SetInitialNLAV(kno3,N_init,L_init)

#######################USER INPUT ENDS#########################################
###############################################################################
#############################INITIALIZATION####################################
t=0
i=0
############################THE RUN############################################
while(t<time):
    CFL=CFL_number
    SetTimeStep(CFL,pfc,water)
    timestep=pfc.dt
    
    SetUBoundary(pfc,flow,zeroflux,noslip,noslip)
    SetVBoundary(pfc,zeroflux,zeroflux,noslip,noslip)
    SetPBoundary(pfc,zeroflux,pressureatm,zeroflux,zeroflux)
    SetTBoundary(pfc,intemp,zeroflux,jackettop,jacketbottom)
    SetCBoundary(pfc,inconc,zeroflux,zeroflux,zeroflux)
    
    GetStarredVelocities(pfc,water)
    SolvePressurePoisson(pfc,water,zeroflux,pressureatm,zeroflux,zeroflux)
    SolveMomentumEquation(pfc,water)
    AdjustUV(pfc)
    
    SolveEnergyEquation(pfc,kno3,water)
    AdjustT(pfc)
    
    SolveSpeciesEquation(pfc,kno3,water)
    AdjustC(pfc)
    
    if(t>crystal_timedelay):
        GetCstar(pfc,kno3)
        GetG(pfc,kno3)
        
        SolvePopulationBalanceEquation(pfc,kno3)
        AdjustNLAV(kno3)
    
    SetCentrePUVTCNLAV(pfc,kno3)
    
    if(t<40):
        file_flag=1
    else:
        file_flag=0
    
    if(file_flag==1):
        WriteToFile(pfc,kno3,i,file_interval)
    
    u=pfc.u
    v=pfc.v
    p=pfc.p
    T=pfc.T
    C=pfc.C
    G=kno3.G
    N=kno3.N
    L=kno3.L
    u_c=pfc.u_c
    v_c=pfc.v_c
    p_c=pfc.p_c
    T_c=pfc.T_c
    C_c=pfc.C_c
    N_c=kno3.N_c
    L_c=kno3.L_c
    A_c=kno3.A_c
    V_c=kno3.V_c
    t+=timestep
    i+=1
    #nan check
    if(sci.isnan(sci.sum(C))):
        break
    print(time-t)

###########################END OF RUN##########################################
###############################################################################
#######################SET ARRAYS FOR PLOTTING#################################
x=sci.linspace(0,length,colpts)
y=sci.linspace(0,breadth,rowpts)
[X,Y]=sci.meshgrid(x,y)

u=pfc.u
v=pfc.v
p=pfc.p
T=pfc.T
C=pfc.C
G=kno3.G
N=kno3.N
L=kno3.L
u_c=pfc.u_c
v_c=pfc.v_c
p_c=pfc.p_c
T_c=pfc.T_c
C_c=pfc.C_c
N_c=kno3.N_c
L_c=kno3.L_c

######################EXTRA PLOTTING CODE BELOW################################
if(plot_flag==1):
    plt.figure(figsize=(16,8))
    plt.contourf(X,Y,p_c,cmap=cm.viridis)
    plt.colorbar()
    plt.quiver(X,Y,u_c,v_c)
    plt.title("Velocity and Pressure Plot")
    
    plt.figure(figsize=(16,8))
    plt.contourf(X,Y,T_c,cmap=cm.inferno)
    plt.colorbar()
    plt.title("Temperature Plot")
    
    plt.figure(figsize=(16,8))
    plt.contourf(X,Y,C_c,cmap=cm.viridis)
    plt.colorbar()
    plt.title("Concentration Plot")
    
    plt.figure(figsize=(16,8))
    plt.contourf(X,Y,L_c/N_c,cmap=cm.viridis)
    plt.colorbar()
    plt.title("Crystal Length Plot")




###########################END OF FILE#########################################