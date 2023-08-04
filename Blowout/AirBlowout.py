# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 10:26:50 2023

@author: Bobke
"""
import numpy as np
import matplotlib.pyplot as plt

class Solver:
    T0 = 288.15         #K, atm temp
    P0 = 1.01325*10**5   #Atmospheric Pressure
    R = 287.17          #Specific Gas constant of air
    
    def __init__(self, Ptank, Vtank, Aexit):
        self.Ptank = Ptank
        self.Vtank = Vtank
        self.Aexit = Aexit
    
    def terminate(self, P):
        run = True
        Perror = (P-self.P0)/self.P0 
        if Perror <= 0.001:
            run = False
        return run
    
    def airTank(self, dt):
        "initial Parameters"
        Pstart = self.Ptank
        m0 = (self.Ptank*self.Vtank)/(self.R*self.T0)
        #print("m_start:", m0, "Ptank: ", self.Ptank)
        t = 0
        
        #Arrays
        P = [Pstart]
        m = [m0]
        rho = [m0/V]
        time = [t]
        mdot = [-(m0*self.Aexit/self.Vtank)*np.sqrt(2*self.R*self.T0 - self.P0*self.Vtank)]
        print("mdot start: ", mdot[0])
        
        run = True
        count = 0
        
        while run:
            
            mdot1 = -(m[count]*self.Aexit/self.Vtank)*np.sqrt(2*self.R*self.T0 - self.P0*self.Vtank)
            mdot.append(mdot1)
            
            m1 = m[count] + (mdot[count+1]+mdot[count])*dt*0.5
            m.append(m1)
            
            P1 = m[count]*self.R*self.T0/self.Vtank
            P.append(P1)
            
            t += dt
            time.append(t)
            
            run = self.terminate(P[count])
            count += 1

        
        print("Time: ",t,"mass: ", m[count], "mdot: ", mdot[count], "Pressure: ",P[count])
        return mdot, P, time
        
    

Pt = 11*10**5        #Pa, initial pressure in air tank
V = 0.00806903       #m^3, volume of air tank
Ae = 5.04112*10**-6

airInstance = Solver(Pt,V,Ae)

mdot, P, t = airInstance.airTank(0.001)


plt.figure(1)
plt.clf()
plt.plot(t, mdot)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Mass Flow Rate [kg/s]")
plt.savefig('Mdot_vs_time.png', dpi=300)



plt.figure(2)
plt.clf()
plt.plot(t, P)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Tank Pressure [Pa]")
plt.savefig('Pressure_vs_time.png', dpi=300)
plt.ion()


