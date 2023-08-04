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
        """Numerically solves air tank blowout, using tank parameters and timestp dt"""
        #initial Parameters
        Pstart = self.Ptank
        m0 = (self.Ptank*self.Vtank)/(self.R*self.T0)
        #print("m_start:", m0, "Ptank: ", self.Ptank)
        t = 0
        
        #Arrays
        P = [Pstart]
        m = [m0]
        rho = [m0/self.Vtank]
        time = [t]
        mdot = [-(m0*self.Aexit/self.Vtank)*np.sqrt(2*self.R*self.T0 - self.P0*self.Vtank)]
        print("mdot start: ", mdot[0])
        
        run = True
        count = 0
        
        while run:
            
            #mass flowrate
            mdot1 = -(m[count]*self.Aexit/self.Vtank)*np.sqrt(2*self.R*self.T0 - self.P0*self.Vtank)
            mdot.append(mdot1)
            
            #mass
            m1 = m[count] + (mdot[count+1]+mdot[count])*dt*0.5
            m.append(m1)
            
            #Pressure
            P1 = m[count]*self.R*self.T0/self.Vtank
            P.append(P1)
            
            t += dt
            time.append(t)
            
            run = self.terminate(P[count])
            count += 1

        
        print("Time: ",t,"mass: ", m[count], "mdot: ", mdot[count], "Pressure: ",P[count])
        return mdot, P, time
    
    def fuelTank(self, dt, VR, m0):
        """Numerically solves fuel tank blowout, using tank parameters and timestp dt"""
        
        #initial Parameters
        rho_f = 1590 #kg/m^3, density of sucrose
        Pstart = self.Ptank
        V0 = rho_f*m0/VR
        
        t = 0
        
        #Arrays
        P = [Pstart]
        V = [V0]
        m = [m0]
        time = [t]
        mdot = [-(rho_f*self.Aexit*np.sqrt(2*(Pstart-self.P0)/rho_f))]
        Vdot = [mdot[0]*rho_f]
        print("mdot start: ", mdot[0])
        
        count = 0
        while m[count] > 0:
            
            #Rate of chanfge of volume of air in tank
            Vdot1 = -mdot[count]*rho_f
            Vdot.append(Vdot1)
            
            #Volume in tank
            V1 = V[count]+ (Vdot[count+1]+Vdot[count])*dt*0.5
            V.append(V1)
            
            #Pressure in tank
            P1 = P[count]*V[count]/V1
            P.append(P1)
                        
            #Mass Flowrate
            mdot1 = -(rho_f*self.Aexit*np.sqrt(2*(P[count+1]-self.P0)/rho_f))
            mdot.append(mdot1)
            
            #Mass of Fuel in tank
            m1 = m[count] + (mdot[count+1]+mdot[count])*dt*0.5
            m.append(m1)
            
            t += dt
            time.append(t)
            
            count += 1
            
        print("Time: ",t,"mass: ", m[count], "mdot: ", mdot[count], "Pressure: ",P[count])
        return mdot, P, time
            
            
        
    
"""Air tank calcs"""
Pt_a = 11*10**5        #Pa, initial pressure in air tank
V_a = 0.00806903       #m^3, volume of air tank
Ae_a = 5.04112*10**-6  #m^2, Exit Area

airInstance = Solver(Pt_a,V_a,Ae_a)

mdot_a, P_a, t_a = airInstance.airTank(0.001)


plt.figure(1)
plt.clf()
plt.plot(t_a, mdot_a)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Mass Flow Rate [kg/s]")
plt.savefig('Mdot_vs_time.png', dpi=300)



plt.figure(2)
plt.clf()
plt.plot(t_a, P_a)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Tank Pressure [Pa]")
plt.savefig('Pressure_vs_time.png', dpi=300)
plt.ion()


"""Fuel tank Calcs"""
Pt_f = 11*10**5         #Pa, initial pressure in air tank
V_f = 0.00806903        #m^3, volume of air tank
Ae_f = 9.77809*10**-8   #m^2, Exit Area
VR = 0.5                #Volume ratio of fuel to air in fuel tank
m_f = 0.015             #kg, inital fuel mass

fuelInstance = Solver(Pt_f,V_f,Ae_f)

mdot_f, P_f, t_f = fuelInstance.fuelTank(0.001, VR, m_f)

plt.figure(3)
plt.clf()
plt.plot(t_f, mdot_f)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Mass Flow Rate [kg/s]")
plt.savefig('Mdot_vs_time.png', dpi=300)



plt.figure(4)
plt.clf()
plt.plot(t_f, P_f)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Tank Pressure [Pa]")
plt.savefig('Pressure_vs_time.png', dpi=300)
plt.ion()


