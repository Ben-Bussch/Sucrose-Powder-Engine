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
    gamma = 1.4         #Ratio of specific heats
    
    def __init__(self, Ptank):
        self.Ptank = Ptank
    
    
    def airTank(self, dt, Vtank, Aexit, tf):
        """Numerically solves air tank blowout, using tank parameters and timestp dt"""
        #initial Parameters
        Pstart = self.Ptank
        m0 = (self.Ptank*Vtank)/(self.R*self.T0)
        print("m_start:", m0, "Ptank: ", self.Ptank)
        t = 0
        
        #Arrays
        P = [Pstart]
        m = [m0]
        T = [self.T0]
        rho =[m0/Vtank]
        time = [t]
        mdot = [-(rho[0]*Aexit*np.sqrt(2*(Pstart-self.P0)/rho[0]))]
        print("mdot start: ", mdot[0])
        
        
        count = 0
        
        while t < tf:
            
            #Initiate temperature variable to use in first loop
            if count == 0:
                T1 = T[count]
            
            #mass flowrate
            mdot1 = -(rho[count]*Aexit*np.sqrt(2*(P[count]-self.P0)/rho[count]))
            mdot.append(mdot1)
                
            #mass
            m1 = m[count] + (mdot[count+1]+mdot[count])*dt*0.5
            m.append(m1)
                        
            #Pressure
            P1 = P[count]*m1*T1/(m[count]*T[count])
            P.append(P1)
            
            #Temperature
            T1 = P1*T[count]/P[count]
            T.append(T1)
            
            #Density
            rho1 = m[count+1]/Vtank
            rho.append(rho1)

           
            
            t += dt
            time.append(t)
            
            count += 1

        
        print("Time: ",t,"mass: ", m[count], "mdot: ", mdot[count], "Pressure: ",P[count])
        return mdot, P,T, time
    
    def fuelTank(self, dt, Aexit, VR, m0):
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
        mdot = [-(rho_f*Aexit*np.sqrt(2*(Pstart-self.P0)/rho_f))]
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
            P1 = P[count]*(V[count]**self.gamma)/(V1**self.gamma)
            P.append(P1)
                        
            #Mass Flowrate
            mdot1 = -(rho_f*Aexit*np.sqrt(2*(P[count+1]-self.P0)/rho_f))
            mdot.append(mdot1)
            
            #Mass of Fuel in tank
            m1 = m[count] + (mdot[count+1]+mdot[count])*dt*0.5
            m.append(m1)
            
            t += dt
            time.append(t)
            
            count += 1
            
        print("Time: ",t,"mass: ", m[count], "mdot: ", mdot[count], "Pressure: ",P[count])
        return mdot, P, time
        
    def OFratio(self, mdot_a,mdot_f):
        """Calculates OF ratio from air and fuel tank blowout"""
        OF = []
        for count in range(len(mdot_f)):
            OFtemp = mdot_a[count]/mdot_f[count]
            OF.append(OFtemp)

        return OF
    
    def OFstats(self, OF, opt):
        """Calculates variance from optimal OF ratio and average OF ratio from tank blowout"""
        varianceSum,avgSum = 0,0
        for count in range(len(OF)):
           varianceSum += (OF[count]-opt)**2
           avgSum += OF[count]
        vari =  varianceSum/len(OF)
        avg = avgSum/len(OF)
        return vari, avg
        
            

"""Fuel tank Calcs"""
Pt = 11*10**5           #Pa, initial pressure in both tanks
Ae_f = 9.77809*10**-8   #m^2, Exit Area
VR = 1.000              #Volume ratio of fuel to air in fuel tank
m_f = 0.01100           #kg, inital fuel mass

dt = 0.001

SolverInstance = Solver(Pt)

mdot_f, P_f, t = SolverInstance.fuelTank(dt, Ae_f, VR, m_f)

"""Air tank calcs"""
V_a = 0.00806903       #m^3, volume of air tank
Ae_a = 5.04112*10**-6  #m^2, Exit Area

tf = t[-1]

mdot_a, P_a, T_a, t_a = SolverInstance.airTank(dt, V_a, Ae_a, tf)

OF = SolverInstance.OFratio(mdot_a,mdot_f)

OF_opt = 4.8

OF_variance, OF_avg = SolverInstance.OFstats(OF, OF_opt)

print("OF variance: ", OF_variance, "OF average: ", OF_avg)

plt.figure(1)
plt.clf()
plt.plot(t, mdot_f, label='Fuel Tank')
plt.plot(t_a, mdot_a, label='Air Tank')

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Mass Flow Rate [kg/s]")
plt.savefig('Mdot_vs_time.png', dpi=300)
plt.legend()


plt.figure(2)
plt.clf()
plt.plot(t, P_f, label='Fuel Tank')
plt.plot(t_a, P_a, label='Air Tank')

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Tank Pressure [Pa]")
plt.savefig('Pressure_vs_time.png', dpi=300)
plt.legend()

plt.figure(3)
plt.clf()
plt.plot(t, OF)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("OF ratio")
plt.savefig('OF_vs_time.png', dpi=300)

plt.figure(4)
plt.clf()
plt.plot(t_a, T_a)

plt.grid(1)
plt.xlabel("Time [s]")
plt.ylabel("Tempurature")
plt.savefig('Temp_vs_time.png', dpi=300)



plt.ion()


