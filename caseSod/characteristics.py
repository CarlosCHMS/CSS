#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 12:27:12 2020

@author: carlos
"""

import matplotlib.pyplot as plt
import numpy as np


class problem():

    def __init__(self, xchange=50.0, p1=1e5, T1=300, p4=1e6, T4=300, R=287,
                 gamma=1.4):

        # Solution of the Riemann unidimensional flow problem using
        # characteristics method
        # Inputs: gas properties in the region 1 and 4;
        #         position of the properties change

        self.p1 = p1
        self.T1 = T1
        self.p4 = p4
        self.T4 = T4

        self.xchange = xchange

        self.R = R
        self.g = gamma

        # Tolerance for bissection method application
        self.tol = 1e-7

        # Set gases
        self.g1 = Gas(self.p1, self.T1, self.R, self.g)
        self.g4 = Gas(self.p4, self.T4, self.R, self.g)

        # A central problem is calculate the shock velocity. It is made by
        # considering contour conditions in all regions
        # The velocity is numerically obtained by applicaiton of bissection
        # method.

        # Final value of shock velocity
        self.vs = self.bisec()

        # Properties of the gas in region 2
        self.g2, v2r = self.g1.calcShock(self.vs)

        # Velcity of the contact region
        self.vc = self.vs - v2r

        # Gas propeties in the region 3
        self.g3 = self.g4.calcIsentropic(0.0, self.vc)

        # Velocity of the expansion fan in the region 3. It is a characteristic
        # curve.
        self.ve3 = self.vc - self.g3.a

        # Velocity of the expansion fan in the region 4. It is a characteristic
        # curve.
        self.ve4 = -self.g4.a

    def calcError(self, vs):

        # Input: gess of shock velocity

        # Gas properties after the shock
        g2, v2 = self.g1.calcShock(vs)

        # Flow velocity before de shock (region 2) relativelly the stagnant gas
        # It is also the contact velocity
        vc = vs-v2

        # Expasion of the gas from region 4 considering the contact velocity
        g3 = self.g4.calcIsentropic(0.0, vc)

        # In the contact the pressure and velocity must be the same in both
        # sides
        # The error consists in the relative difference between the two
        # pressure values.

        #print('vs: ', vs, 'vc: ', vc, ', p2: ', g2.p, ', p3:', g3.p)
        return g2.p/g3.p - 1

    def bisec(self):

        # Bissection method. The objective is obtain the zero of calError(vs).
        # The result is the shock velocity.
        dv = self.g1.a*0.1

        v0 = self.g1.a
        v1 = v0 + dv

        e0 = self.calcError(v0)
        e1 = self.calcError(v1)

        ii = 0
        while abs(e1) > self.tol:

            step = (v1-v0)*e1/(e1-e0)

            if step > dv:
                step = dv
            elif step < -dv:
                step = -dv

            v2 = v1 - step
            e2 = self.calcError(v2)

            v0 = v1
            e0 = e1
            v1 = v2
            e1 = e2

            ii += 1

        return v1

    def calcVfan(self, t, x):

        # Calculate the gas velocity in the expansion fan.
        # Inputs: instant t and the position x.

        # Negative characteristic velocity
        # av = u - a
        av = x/t

        # Positive characteristic invariant
        # Cpp = (g-1)*u/2 + a
        # Solving the system:
        v = 2*(self.g4.a + av)/(self.g + 1)

        return v

    def plotXT(self, t):

        plt.figure()
        tt = np.array([0, t])
        plt.plot(self.vs*tt + self.xchange, tt)
        plt.plot(self.vc*tt + self.xchange, tt)
        plt.plot(self.ve3*tt + self.xchange, tt)
        plt.plot(self.ve4*tt + self.xchange, tt)
        plt.grid(True)
        plt.xlabel('x')
        plt.ylabel('t')
        plt.show()

        return None

    def calcVar(self, t, xx):

        xs = self.vs*t + self.xchange
        xc = self.vc*t + self.xchange
        xe3 = self.ve3*t + self.xchange
        xe4 = self.ve4*t + self.xchange

        var = dict()

        # print(xe4, xe3, xc, xs)

        var['t'] = t
        var['x'] = xx
        var['rho'] = xx*0
        var['p'] = xx*0
        var['T'] = xx*0
        var['u'] = xx*0
        var['hs'] = xx*0

        # Claculation by regions
        for ii in range(0, len(xx)):

            if xx[ii] < xe4:
                # Region 4: High pressure reservoir
                var['p'][ii] = self.g4.p
                var['rho'][ii] = self.g4.rho
                var['T'][ii] = self.g4.T
                var['u'][ii] = 0.0
                var['hs'][ii] = self.g4.h

            elif xx[ii] < xe3:
                # Region 4-3: expansion fan
                v = self.calcVfan(t, xx[ii] - self.xchange)

                gg = self.g4.calcIsentropic(0.0, v)

                var['p'][ii] = gg.p
                var['rho'][ii] = gg.rho
                var['T'][ii] = gg.T
                var['u'][ii] = v
                var['hs'][ii] = gg.h

            elif xx[ii] < xc:
                # Region 3: between the expansion fan and the constact
                var['p'][ii] = self.g3.p
                var['rho'][ii] = self.g3.rho
                var['T'][ii] = self.g3.T
                var['u'][ii] = self.vc
                var['hs'][ii] = self.g3.h

            elif xx[ii] < xs:
                # Region 2: between the constact and the shock
                var['p'][ii] = self.g2.p
                var['rho'][ii] = self.g2.rho
                var['T'][ii] = self.g2.T
                var['u'][ii] = self.vc
                var['hs'][ii] = self.g2.h

            else:
                # Region 1: low pressure reservoir
                var['p'][ii] = self.g1.p
                var['rho'][ii] = self.g1.rho
                var['T'][ii] = self.g1.T
                var['u'][ii] = 0.0
                var['hs'][ii] = self.g1.h

        # Total entalpy massic density
        var['h'] = var['hs'] + (var['u']*var['u'])/2
        # Total internal energy volumetric density
        var['e'] = var['h']*var['rho'] - var['p']
        # Mass flux
        var['m'] = var['rho']*var['u']

        return var


class Gas():

    def __init__(self, p, T, R, g):

        # Gas state class
        # Inputs: pressure [Pa], temperature [K]

        # Gas constants
        self.g = g  # Ratio of specific heats
        self.R = R  # Gas constant in therms of mass

        self.p = p
        self.T = T

        # Density calculation
        self.rho = self.p/(self.R*self.T)

        # Sound velocity calculation
        self.a = np.sqrt(self.g*self.R*self.T)

        # Static entalpy massic density calculation
        self.h = self.g*self.R*self.T/(self.g - 1)

    def calcShock(self, v0):

        # Calculate the resulting gas state after a shock.
        # Input: gas velocity relative to the shock wave
        # Output: gas state after the shock, velocity after the shock relative
        # to the shock wave

        # Calculations made by using typical normal shock relations

        # Initial mach number
        M = v0/self.a

        # Auxiliary values:
        b = ((self.g - 1)*M*M + 2)
        c = (2*self.g*M*M - (self.g - 1))

        # Mach number after the shock
        M1 = np.sqrt(b/c)

        # Pressure after the shock
        p1 = self.p*c/(self.g + 1)

        # Temperature after the shock
        T1 = self.T*b*c/((self.g+1)*M*(self.g+1)*M)

        # Sound velocity after the shock
        a1 = self.a*np.sqrt(T1/self.T)

        # After shock velocity
        v1 = M1*a1

        return Gas(p1, T1, self.R, self.g), v1

    def calcIsentropic(self, v0, v1):

        # Calculate the resulting gas state after a expansion.
        # Input: initial gas velocity and final gás velocity
        # Output: gas state after the expansion

        # Positive characteristic invariant
        Cp = self.a + v0*(self.g - 1)/2

        # Sound velocity after the expansion
        a1 = Cp - v1*(self.g - 1)/2

        # Temperature calculated from sound velocity definition
        T1 = ((a1/self.a)**2)*self.T

        # Isentropic relation
        p1 = self.p*((T1/self.T)**(self.g/(self.g-1)))

        return Gas(p1, T1, self.R, self.g)

    def printVar(self):

        # Print gas properties
        print("\np = ", self.p)
        print("T = ", self.T)
        print("rho = ", self.rho)
        print("a = ", self.a)


if __name__ == "__main__":

    prob1 = problem(xchange=50.0, p1=1e5, T1=300, p4=1e6, T4=300)

    xx = np.linspace(0, 100, num=100)
    var = prob1.calcVar(0.05, xx)

    plt.close('all')
    plt.figure()

    v = 'rho'
    # ploting characteristics results
    plt.plot(var['x'], var[v])

    plt.xlabel("%s" % ('x'))
    plt.ylabel("%s" % (v))

    title = "t: %f [s]" % var['t']
    plt.title(title)

    plt.grid(True)
    plt.show()

    prob1.plotXT(0.05)
