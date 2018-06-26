#!/usr/bin/env python 

## \file scipy_tools.py
#  \brief tools for interfacing with scipy
#  \author T. Lukaczyk, F. Palacios
#  \version 4.3.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy

from .. import eval as su2eval
from numpy import array, zeros
from numpy.linalg import norm
from numpy import *
from numpy.polynomial.hermite import hermgauss as hermite
from numpy.random import normal as samplenormal
from numpy import exp as exp


# -------------------------------------------------------------------
#  Scipy SLSQP
# -------------------------------------------------------------------

def calcdistance(projectx,config,state,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    project=projectx
    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []
        
    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv
    
    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
#    nquad=12
    nquad=4
    sigma=sqrt(0.0001)
    mu=0.8

    print "Distance"
    distance=distancefunc(x0,project,nquad,mu,sigma)
    print repr(distance)

    return "finished"

def calclosses(projectx,config,state,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    project=projectx
    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []
        
    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv
    
    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
#    nquad=12
    nquad=4
    sigma=sqrt(0.0001)
    mu=0.8

    print "Losses"
    losses=lossesfunc(x0,project,nquad,mu,sigma)
    print repr(losses)

    return "finished"

def stochastic(projectx,config,state,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_slsqp(project,x0=[],xb=[],its=100,accu=1e-10)
    
        Runs the Scipy implementation of SLSQP with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10
        
        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """
    project=projectx
    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []
        
    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv
    
    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
    nquad=12
    sigma=sqrt(0.0001)
    mu=0.8

    obj=obj_f(x0,project)
    print "Function Value"
    print repr(obj)

    print "Expected Value"
    exp=exp_f(x0,project,nquad,mu,sigma)
    print repr(exp) 
    
    print "Variance"
    var=var_f(x0,project,nquad,mu,sigma)
    print repr(var)
 
    obj=obj_g(x0,project)
    print "Function Value g"
    print repr(obj)

    print "Expected Value g"
    exp=exp_g(x0,project,nquad,mu,sigma)
    print repr(exp) 
    
#    print "Losses"
#    losses=lossescalc(x0,project,nquad,mu,sigma)
#    print repr(losses)

    print "Variance g"
    var=var_g(x0,project,nquad,mu,sigma)
    print repr(var)

    print "Sample"
    sampling_number=200
    s = samplenormal(mu, sigma, sampling_number)
    out=sampling_f(x0,project,mu,sigma,sampling_number, s)

    print "Sample g"
    out=sampling_g(x0,project,mu,sigma,sampling_number, s)

    #print "Scatter"
    #out=scatter_f(x0,project)

    #print "Scatter g"
    #out=scatter_g(x0,project)
    # Done
    return "finished"

def test_example(projectx,config,state,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    #same normal distribution, but sample function, where I now the result 
    #sigma=sqrt(0.0001)
    #mu=0.8
    #mu=0.0
    #sigma=1.0
    mu=1.0
    sigma=0.25
    tin=1.0 
    beta=1.0

    print "f(xmean)"
    f=exp(-mu*tin)
    print repr(f)

    print "exp(f)"
    expf=beta*exp(-mu*tin+sigma*sigma*tin*tin*0.5)
    print repr(expf)
    expfreal=expf
    for nquad in range(2,20):
        points, weights = hermite(nquad)
        sumobj=0
        for count in range(0,nquad):
            inp=sqrt(2)*sigma*points[count]+mu
            #sumobj=sumobj+weights[count]*beta*exp(-(mu+sigma*inp)*tin)
            sumobj=sumobj+weights[count]*beta*exp((-inp)*tin)
        expf=sumobj/sqrt(pi)
        print repr(nquad)+" "+repr(expf)
    #expfreal=expf
    mcf=[]
    sumobj=0
    s = samplenormal(mu, sigma, 10000)
    for count in range(0,10000):
        inp=s[count]
        #sumobj=sumobj+beta*exp(-(mu+sigma*inp)*tin)
        sumobj=sumobj+beta*exp(-inp*tin)
    expf=sumobj/10000
    print "mc: "+repr(expf)

    print "var(f)"
    varf=exp(-2*mu*tin)*beta*beta*(exp(2*sigma*sigma*tin*tin)-exp(sigma*sigma*tin*tin))
    print repr(varf)
    varf=0
    for degree in range(1,5):
        nquad=20
        points, weights = hermite(nquad)
        sumobj=0
        for count in range(0,nquad):
            inp=sqrt(2)*sigma*points[count]+mu
            inp1=sqrt(2)*points[count]
            if(degree==1): 
                #sumobj=sumobj+weights[count]*beta*exp(-(mu+sigma*inp)*tin)*inp
                sumobj=sumobj+weights[count]*beta*exp(-inp*tin)*inp1
            if(degree==2):
                #sumobj=sumobj+weights[count]*beta*exp(-(mu+sigma*inp)*tin)*(inp*inp-1)
                sumobj=sumobj+weights[count]*beta*exp(-inp*tin)*(inp1*inp1-1)
            if(degree==3):
                #sumobj=sumobj+weights[count]*beta*exp(-(mu+sigma*inp)*tin)*(inp*inp*inp-3*inp)
                sumobj=sumobj+weights[count]*beta*exp(-inp*tin)*(inp1*inp1*inp1-3*inp1)
            if(degree==4):
                #sumobj=sumobj+weights[count]*beta*exp(-(mu+sigma*inp)*tin)*(inp*inp*inp*inp-6*inp*inp+3)
                sumobj=sumobj+weights[count]*beta*exp(-inp*tin)*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
        varf+=sumobj*sumobj/pi*(1./factorial(degree))
        print repr(degree)+" "+repr(varf)
    mcf=[]
    number=1000
    for count in range(1,10):
        number=number*2
        sumobj=0
        sumobj2=0
        expfreal=0
        s = samplenormal(mu, sigma, number)
        for count in range(0,number):
            inp=s[count]
            expfreal=expfreal+beta*exp(-inp*tin)
            #expfreal=expfreal+beta*exp(-(mu+sigma*inp)*tin)
        expfreal=expfreal/number
        for count in range(0,number):
            inp=s[count]
            sumobj=sumobj+(expfreal-beta*exp(-inp*tin))*(expfreal-beta*exp(-inp*tin))
            sumobj2=sumobj2+beta*exp(-inp*tin)*beta*exp(-inp*tin)
            #sumobj=sumobj+(expfreal-beta*exp(-(mu+sigma*inp)*tin))*(expfreal-beta*exp(-(mu+sigma*inp)*tin))
            #sumobj2=sumobj2+beta*exp(-(mu+sigma*inp)*tin)*beta*exp(-(mu+sigma*inp)*tin)
        sumobjnew=sumobj2/(number)-expfreal*expfreal
        sumobj=sumobj/(number-1)
        print "mc: "+repr(sumobj)+" "+repr(sumobjnew)

def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

def convergence(projectx,config,state,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    project=projectx
    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []

    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv

    # Initial guess
    if not x0: x0 = [0.0]*n_dv

    sigma=sqrt(0.0001)
    mu=0.8

    for nquad in range(2,20):
        print "Expected Value"
        exp=exp_f(x0,project,nquad*2,mu,sigma)
        print repr(exp)
        print "Variance"
        var=var_f(x0,project,nquad*2,mu,sigma)
        print repr(var)
        print "Expected Value g"
        exp=exp_g(x0,project,nquad*2,mu,sigma)
        print repr(exp) 
        print "Variance g"
        var=var_g(x0,project,nquad*2,mu,sigma)
        print repr(var)
    
def obj_f(x,project):
    """ obj = obj_f(x,project)
        
        Objective Function
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         minimize f(x), list[nobj]
        scipy_slsqp: minimize f(x), float
    """
    project.config.MACH_NUMBER=0.8    
    obj = project.obj_f(x)
    return obj

def exp_f(x,project,nquad, mu,sigma):
    points, weights = hermite(nquad)
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        obj = project.obj_f(x)
        #print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(obj)
        sumobj=sumobj+weights[count]*obj[0]
    exp=sumobj/sqrt(pi)
    return exp

def var_f(x,project,nquad, mu, sigma):
    points, weights = hermite(nquad)
    varf=0
    for degree in range(1,8):
        sumobj=0
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
            obj = project.obj_f(x) 
            inp1=sqrt(2)*points[count]
            if(degree==1):
                sumobj=sumobj+weights[count]*obj[0]*inp1
            if(degree==2):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1-1)
            if(degree==3):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1-3*inp1)
            if(degree==4):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
            if(degree==5):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
            if(degree==6):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
            if(degree==7):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
            if(degree==8):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
            if(degree==9):
                sumobj=sumobj+weights[count]*obj[0]*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
        varf+=sumobj*sumobj/pi*(1./factorial(degree))
        #print repr(degree)+" "+repr(varf)
    return varf


def obj_g(x,project):     
    project.config.MACH_NUMBER=0.8
    econs = project.con_ceq(x)
    iecons = project.con_cieq(x)
    cons = array(econs+iecons, float_)
    return cons

def exp_g(x,project,nquad,mu,sigma): 
    points, weights = hermite(nquad)
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        econs = project.con_ceq(x)
        iecons = project.con_cieq(x)
        cons = array(econs+iecons, float_)
        #print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(cons[0])
        sumobj=sumobj+weights[count]*cons
    exp=sumobj/sqrt(pi)
    return exp

def var_g(x,project,nquad, mu, sigma):  
    points,weights=hermite(nquad)
    varf=0
    for degree in range(1,8):
        sumobj=0
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu 
            econs = project.con_ceq(x)
            iecons = project.con_cieq(x)
            cons = array(econs+iecons, float_)
            inp1=sqrt(2)*points[count]
            if(degree==1):
                sumobj=sumobj+weights[count]*cons*inp1
            if(degree==2):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1-1)
            if(degree==3):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1-3*inp1)
            if(degree==4):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
            if(degree==5):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
            if(degree==6):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
            if(degree==7):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
            if(degree==8):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
            if(degree==9):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
        varf+=sumobj*sumobj/pi*(1./factorial(degree))
        #print repr(degree)+" "+repr(varf)
    return varf

def scatter_f(x,project):
    mach_start=0.75
    for count in range(0,51):
        project.config.MACH_NUMBER=mach_start+count*(0.1/50)
        obj = project.obj_f(x)
        print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(obj)
    return 0.0

def scatter_g(x,project):
    mach_start=0.75
    for count in range(0,51):
        project.config.MACH_NUMBER=mach_start+count*(0.1/50)
        econs = project.con_ceq(x)
        iecons = project.con_cieq(x)
        cons = array(econs+iecons, float_)
        print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(cons[0])
    return 0.0

def sampling_f(x,project, mu, sigma, number, s):
    for count in range(0,number):
        project.config.MACH_NUMBER=s[count]
        obj = project.obj_f(x)
        print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(obj)
    return 0.0

def sampling_g(x,project, mu, sigma,number, s):
    for count in range(0,number):
        project.config.MACH_NUMBER=s[count]
        econs = project.con_ceq(x)
        iecons = project.con_cieq(x)
        cons = array(econs+iecons, float_)
        print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(cons[0])
    return 0.0

def lossesfunc(x,project,nquad,mu,sigma):
    points, weights = hermite(nquad)
    #calculate fmean
    project.config.MACH_NUMBER=0.8
    obj=project.obj_f(x)
    objmean=obj[0]

    #calculate E(max(0,f-fmean))
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        obj=project.obj_f(x)
        objnew = max(0,obj[0]-objmean)
        #print "f: "+repr(obj[0])+" "+repr(objmean)+" "+repr(objnew)
        sumobj=sumobj+weights[count]*objnew
    sumobj=sumobj/(sqrt(pi))
    #obj_scale = obj[obj.keys()[0]]['SCALE']
    #print "sumobj: "+repr(sumobj)
    #print "objmean: "+repr(objmean)
    sumobj=sumobj/abs(objmean) #divide by objmean?
    print "fvalue: "+repr(sumobj)

    #calculate gmean    
    project.config.MACH_NUMBER=0.8
    econs = project.con_ceq(x)
    iecons=project.con_cieq(x)
    consres = array(econs+iecons, float_)
    gmean=consres[0]

    #calculate E(max(0,g-gmean))
    sumcons=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        econs = project.con_ceq(x)
        iecons = project.con_cieq(x)
        cons = array(econs+iecons, float_)
        #print "cons"
        #print cons[0]
        #print gmean
        consnew = max(0,cons[0]-gmean)
        print "g: "+repr(cons[0])+" "+repr(gmean)+" "+repr(consnew)
        sumcons=sumcons+weights[count]*consnew
    sumcons=sumcons/sqrt(pi)
    #cons_scale = cons[cons.keys()[0]]['SCALE']
    #print "gmean"+repr(gmean)
    #print "sumcons: "+repr(sumcons)
    sumcons=sumcons/abs(gmean)
    print "gvalue: "+repr(sumcons)
    return sumcons+sumobj

def distancefunc(x,project,nquad,mu,sigma):
    points, weights = hermite(nquad)
    #calculate E(max(dist(p,paretop),0)
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        obj=project.obj_f(x)
        econs = project.con_ceq(x)
        iecons = project.con_cieq(x)
        cons = array(econs+iecons, float_)
        objnew = max(0,distance_to_front(obj[0],cons[0])) #nomax before
        sumobj=sumobj+weights[count]*objnew
    sumobj=sumobj/(sqrt(pi))
    return sumobj

def distance_to_front(x_val,y_val):
    scalef=0.00907643442
    scaleg=0.479704775
    x_front=[0.000766592145/scalef,0.00123225475/scalef,0.0022372218/scalef,0.00489305698/scalef,0.00502461309/scalef,0.00885973466/scalef,0.0132521623/scalef,0.0173862767/scalef]
    y_front=[-0.15919563/scaleg,-0.31068785/scaleg,-0.47404976/scaleg,-0.63391492/scaleg,-0.63638229/scaleg,-0.71876571/scaleg,-0.76617848/scaleg,-0.80021392/scaleg]
    distance=100000
#    print x_val, y_val, x_val/scalef, y_val/scaleg
    #print x_front
    #print y_front
    for i in range(0,7):
        normal_vector_x=-(y_front[i+1]-y_front[i])
        normal_vector_y=(x_front[i+1]-x_front[i])
        length=sqrt(normal_vector_x*normal_vector_x+normal_vector_y*normal_vector_y)
        dist=1./length*(normal_vector_x*(x_val/scalef-x_front[i])+normal_vector_y*(y_val/scaleg-y_front[i]))
        #print "dist "+repr(distance)
        if(dist<distance):
            distance=dist
    #print "Distance "+repr(distance)
    return distance
