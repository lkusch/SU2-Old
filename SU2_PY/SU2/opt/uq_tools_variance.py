## \file scipy_tools.py
#  \brief tools for interfacing with scipy
#  \author T. Lukaczyk, F. Palacios
#  \version 3.2.9 "eagle"
#
# SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

import SU2
from .. import eval as su2eval
from numpy import *
from numpy.linalg import norm
from numpy.polynomial.hermite import hermgauss as hermite


# -------------------------------------------------------------------
#  CoinIpopt (BFGS)
# -------------------------------------------------------------------
global project,ncon, points, weights, nquad, sigma, mu, liftdrag, constraint

#TODO project als globale Variable / Klassenvariable?    
def eval_f(x,user_data=None):
    print "eval_f"
    obj = exp_f(x,project,nquad, mu,sigma)
    print obj
    return obj

def eval_grad_f(x,user_data=None):
    """ dobj = eval_grad_f(x)
        
        Objective Function Gradients
        SU2 Project interface to Ipopt
        
        su2:         df(x), list[nobj x dim]
        ipopt: 	     df(x), ndarray[dim]
    """    
    assert len(x)==project.n_dv
    #print "eval_grad_f"
    dobj = dexp_f(x,project,nquad, mu,sigma)
    print dobj 
    #dobj = array( dobj[0], float_)
    #print dobj
    return dobj

def eval_g(x,user_data=None):
    print "eval_g"
    """ cons = eval_g(x)
        
        Equality Constraint and Inequality Constraint Functions
        SU2 Project interface to Ipopt
        
        su2:         ceq(x),cieq(x) =/< 0.0, list[nceq]
        ipopt:       cons(x) = 0.0, ndarray[nceq]
    """
    assert len(x)==project.n_dv
    cons=exp_g(x,project,nquad,mu,sigma)   
    if(liftdrag==1):
        consf=(var_f(x,project,nquad,mu,sigma)*100-0.001)*1E-3
        consg=(var_g(x,project,nquad,mu,sigma)-0.001)*1E-3
    else:
        consg=(var_g(x,project,nquad,mu,sigma)*100-0.001)*1E-3
        consf=(var_f(x,project,nquad,mu,sigma)-0.001)*1E-3
    consbetw=append(cons, consf)
    consend=append(consbetw,consg)
    print consend
    return consend

def eval_jac_g(x,flag, user_data=None):
    """ dcons = eval_jac_g(x,flag)
        
        Equality Constraint Gradients
        SU2 Project interface to Ipopt
        
        su2:         dceq(x), dcieq(x), list[nceq x dim]
        ipopt: 	     dcons(x), ndarray[nceq*dim]
    """
    if flag:
	iRow=[]
	jCol=[]
  	for i in range(0, ncon):
	    for j in range(0, project.n_dv):
		iRow.append(i);
		jCol.append(j);
	return (array(iRow), array(jCol))
    else:
	assert len(x)==project.n_dv
        print "eval_jac_g"
#        dconsnew=(sumdobj+sumdcons)*1E-3
        dcons=dexp_g(x,project,nquad,mu,sigma)
        if(liftdrag==1):   
            dconsf=dvar_f(x,project,nquad,mu,sigma)*100*1E-3
            dconsg=dvar_g(x,project,nquad,mu,sigma)*1E-3
        else:   
            dconsg=dvar_g(x,project,nquad,mu,sigma)*100*1E-3
            dconsf=dvar_f(x,project,nquad,mu,sigma)*1E-3
        #print dconsf
        #print "dconsg"
        #print dconsg
        #print "dcons"
        #print dcons
        dconsbetw=append(dcons, dconsf)
        dconsend=append(dconsbetw, dconsg)
        #print "all"
        #print dconsend
	values=[]
	for j in range(0, ncon):
	    for k in range(0, project.n_dv):
 		values.append(dconsend[j*project.n_dv+k])
        print array(values)
    	return array(values)

def eval_h(x, lagrange, obj_factor, flag, user_data = None):
    if flag:
        return 0
    else:
        return 0


def apply_new(x):
    return True

def ipopt_run(projectx, config, state, x0=None,xb_low=None, xb_up=None, its=100,accun=1e-10,used_constraint=0,grads=True):
    """ result = ipopt_run(project,x0=[],xb_low=[], xb_up=[], its=100,accu=1e-10)
    
        Runs the BFGS method of Ipopt with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10
        
        Outputs:
           result - the outputs from Ipopt
     """
    global project,ncon, points, weights, nquad, sigma, mu, liftdrag, constraint
    import pyipopt
    project =projectx
    constraint=used_constraint
    
    nquad=10
    sigma=sqrt(0.0001)
    mu=0.8
    points, weights = hermite(nquad)

    liftdrag=0

    if x0 is None: x0 = []
    if xb_low is None: xb_low = []
    if xb_up is None: xb_up = []
    
	# gradient handles
	# TODO no gradients given -> ipopt
    	#if project.config.get('GRADIENT_METHOD','NONE') == 'NONE': 
    	#    fprime         = None
    	#    fprime_eqcons  = None
    	#    fprime_ieqcons = None
    	#else:
    	#    fprime         = obj_df
    	#    fprime_eqcons  = con_dceq
   	#    fprime_ieqcons = con_dcieq        
    
    	# number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    nvar = sum( dv_size)
    project.n_dv = nvar
	
    x_L=array(xb_low)
    x_U=array(xb_up)

    nceq=len(project.config['OPT_CONSTRAINT']['EQUALITY'])
    ncieq=len(project.config['OPT_CONSTRAINT']['INEQUALITY'])+2
    ncon=nceq+ncieq

    con_low=[]
    con_up=[]
    for i in range(0,ncieq):
	con_low.append(-2.0*pow(10.0,19))
	con_up.append(0.0)
    for i in range(0,nceq):
	con_low.append(0.0)
	con_up.append(0.0)
    g_L=array(con_low)
    g_U=array(con_up)

    nnzj=nvar*ncon
	
    nnzh= 0;
    for i in range(1,ncon+1):
  	nnzh = nnzh+i;
	#nnzh=? TODO pi0...
    
    nlp = pyipopt.create(nvar, x_L, x_U, ncon, g_L, g_U, nnzj, nnzh, eval_f, eval_grad_f, eval_g, eval_jac_g)

    # Initial guess

    if not x0: x0 = [0.0]*nvar
    
    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    k = 0
    for i, dv_scl in enumerate(dv_scales):
        for j in range(dv_size[i]):
            x0[k] =x0[k]/dv_scl;
            k = k + 1
    x0P=array(x0)  
   
    accu=1e-2 
    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    accu = accu*obj_scale
	
    if(accu<0):
	accu=accu*(-1)

    # scale accuracy
    eps = 1.0e-04

    nlp.num_option('tol', accun)
    #nlp.num_option('constr_viol_tol',0.0001)
    #nlp.num_option('compl_inf_tol',0.0001)
    #nlp.num_option('tol', 7e-4)
    nlp.num_option('constr_viol_tol',1.0e-12)
    nlp.num_option('compl_inf_tol',accun)
    nlp.int_option('max_iter', its)
    nlp.str_option('output_file', 'ipopt.out')
    nlp.str_option('hessian_approximation', 'limited-memory')

    
    	# optimizer summary
    sys.stdout.write('Ipopt parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(project.n_dv) + '\n')
    sys.stdout.write('Number of constraints: ' + str(ncon) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accun) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower bound for each independent variable: ' + str(xb_low) + '\n')
    sys.stdout.write('Upper bound for each independent variable: ' + str(xb_up) + '\n')

    # Run Optimizer
    x, zl, zu, constraint_multipliers, objfun, status = nlp.solve(x0P)
    nlp.close()
    
    constr=eval_g(x)

    def print_variable(variable_name, value):
  	for i in xrange(len(value)):
            print variable_name + "["+str(i)+"] =", value[i]
		
    #print
    #print "Solution of the primal variables, x"
    #print_variable("x", x)
    #print
    #print "Solution of the bound multipliers, z_L and z_U"
    #print_variable("z_L", zl)
    #print_variable("z_U", zu)
    #print
    #print "Solution of the constraint multipliers, lambda"
    #print_variable("lambda", constraint_multipliers)
    print
    print "Objective value"
    print "f(x*) = ", objfun
    print "Constraint value"
    print "Con = ", constr
    
    # Done
    return objfun, constr, status, project
    

    
def objnormal_f(x,project):
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
        print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(obj)
        sumobj=sumobj+weights[count]*obj[0]
    exp=sumobj/sqrt(pi)
    return exp

def dexp_f(x,project,nquad, mu,sigma):
    points, weights = hermite(nquad)
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        dobj = project.obj_df(x)
        dobj = array( dobj[0], float_)
        #print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(obj)
        sumobj=sumobj+weights[count]*dobj
    exp=sumobj/sqrt(pi)
    return exp

def var_f(x,project,nquad, mu, sigma):
    points, weights = hermite(nquad)
    varf=0
    for degree in range(1,10):
        sumobj=0
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
            obj = project.obj_f(x)[0]*1E4 
            inp1=sqrt(2)*points[count]
            if(degree==1):
                sumobj=sumobj+weights[count]*obj*inp1
            if(degree==2):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1-1)
            if(degree==3):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1-3*inp1)
            if(degree==4):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
            if(degree==5):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
            if(degree==6):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
            if(degree==7):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
            if(degree==8):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
            if(degree==9):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
        varf+=sumobj*sumobj/pi*(1./factorial(degree))
    print "variance g"+repr(varf)
    return varf

def dvar_f(x,project,nquad, mu, sigma):
    points, weights = hermite(nquad)
    varf=0
    for degree in range(1,10):
        sumobj=0
        dsumobj=0
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
            obj = project.obj_f(x)[0]*1E4 
            dobj= project.obj_df(x)
            dobj = array( dobj[0], float_)*1E4
            inp1=sqrt(2)*points[count]
            if(degree==1):
                sumobj=sumobj+weights[count]*obj*inp1
                dsumobj=dsumobj+weights[count]*dobj*inp1
            if(degree==2):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1-1)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1-1)
            if(degree==3):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1-3*inp1)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1-3*inp1)
            if(degree==4):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
            if(degree==5):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
            if(degree==6):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
            if(degree==7):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
            if(degree==8):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
            if(degree==9):
                sumobj=sumobj+weights[count]*obj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
                dsumobj=dsumobj+weights[count]*dobj*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
        varf+=2*sumobj*dsumobj/pi*(1./factorial(degree))
        #print repr(degree)+" "+repr(varf)
    return varf

def objnormal_g(x,project):     
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

def dexp_g(x,project,nquad,mu,sigma): 
    points, weights = hermite(nquad)
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        decons = project.con_dceq(x)
        diecons = project.con_dcieq(x)
        dcons = array(decons+diecons, float_)
        #print str(count)+" "+str(project.config.MACH_NUMBER)+" "+str(cons[0])
        sumobj=sumobj+weights[count]*dcons
    dexp=sumobj/sqrt(pi)
    return dexp

def var_g(x,project,nquad, mu, sigma):  
    points,weights=hermite(nquad)
    varf=0
    for degree in range(1,10):
        sumobj=0
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu 
            econs = project.con_ceq(x)
            iecons = project.con_cieq(x)
            consn = array(econs+iecons, float_)
            if(liftdrag==1):
                cons = consn[0]*1E3-constraint
            else:
	        cons= consn[0]*1E3+constraint 
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
    print "variance g: "+repr(varf)
    return varf

def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

def dvar_g(x,project,nquad, mu, sigma):  
    points,weights=hermite(nquad)
    varf=0
    for degree in range(1,10):
        sumobj=0
        dsumobj=0
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu 
            econs = project.con_ceq(x)
            iecons = project.con_cieq(x)
            consn = array(econs+iecons, float_)
            decons = project.con_dceq(x)
            diecons = project.con_dcieq(x)
            dconsn = array(decons+diecons, float_)
            if(liftdrag==1):
                cons = consn[0]*1E3-constraint
            else:
	        cons= consn[0]*1E3+constraint 
            dcons=dconsn[0]*1E3
            inp1=sqrt(2)*points[count]
            if(degree==1):
                sumobj=sumobj+weights[count]*cons*inp1
                dsumobj=dsumobj+weights[count]*dcons*inp1
            if(degree==2):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1-1)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1-1)
            if(degree==3):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1-3*inp1)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1-3*inp1)
            if(degree==4):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1*inp1-6*inp1*inp1+3)
            if(degree==5):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1*inp1*inp1-10*inp1*inp1*inp1+15*inp1)
            if(degree==6):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1*inp1*inp1*inp1-15*inp1*inp1*inp1*inp1+45*inp1*inp1-15)
            if(degree==7):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1-21*inp1*inp1*inp1*inp1*inp1+105*inp1*inp1*inp1-105*inp1)
            if(degree==8):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-28*inp1*inp1*inp1*inp1*inp1*inp1+210*inp1*inp1*inp1*inp1-420*inp1*inp1+105)
            if(degree==9):
                sumobj=sumobj+weights[count]*cons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
                dsumobj=dsumobj+weights[count]*dcons*(inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1*inp1-36*inp1*inp1*inp1*inp1*inp1*inp1*inp1+378*inp1*inp1*inp1*inp1*inp1-1260*inp1*inp1*inp1+945*inp1)
        #varf+=sumobj*sumobj/pi*(1./factorial(degree))
        varf+=2*sumobj*dsumobj/pi*(1./factorial(degree))
        #print repr(degree)+" "+repr(varf)
    return varf
