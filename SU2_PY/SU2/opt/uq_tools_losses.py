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
    """ obj = eval_f(x)
        
        Objective Function
        SU2 Project interface to Ipopt
        
        su2:         minimize f(x), list[nobj]
        ipopt: 	     minimize f(x), float
    """
    assert len(x)==project.n_dv
  
    project.config.MACH_NUMBER=0.8
    obj = project.obj_f(x)
    
    obj = obj[0]
    return obj

def eval_grad_f(x,user_data=None):
    """ dobj = eval_grad_f(x)
        
        Objective Function Gradients
        SU2 Project interface to Ipopt
        
        su2:         df(x), list[nobj x dim]
        ipopt: 	     df(x), ndarray[dim]
    """    
    assert len(x)==project.n_dv

    project.config.MACH_NUMBER=0.8
    dobj = project.obj_df(x)
    
    dobj = array( dobj[0], float_)
    return dobj

def eval_g(x,user_data=None):
    """ cons = eval_g(x)
        
        Equality Constraint and Inequality Constraint Functions
        SU2 Project interface to Ipopt
        
        su2:         ceq(x),cieq(x) =/< 0.0, list[nceq]
        ipopt:       cons(x) = 0.0, ndarray[nceq]
    """
    assert len(x)==project.n_dv
    
    #calculate fmean
    project.config.MACH_NUMBER=0.8
    obj=project.obj_f(x)
    objmean=obj[0]
    
    #calculate E(max(0,f-fmean))
    sumobj=0
    for count in range(0,nquad):
        project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
        obj=project.obj_f(x)
        #print "obj"
        #print obj[0]
        #print objmean
        objnew = max(0,obj[0]-objmean)
        print "f: "+repr(obj[0])+" "+repr(objmean)+" "+repr(objnew)
        sumobj=sumobj+weights[count]*objnew
    sumobj=sumobj/(sqrt(pi)) 
    #obj_scale = obj[obj.keys()[0]]['SCALE']
    print "sumobj: "+repr(sumobj)
    print "objmean: "+repr(objmean)
    sumobj=sumobj/abs(objmean) #divide by objmean?

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
    if(liftdrag==1):
    	gmean=gmean-1E-3*constraint
    else:    
    	gmean=gmean+1E-3*constraint
    #cons_scale = cons[cons.keys()[0]]['SCALE']
    #print "gmean"+repr(gmean)
    print "sumcons: "+repr(sumcons)
    sumcons=sumcons/abs(gmean)
    print "gmean: "+repr(gmean)
    consnew=((sumobj+sumcons)-0.15)*1E-3
    consend=append(consres, consnew)
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
        #calculate dfmean
        project.config.MACH_NUMBER=0.8
        dobj=project.obj_df(x)
        dobjmean=array(dobj[0], float_)
        #calculate fmean
        obj=project.obj_f(x)
        objmean=obj[0]
        #calculate derivative
	sumdobj=0
    	for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
            dobj = project.obj_df(x)
            obj=project.obj_f(x)
            if(obj[0]-objmean>0):
                dobjnew=array(dobj[0],float_)-dobjmean
            else:
                dobjnew=0.0
            sumdobj=sumdobj+weights[count]*dobjnew
        sumdobj=sumdobj/sqrt(pi)
        sumdobj=sumdobj/(objmean)
        
        #calculate dgmean
	project.config.MACH_NUMBER=0.8	
        decons = project.con_dceq(x)
    	diecons = project.con_dcieq(x)
    	dconsres = array(decons+diecons, float_)
	dgmean=dconsres[0] 
        #calculate gmean    
        econs = project.con_ceq(x)
        iecons=project.con_cieq(x)
        consres = array(econs+iecons, float_)
        gmean=consres[0]
        #calculate derivative
        sumdcons=0 
        for count in range(0,nquad):
            project.config.MACH_NUMBER=sqrt(2)*sigma*points[count]+mu
            decons = project.con_dceq(x)
            diecons = project.con_dcieq(x)
            dcons = array(decons+diecons, float_)
            econs = project.con_ceq(x)
            iecons=project.con_cieq(x)
            cons = array(econs+iecons, float_)
            if(cons[0]-gmean>0):
                dconsnew = dcons[0]-dgmean
            else:
		dconsnew=0.0
            sumdcons=sumdcons+weights[count]*dconsnew
        sumdcons=sumdcons/sqrt(pi)
        if(liftdrag==1):
    	    gmean=gmean-1E-3*constraint
        else:    
    	    gmean=gmean+1E-3*constraint
        #cons_scale = cons[cons.keys()[0]]['SCALE']
        sumdcons=sumdcons/abs(gmean)
        dconsnew=(sumdobj+sumdcons)*1E-3

        dconsend=append(dconsres, dconsnew)
	values=[]
	for j in range(0, ncon):
	    for k in range(0, project.n_dv):
 		values.append(dconsend[j*project.n_dv+k])
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
    
    nquad=4
    sigma=sqrt(0.0001)
    mu=0.8
    points, weights = hermite(nquad)
    print repr(points)
    print repr(weights)
    liftdrag=1

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
    ncieq=len(project.config['OPT_CONSTRAINT']['INEQUALITY'])+1
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
    nlp.int_option('max_iter', 100)
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
    



 
