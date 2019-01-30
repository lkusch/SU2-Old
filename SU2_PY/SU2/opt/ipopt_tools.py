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


# -------------------------------------------------------------------
#  CoinIpopt (BFGS)
# -------------------------------------------------------------------
global project,ncon

#TODO project als globale Variable / Klassenvariable?    
def eval_f(x,user_data=None):
    """ obj = eval_f(x)
        
        Objective Function
        SU2 Project interface to Ipopt
        
        su2:         minimize f(x), list[nobj]
        ipopt: 	     minimize f(x), float
    """
    assert len(x)==project.n_dv
  
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

    econs = project.con_ceq(x)
    iecons = project.con_cieq(x)
    cons = array(econs+iecons, float_)   
    return cons

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

    	decons = project.con_dceq(x)
    	diecons = project.con_dcieq(x)

    	dcons= array(decons+diecons, float_)
	values=[]
	for j in range(0, ncon):
	    for k in range(0, project.n_dv):
 		values.append(dcons[j][k]);
    	return array(values)

def eval_h(x, lagrange, obj_factor, flag, user_data = None):
    if flag:
        return 0
    else:
        return 0

def apply_new(x):
    return True

def ipopt_run(projectx, x0=None,xb_low=None, xb_up=None, its=100,accu=1e-10,grads=True):
    """ result = ipopt_run(project,x0=[],xb_low=[], xb_up=[], its=100,accu=1e-10)
    
        Runs the BFGS method of Ipopt with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb_low,xb_up - optional, design variable bounds
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10
        
        Outputs:
           result - the outputs from Ipopt
     """
    global project,ncon
    import pyipopt
    project =projectx
    if x0 is None: x0 = []
    if xb_low is None: xb_low = []
    if xb_up is None: xb_up = []
    
    # number of design variables
    n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    project.n_dv = n_dv
	
    x_L=array(xb_low)
    x_U=array(xb_up)

    nceq=len(project.config['OPT_CONSTRAINT']['EQUALITY'])
    ncieq=len(project.config['OPT_CONSTRAINT']['INEQUALITY'])
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

    nnzj=n_dv*ncon
	
    nnzh= 0;
    for i in range(1,ncon+1):
  	nnzh = nnzh+i;
    
    nlp = pyipopt.create(n_dv, x_L, x_U, ncon, g_L, g_U, nnzj, nnzh, eval_f, eval_grad_f, eval_g, eval_jac_g)

    # Initial guess

    if not x0: x0 = [0.0]*n_dv
    
    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]
    x0P=array(x0)  
   
    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = obj[obj.keys()[0]]['SCALE']
    accu = accu*obj_scale

    # scale accuracy
    eps = 1.0e-04

    nlp.num_option('tol', accu)
    nlp.num_option('constr_viol_tol',accu)
    nlp.num_option('compl_inf_tol',accu)
    nlp.int_option('max_iter', its)
    nlp.str_option('output_file', 'ipopt.out')
    nlp.str_option('hessian_approximation', 'limited-memory')

    
    # optimizer summary
    sys.stdout.write('Ipopt parameters:\n')
    sys.stdout.write('Number of design variables: ' + str(project.n_dv) + '\n')
    sys.stdout.write('Number of constraints: ' + str(ncon) + '\n')
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower bound for each independent variable: ' + str(xb_low) + '\n')
    sys.stdout.write('Upper bound for each independent variable: ' + str(xb_up) + '\n')

    # Run Optimizer
    x, zl, zu, constraint_multipliers, objfun, status = nlp.solve(x0P)
    nlp.close()
    
    # Done
    return status, project
    



 
