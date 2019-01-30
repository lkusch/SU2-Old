#!/usr/bin/env python 

## \file multi_objective_optimization.py
#  \brief Python script for performing shape optimization with multiple objective (scalarization).
#  \author L. Kusch (derived from shape_optimization.py)
#  \version 6.1.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

import os, sys, shutil, copy, string
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2
from SU2 import io as su2io #needed to get sign of constraints
from SU2.util.ordered_dict import OrderedDict

# List of possible code optimizations:
#
# - Changing objectives and constraints directly in config:
#   This will need a config updating routine and an extended config comparison when calculating new results, but avoids generation of a new 
#   project for every optimization. Another nice aspect is that old results of former optimizations can be used (for restart or to avoid recalculation).
#   This is a possible way to get information on scaling etc. to change the constraint defiition  
#    def_cons = project.config['OPT_CONSTRAINT']['INEQUALITY']
#    constraints = def_cons.keys()
#    scale=[]
#    value=[]
#    operator = []
#    for i_obj,this_con in enumerate(constraints):
#        scale.append(def_cons[this_con]['SCALE'])
#        value.append(def_cons[this_con]['VALUE'])
#	 operator.append(def_cons[this_con]['SIGN']) #su2io.get_constraintSign would give sign
#   cons['INEQUALITY'][this_con] = {'SIGN':def_cons[this_con]['SIGN'], 'VALUE':def_cons[this_con]['VALUE'], 'SCALE':def_cons[this_con]['SCALE']}
#
# - Class for optimization, objective and constraints to avoid global variables and function parameters
# 
# - Initial condition for optimizer is result of nearest optimization

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-r", "--name", dest="projectname", default='',
                      help="try to restart from project file NAME", metavar="NAME")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-g", "--gradient", dest="gradient", default="DISCRETE_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_option("-o", "--optimization", dest="optimization", default="IPOPT",
                      help="OPTIMIZATION technique (SLSQP, CG, BFGS, POWELL, IPOPT)", metavar="OPTIMIZATION")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")
    parser.add_option("-m", "--multiobjective", dest="moofile", default='file.moo',
                      help="read information on multi-objective problem from file", metavar="MOOFILE")
    
    (options, args)=parser.parse_args()
    
    # process inputs
    options.partitions  = int( options.partitions )
    options.quiet       = options.quiet.upper() == 'TRUE'
    options.gradient    = options.gradient.upper()
    options.nzones      = int( options.nzones )
    
    sys.stdout.write('\n-------------------------------------------------------------------------\n')
    sys.stdout.write('|    ___ _   _ ___                                                      |\n')
    sys.stdout.write('|   / __| | | |_  )                                			|\n')
    sys.stdout.write('|   \\__ \\ |_| |/ /                                                      |\n')
    sys.stdout.write('|   |___/\\___//___|   Multi-Objective Optimization Script             	|\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')

    multi_objective_optimization( options.filename    , options.moofile,
                        options.projectname ,
                        options.partitions  ,
                        options.gradient    ,
                        options.optimization ,
                        options.quiet       ,
                        options.nzones)
    
#: main()

def calculate_pareto_point(outer, inner, obj, project, x0, xb_low, xb_up, its, accu, optimization):
    #outer defines the current objective function
    #inner defines the current step (0 stands for no additional constraint)
    global config, state

    #define objective and constraints for optimization problem
    objec = {}
    cons= {}
    cons['EQUALITY'] = OrderedDict()
    cons['INEQUALITY'] = OrderedDict()
    objec[nameObj[outer]] = {'OBJTYPE': 'DEFAULT', 'MARKER': markerObj[outer], 'SCALE':scaleObj[outer], 'VALUE': 0.0}
    if (inner!=0):
        for ihelp in range(0,numObj):
            if (ihelp != outer):
                #different treatment for maximization/minimization
                if (signObj[ihelp]==1.0):
                    cons['INEQUALITY'][nameObj[ihelp]] = {'SIGN':'<', 'VALUE':obj[ihelp], 'SCALE':scaleObjCons[ihelp]}
                else:
                    cons['INEQUALITY'][nameObj[ihelp]] = {'SIGN':'>', 'VALUE':obj[ihelp], 'SCALE':scaleObjCons[ihelp]}
    for ihelp in range (0,numIneqCons+numEqCons):
        cons[typeCons[ihelp]][nameCons[ihelp]] = {'SIGN': operatorCons[ihelp], 'VALUE':valueCons[ihelp], 'SCALE':scaleCons[ihelp]}

    config.OPT_OBJECTIVE=objec
    config.OPT_CONSTRAINT=cons

    project = SU2.opt.Project(config,state)
    
    #run optimization problem
    project=shape_optimization_moo(project, x0, xb_low, xb_up, its, accu, optimization)

    #scan optimization results for best solution
    optimum = 9999.0*signObj[outer]
    storeoptimum = 1
    for ihelp in range (0, len(project.results.FUNCTIONS[nameObj[outer]])):
        betterResult = 0
        #check for better objective function
        if(signObj[outer]==1.0):
            betterResult = (project.results.FUNCTIONS[nameObj[outer]][ihelp]<optimum)
        else:
            betterResult = (project.results.FUNCTIONS[nameObj[outer]][ihelp]>optimum)
        #check for constraint (objective functions)
        if (betterResult == 1 and inner!=0):
            for jhelp in range(0,numObj):
                if(jhelp != outer):
	            if(signObj[jhelp]==1.0):
                        betterResult=betterResult*(project.results.FUNCTIONS[nameObj[jhelp]][ihelp]-accu<obj[jhelp]) 
                    else: 
                        betterResult=betterResult*(project.results.FUNCTIONS[nameObj[jhelp]][ihelp]+accu>obj[jhelp]) 
        #check for additional constraints
        if (betterResult == 1):
            for jhelp in range(0,numIneqCons+numEqCons):
                if(operatorCons[jhelp]=="<"):
                    betterResult = betterResult*(project.results.FUNCTIONS[nameCons[jhelp]][ihelp]-accu<valueCons[jhelp])
                if(operatorCons[jhelp]==">"):
                    betterResult = betterResult*(project.results.FUNCTIONS[nameCons[jhelp]][ihelp]+accu>valueCons[jhelp])
                if(operatorCons[jhelp]=="="):
                    betterResult = betterResult*(abs(valueCons[jhelp]-project.results.FUNCTIONS[nameCons[jhelp]][ihelp])<accu)
        #store result    
        if(betterResult == 1):
            optimum = project.results.FUNCTIONS[nameObj[outer]][ihelp]
            storeoptimum = ihelp
    print "Best solution is design number: ",storeoptimum+1, ", Value: ",optimum 

    objective=[]
    for ihelp in range(0,numObj):
        objective.append(project.results.FUNCTIONS[nameObj[ihelp]][storeoptimum])
            
    print "Pareto point: ", objective+[outer,inner]

    os.system("mkdir Pareto"+str(outer)+"_"+str(inner))
#   os.system("cp DESIGNS/ Pareto"+str(outer)+"_"+str(inner)+" -rf")  #use if all designs produced during optimization shall be stored
    os.system("cp history_project.dat Pareto"+str(outer)+"_"+str(inner)+" -rf")
    os.system("cp ipopt.out Pareto"+str(outer)+"_"+str(inner)+" -rf")
    if (storeoptimum+1<10):
        os.system("mkdir Pareto"+str(outer)+"_"+str(inner)+"/DESIGN_00"+str(storeoptimum+1))
        os.system("cp DESIGNS/DSN_00"+str(storeoptimum+1)+" Pareto"+str(outer)+"_"+str(inner)+"/DESIGN_00"+str(storeoptimum+1)+" -rf")
    elif (storeoptimum+1<100):
        os.system("mkdir Pareto"+str(outer)+"_"+str(inner)+"/DESIGN_0"+str(storeoptimum+1))
        os.system("cp DESIGNS/DSN_0"+str(storeoptimum+1)+" Pareto"+str(outer)+"_"+str(inner)+"/DESIGN_0"+str(storeoptimum+1)+" -rf")
    else:
        os.system("mkdir Pareto"+str(outer)+"_"+str(inner)+"/DESIGN_"+str(storeoptimum+1))
        os.system("cp DESIGNS/DSN_"+str(storeoptimum+1)+" Pareto"+str(outer)+"_"+str(inner)+"/DESIGN_"+str(storeoptimum+1)+" -rf")

    return objective
#: calculate_pareto_point()

def loop(s, itbegin, itend, border, obj, deltaObj, project, x0, xb_low, xb_up, its, accu, optimization):
    global pareto, itCounter, firstiteration, stepNumber, paretofile
    
    #iterate only constraints and skip current objective function (s)
    if(itend-itbegin == numObj-1): 
        if(s==itend): 
	    itend=itend-1
    if(itbegin==s): 
        itbegin=itbegin+1

    #reset border for new loop
    obj[itbegin]=border[s][itbegin] 

    #find constraints for each step
    for k in range(0,stepNumber):
        #recursively find all possible combinations
	if(itbegin!=itend):
            loop(s, itbegin+1, itend, border, obj, deltaObj, project, x0, xb_low, xb_up, its, accu, optimization)
        #solve optimization problem for current combination
	else:
            if(firstiteration):
                #the first combination is already the Pareto point at the border and does not need to be calculated
                firstiteration = 0
            else:
            	store=calculate_pareto_point((s),(itCounter), obj,project, x0, xb_low, xb_up, its, accu, optimization);
		pareto.append(store+[(s),(itCounter)])
                output = ""
                for l in range(0,numObj):
                    output = output+str(store[l])+" "
                paretofile.write(output+"\n")
                itCounter=itCounter+1;
        #update constraint
        obj[itbegin]=obj[itbegin]-deltaObj[itbegin]
    return 0
#: loop()

def multi_objective_optimization( filename, moofile, 
                        projectname = ''        ,
                        partitions  = 0         , 
                        gradient    = 'CONTINUOUS_ADJOINT' ,
                        optimization = 'IPOPT'             ,
                        quiet       = False                ,
                        nzones      = 1                    ):

    global paretofile, config, state, itCounter, firstiteration, scaleObj, nameObj, signObj, scaleObjCons, markerObj, typeCons, nameCons, operatorCons, valueCons, scaleCons, numIneqCons, numEqCons, numObj, stepNumber, pareto
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int( nzones )
    if quiet: config.CONSOLE = 'CONCISE'
    config.GRADIENT_METHOD = gradient

    #read information on design space and optimization parameters
    
    its              = int ( config.OPT_ITERATIONS )                      # number of opt iterations
    bound_upper      = float ( config.OPT_BOUND_UPPER )                   # variable bound to be scaled by the line search
    bound_lower      = float ( config.OPT_BOUND_LOWER )                   # variable bound to be scaled by the line search
    relax_factor     = float ( config.OPT_RELAX_FACTOR )                  # line search scale
    gradient_factor  = float ( config.OPT_GRADIENT_FACTOR )               # objective function and gradient scale
    def_dv           = config.DEFINITION_DV                               # complete definition of the design variable
    n_dv             = sum(def_dv['SIZE'])                                # number of design variables
    accu             = float ( config.OPT_ACCURACY ) * gradient_factor    # optimizer accuracy
    x0          = [0.0]*n_dv # initial design
    xb_low           = [float(bound_lower)/float(relax_factor)]*n_dv      # lower dv bound it includes the line search acceleration factor
    xb_up            = [float(bound_upper)/float(relax_factor)]*n_dv      # upper dv bound it includes the line search acceleration fa
    xb          = list(zip(xb_low, xb_up)) # design bounds
    # State
    state = SU2.io.State()
    state.find_files(config)
    # Project
    if os.path.exists(projectname):
        project = SU2.io.load_data(projectname)
        project.config = config
    else:
        project = SU2.opt.Project(config,state)

    #additional variables
    obj=[] 	#store current objective function values (for constraints)	
    deltaObj=[] #store stepsize for scanning
    border=[] 	#store borders of Pareto front
    pareto=[]   #store all Pareto optimal points 
    objec=[]	#store solution of single-objective optimization
    calcBorder = 1

    numObj = 0
    numIneqCons = 0
    numEqCons = 0
    stepNumber = 5 #default
    nameObj = []
    signObj = []
    scaleObj = []
    scaleObjCons = []
    markerObj = []
    typeCons = []
    nameCons = []
    operatorCons = []
    valueCons = []
    scaleCons = []
    
    #read information on multi-objective problem setup
    if os.path.exists(moofile):
        #readfile = open(moofile, "r")
        with open(moofile,"r") as moo:
            for line in moo:
        #line = readfile.readline()
                if(string.find(line,"NUMOBJ")>=0):
                    words =line.split("=")
                    numObj = string.atoi(words[1])
                if(string.find(line,"NUMSTEPS")>=0):
                    words = line.split("=")
                    stepNumber = string.atoi(words[1])
                if(string.find(line,"BORDER")>=0):
                    calcBorder = 0
                    words = line.split("=")
                    values = words[1].split(",")
                    for i in range(0,numObj):
                        helper = []
                        for j in range(0,numObj):
                            helper.append(string.atof(values[i*numObj+j]))
                        border.append(helper)
                if(string.find(line, "OBJEC")>=0):
                    words = line.split("=")
                    names = words[1].split(",")
                    for i in range(0,numObj):
                        nameObj.append(names[i*5].strip())
                        signObj.append(string.atof(names[i*5+1]))
                        scaleObj.append(string.atof(names[i*5+2]))
                        scaleObjCons.append(string.atof(names[i*5+3]))
                        markerObj.append(names[i*5+4])
                if(string.find(line,"NUMINCONS")>=0):
                    words = line.split("=")
                    numIneqCons = string.atoi(words[1])
                if(string.find(line,"NUMCONS")>=0):
                    words = line.split("=")
                    numEqCons = string.atoi(words[1])
                if(string.find(line, "INEQ")>=0):
                    words = line.split("=")
                    names = words[1].split(",")
                    for i in range(0,numIneqCons):
                        typeCons.append("INEQUALITY")
                        nameCons.append(names[i*4].strip())
                        operatorCons.append(names[i*4+1].strip())
                        valueCons.append(string.atof(names[i*4+2]))
                        scaleCons.append(string.atof(names[i*4+3]))
                if(string.find(line, "EQUAL")>=0):
                    words = line.split("=")
                    names = words[1].split(",")
                    for i in range(0,numEqCons):
                        typeCons.append("EQUALITY")
                        nameCons.append(names[i*4].strip())
                        operatorCons.append(names[i*4+1].strip())
                        valueCons.append(string.atof(names[i*4+2]))
                        scaleCons.append(string.atof(names[i*4+3]))
    else:
        print "Error: No information on MOO problem found!" 
        return project

    for i in range(0,numObj):
        obj.append(0.0)		
        deltaObj.append(0.0)
    
    #output
    paretofile=open("ParetoPoints",'w')
    
    if calcBorder:
        #calculate SOO problems to find border of Pareto front
        for i in range(0,numObj):
            objec=calculate_pareto_point(i,0, obj,project, x0, xb_low, xb_up, its, accu, optimization)
            border.append(objec)

    for i in range(0,numObj):
        pareto.append(border[i]+[i,0])
        output = ""
        for j in range(0,numObj):
            output = output+str(border[i][j])+" "
        paretofile.write(output+"\n")

    #perform scanning of front
    if (numObj>1):
	for i in range(0,numObj):
            firstiteration=1
            itCounter=1
	    for j in range(0,numObj):
		deltaObj[j]=float(border[i][j]-border[j][j])/float(stepNumber)
            	obj[j]=border[i][j]
            print "Delta: ", deltaObj
	    loop(i,0,(numObj-1), border, obj, deltaObj,project, x0, xb_low, xb_up, its, accu, optimization)

    paretofile.close()
    
    # rename project file
    if projectname:
        shutil.move('project.pkl',projectname)
    
    return project
#: multi_objective_optimization()

def shape_optimization_moo(project, x0, xb_low, xb_up, its, accu, optimization):
    sys.stdout.write('DEBL shape_optimization()\n')
    function_call = "SU2.opt."+optimization
    exitcode, project=eval(function_call)(project,x0,xb_low,xb_up,its,accu)
    return project
#: shape_optimization_moo

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

