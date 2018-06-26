#!/usr/bin/env python 

## \file shape_optimization.py
#  \brief Python script for performing the shape optimization.
#  \author T. Economon, T. Lukaczyk, F. Palacios
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

import os, sys, shutil, copy, string
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

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
    parser.add_option("-g", "--gradient", dest="gradient", default="CONTINUOUS_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
    
    (options, args)=parser.parse_args()
    
    # process inputs
    options.partitions  = int( options.partitions )
    options.quiet       = options.quiet.upper() == 'TRUE'
    options.gradient    = options.gradient.upper()
    
    sys.stdout.write('\n-------------------------------------------------------------------------\n')
    sys.stdout.write('|    ___ _   _ ___                                                      |\n')
    sys.stdout.write('|   / __| | | |_  )   Release 4.3.0 \"Cardinal\"                          |\n')
    sys.stdout.write('|   \\__ \\ |_| |/ /                                                      |\n')
    sys.stdout.write('|   |___/\\___//___|   Aerodynamic Shape Optimization Script             |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Lead Dev.: Dr. Francisco Palacios, Francisco.D.Palacios@boeing.com|\n')
    sys.stdout.write('|                Dr. Thomas D. Economon, economon@stanford.edu          |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Developers:                                                       |\n')
    sys.stdout.write('| - Prof. Juan J. Alonso\'s group at Stanford University.                |\n')
    sys.stdout.write('| - Prof. Piero Colonna\'s group at Delft University of Technology.      |\n')
    sys.stdout.write('| - Prof. Nicolas R. Gauger\'s group at Kaiserslautern U. of Technology. |\n')
    sys.stdout.write('| - Prof. Alberto Guardone\'s group at Polytechnic University of Milan.  |\n')
    sys.stdout.write('| - Prof. Rafael Palacios\' group at Imperial College London.            |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| Copyright (C) 2012-2016 SU2, the open-source CFD code.                |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is free software; you can redistribute it and/or                  |\n')
    sys.stdout.write('| modify it under the terms of the GNU Lesser General Public            |\n')
    sys.stdout.write('| License as published by the Free Software Foundation; either          |\n')
    sys.stdout.write('| version 2.1 of the License, or (at your option) any later version.    |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is distributed in the hope that it will be useful,                |\n')
    sys.stdout.write('| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n')
    sys.stdout.write('| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n')
    sys.stdout.write('| Lesser General Public License for more details.                       |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| You should have received a copy of the GNU Lesser General Public      |\n')
    sys.stdout.write('| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')

    shape_optimization( options.filename    ,
                        options.projectname ,
                        options.partitions  ,
                        options.gradient    ,
                        options.quiet        )
    
#: main()

def shape_optimization( filename                , 
                        projectname = ''        ,
                        partitions  = 0         , 
                        gradient    = 'CONTINUOUS_ADJOINT' ,
                        quiet       = False      ):
  
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    if quiet: config.CONSOLE = 'CONCISE'
    config.GRADIENT_METHOD = gradient
    
    its         = int ( config.OPT_ITERATIONS )
    accu        = float ( config.OPT_ACCURACY )
    bound_upper = float ( config.OPT_BOUND_UPPER )
    bound_lower = float ( config.OPT_BOUND_LOWER )
    def_dv      = config.DEFINITION_DV
    n_dv        = sum(def_dv['SIZE'])
    x0          = [0.0]*n_dv # initial design
    xb_low      = [float(-0.005)]*n_dv # lower dv bound
    xb_up       = [float(0.005)]*n_dv # upper dv bound
    xb          = zip(xb_low,xb_up) # design bounds
    
    # State
    state = SU2.io.State()
    state.find_files(config)
    
    # Project
    if os.path.exists(projectname):
        project = SU2.io.load_data(projectname)
        project.config = config
    else:
        project = SU2.opt.Project(config,state)

    # Optimize

    file=open("design.txt", "r")
    line=file.readline()
    splitline=string.split(line, ",")
    for i in range(0,n_dv):
        x0[i]=string.atof(splitline[i])
    print x0 
    #x0=[0.00332610093355, 0.00335170347623, 0.00374332400682, 0.00462875711942, 0.00499992446719, -0.000385232335349, -0.00354713058659, 0.0031695003685, 0.00476152845504, 0.0048608082812, 0.00481119018361, 0.00082829744884, -0.00492758889176, -0.00496110324433, -0.00496982446174, -0.00497018471013, -0.00496556435779, -0.00494462825149, -0.00451347765238, -0.00451765940531, -0.0049609365035, -0.00495359010145, -0.00491954328652, -0.0048356882413, -0.00347857840348, -0.00357134394465, 8.22761431737e-05, 0.00151612355899, 0.00205112048417, -2.6923704913e-05, -0.00263924753369, -0.00413425590605, 0.00232274325063, 0.00493909641849, 0.00498616209446, 0.0049942556893, 0.00499494887899, 0.00494231364797]
    #x0=[0.00480829559943, 0.00477076318634, 0.00187690108304, -0.000478760998222, -0.0013415668497, -0.00104604612913, 2.77761185838e-05, 0.00135969768414, 0.00230111761532, 0.00262694728625, 0.00224540312538, 0.001856320273, 0.00222343875774, 0.00328119275018, 0.00383825479819, 0.00353510621017, 0.00243775609505, 0.00242570928677, 0.00312860675774, 0.000658736368508, -0.00472471552139, -0.00471849236417, -0.0041232180627, -0.00162819108919, -0.000359388528859, -0.000345170149444, -0.000730457905124, -0.000736499006622, 0.000168129486131, 0.0014729580001, 0.00241012457281, 0.00234616659423, 0.00178723977252, 0.00265726523811, 0.00495738763572, 0.00494541520429, 0.00112795567129, -0.000971063992412]
#    SU2.opt.STOCHLOSSES(project,config,state,x0,xb,its,accu)
    SU2.opt.STOCHASTIC(project,config,state,x0,xb,its,accu)
#    SU2.opt.LOSSES(project,config,state,x0,xb_low,xb_up,its,accu,0.64)
    #SU2.opt.CONV(project,config,state,x0,xb,its,accu)
    #SU2.opt.TESTFUNC(project,config,state,x0,xb,its,accu)
    
    # rename project file
    if projectname:
        shutil.move('project.pkl',projectname)
    
    return project

#: shape_optimization()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

