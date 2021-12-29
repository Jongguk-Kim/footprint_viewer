from abaqus import *
from abaqusConstants import *
from viewerModules import *
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

import os, sys 
import numpy as np 
def extractFootdata(step=3):

    cont = 'CONT'
    step =  3 
    dir = 'X+'
    odbfile = 'remote.odb'

    if os.path.isfile('abq'): 
        with open('abq') as AB: 
            lines = AB.readlines()

        cont = lines[0].strip()
        step = lines[1].strip()
        dir = lines[2].strip()
        if cont =='': cont = 'CONT'
        if step =='': step = 3
        if dir =='': dir = 'X+'
        else: step = int(step)
        if lines[3].strip() =="": 
            odbfile = 'remote.odb'
        else: 
            odbfile = lines[3].strip()

    outputfilename = 'cpress.txt'
    fp = open(outputfilename, 'w')
    
    # fp.write(odbfile)
    # fp.write("\n")
    # fp.write (cont)
    # fp.write("\n")
    
    o1 = session.openOdb(name=odbfile)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    leaf = dgo.LeafFromSurfaceSets(surfaceSets=(cont, ))
    session.viewports['Viewport: 1'].odbDisplay.displayGroup.replace(leaf=leaf)
    odb = session.odbs[odbfile]
    targetStep=odb.steps['Step-%d'%(step)]
    frames = len(targetStep.frames) 

    try: 
        session.writeFieldReport(fileName=outputfilename, append=OFF, 
            sortItem='Node Label', odb=odb, step=step-1, frame=frames-1, outputPosition=NODAL,
            variable=(('COORD', NODAL), ('CPRESS', ELEMENT_NODAL), ))
        ## NODE COORDINATES :  646     499.014E-03     477.909E-03    -143.589E-03     2.13068E-06
        ## PRESSURE  :    Element Label           Facet      Node Label          CPRESS
        #                     302947               5          402715     348.377E+03

    except: 
        print ("*DYNAMIC SIMULATION")
        nset = odb.rootAssembly.instances['PART-1-1'].nodes
        tfp = open('original_coordinates.tmp', 'w')
        for nd in nset: 
            if nd.label > 10**7: 
                tfp.write("%d, %.7E, %.7E, %.7E\n"%(nd.label, nd.coordinates[0], nd.coordinates[1], nd.coordinates[2]))
        tfp.close()
    
        session.writeFieldReport(fileName=outputfilename, append=ON, 
            sortItem='Node Label', odb=odb, step=step-1, frame=frames-1, outputPosition=NODAL,
            variable=( ('U', NODAL), ('CPRESS', ELEMENT_NODAL), )) ## dynamic 


    odb.close()
    

if __name__ == "__main__":

    extractFootdata(step=3)

    
    


