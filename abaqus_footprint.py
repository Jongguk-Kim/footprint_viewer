from abaqus import *
from abaqusConstants import *
from viewerModules import *
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

import os, sys 
def extractFootdata(step=3):

    cont = 'CONT'
    step =  3 
    odbfile = 'remote.odb'

    if os.path.isfile('abq'): 
        with open('abq') as AB: 
            lines = AB.readlines()

        cont = lines[0].strip()
        step = lines[1].strip()
        if cont =='': cont = 'CONT'
        if step =='': step = 3
        else: step = int(step)
        if lines[2].strip() =="": 
            odbfile = 'remote.odb'
        else: 
            odbfile = lines[2].strip()

    outputfilename = 'cpress.txt'

    print(odbfile)
    o1 = session.openOdb(name=odbfile)
    
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)

    leaf = dgo.LeafFromSurfaceSets(surfaceSets=(cont, ))
    session.viewports['Viewport: 1'].odbDisplay.displayGroup.replace(leaf=leaf)

    odb = session.odbs[odbfile]
    session.writeFieldReport(fileName=outputfilename, append=OFF, 
        sortItem='Node Label', odb=odb, step=step-1, frame=1, outputPosition=NODAL, 
        variable=(('COORD', NODAL), ('CPRESS', ELEMENT_NODAL), ))
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=step-1, frame=1)

    odb.close()

    

if __name__ == "__main__":

    extractFootdata(step=3)

    
    


