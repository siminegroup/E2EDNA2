# scripts for aptamer-analyte docking
import os




def lightDock(aptamer, analyte):
    dockparams = {}
    dockparams['setup path'] = 'lightdock3_setup.py'
    dockparams['lightdock path'] = 'lightdock3.py'
    dockparams['lgd path'] = 'lgd_generate_conformations.py'

    dockparams['swarms'] = 1 # number of glowworm swarms
    dockparams['glowworms'] = 10 # number of glowworms per swarm
    dockparams['steps'] = 10 # number of steps per docking run

    # run setup
    os.system(dockparams['setup path'] + ' ' + aptamer + ' ' + analyte + ' ' + str(dockparams['swarms']) + ' ' + str(dockparams['glowworms']))

    # run docking
    os.system(dockparams['lightdock path'] + ' setup.json ' + str(dockparams['steps']))

    # generate docked structures
    for i in range(dockparams['swarms']):
        os.chdir('swarm_%d'%i)
        os.system(dockparams['lgd path'] + ' ../' + aptamer + ' ../' + analyte + ' gso_%d'%i + ' %d'%i)
        os.chdir('../')