import rosetta
from rosetta import *
from rosetta import protocols
from rosetta.protocols import *
from rosetta.protocols.hybridization import *
from numpy import array as nparray 
from numpy import mean as npmean
from numpy import float64 as npfloat64
from numpy import tile as nptile
from numpy import dot as npdot
from numpy import transpose as nptranspose
from numpy.linalg import svd as nplinalgsvd
from numpy.linalg import det as nplinalgdet

def tmalign( pose, ref_pose ):
    '''
    Takes a pose (with or without a ligand) and overlays it with the ref pose
    Assumes the ligand (if any) is in the pose to be moved
    input: pose, reference_pose
    output: atommap, and tm object
    use: atommap , tm = tmalign(model,crystal)
    
    '''

    print 'Running tmalign on poses'
    print 'Starting xyz coords for res1 CA pose and reference pose'
    print pose.residue(1).xyz('CA')
    print ref_pose.residue(1).xyz('CA')
    
    ligs = False
    listlig = []    
    for i in range(1,pose.total_residue()+1):
        if pose.residue(i).is_ligand():
            listlig.append(i)

    if len(listlig) > 0:
        print 'Ligand detected, will also overlay'
        ligs = True
        copy = Pose()
        copy.assign(pose)

        offset = len(listlig)

    # construct the TMalign object
    tm = rosetta.protocols.hybridization.TMalign()
    tm.apply(pose, ref_pose)
    longest = max(pose.n_residue()+1, ref_pose.n_residue()+1)
    
    print 'TMScore = %s ' %tm.TMscore(longest)
    #print tm.TMscore(longest)
    
    # Now pull the atom mapping from tmalign
    # tmalign makes it's own alignment so we use that to do the ''partial'' threading

    # some setup for alignment2AtomMap method
    atommap =  rosetta.core.id.AtomID_Map_T_core_id_AtomID_T()
#    atommap =  rosetta.core.id.AtomID_Map_core_id_AtomID_t()
    rosetta.core.pose.initialize_atomid_map( atommap, pose )
    tm.alignment2AtomMap( pose, ref_pose, atommap )
    
    # some setup for partial thread
    aln_cutoff = rosetta.utility.vector1_Real()
#    aln_cutoff = rosetta.utility.vector1_double()
    for i in [2,1.5,1.0,.5]:
        aln_cutoff.append(i)
#    min_coverage = .2
    min_coverage = .5
    rosetta.protocols.hybridization.partial_align(pose,ref_pose, atommap, True, aln_cutoff, min_coverage)
    
    print 'Hopefully these coordinates have changed, use the PyMolMover / Observer to watch in realtime'
    print pose.residue(1).xyz('CA')
    print ref_pose.residue(1).xyz('CA')

    print "Checking for ligand, will move it too!"
    
    if ligs:
        precoord = []
        postcoord = []
        for i in range(1,pose.total_residue()+1-offset):
            try:
                pre_i = copy.residue(i).atom('CA').xyz()
                prexyz = [ pre_i.x, pre_i.y, pre_i.z ]
                precoord.append( prexyz )
                post_i = pose.residue(i).atom('CA').xyz()
                postxyz = [ post_i.x, post_i.y, post_i.z ]
                postcoord.append(postxyz)
            except:
                pass
        premat = nparray(precoord,dtype=npfloat64)
        postmat = nparray(postcoord,dtype=npfloat64)
        R,t = rigid_transform_3D( premat, postmat)
        for j in listlig:
            for k in range(1,pose.residue(j).natoms()+1):
                ligatm = pose.residue(j).atom(k)
                prevec = nparray([ ligatm.xyz().x, ligatm.xyz().y, ligatm.xyz().z ], dtype = npfloat64 )
                postvec = npdot(R,prevec) + t
                newrosettavec = rosetta.numeric.xyzVector_Real( postvec[0], postvec[1], postvec[2] )
                pose.set_xyz(rosetta.core.id.AtomID(k,j), newrosettavec)

    return (atommap,tm)


def active_site_rmsd( crystal, model ):
    crystalpdbinfo = crystal.pdb_info()
    print crystalpdbinfo.get_num_unrecognized_atoms()
    ua = crystalpdbinfo.get_unrecognized_atoms()

    activesiteresnumberset = set()
    cutoff = 8
#    cutoff = 16

    for i in ua:
        ligcoords = i.coords()
        for j in range(1,crystal.n_residue()+1):
            resj = crystal.residue( j )
            for k in resj.atoms():
                disttolig = ligcoords.distance( k.xyz() )  ## Note, type information and methods are in the .hh files in rosetta
                if disttolig < cutoff:
                    activesiteresnumberset.add(int(j) )  
    print activesiteresnumberset


    mycrystalpdbnumberset = set()
    for i in list(activesiteresnumberset):
        mycrystalpdbnumberset.add(crystalpdbinfo.pose2pdb(i) )
        print "Pose Number %s <-----> Crystal Number %s" %(i,crystalpdbinfo.pose2pdb(i))

    from sjb_util import tmalign
    atommap, tm = tmalign(model,crystal)
    for i in range(1,model.n_residue()+1):
        if model.residue(i).is_protein():
            atom_id_CA_resi = rosetta.core.id.AtomID( model.residue(i).atom_index("CA"),i)
#        print "Model CA %s  maps to Crystal CA %s " %(atom_id_CA_resi.rsd(), atommap.get(atom_id_CA_resi).rsd())
    print atommap
    model_map_to_crystal = dict()
    for i in list(activesiteresnumberset):
        try:
            if model.residue(i).is_protein():
        #print i
        #print "CA atom index = %s" %model.residue(i).atom_index("CA")
                atom_id_CA_resi = rosetta.core.id.AtomID( model.residue(i).atom_index("CA"),i)
                print "Model CA %s  maps to Crystal CA %s " %(atom_id_CA_resi.rsd(), atommap.get(atom_id_CA_resi).rsd())
                model_map_to_crystal[atom_id_CA_resi.rsd()] = atommap.get(atom_id_CA_resi).rsd()
        except:
            print "ERROR %s" %i # Dies on the metal ions
            pass


    running_activesite_rms = 0.0
    natoms = 0.0

    running_activesite_rms_lig = 0.0
    nligatoms = 0.0

    running_activesite_rms_ca = 0.0
    ncaatoms = 0.0

    for i in model_map_to_crystal.iterkeys():   ## iter over the active site residues
        if model_map_to_crystal[i] != 0:
          modelresi = model.residue(i)
          crystalresi = crystal.residue(model_map_to_crystal[i]) # Get the crystal mapped residue
          assert modelresi.name3() == crystalresi.name3(), "These are different residues~!"
    #print i
          for j in range(1,crystalresi.natoms()):     # loop over atoms in crystal
            if not crystalresi.atom_is_hydrogen( j ):  #skip the hydrogens
                atomcrystal = crystalresi.atom(j)   

            #print crystalresi.atom_name(1)
                for l in range(1,modelresi.natoms()): #loop over atoms in model
                    if not modelresi.atom_is_hydrogen( l ): #skip hydrogens
                        atommodel = modelresi.atom(l)
                    
                    
                        if (crystalresi.atom_name(j) == modelresi.atom_name(l)): # same type of atom
                        #print "Match"
                        #print "|%s|" %crystalresi.atom_name(j)
                        #print modelresi.atom_name(l)
                        
                        #print atomcrystal.type()
                        #print type( atomcrystal.type())
                        #atmtyp =  crystalresi.atom_type( atomcrystal.type() )
                        # All atom active site rmsd counters
                        #print "Match"
                        #print atommodel.xyz()
                        #print atomcrystal.xyz()
                             rms_j_l = atommodel.xyz().distance( atomcrystal.xyz() )
                   
                             running_activesite_rms += rms_j_l
                             natoms += 1.0
                        
                        # Ca active site rsmd counters here
                             if modelresi.atom_name(j) == ' CA ':
                        
                        #print "This is the CA"
                                  rms_jca_lca = atommodel.xyz().distance( atomcrystal.xyz() )
                                  running_activesite_rms_ca += rms_jca_lca
                                  ncaatoms +=1
                        
    #print "Over %s atoms in active site" %natoms
    active_rms = (running_activesite_rms)/(natoms)

    #print "Over %s CA atoms  in active site " %ncaatoms #this should be 55 - 3 metals
    ca_active_rms = (running_activesite_rms_ca)/(ncaatoms)



    import pandas as pd
    from rosetta.core.scoring.methods import EnergyMethodOptions

    sfxn = get_fa_scorefxn()  #get a scorefunction (default talaris2013)
    talaris2013_energy_methods = sfxn.energy_method_options() #have to copy the default energy methods from talaris first
    emo = EnergyMethodOptions( talaris2013_energy_methods)    #must do this to get per res hbond_bb terms in breakdown
    emo.hbond_options().decompose_bb_hb_into_pair_energies( True )  # set to true, defaults False
    sfxn.set_energy_method_options( emo ) #set the sfxn up with the energy method options
    print sfxn(model)

    score_types = []
    for i in range(1, rosetta.core.scoring.end_of_score_type_enumeration+1):
        ii = rosetta.core.scoring.ScoreType(i)
        if model.energies().weights()[ii] != 0: score_types.append(ii)
        
    listofseries = []
    for j in range(1,model.total_residue()+1):
        mydict = {}
        for i in score_types:
            myweight = model.energies().weights()[i]
            mydict['%s' %core.scoring.ScoreTypeManager.name_from_score_type(i)] = myweight*model.energies().residue_total_energies(j)[i]

        listofseries.append( pd.Series(mydict))

    df = pd.DataFrame(listofseries)
    df.index +=1 #makes index start at 1, not 0. Now, each row refers to its proper residue number (ie resi =1 -> row1)
    df = df.T # Need to do this to be able to access per resi energy with df[[1]] indexing
    print model.n_residue() # Just to check the lengths are the same
    active_energy = df[ [x for x in model_map_to_crystal.iterkeys()] ].sum().sum()
    
    return (active_energy, active_rms, ca_active_rms)


def rigid_transform_3D( MatA, MatB ):
    '''
    Pass in 2 numpy arrays to get the R and t
    '''
    assert len(MatA) == len(MatB)
    N = MatA.shape[0]
    comA = npmean(MatA, axis = 0 )
    comB = npmean(MatB, axis = 0 )
    A = MatA - nptile(comA, (N,1))
    B = MatB - nptile(comB, (N,1))
    H = npdot( nptranspose(A), B )
    U,S,V = nplinalgsvd(H)
    R = npdot(nptranspose(V),nptranspose(U))
    if nplinalgdet(R) < 0:
        V[2,:] *= -1
        R = npdot(nptranspose(V),nptranspose(U))
    t = -npdot(R,nptranspose(comA))+nptranspose(comB)

    return R,t
