from pylinus_2_4 import linusSim, pdbio, construct, linusHbdSat, linusScore, ran #only need ran if need gen seeds
from pylinus_2_4.Atom3d import ztox,angle,torsion,commit_internal_coords, restore_internal_coords
from pylinus_2_4.bumpcheck import bump_check, bump_check_intra
from math import sqrt, pi, sin, cos, acos, exp
from random import randint,choice,uniform
import sys, os, string
#import Numeric
from Numeric import *

DEGRAD = pi/180.0

#do i need to recalc Ome if ome = 185? (i.e. -175)
## should these energies be summed?

try:
    from linusc_24 import asa_evaluate
except:
    print 'Accessible surface area calculations disabled'
    asa_evaluate = lambda a, b, c, d, e, f, g, h, i, j: 0.0

CHI_VALUES = {
    'ARG': (((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0),(  170.0,  190.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0),(  -95.0,  -75.0)),
             ((  -77.0,  -57.0),( -177.0, -157.0),(  -75.0,  -55.0),(  -95.0,  -75.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(   55.0,   75.0),( -185.0, -165.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  -75.0,  -55.0),(  165.0,  185.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0),(   75.0,   95.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(   55.0,   75.0),(   75.0,   95.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(   55.0,   75.0),( -185.0, -165.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0),(  -95.0,  -75.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  -75.0,  -55.0),(  -95.0,  -75.0)),
             ((   52.0,   72.0),(  170.0,  190.0),(  170.0,  190.0),(   75.0,   95.0)),
             ((   52.0,   72.0),(  170.0,  190.0),(  170.0,  190.0),(  170.0,  190.0)),
             ((   52.0,   72.0),(  170.0,  190.0),(  170.0,  190.0),(  -95.0,  -75.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(   55.0,   75.0),(   75.0,   95.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  -75.0,  -55.0),(   95.0,  115.0)),
             ((  -72.0,  -52.0),(  -78.0,  -58.0),(  170.0,  190.0),(  170.0,  190.0)),
             ((  -72.0,  -52.0),(  -78.0,  -58.0),(  170.0,  190.0),(  -95.0,  -75.0)),
             ((  -72.0,  -52.0),(  -78.0,  -58.0),(  -75.0,  -55.0),(  -95.0,  -75.0)),
             (( -187.0, -167.0),(   55.0,   75.0),(  170.0,  190.0),(   75.0,   95.0)),
             (( -187.0, -167.0),(   55.0,   75.0),(  170.0,  190.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0),(   75.0,   95.0)),
           ),

    'LYS': ( ((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0),(  170.0,  190.0)),
             ((  -72.0,  -52.0),(  -78.0,  -58.0),(  170.0,  190.0),(  170.0,  190.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0),(  -75.0,  -55.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0),(   55.0,   75.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(   58.0,   78.0),(  170.0,  190.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0),(   55.0,   75.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  -78.0,  -58.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(   58.0,   78.0),(  170.0,  190.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0),(  -75.0,  -55.0)),
             ((   52.0,   72.0),(  170.0,  190.0),(  170.0,  190.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(   58.0,   78.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  -78.0,  -58.0),(  170.0,  190.0)),
           ),

    'MET': ( ((  -75.0,  -55.0),(  -75.0,  -55.0),(  -80.0,  -60.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(   65.0,   85.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  -85.0,  -65.0)),
             ((  -77.0,  -57.0),(  170.0,  190.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  -85.0,  -65.0)),
             (( -187.0, -167.0),(   55.0,   75.0),(   65.0,   85.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(   65.0,   85.0)),
             ((   52.0,   72.0),(  170.0,  190.0),(  -85.0,  -65.0)),
             ((  -75.0,  -55.0),(  -75.0,  -55.0),(   93.0,  113.0)),
             (( -187.0, -167.0),(  170.0,  190.0),(  170.0,  190.0)),
             ((   52.0,   72.0),(  170.0,  190.0),(   65.0,   85.0)),
             ((  -75.0,  -55.0),(  -75.0,  -55.0),(  170.0,  190.0)),
             (( -187.0, -167.0),(   55.0,   75.0),(  170.0,  190.0)),
           ),

    'GLN': (((  -77.0,  -57.0),(  170.0,  190.0),(  -35.0,  -15.0)),
            (( -187.0, -167.0),(  170.0,  190.0),(  -10.0,   10.0)),
            ((  -75.0,  -55.0),(  -75.0,  -55.0),(  -50.0,  -30.0)),
            (( -187.0, -167.0),(   55.0,   75.0),(   50.0,   70.0)),
            ((   52.0,   72.0),(  170.0,  190.0),(   10.0,   30.0)),
            ((  -75.0,  -55.0),(   75.0,   95.0),(  -10.0,   10.0)),
            ((  -75.0,  -55.0),(  -75.0,  -55.0),(   90.0,  110.0)),
            ((   60.0,   80.0),(  -85.0,  -65.0),(  -10.0,   10.0)),
            (( -187.0, -167.0),(   55.0,   75.0),( -110.0,  -90.0)),
           ),

    'GLU': (((  -77.0,  -57.0),(  170.0,  190.0),(  -20.0,    0.0)),
            (( -187.0, -167.0),(  170.0,  190.0),(  -10.0,   10.0)),
            ((  -75.0,  -55.0),(  -75.0,  -55.0),(  -50.0,  -30.0)),
            ((  -75.0,  -55.0),(   75.0,   95.0),(  -10.0,   10.0)),
            (( -187.0, -167.0),(   55.0,   75.0),(    0.0,   20.0)),
            ((   52.0,   72.0),(  170.0,  190.0),(  -30.0,  -10.0)),
            ((   60.0,   80.0),(  -90.0,  -70.0),(  -10.0,   10.0)),
           ),

    'PHE': (((  -75.0,  -55.0),(  -95.0,  -75.0)),
            (( -187.0, -167.0),(   70.0,   90.0)),
            ((   52.0,   72.0),(   80.0,  100.0)),
            ((  -75.0,  -55.0),(  -40.0,  -20.0)),
           ),

    'TYR': (((  -75.0,  -55.0),(  -95.0,  -75.0)),
            (( -187.0, -167.0),(   70.0,   90.0)),
            ((   52.0,   72.0),(   80.0,  100.0)),
            ((  -75.0,  -55.0),(  -40.0,  -20.0)),
           ),

    'ASP': (((  -80.0,  -60.0),(  -25.0,   -5.0)),
            (( -187.0, -167.0),(  -10.0,   10.0)),
            ((   52.0,   72.0),(  -20.0,    0.0)),
            ((   52.0,   72.0),(   20.0,   40.0)),
            (( -187.0, -167.0),(   55.0,   75.0)),
           ),

    'TRP': (((  -75.0,  -55.0),(   85.0,  105.0)),
            (( -187.0, -167.0),(   80.0,  100.0)),
            (( -187.0, -167.0),( -115.0,  -95.0)),
            ((   52.0,   72.0),( -100.0,  -80.0)),
            ((  -75.0,  -55.0),(  -15.0,    5.0)),
            ((   52.0,   72.0),(   80.0,  100.0)),
            ((  -75.0,  -55.0),( -100.0,  -80.0)),
           ),

    'ASN': (((  -75.0,  -55.0),(  -30.0,  -10.0)),
            (( -187.0, -167.0),(   20.0,   40.0)),
            (( -184.0, -164.0),(  -30.0,  -10.0)),
            ((   52.0,   72.0),(   20.0,   40.0)),
            ((  -75.0,  -55.0),(  -85.0,  -65.0)),
            ((   52.0,   72.0),(  -20.0,    0.0)),
            ((  -75.0,  -55.0),(  110.0,  130.0)),
           ),

    'HIS': (((  -75.0,  -55.0),(  -80.0,  -60.0)),
            (( -187.0, -167.0),(   50.0,   70.0)),
            ((  -75.0,  -55.0),(   70.0,   90.0)),
            (( -187.0, -167.0),(  -90.0,  -70.0)),
            ((   52.0,   72.0),(  -85.0,  -65.0)),
            ((  -75.0,  -55.0),(  155.0,  175.0)),
            (( -187.0, -167.0),( -175.0, -155.0)),
            ((   52.0,   72.0),(   70.0,   90.0)),
           ),

    'CYS': (((  -75.0,  -55.0),),
            (( -187.0, -167.0),),
            ((   52.0,   72.0),),
           ),

    'SER': (((   52.0,   72.0),),
            ((  -75.0,  -55.0),),
            (( -187.0, -167.0),),
           ),

    'THR': (((   52.0,   72.0),),
           ((  -75.0,  -55.0),),
           (( -185.0, -165.0),),
           ),

    'VAL': (((  165.0,  185.0),),
           ((  -70.0,  -50.0),),
           ((   53.0,   73.0),),
           ),

    'ILE': (((  -75.0,  -55.0),(  160.0,  180.0)),
           ((  -67.0,  -47.0),(  -70.0,  -50.0)),
           ((   52.0,   72.0),(  160.0,  180.0)),
           (( -187.0, -167.0),(  155.0,  175.0)),
           (( -187.0, -167.0),(   56.0,   76.0)),
           ),

    'LEU': (((  -75.0,  -55.0),(  165.0,  185.0)),
           (( -187.0, -167.0),(   55.0,   75.0)),
           ((  -95.0,  -75.0),(   55.0,   75.0)),
           (( -182.0, -162.0),(  135.0,  155.0)),
           ),
    'ALA': (1,),
    'GLY': (1,),
    'PRO': (1,),
    
}

def gen_seeds(): #do i need to use gen seeds? 
    from random import uniform
    s1 = long(uniform(1.0, ran.M1))
    s2 = long(uniform(1.0, ran.M1))
    s3 = long(uniform(1.0, ran.M1))
    s4 = long(uniform(1.0, ran.M2))
    s5 = long(uniform(1.0, ran.M2))
    s6 = long(uniform(1.0, ran.M2))
    return s1, s2, s3, s4, s5, s6


def set_sy(sy):
    sy = sy - 180.0
    if ((sy > 180.0) and (sy < 900.0)):
        sy = sy - 360.0
    elif sy < -180.0:
        sy = sy + 360.0
    return sy

def set_conf(line,prot):
    ##line format: identnum, phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3,  phi4, psi4, ome4, energy, hexostate, dh1, dh2, dh3 
    prot.phi_atoms[1].tp_torsion =  float(line[1])
    prot.psi_atoms[1].tp_torsion = set_sy(float(line[2]))
    prot.omega_atoms[1].tp_torsion = float(line[3])
    
    prot.phi_atoms[2].tp_torsion =  float(line[4])
    prot.psi_atoms[2].tp_torsion = set_sy(float(line[5]))
    prot.omega_atoms[2].tp_torsion = float(line[6])

    prot.phi_atoms[3].tp_torsion = float(line[7])
    prot.psi_atoms[3].tp_torsion = set_sy(float(line[8]))
    prot.omega_atoms[3].tp_torsion = float(line[9])

    prot.phi_atoms[4].tp_torsion = float(line[10])
    prot.psi_atoms[4].tp_torsion = set_sy(float(line[11]))
    prot.omega_atoms[4].tp_torsion = float(line[12])

def write_pdb(prot, name, title):
    #write file
    pfile = open(name, 'wb')
    pdbio.write_pdb(prot, pfile, title)
    pfile.close()
    
def calc_chasa_solv2(mol, atoms, fai, residue_names, numint_loos, numvirt_loos, bbtot, numbb,ext_atoms, solv_list):
    """Return the CHASA solvation energy for use in simulations
    sim.set_chasa_parameters(use_chasa=1)
    # if chasa is used must use following
    sim.set_hbond_parameters(use_hbond=1,
#                            hbond_score_short=2.5,
#                            hbond_score_long=2.5,
                             hbond_winmax=NUMRES,
                             use_sidechain_hbond=1)

    """
    #atoms = mol.atoms
    #fai = mol.residue_first_atom_indices
    #residue_names = mol.residue_names
    minres, maxres = construct.get_res_extents(mol)

    use_ext = 0
    ext_coords = None
    use_data = 1
    data = zeros(len(atoms), 'd')
    
    if ext_atoms:
        ext_coords = []
        map(lambda x: map(lambda y: ext_coords.append(y), x), ext_atoms)
        ext_coords = array(ext_coords, 'd')
        use_ext = len(ext_atoms)

    flags = construct.make_asa_list(mol)
    probe = 1.4
    ndiv = 3
    ext_radius = 1.4
    tot_asa =  asa_evaluate(atoms, data, ext_coords, flags, probe,
                            ext_radius, use_data, use_ext, ndiv)
    p_solv_nrg = 0.0
    ap_solv_nrg = 0.0
    Gamma_p = 3.0/5.0
    Gamma_hb_oxy = 0.6
    Gamma_ap = 0.03
    CHASA = 0.0
    for i in xrange(minres,maxres):
        rname = residue_names[i]
        start = fai[i]
        end = fai[i+1]
        occ = 0.0
        for j in range(start, end):
            atom = atoms[j]
            residue_num = int(mol.res_pdb_number[atom.resnum])
            if atom.name == ' N  ':
                if solv_list[i][0][2] > 0:
                    p_solv_nrg = p_solv_nrg - (Gamma_p *(solv_list[i][0][2]))
#               elif solv_list[i][0][2] < 0:
#                   p_solv_nrg = p_solv_nrg + non_hbd_score

            elif atom.name == ' O  ':
                if solv_list[i][1][2] > 0:
                    if solv_list[i][1][3] == 0:
                        p_solv_nrg = p_solv_nrg - (Gamma_p *(solv_list[i][1][2]))
                    elif solv_list[i][1][3] > 0:
                        p_solv_nrg = p_solv_nrg - (Gamma_hb_oxy)
#               elif solv_list[i][1][2] < 0:
#                   p_solv_nrg = p_solv_nrg + non_hbd_score

            elif 'C' in atom.name:
                ap_solv_nrg = ap_solv_nrg + (Gamma_ap * data[j])
#               CHASA = CHASA + data[j]

    tot_solv_nrg = ap_solv_nrg + p_solv_nrg
#   print ap_solv_nrg,  p_solv_nrg

    return tot_solv_nrg


def run(pdbfile, allowfile):
    print "Started"
    gen = apply(ran.generator, gen_seeds())   ##do i need to use gen_seeds here?
    prot = linusSim.Linus(pdbfile, gen).protein
    numres = prot.num_residues
    numatm = prot.num_atoms
    atoms = prot.atoms
    chi_atoms = prot.chi_atoms
    rn = prot.residue_names
    fai = prot.residue_first_atom_indices
    r_ids = range(0,numres)
    hbparms = linusSim.Linus(pdbfile, None).hbdpar
    hbparms['use_hbond'] = 1
    hbparms['use_sidechain_hbond'] = 1
    wmin = 2
    wmax = numres-1 #wmax=prot.num_residues-2
    BETA=2.0
    
    hblist = construct.make_hbond_list(prot, wmin, wmax, hbparms)    
    #qlist = construct.make_coul_list(prot)   #not using coulombic list
    vlist = construct.make_vdw_list(prot)

    soucedir='/raid/npanasik/PolyAboxes_random/'
    #soucedir='/scratch/npanasik/all_atom_sim/box_of_tetra_alas/randomization/'
    #soucedir=''

    afile = open(soucedir+allowfile, 'rb')
    ofile = open(sys.argv[1][:-4]+"_"+sys.argv[2][7:]+'_succesful_decorations', 'wb')
    ofile2 = open(sys.argv[1][:-4]+"_"+sys.argv[2][7:]+'_succesful_decorations_summary', 'wb')
    print "side chain rotamer set sizes", rn[1], len(CHI_VALUES[rn[1]]), rn[2], len(CHI_VALUES[rn[2]]),  rn[3], len(CHI_VALUES[rn[3]]),rn[4], len(CHI_VALUES[rn[4]]), "  Total possible combinations:",len(CHI_VALUES[rn[1]])* len(CHI_VALUES[rn[2]])*len(CHI_VALUES[rn[3]])*len(CHI_VALUES[rn[4]])
    
    accpt = 0
    num=0
    all_energy_for_set=0
    all_Henergy_for_set=0
    total_attmp=0

    

    #for writing outfile lines in batches rather than 1 at a time (for speed)
    outnum=0
    outlist=[]

    ##for debugging
    failedintra=0
    bump=0
    set=1
    startnum=0
    line2="default"

    ok=1
    nextline = afile.readline
        
    line1 = nextline()
    line=string.split(line1)
    set_conf(line, prot)
    ljnum=0
    no_bump=0
    while ok:
        chi_acpt=0
        total_e=0
        total_hbonde=0
        lowest_e=100000
        num=num+1 #keep track of total number of backbone confs tried
        for chi1 in CHI_VALUES[rn[1]]:      
            for chi2 in CHI_VALUES[rn[2]]:                
                for chi3 in CHI_VALUES[rn[3]]:
                    for chi4 in CHI_VALUES[rn[4]]:
                        
                        chi_atom = chi_atoms[1]
                        if chi_atom:
                            pos = 0
                            for atom in chi_atom:
                                chi = apply(uniform, chi1[pos])
                                if chi > 180.0: chi = chi - 360.0
                                elif chi < -180.0: chi = chi + 360.0
                                atom.tp_torsion = chi
                                pos = pos + 1

                        chi_atom = chi_atoms[2]
                        if chi_atom:
                            pos = 0
                            for atom in chi_atom:
                                chi = apply(uniform, chi2[pos])
                                if chi > 180.0: chi = chi - 360.0
                                elif chi < -180.0: chi = chi + 360.0
                                atom.tp_torsion = chi
                                pos = pos + 1

                        chi_atom = chi_atoms[3]
                        if chi_atom:
                            pos = 0
                            for atom in chi_atom:
                                chi = apply(uniform, chi3[pos])
                                if chi > 180.0: chi = chi - 360.0
                                elif chi < -180.0: chi = chi + 360.0
                                atom.tp_torsion = chi
                                pos = pos + 1


                        chi_atom = chi_atoms[4]
                        if chi_atom:
                            pos = 0
                            for atom in chi_atom:
                                chi = apply(uniform, chi4[pos])
                                if chi > 180.0: chi = chi - 360.0
                                elif chi < -180.0: chi = chi + 360.0
                                atom.tp_torsion = chi
                                pos = pos + 1
                        total_attmp=total_attmp+1
                        
                        ztox(atoms, 0, numatm) # convert for bumpcheck and ecalc
                        int_ok = 1
                        for r_id in r_ids:  #nescessary when have side chains
                            if bump_check_intra(atoms, fai, r_id):
                                int_ok = 0
                                break
                        if int_ok and not bump_check(atoms, fai, 1, numres-2, 2, numres-1): ##in linusc_24.c
                            numint, numvirt, bbtot, numbb, wat_list, cntmult_list = linusHbdSat.mk_virt_bb_cntmult_loosHbd(prot,hblist) ##not in linusc_24.c
                            hbs = numbb
                            if hbs > 7: #max should be 8
                                hbondE = linusScore.hbond_score(hblist)
                                solvE = calc_chasa_solv2(prot, atoms, fai, rn, numint, numvirt, bbtot, numbb, wat_list, cntmult_list)  ##not in linusc_24.c
                                LJE = linusScore.LJ_score(prot,vlist)            ##not in linusc_24.c
                                energyV = hbondE +solvE  + LJE  #+ coulE      #find energy for this single rotamer conf
                                energyH = exp( -1*hbondE*BETA)
                                energy = exp( -1*energyV*BETA)
                                if lowest_e < energy:
                                    lowest_e = energy
                                total_hbonde=total_hbonde+energyH
                                total_e=total_e+energy   #add energy of this rotamer conf to all rotamer confs of this backbone set
                                accpt = accpt + 1  #number of acceptances for rotamer confs of all backbone sets (total numstates)
                                chi_acpt=chi_acpt+1  #keep track of number of acceptances for rotamer confs of this backbone set (numstates)


        outlist.append(line1[:-1]+" "+str(chi_acpt)+" "+str(total_e)+" "+str(total_hbonde)+" "+str(lowest_e)+" \n")#need to have heat capacity, mesostate?
        all_energy_for_set=all_energy_for_set+total_e #keep running total of total energy of all backbone sets
        all_Henergy_for_set=all_Henergy_for_set+total_hbonde #keep running total of total energy of all backbone sets
        outnum=outnum+1
        if outnum==10000:
            for t in xrange(len(outlist)):
                ofile.write(outlist[t])
            outlist=[]
            outnum=0
        line1 = nextline()
        if not line1:
            ok=0
        else:
            line=string.split(line1)
            set_conf(line, prot)

    for t in xrange(len(outlist)):
        ofile.write(outlist[t])

    print "number of backbone struct = ", num,
    print 'Total attempts to decorate = ',total_attmp
    print "total successes of decoration = ",accpt
    print "sum of all enegy of complete set = ", all_energy_for_set
    ofile2.write("num backbones: "+str(num)+"  total attempts to decorate: "+str(total_attmp)+"  success decorations: "+str(accpt)+"  Total Energy: "+str(all_energy_for_set)+"  Total Hbond Energy: "+str(all_Henergy_for_set)+" \n")
    ofile.close()                  
    print 'Finished'
    
if __name__ == "__main__":
    # sysargv1 = seq pdbfile, sysargv2 = which box of polyala
    run(sys.argv[1],sys.argv[2])
    #profile.run("run(sys.argv[1],sys.argv[2])")
