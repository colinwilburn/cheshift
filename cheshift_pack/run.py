'''
Runnable script version of: 
cheshift plugin-in: Validate your protein model with PyMOL
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/cheshift

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: September 2014
License: GNU General Public License
Version 3.6
'''

import Tkinter
from Tkinter import *
import tkFileDialog
import Pmw
from pymol import cmd, stored
import sys
import re
import os
import numpy as np
from scipy.interpolate import griddata

import traceback


def run():
    """Checks if files were provided and calls the prediction routine"""
    path = os.path.dirname(os.path.abspath(__file__)) 
    print path

    try:
        pdb_filename = cmd.get_names('all')[0]
        print pdb_filename

    except:
        Pmw.MessageDialog(title = 'Error',message_text = 'Please choose a\n PDB file')

    prediction(pdb_filename, path)

def prediction(pdb_filename, path):
    """Run the CheShift CS prediction routine"""
    cmd.set('suspend_updates', 'on')
    pose, res_label_list, res_num_list, states = pose_from_pdb(pdb_filename)
    res_num_list = res_num_list[:len(res_label_list)] 
    Db = load(path)
    raw(pose, res_label_list, res_num_list, states, Db)
    print '<'*80 + '\nYou didn`t provide a file with chemical Shifts, hence CheShift-2 assumed you\n only wanted the predicted CS. The predicted chemical shifts can be found in the file %s.txt\n' % pose + '>'*80
    for sel in ['alpha', 'beta', 'gamma', 'delta']:
        cmd.delete(sel)
    cmd.set('suspend_updates', 'off')
    '''
    raw(pose, residues, total_residues, states, Db) 
    print '<'*80 + '\nYou didn`t provide a file with chemical Shifts, hence CheShift-2 assumed you\n only wanted the predicted CS. The predicted chemical shifts can be found in the file %s.txt\n' % pose + '>'*80
    for sel in ['A', 'B', 'C', 'D']:
        cmd.delete(sel)
    cmd.set('suspend_updates', 'off')'''



def pose_from_pdb(pdb_file_name):
    """Gets information from the pdb like the number of residues, the sequence,
    the number of states and the name of the object"""
    pose = pdb_file_name
    remStr = "all and not (alt ''+A)"
    cmd.remove(remStr)
    cmd.alter('all', "alt=''")
    stored.residues = []
    stored.ResiduesNumber = []
    cmd.iterate('(name ca)','stored.residues.append((resn))')
    cmd.iterate('all','stored.ResiduesNumber.append((resi))')
    first = int(stored.ResiduesNumber[0])
    #cmd.alter(pose, 'resi=str(int(resi)-%s)' % (first))
    cmd.sort(pose)
    states = cmd.count_states('all') + 1
    
    res_num_set = set(stored.ResiduesNumber)

    res_num_list = list(res_num_set)
    res_num_list = [int(_) for _ in res_num_list]
    res_num_list.sort()
    res_label_list = stored.residues[:len(res_num_list)]

    return pose, res_label_list, res_num_list, states




def get_phi(res_num, state):
    """Computes the phi torsional angle"""
    if res_num != 0:
        cmd.select('alpha', 'resi %s and name C' % (res_num-1))
        cmd.select('beta', 'resi %s and name N' % res_num)
        cmd.select('gamma', 'resi %s and name CA' % res_num)
        cmd.select('delta', 'resi %s and name C' % res_num)
        phi = cmd.get_dihedral('alpha', 'beta', 'gamma', 'delta', state)
        return phi
    else:
        return float('nan')


def get_psi(res_num, state, res_num_list):
    """Computes the psi torsional angle"""
    if res_num != res_num_list[-1]:
        cmd.select('alpha', 'resi %s and name N and chain A' % res_num)
        cmd.select('beta', 'resi %s and name CA and chain A' % res_num)
        cmd.select('gamma', 'resi %s and name C and chain A' % res_num)
        cmd.select('delta', 'resi %s and name N and chain A' % (res_num+1))
        psi = cmd.get_dihedral('alpha', 'beta', 'gamma', 'delta', state)
        return psi
    else:
        return float('nan')


def get_omega(res_num, state, res_num_list):
    """Computes the omega torsional angle"""
    if res_num != res_num_list[-1]:
        cmd.select('alpha', 'resi %s and name CA' % res_num)
        cmd.select('beta', 'resi %s and name C' % res_num)
        cmd.select('gamma', 'resi %s and name N' % (res_num+1))
        cmd.select('delta', 'resi %s and name CA' % (res_num+1))
        omega = cmd.get_dihedral('alpha', 'beta', 'gamma', 'delta', state)
        return omega
    else:
        return float('nan')


def get_chi1(res_num, res_name, state):
    """Computes the chi1 torsional angle"""
    if res_name not in ['ALA', 'GLY', 'PRO']:
        cmd.select('alpha', 'resi %s and name N' % res_num)
        cmd.select('beta', 'resi %s and name CA' % res_num)
        cmd.select('gamma', 'resi %s and name CB' % res_num)
        cmd.select('delta', 'resi %s and (name CG or name CG1 or name OG1 or name OG or name SG)' % res_num)
        chi1 = cmd.get_dihedral('alpha', 'beta', 'gamma', 'delta', state)
        return chi1
    else:
        return float('nan')


def get_chi2(res_num, res_name, state):
    """Computes the chi2 torsional angle"""
    if res_name not in ['ALA', 'GLY', 'PRO', 'SER', 'THR', 'VAL', 'CYS']:
        cmd.select('alpha', 'resi %s and name CA' % res_num)
        cmd.select('beta', 'resi %s and name CB' % res_num)
        cmd.select('gamma', 'resi %s and (name CG or name CG1 or name OG1 or name OG)' % res_num)
        cmd.select('delta', 'resi %s and (name CD or name CD1 or name OD1 or name ND1 or name SD)' % res_num)
        chi2 = cmd.get_dihedral('alpha', 'beta', 'gamma', 'delta', state)
        return chi2
    else:
        return float('nan')


def load(path):
    """Load the files containing the theoretical chemical shifts. Creates a 
    dictionary to store the data."""
    aminoacids = ['ALA','ARG','ASN','ASP','GLU','GLN','GLY','HIS','ILE','LEU',
    'LYS','MET','PHE','PRO','SER','THR','TYR','TRP','VAL']
    Db = {}
    for aminoacid in aminoacids:
        vector = []
        print(path)
        for line in open(os.path.join(path, 'CS_DB', 'CS_db_%s' % aminoacid)).readlines():
            vector.append(map(float, line.split()))
        Db[aminoacid] = vector
    return Db


def near_pro(omega, psi, Db):
    """Computes the chemical shifts from the psi and omega torsional angles
    by linear interpolation from the theoretical values stored in Db. 
    This funcion works only for proline"""
    points = []
    values_Ca = []
    values_Cb = []
    if omega  <= -90: # torsional angles are circular, for example -180=180. angles smaller than -90
        omega = 180   # are closer to 180 than to 0
    lista = np.array([0., 180.])# cheshift databse has only two values for proline omega angle, this
    index = (np.abs(lista-omega)).argmin() # two lines calculate the nearest omega angle, in the datbase,
    nearestOmega = lista[index] 
    for line in Db['PRO']: # PRO database is small. hence to calculate the theoretical CS i just take
        if line[0] == nearestOmega: # all the values with the nearestOmega
            vi, yi, zi = line[1], line[4], line[5]
            points.append(vi)
            values_Ca.append(yi), values_Cb.append(zi)
    points = np.array(points)
    values_Ca = np.array(values_Ca)
    values_Cb = np.array(values_Cb)
    values_Ca_New = griddata(points, values_Ca, (psi), method='linear') #linear interpolation
    values_Cb_New = griddata(points, values_Cb, (psi), method='linear') #linear interpolation
    return values_Ca_New, values_Cb_New


def near(phi, psi, chi1, chi2, res_name, Db):
    """Computes the chemical shifts from the torsional angles by linear 
    interpolation from the theoretical values stored in Db. 
    This funcion works for non-proline residues"""
    points = []
    values_Ca = []
    values_Cb = []
    phi_list = []
    phi_round = int(round(phi, -1)) # round to the nearest values in the database
    psi_round = int(round(psi, -1))
    chi1_rotamers =  np.array([-180., -150., -120., -90., -60., -30., 0., 30., 60., 90., 120., 150., 180.])
    index = (np.abs(chi1_rotamers-chi1)).argmin()
    nearestChi1_A = chi1_rotamers[index]
    chi1_rotamers_new = np.delete(chi1_rotamers, index)
    index = (np.abs(chi1_rotamers_new-chi1)).argmin()
    nearestChi1_B = chi1_rotamers_new[index]
    if phi > phi_round: # for phi and psi angles get the two nearest values in the database
        phi_range = range(phi_round, phi_round+20, 10)
    else:
        if phi_round == -180:
            phi_round = 180
        phi_range = range(phi_round-10, phi_round+10, 10)
    if psi > psi_round:
        psi_range = range(psi_round, psi_round+20, 10)
    else:
        if psi_round == -180:
            psi_round = 180
        psi_range = range(psi_round-10, psi_round+10, 10)
    for phi_value in phi_range: # A trick to avoid reading the whole list. The indexes where
        y = int(phi_value * 0.1 + 19)#  the necesarry values are stored (start and end) are calculated.
        if y > 37: 
            y = -(37-y)
        lenght = (len(Db[res_name])/37)
        end = lenght * y
        start = end-lenght
        for i in range(start, end):
            phi_list.append(Db[res_name][i])
    if res_name in ['ALA', 'GLY']:
        for line in phi_list:
            for psi_value in psi_range:
                if psi_value == line[1]:
                    ui, vi, yi, zi = line[0], line[1], line[4], line[5] 
                    vector = ui, vi
                    points.append(vector)
                    values_Ca.append(yi), values_Cb.append(zi)
        points = np.array(points)
        values_Ca = np.array(values_Ca)
        values_Cb = np.array(values_Cb)
        values_Ca_New = griddata(points, values_Ca, (phi, psi), method='linear')
        values_Cb_New = griddata(points, values_Cb, (phi, psi), method='linear')
        return  values_Ca_New, values_Cb_New
    elif res_name in ['SER', 'THR', 'VAL', 'CYS']:
        for line in phi_list:
            for psi_value in psi_range:
                if psi_value == line[1] and (line[2] == nearestChi1_A or line[2] == nearestChi1_B):
                    ui, vi, wi, yi, zi = line[0], line[1], line[2], line[4], line[5] 
                    vector = ui, vi, wi
                    points.append(vector)
                    values_Ca.append(yi), values_Cb.append(zi)
        points = np.array(points)
        points = np.array(points)
        values_Ca = np.array(values_Ca)
        values_Cb = np.array(values_Cb)
        values_Ca_New = griddata(points, values_Ca, (phi, psi, chi1), method='linear')
        values_Cb_New = griddata(points, values_Cb, (phi, psi, chi1), method='linear')
        return values_Ca_New, values_Cb_New
    else:
        lista = []
        for i in range (0, 3):
            rotamer = phi_list[i][3]
            if rotamer < 0:
                rotamer = rotamer + 360
            lista.append(rotamer)
        if 0. in lista:
            lista.append(360)
        lista = np.array(lista)
        if chi2 < 0:
            chi2 = chi2 + 360
        index = (np.abs(lista-chi2)).argmin()
        nearestChi2 = lista[index] 
        if nearestChi2 > 180:
            nearestChi2 = nearestChi2 - 360
        for line in phi_list:
            for psi_value in psi_range:
                if psi_value == line[1] and line[3] == nearestChi2 and (line[2] == nearestChi1_A or line[2] == nearestChi1_B):
                    ui, vi, wi, yi, zi = line[0], line[1], line[2], line[4], line[5] 
                    vector = ui, vi, wi
                    points.append(vector)
                    values_Ca.append(yi), values_Cb.append(zi)
        points = np.array(points)
        values_Ca = np.array(values_Ca)
        values_Cb = np.array(values_Cb)
        values_Ca_New = griddata(points, values_Ca, (phi, psi, chi1), method='linear')
        values_Cb_New = griddata(points, values_Cb, (phi, psi, chi1), method='linear')
        return values_Ca_New, values_Cb_New



def get_chemical_shifts_raw(res_label_list, res_num_list, state, Db):
    """Call the near and near_pro function for all the residues. This function 
    is exclusive of the prediction routines"""
    chemical_shifts = []
    for i, res_num in enumerate(res_num_list): 
        print(i, len(res_num_list), len(res_label_list))
        res_name = res_label_list[i]
        try:
            res_name_next = residues[i+1]
        except:
            res_name_next = 'GLY'
        if res_name != 'PRO' and res_name != 'CYS':
            try:
                phi = get_phi(res_num, state)
                psi = get_psi(res_num, state, res_num_list)
                chi1 = get_chi1(res_num, res_name, state)
                chi2  = get_chi2(res_num, res_name, state)
                values_Ca_New, values_Cb_New = near(phi, psi, chi1, chi2, res_name, Db)
            except Exception as err:
                traceback.print_exc()
                values_Ca_New = 999.00
                values_Cb_New = 999.00
        elif res_name == 'CYS':
            values_Ca_New = 999.00
            values_Cb_New = 999.00
        else:
            try:
                omega = get_omega(res_num-1, state, res_num_list)
                values_Ca_New, values_Cb_New = near_pro(omega, psi, Db)
            except:
                values_Ca_New, values_Cb_New = 999.00, 999.00 
        if res_name_next == 'PRO':
            a = [res_name, res_num, round((values_Ca_New -1.95), 2), round((values_Cb_New),2)]
            chemical_shifts.append(a)
        else:
           a = [res_name, res_num, round((values_Ca_New), 2), round((values_Cb_New),2)]
           chemical_shifts.append(a)
    return chemical_shifts


def raw(pose, res_label_list, res_num_list, states, Db):
    """Calculates the theoretical chemical shifts. This function is exclusive 
    of the prediction routine"""
    cs_theo_list = []
    for state in range(1, states):
        cs_list = get_chemical_shifts_raw(res_label_list, res_num_list, state, Db)
        cs_theo_list.append(cs_list)
    fd = open('%s.txt' % pose, 'w')
    fd.write('Ca Chemical Shifts\n')
    for i, res_num in enumerate(res_num_list): 
        cs_theo_line = []
        for conformation in range(0, len(cs_theo_list)):
            res_name, res_num, Ca_shift, Cb_shift = cs_theo_list[conformation][i]
            if float(Ca_shift) > 100.:
                Ca_shift = 999.00
            cs_theo_line.append('%6.2f' % (Ca_shift))
        res_line = "\t".join(cs_theo_line)
        fd.write('%s\t %s\t %s\n' % (res_num, res_label_list[i], res_line))
    fd.write('\nCb Chemical Shifts\n')
    for i, res_num in enumerate(res_num_list): 
        cs_theo_line=[]
        for conformation in range(0, len(cs_theo_list)):
            res_name, res_num, Ca_shift, Cb_shift = cs_theo_list[conformation][i]
            if float(Cb_shift) > 100.:
                Cb_shift = 999.00
            cs_theo_line.append('%6.2f' % (Cb_shift))
        res_line = "\t".join(cs_theo_line)
        fd.write('%s\t %s\t %s\n' % (res_num, res_label_list[i], res_line))
    fd.close()

cmd.extend("cheshift_prediction", run)
