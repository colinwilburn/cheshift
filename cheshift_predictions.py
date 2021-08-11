import subprocess

from convert_file import convert_file
from cheshift_pack import run


def predict(filename_pdb, filename_out):
    args = [
        '/usr/software/bin/pymol', '-c', filename_pdb, '-d', 
        'run cheshift_predictions.py', '-d', 'cheshift_prediction'
    ]
    subprocess.call(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE) # get CS predictions from cheshift
    filename_cheshift = filename_pdb[:-3] + 'txt' # take off the pdb and append txt; this is where the cheshift predictions were written
    convert_file(filename_cheshift, filename_out) # convert to CSV and write
