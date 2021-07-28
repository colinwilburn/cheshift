import subprocess
from cheshift_pack import run


def predict(pdb_filename):
    args = [
        'pymol', '-c', pdb_filename, '-d', 'run cheshift_predictions.py', '-d', 
        'cheshift_prediction'
    ]
    subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
