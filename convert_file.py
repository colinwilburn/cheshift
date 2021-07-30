import csv


def convert_file(filename_old, filename_new):
    residues_dict = read_from_cheshift(filename_old)
    write_to_out(filename_new, residues_dict)

def read_from_cheshift(filename_old):

    atom_names = ['CA', 'CB'] # the files lists CA shifts first, and then CB
    half_counter = -1

    residues_dict = {}

    with open(filename_old, 'r') as infile:
        filereader = csv.reader(infile, delimiter='\t')
        for row in filereader:
            if len(row) == 1: # we've hit one of the header files, so transition to CA or CB
                half_counter += 1
                atom_name = atom_names[half_counter]
                res_num = 0
            elif len(row) > 1: # not one of the specific rows
                res_num += 1
                res_label, cs = convert_row(row)
                if res_num not in residues_dict:
                    residues_dict[res_num] = {}
                    residues_dict[res_num]['res_label'] = res_label
                residues_dict[res_num][atom_name] = cs

    return residues_dict
    
def convert_row(row):
    res_label = row[0]

    if float(row[1]) == 999.0: # this means that chemshift couldn't calculate the chemical shift
        cs = "None"
    else:
        cs = 0 # and we will take the average
        for cs_model in row[1:]:
            cs += float(cs_model)
        cs = cs / len(row[1:])

    return res_label, cs

def write_to_out(filename_new, residues_dict):
    
    rows = [] # to write
    res_num_list = list(residues_dict.keys())
    res_num_list.sort()
    for res_num in res_num_list:
        info_dict = residues_dict[res_num]
        res_label = info_dict['res_label']
        cs_ca = info_dict['CA']
        cs_cb = info_dict['CB']
        row = [res_num, res_label, cs_ca, cs_cb]
        rows.append(row)

    header = ['res_num', 'res_label', 'cs_ca', 'cs_cb']
    with open(filename_new, 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    