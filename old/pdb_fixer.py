xyz = open('Pdplussixprotein.pdb')
'''
toggle = 0
lines = ''
for line in xyz:
    if len(line.split()) > 1:
        atom_num = int(line.split()[1])

        if (toggle == 0) and (atom_num != 1963):
            new_line = line
            prev_num = atom_num * 1

        if toggle == 1:
            replacement = str(atom_num + 1)
            new_line = line.replace(str(atom_num),replacement)
            #if len(atom_num) != len(prev_num): do this manually
            prev_num = atom_num * 1

        if atom_num == 1963:
            toggle = 1 # everything after this gets indexed by 1
            replacement = str(prev_num + 1)
            new_line = line.replace(str(atom_num),replacement)
            prev_num = 98

        lines += new_line
'''

#for lena's Pd + protein thinge
toggle = 0
lines = ''
class_num = 1
for line in xyz:
    try:
        if line.split()[2] == 'PD1':
            new_line = line

            replacement = str(class_num)
            n_digits = len(replacement)
            if n_digits == 1:
                new_line = line.replace('   2 ','   ' +replacement+' ')
            elif n_digits == 2:
                new_line = line.replace('   2 ','  ' +replacement+' ')
            elif n_digits == 3:
                new_line = line.replace('   2 ',' ' +replacement+' ')
            elif n_digits == 4:
                new_line = line.replace('   2 ','' + replacement+' ')


            class_num += 1

            lines += new_line
    except:
        aa = 1