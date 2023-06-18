#!/usr/bin/env python
import os
import argparse 

def write(nodes, cores, time, out, alloc, script, directory, write, unparsable, anions):
    writelines = '#!/bin/bash' + '\n'
    writelines += '#SBATCH -J ' + out + '\n'
    writelines += '#SBATCH --time=' + str(time) + ':00:00' + '\n'
    writelines += '#SBATCH -N ' + str(nodes) + '\n'
    writelines += '#SBATCH --tasks ' + str(cores) + '\n'
    writelines += '#SBATCH -o ' + out + '-%j.out' + '\n'
    writelines += '#SBATCH -e ' + out + '-%j.err' + '\n'
    writelines += '#SBATCH --account=' + alloc + '\n'
    if time == 1:
        writelines += '#SBATCH --partition=debug\n'
    elif time <= 4:
        writelines += '#SBATCH --partition=short\n'
    elif time >= 48:
        writelines += '#SBATCH --partition=long\n'
    else:
        writelines += '#SBATCH --partition=standard\n'
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    script_location = os.path.join(__location__, script)
    if os.path.exists(script_location):
        writelines +='python ' + script_location + ' -d ' + directory + ' -w ' + write
        if unparsable is True and anions is True:
            writelines += ' ' + '-c' + ' ' + '-a' + '\n'
        elif unparsable is False and anions is True:
            writelines += ' ' + '-a' + '\n'
        elif unparsable is True and anions is False:
            writelines += ' ' + '-c' + '\n'
        else:
            writelines += '\n'
        writelines +='exit 0'+'\n'

        with open('submit.sh', 'w') as f:
            f.write(writelines)
    else:
        print('%s not found at %s' % (script, __location__))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nodes', help = 'Number of nodes',
                        type = int, default = 1)
    parser.add_argument('-c', '--cores', help = 'Number of cores',
                        type = int, default = 36)
    parser.add_argument('-t', '--time', help = 'Time limit',
                        type = int, default = 4)
    parser.add_argument('-o', '--outfile', help = 'Outfile name',
                        type = str, required = True)
    parser.add_argument('-a', '--allocation', help = 'Allocation',
                        type = str, default = 'custws')
    parser.add_argument('-s', '--script', help = 'Python script to submit.',
                        type = str, required = True)
    parser.add_argument('-d', '--directory', help = 'directory to read from',
                        type = str, required = True)
    parser.add_argument('-w', '--write_file_name', help = '.json written to',
                        type = str, required = True)
    parser.add_argument('-cu', '--consider_unparsable', help='Consider unparsable files; default = False', action='store_true')
    parser.add_argument('-aab', '--anions_AB', help='Anions assigned to A/B sites; default = False', action='store_true')
    args = parser.parse_args()

    write(args.nodes, args.cores, args.time, args.outfile, args.allocation,
          args.script, args.directory, args.write_file_name, args.consider_unparsable, args.anions_AB)
    os.system('sbatch submit.sh')
