#! /usr/bin/python
import sys
import getopt

def matrix(filename):
    with open(filename, 'r') as infile:
        while True:
            line = infile.readline()
            if not line:
                break
            splitline = line.split()
            count = 0
            if line.startswith('MEME version'):
                meme_line = line
            if line.startswith('ALPHABET'):
                alphabet_line = line
            if line.startswith('strands:'):
                strands_line = line
            if line.startswith('Background letter frequencies'):
                next_line = infile.readline()
                if next_line:  # Ensure there is a next line
                	bkg_freq = next_line.strip()          
            if line.startswith('MOTIF'):
                if len(splitline) == 3:
                    outfilename_prefix = splitline[2]
                else:
                    outfilename_prefix = splitline[1]
                with open(outfilename_prefix + '_meme.txt', 'w') as outfile:
                    outfile.write(meme_line + '\n')
                    outfile.write(alphabet_line + '\n')
                    outfile.write(strands_line + '\n')
                    outfile.write(bkg_freq + '\n')
                    outfile.write('\n')
                    outfile.write('{}\t{}\n'.format(splitline[0], splitline[1]))
            if line.startswith('letter-probability matrix:'):# or ('log-odds matrix'):
                with open(outfilename_prefix + '_meme.txt', 'a') as outfile:
                    outfile.write(line)
                    while True:
                        next_line = infile.readline()
                        if len(next_line.split()) != 4:
                            break
                        outfile.write(next_line)
                    outfile.write('\n')

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:", ["help", "input="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)
    name = False
    for opt, arg in opts:
        if opt in ('-i', '--input'):
            name = arg
        elif opt in ('-h', '--help'):
            print('python ~/MEME_individual_from_db.py -i combined_motif_db.txt')
            sys.exit()
    if name:
        matrix(name)
    else:
        print('python ~/MEME_individual_from_db.py -i combined_motif_db.txt')

if __name__ == "__main__":
    main(sys.argv[1:])

