#!/usr/bin/env python

from sys import argv, exit
from itertools import combinations

start_seq = ['A', 'U', 'G']
end_seqs  = [
             ['U', 'A', 'A'],
             ['U', 'A', 'G'],
             ['U', 'G', 'A'],
            ]

matrix = {"A": "U",
          "T": "A",
          "G": "C",
          "C": "G"}

class organism():
    def __init__(self, name, dna):
        self.dna = dna
        self.name = name
    def translate(self):
        self.rna = [matrix[x] for x in self.dna]
    def find_start(self):
        try:
            for i, x in enumerate(self.rna):
                if [self.rna[i], self.rna[i+1], self.rna[i+2]] == start_seq:
                    self.start = i
                    print("Found a start sequence at %i" % i)
                    break
        except IndexError:
            print("No start sequence found, cannot proceed")
            exit(1)

    def sep_rna(self):
        self.rna_sep = []
        for y in range(len(self.rna)//3):
            self.rna_sep.append(self.rna[3*y+self.start:3*y+3+self.start])

    def trim_end(self):
        for pos, seq in enumerate(self.rna_sep):
            if seq in end_seqs:
                self.end = pos
                print("Found end sequence %s\n" % ''.join(seq))
                self.rna_sep = self.rna_sep[:pos+1]
                break

def compare_rna(orgs):
    for x in combinations(orgs, 2):
        pt_muts = 0
        print("Comparison between {} and {}:".format(x[0].name, x[1].name))
        for pos, codon in enumerate(x[0].rna_sep):
            if codon != x[1].rna_sep[pos]:
                print("The strands differ at codon {}: {} vs {}".format((pos+1), ''.join(x[0].rna_sep[pos]), ''.join(x[1].rna_sep[pos])))
                for acid in range(0,3):
                    if x[0].rna_sep[pos][acid] != x[1].rna_sep[pos][acid]:
                        pt_muts = pt_muts + 1
        print("\tNumber of point mutations: {}\n".format(pt_muts))

def main():
    if argv[1] == "-ni":
        dna = ["CTACGTTCATCTGGTCAGAACTGGTTA",
               "TCACCTACGTCGATCTGGTCAGGACTT",
               "CACCTACGTTGATCCGGTCAGGACTGGTTA"]
        names = ["African", "Andes", "Antarctic"]
    else:
        dna = argv[1:]
        names = []
        for x in range(len(dna)):
            name = input("What is strand {}? ".format(x + 1))
            names.append(name)
        print("-" * 20)
    organisms = []
    #rna_list = []
    for num, strand in enumerate(dna):
        organisms.append(organism(names[num], strand))
    
    for pos, org in enumerate(organisms):
        org.translate()
        print("\nRNA sequence %d: %s\n" % ((pos + 1), ''.join(org.rna)))
        org.find_start()
        org.sep_rna()
        org.trim_end()
        for x in org.rna_sep:
            print(''.join(x), end=' ')
        print("\n\n" + "-" * 20)
        
        if len(dna) == 1:
            exit(0)
    compare_rna(organisms)


if __name__ == "__main__":
    main()
            
