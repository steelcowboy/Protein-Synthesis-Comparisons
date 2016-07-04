from sys import argv, exit
from itertools import combinations

VERSION = 0.6

class Organism():
    rna_matrix = { "A": "U",
                   "T": "A",
                   "G": "C",
                   "C": "G" }

    def __init__(self, name, number, dna):
        self.dna_ = dna
        self.organism_number_ = number 
        self.name_ = name
        self.found_start_ = False
        self.is_sequenced_ = False
        self.separated_rna_ = []
    
    def get_name(self):
        return self.name_
    def get_rna(self):
        return self.rna_
    
    def get_output(self):
        oput = ""
        oput += "RNA sequence {}: {}\n\n".format((self.organism_number_ + 1), self.rna_)
        if not self.is_sequenced_:
            if not self.found_start_:
                oput += "No start sequence found, cannot sequence\n"
                oput += "-"*20 + "\n\n" 
                return oput
            else:
                oput += "No end sequence found, cannot sequence\n"
                oput += "-"*20 + "\n\n"
                return oput
        oput += "Found start sequence at {}\n".format(self.start_position_)
        oput += "Found end sequence {}\n\n".format(self.end_codon_)
        oput += " ".join([x for x in self.separated_rna_])
        oput += "\n\n" + "-"*20 + "\n\n"
        return oput

    def get_separated_rna(self):
        return self.separated_rna_
    def is_sequenced(self):
        return self.is_sequenced_
    def given_rna(self):
        self.rna_ = self.dna_
    def translate(self):
        self.rna_ = ''.join([self.rna_matrix[x] for x in self.dna_])
    
    def sequence_protein(self):
        start_seq = "AUG" 
        end_seqs  = ["UAA", "UAG", "UGA"]
        self.start_position_ = self.rna_.find(start_seq)
        if self.start_position_ == -1:
            self.found_start_ = False
            return
        else:
            self.found_start_ = True
        rna_tmp = self.rna_[self.start_position_:]
        length = len(rna_tmp)
        codon_count = length//3
        rna_tmp = rna_tmp[0:codon_count*3]
        self.separated_rna_.append("AUG")
        for i in range(1,codon_count):
            codon = rna_tmp[ 3*i : 3*i+3 ]
            if codon in end_seqs:
                self.separated_rna_.append(codon)
                self.end_codon_ = codon
                self.is_sequenced_ = True
                break
            else:
                self.separated_rna_.append(codon)

def compare_rna(orgs):
    oput = ""
    for x in combinations(orgs, 2):
        pt_muts = 0
        oput += "Comparison between {} and {}:\n".format(x[0].get_name(), x[1].get_name())
        for pos, codon in enumerate(x[0].get_separated_rna()):
            if codon != x[1].get_separated_rna()[pos]:
                oput += "The strands differ at codon {}: {} vs {}\n".format(
                        (pos+1), 
                        ''.join(x[0].get_separated_rna()[pos]), 
                        ''.join(x[1].get_separated_rna()[pos]))
                for acid in range(0,3):
                    if x[0].get_separated_rna()[pos][acid] != x[1].get_separated_rna()[pos][acid]:
                        pt_muts = pt_muts + 1
        oput += "\tNumber of point mutations: {}\n\n".format(pt_muts)
    return oput

def calculate_orgs(orgs, rna_given=False):
    new_orgs = []
    output = ""
    for x in range(0,len(orgs)):
        if rna_given:
            orgs[x].given_rna()
        else:
            orgs[x].translate()
        orgs[x].sequence_protein()
        output += orgs[x].get_output()
        if (orgs[x].is_sequenced()):
            new_orgs.append(orgs[x])
    if len(new_orgs) > 1:
        output += compare_rna(new_orgs)
    return output 

def output_info(output_string, quiet, output_file, filename):
    if output_file:
        outfile = open(filename, 'w')
        outfile.write(output_string)
        outfile.close()
    if quiet:
        return
    print(output_string)

def default_process(test_mode, quiet, output_file, filename):
    x = 0;
    times = 10000 if test_mode else 1 
    
    while True:
        orgs = []
        afri = Organism("African", 0, "CTACGTTCATCTGGTCAGAACTGGTTA")
        ande = Organism("Andes", 1, "TCACCTACGTCGATCTGGTCAGGACTT")
        anta = Organism("Antarctic", 2, "CACCTACGTTGATCCGGTCAGGACTGGTTA")
        
        orgs.append(afri)
        orgs.append(ande)
        orgs.append(anta)

        outinfo = calculate_orgs(orgs)
        output_info(outinfo, quiet, output_file, filename)
        
        x = x + 1
        if x == times:
            break

def showhelpinfo(script_name):
    print("Usage:   " + script_name + " [-option] [argument]\n" 
         +"option:  " + "-h  show help information\n" 
         +"         " + "-d  use default values, useful in testing\n" 
         +"         " + "-q  quiet\n" 
         +"         " + "-o  output file\n" 
         +"         " + "-r  input is rna\n" 
         +"         " + "-t  enable test mode\n" 
         +"         " + "-v  show version infomation\n" 
         +"example: " + script_name + " -o output.txt <DNA1> <DNA2> ...\n")
