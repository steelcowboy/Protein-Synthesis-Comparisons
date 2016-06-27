#include "synthesis.h"
#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>

std::map<char, char> Organism::rna_matrix = {
    {'A', 'U'},
    {'T', 'A'},
    {'G', 'C'},
    {'C', 'G'}
};

std::string Organism::get_name() {
    return this->name;
}

std::string Organism::get_rna() {
    return this->rna;
}

std::string Organism::get_output() {
    std::stringstream oput;
    oput << "RNA sequence " << this->organism_number + 1 << ": " << this->rna << "\n\n";
    if (!this->is_sequenced_) {
        oput << "No start sequence found\n";
        oput << std::string(20, '-') + "\n";
        return oput.str();
    }
    oput << "Found a start sequence at " << this->start_position << std::endl;
    oput << "Found end sequence " << this->end_codon << "\n\n";
    for (int x = 0; x < this->rna_sep.size(); x++) 
        oput << rna_sep[x] << " ";
    oput << "\n\n";
    oput << std::string(20, '-') + "\n";
    return oput.str();
}

std::vector<std::string> Organism::get_rna_sep() {
    return this->rna_sep;
}

void Organism::given_rna() {
    this->rna = this->dna;
}

void Organism::translate() {
    for (int i=0;  i<this->dna.length(); i++) 
        this->rna += this->rna_matrix[this->dna[i]];
}

void Organism::sequence_protein() {
    char start_sequence[] = "AUG";
    this->start_position = this->rna.find(start_sequence);
    if (this->start_position == std::string::npos) {
        this->is_sequenced_ = 0;
     }
    std::string rna_tmp = this->rna.substr(this->start_position);
    int len = rna_tmp.length();
    int codon_count = len/3;
    rna_tmp = rna_tmp.substr(0, codon_count*3);
    rna_sep.push_back("AUG");
    for (int i = 1; i  < codon_count; i++) {
        std::string codon = rna_tmp.substr(3*i, 3);
        if (codon == "UAA" || codon  == "UAG" || codon == "UGA") {
            rna_sep.push_back(codon);
            this->end_codon = codon;
            break;
        }
        else { 
            rna_sep.push_back(codon);
         }
    }
    this->is_sequenced_ = 1;
} 

bool Organism::is_sequenced() { return this->is_sequenced_; }

std::string compare_rna(std::vector<Organism> &orgs) {
    int n = orgs.size();
    
    std::vector<bool> v(n);
    fill(v.begin(), v.end() - n + 2, true);

    std::vector<std::vector<Organism>> permutations;

    do { 
        std::vector<Organism> perm;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                perm.push_back(orgs[i]);
             }
        } 
        permutations.push_back(perm);
    } while (prev_permutation(v.begin(), v.end()));

    std::stringstream output;

    for (int x = 0; x < permutations.size(); x++) {
        int pt_muts = 0;
        output << "Comparison between " << 
          permutations[x][0].get_name() << " and " <<
          permutations[x][1].get_name() << ":" << std::endl;
        std::vector<std::string> rna1 = permutations[x][0].get_rna_sep();
        std::vector<std::string> rna2 = permutations[x][1].get_rna_sep();
        for (int i = 0; i < rna1.size(); i++) {
            if (rna1[i] != rna2[i]) {
                output << "The strands differ at codon " <<
                  i + 1 << ": " << rna1[i] << 
                  " vs " << rna2[i] << std::endl;
                for (int j =0; j < 3; j++) {
                    if (rna1[i][j] != rna2[i][j]) {
                        pt_muts++;
                    }
                }
            }
        }
        output << "\tNumber of point mutations: " << pt_muts << "\n\n";
    }
    return output.str();
 } 

std::string calculate_orgs(std::vector<Organism> &orgs, bool rna_given=false) {
    std::vector<Organism> new_orgs;
    std::string output;
    for (int x=0; x<orgs.size(); x++) {
        if (rna_given) 
            orgs[x].given_rna();
        else 
            orgs[x].translate();
        orgs[x].sequence_protein();
        output += orgs[x].get_output();
        if (orgs[x].is_sequenced()) {
            new_orgs.push_back(orgs[x]);
         }
      }
    if (new_orgs.size() > 1) 
        output += compare_rna(new_orgs);
    return output;
 }

void output_info(std::string &output, bool q, bool oput, std::string fn) {
    if (oput) {
        std::ofstream outfile;
        outfile.open(fn);
        outfile << output;
        outfile.close();
    }
    if (q) {
        return;
    }
    std::cout << output;

}

void default_process(bool t, bool q, bool otf, std::string fn) {
    int times;
    int x=0;
    
    if (t) 
        times = 10000;
    else
        times = 1;
    
    do {
        std::vector<Organism> orgs;
        Organism afri("African", 0, "CTACGTTCATCTGGTCAGAACTGGTTA");
        Organism ande("Andes", 1, "TCACCTACGTCGATCTGGTCAGGACTT");
        Organism anta("Antarctic", 2, "CACCTACGTTGATCCGGTCAGGACTGGTTA");
        
        orgs.push_back(afri);
        orgs.push_back(ande);
        orgs.push_back(anta);

        std::string outinfo = calculate_orgs(orgs);
        output_info(outinfo, q, otf, fn);
        x++;
     } while (x < times);
} 

void showhelpinfo(char *s) {
  std::cout<<"Usage:   "<< s <<" [-option] [argument]"<<std::endl;
  std::cout<<"option:  "<<"-h  show help information"<<std::endl;
  std::cout<<"         "<<"-d  use default values, useful in testing"<<std::endl;
  std::cout<<"         "<<"-q  quiet"<<std::endl;
  std::cout<<"         "<<"-o  output file"<<std::endl;
  std::cout<<"         "<<"-r  input is rna"<<std::endl;
  std::cout<<"         "<<"-t  enable test mode"<<std::endl;
  std::cout<<"         "<<"-v  show version infomation"<<std::endl;
  std::cout<<"example: "<< s <<" -o output.txt <DNA1> <DNA2> ..."<<std::endl;
}
