#include "synthesis.h"
#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>

        std::map<char, char> rna_matrix = {
            {'A', 'U'},
            {'T', 'A'},
            {'G', 'C'},
            {'C', 'G'}
        };

std::string organism::get_name() {
    return this->name;
}

std::string organism::get_rna() {
    return this->rna;
}

std::string organism::get_output() {
    return this->output;
}

std::vector<std::string> organism::get_rna_sep() {
    return this->rna_sep;
}

void organism::add_separator() {
    output += std::string(20, '-') + "\n";
}

void organism::translate(int num, bool r) {
    std::stringstream oput;
    if (!r) {
        int x = this->dna.length();
        for (int i=0;  i<x; i++) {
            this->rna += rna_matrix[this->dna[i]];
        }
    } else {
        this->rna = this->dna;
    }
    oput << "RNA sequence " << num + 1 << ": " << this->rna << "\n\n";
    this->output += oput.str();
}

int organism::sequence_protein() {
    std::stringstream oput;
    int start = 0;
    int i;
    for (i = 0; i < this->rna.length()-3; i++) {
        if (this->rna.substr(i, 3) == "AUG") {
            start = i;
            oput << "Found a start sequence at " << i << std::endl;
            break;
         }

     }
     if (start == 0) {
        this->output += "No start sequence found\n";
        return 0;
     }
    this->rna = this->rna.substr(start);
    int len = this->rna.length();
    int codon_count = len/3;
    this->rna = this->rna.substr(0, codon_count*3);
    rna_sep.push_back("AUG");
    for (i = 1; i < codon_count; i++) {
        std::string codon = this->rna.substr(3*i, 3);
        if (codon == "UAA" || codon == "UAG" || codon == "UGA") {
            rna_sep.push_back(codon);
            oput << "Found end sequence " << codon << "\n\n";
            break;
        }
        else { 
            rna_sep.push_back(codon);
        }
    }
    for (int x = 0; x < rna_sep.size(); x++) {
        oput << rna_sep[x] << " ";
    }
    oput << "\n\n";
    this->output += oput.str();
    return 1;
}

std::string compare_rna(std::vector<organism> &orgs) {
    int n = orgs.size();
    
    std::vector<bool> v(n);
    fill(v.begin(), v.end() - n + 2, true);

    std::vector<std::vector<organism>> permutations;

    do { 
        std::vector<organism> perm;
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

std::string calculate_orgs(std::vector<organism> &orgs, bool r) {
    std::vector<organism> new_orgs;
    std::stringstream output;
    for (int x=0; x<orgs.size(); x++) {
        orgs[x].translate(x, r);
        int sequenced = orgs[x].sequence_protein();
        orgs[x].add_separator();
        output << orgs[x].get_output();
        if (sequenced) {
            new_orgs.push_back(orgs[x]);
        }
     }
    if (new_orgs.size() > 1) 
        output << compare_rna(new_orgs);
    return output.str();
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

void default_process(bool t, bool q, bool r, bool otf, std::string fn) {
    int times;
    int x=0;
    
    if (t) 
        times = 10000;
    else
        times = 1;
    
    do {
        std::vector<organism> orgs;
        organism afri("African", "CTACGTTCATCTGGTCAGAACTGGTTA");
        organism ande("Andes", "TCACCTACGTCGATCTGGTCAGGACTT");
        organism anta("Antarctic", "CACCTACGTTGATCCGGTCAGGACTGGTTA");
        
        orgs.push_back(afri);
        orgs.push_back(ande);
        orgs.push_back(anta);

        std::string outinfo = calculate_orgs(orgs, r);
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
