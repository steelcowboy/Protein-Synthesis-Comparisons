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

std::string Organism::get_name() const {
    return name_;
}

std::string Organism::get_rna() const {
    return rna_;
}

std::string Organism::get_output() const {
    std::stringstream oput;
    oput << "RNA sequence " << organism_number_ + 1 << ": " << rna_ << "\n\n";
    if (!is_sequenced_) {
        if (!found_start_) {
            oput << "No start sequence found, cannot sequence\n";
            oput << std::string(20, '-') + "\n";
            return oput.str();
        }
        else {
            oput << "No end sequence found, cannot sequence\n";
            oput << std::string(20, '-') + "\n";
            return oput.str();
        }

    }
    oput << "Found start sequence at " << start_position_ << std::endl;
    oput << "Found end sequence " << end_codon_ << "\n\n";
    for (size_t x = 0; x < separated_rna_.size(); x++) 
        oput << separated_rna_[x] << " ";
    oput << "\n\n" << std::string(20, '-') + "\n\n";
    return oput.str();
 }

std::vector<std::string> Organism::get_separated_rna() const {
    return separated_rna_;
} 

void Organism::given_rna() {
    rna_ = dna_;
} 

void Organism::translate() {
    for (size_t i=0;  i<dna_.length(); i++) 
        rna_ += rna_matrix[dna_[i]];
}

void Organism::sequence_protein() {
    char start_sequence[] = "AUG"; 
    
    start_position_ = rna_.find(start_sequence);
    if (start_position_ == std::string::npos) {
        return;
    } else { found_start_ = 1; }
    std::string rna_tmp = rna_.substr(start_position_);
    int len = rna_tmp.length();
    int codon_count = len/3;
    rna_tmp = rna_tmp.substr(0, codon_count*3);
    separated_rna_.push_back("AUG");
     for (int i = 1; i  < codon_count; i++) {
        std::string codon = rna_tmp.substr(3*i, 3);
        if  (codon == "UAA" || codon  == "UAG" || codon == "UGA") {
            separated_rna_.push_back(codon);
            end_codon_ = codon;
            is_sequenced_ = 1;
            break;
        }
        else { 
            separated_rna_.push_back(codon);
         }
    }
} 

bool Organism::is_sequenced() { return is_sequenced_; }

std::string compare_rna(const std::vector<Organism> &orgs) {
    auto n = orgs.size();
    
    std::vector<bool> v(n);
    fill(v.begin(), v.end() - n + 2, true);

    std::vector<std::vector<Organism>> permutations;

    do { 
        std::vector<Organism> perm;
        for (size_t i = 0; i < n; ++i) {
            if (v[i]) {
                perm.push_back(orgs[i]);
             }
        } 
        permutations.push_back(perm);
    } while (prev_permutation(v.begin(), v.end()));

    std::stringstream output;

    for (size_t x = 0; x < permutations.size(); x++) {
        int pt_muts = 0;
        output << "Comparison between " << 
          permutations[x][0].get_name() << " and " <<
          permutations[x][1].get_name() << ":" << std::endl;
        std::vector<std::string> rna1 = permutations[x][0].get_separated_rna();
        std::vector<std::string> rna2 = permutations[x][1].get_separated_rna();
        for (size_t i = 0; i < rna1.size(); i++) {
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
    for (size_t x = 0; x<orgs.size(); x++) {
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

void output_info(const std::string &output_string, bool quiet, bool output_file, std::string filename) {
    if (output_file) {
        std::ofstream outfile;
        outfile.open(filename);
        outfile << output_string;
        outfile.close();
    }
    if (quiet) {
        return;
    }
    std::cout << output_string;

}

void default_process(bool test_mode, bool quiet, bool output_file, std::string filename) {
    int x = 0;
    auto times = test_mode ? 10000 : 1; 
    
    do { 
        std::vector<Organism> orgs;
        Organism afri("African", 0, "CTACGTTCATCTGGTCAGAACTGGTTA");
        Organism ande("Andes", 1, "TCACCTACGTCGATCTGGTCAGGACTT");
        Organism anta("Antarctic", 2, "CACCTACGTTGATCCGGTCAGGACTGGTTA");
        
        orgs.push_back(afri);
        orgs.push_back(ande);
        orgs.push_back(anta);

        std::string outinfo = calculate_orgs(orgs);
        output_info(outinfo, quiet, output_file, filename);
        x++;
     } while (x < times); 
} 

void showhelpinfo(char *s) {
  std::cout<<"Usage:   "<< s <<" [-option] [argument]\n"
           <<"option:  "<<"-h  show help information\n"
           <<"         "<<"-d  use default values, useful in testing\n"
           <<"         "<<"-q  quiet\n"
           <<"         "<<"-o  output file\n"
           <<"         "<<"-r  input is rna\n"
           <<"         "<<"-t  enable test mode\n"
           <<"         "<<"-v  show version infomation\n"
           <<"example: "<< s <<" -o output.txt <DNA1> <DNA2> ...\n"<<std::flush;
}
