#ifndef SYNTHESIS_H
#define SYNTHESIS_H

#include<string>
#include<vector>
#include<map>

#define VERSION 0.6

class Organism {
    private:
        std::string name_;
        int organism_number_;
        std::string dna_;
        static std::map<char, char> rna_matrix;
        std::string rna_;
        std::size_t start_position_ = 0;
        bool found_start_ = 0;
        bool is_sequenced_ = 0;
        std::string end_codon_;
        std::vector<std::string> separated_rna_;
   
    public:
        Organism(std::string organism_name, int organism_number, std::string organism_dna) : 
            name_(organism_name), organism_number_(organism_number), dna_(organism_dna) {}
        std::string get_name() const;
        std::string get_rna() const;
        std::string get_output() const;
        std::vector<std::string> get_separated_rna() const; 
        bool is_sequenced();
        void given_rna();
        void translate();
        void sequence_protein(); 
};
std::string compare_rna(const std::vector<Organism> &orgs);
std::string calculate_orgs(std::vector<Organism> &orgs, bool rna_given);
void output_info(const std::string &output_string, bool quiet, bool output_file, std::string filename);
void default_process(bool test_mode, bool quiet, bool output_file, std::string filename);
void showhelpinfo(char *s);

#endif
