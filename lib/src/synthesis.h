#ifndef SYNTHESIS_H
#define SYNTHESIS_H

#include<string>
#include<vector>
#include<map>

#define VERSION 0.4

class Organism {
    private:
        std::string name;
        int organism_number;
        std::string dna;
        static std::map<char, char> rna_matrix;
        std::string rna;
        std::size_t start_position = 0;
        bool is_sequenced_ = 0;
        std::string end_codon;
        std::vector<std::string> rna_sep;
   
    public:
        Organism(std::string na, int nu, std::string d) : name(na), organism_number(nu), dna(d) {}
        std::string get_name();
        std::string get_rna();
        std::string get_output();
        std::vector<std::string> get_rna_sep(); 
        void given_rna();
        void translate();
        void sequence_protein(); 
        bool is_sequenced();
};
std::string compare_rna(std::vector<Organism> &orgs);
std::string calculate_orgs(std::vector<Organism> &orgs, bool r);
void output_info(std::string &output, bool q, bool oput, std::string fn);
void default_process(bool t, bool q, bool otf, std::string fn);
void showhelpinfo(char *s);

#endif
