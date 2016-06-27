#ifndef SYNTHESIS_H
#define SYNTHESIS_H

#include<string>
#include<vector>

#define VERSION 0.4

class organism {
    private:
        std::string name;
        std::string dna;
        std::string rna;
        std::vector<std::string> rna_sep;
        std::string output;
   
    public:
        organism(std::string n, std::string d) : name(n), dna(d) {}
        std::string get_name();
        std::string get_rna();
        std::string get_output();
        std::vector<std::string> get_rna_sep(); 
        void add_separator(); 
        void translate(int num, bool r);
        int sequence_protein(); 
};
std::string compare_rna(std::vector<organism> &orgs);
std::string calculate_orgs(std::vector<organism> &orgs, bool r);
void output_info(std::string &output, bool q, bool oput, std::string fn);
void default_process(bool t, bool q, bool r, bool otf, std::string fn);
void showhelpinfo(char *s);

#endif
