#include<unistd.h> 
#include<iostream> 
#include<string> 
#include<vector> 
#include "synthesis.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        showhelpinfo(argv[0]);
        return 1;
    }
    
    char tmp;
    bool output_to_file = false;
    bool test_mode = false;
    bool quiet = false;
    bool def = false;
    bool given_rna = false;
    std::string filename = "false";

    while((tmp = getopt(argc,argv,"ho:vtqdr")) != -1) {
        switch(tmp) {
            case 'h':
                showhelpinfo(argv[0]);
                return 0;
            case 'v':
                std::cout << "Program version: " << VERSION << std::endl;
                return 0;
            case 'o':
                output_to_file = true;
                filename = std::string(optarg);
                break;
            case 't':
                test_mode = true;
                break;
            case 'q':
                quiet = true;
                break;
            case 'd':
                def = true;
                break;
            case 'r':
                given_rna = true;
                break;
            default:
                showhelpinfo(argv[0]);
                return 1;
        }
    }           
    
    if (def) {  
        default_process(test_mode, quiet, given_rna, output_to_file, filename);    
    }

    else {
        std::vector<organism> orgs;
        int x;
        std::vector<std::string> dna;
        
        for (x=optind; x < argc; x++) {
            dna.push_back(std::string(argv[x]));
        }
        
        int numorgs = dna.size();
        
        if (numorgs == 1) {
            organism idk("idk", dna[0]);
            orgs.push_back(idk);
        }
        else {
            for (x=0; x < numorgs; x++) {
                std::string name;
                std::cout << "What is strand " << x + 1 << "? ";
                getline(std::cin, name);
                organism tmp(name, dna[x]);
                orgs.push_back(tmp);
            }
            std::cout << std::string(20, '-') << std::endl;

        }
        std::string outinfo = calculate_orgs(orgs, given_rna);
        output_info(outinfo, quiet, output_to_file, filename);
    }
    return 0;
}



