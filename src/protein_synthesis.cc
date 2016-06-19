#include<unistd.h>
#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>

#define VERSION 0.3

using namespace std;

map<char, char> rna_matrix = {
    {'A', 'U'},
    {'T', 'A'},
    {'G', 'C'},
    {'C', 'G'}
};

class organism {
    private:
        string name;
        string dna;
        string rna;
        vector<string> rna_sep;
        string output;
   
    public:
        organism(string n, string d) : name(n), dna(d) {}

        string get_name() {
            return this->name;
        }

        string get_rna() {
            return this->rna;
        }

        string get_output() {
            return this->output;
        }

        vector<string> get_rna_sep() {
            return this->rna_sep;
        }

        void add_separator() {
            output += string(20, '-') + "\n";
        }
         void translate(int num, bool r) {
            stringstream oput;
            int x = this->dna.length();
            if (!r) {
                for (int i=0; i<x; i++) {
                    this->rna += rna_matrix[this->dna[i]];
                }
            } else {
                this->rna = this->dna;
            }
            oput << "RNA sequence " << num + 1 << ": " << this->rna << "\n\n";
            this->output += oput.str();
        }

        int sequence_protein() {
            stringstream oput;
            int start = 0;
            int i;
            for (i = 0; i < this->rna.length()-3; i++) {
                if (this->rna.substr(i, 3) == "AUG") {
                    start = i;
                    oput << "Found a start sequence at " << i << endl;
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
                string codon = this->rna.substr(3*i, 3);
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
};

string compare_rna(vector<organism> &orgs);
string calculate_orgs(vector<organism> &orgs, bool r);
void output_info(string &output, bool q, bool oput, string fn);
void default_process(bool t, bool q, bool r, bool otf, string fn);
void showhelpinfo(char *s);

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
    string filename = "false";

    while((tmp = getopt(argc,argv,"ho:vtqdr")) != -1) {
        switch(tmp) {
            case 'h':
                showhelpinfo(argv[0]);
                return 0;
            case 'v':
                cout << "Program version: " << VERSION << endl;
                return 0;
            case 'o':
                output_to_file = true;
                filename = string(optarg);
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
        vector<organism> orgs;
        int x;
        vector<string> dna;
        
        for (x=optind; x < argc; x++) {
            dna.push_back(string(argv[x]));
        }
        
        int numorgs = dna.size();
        
        if (numorgs == 1) {
            organism idk("idk", dna[0]);
            orgs.push_back(idk);
        }
        else {
            for (x=0; x < numorgs; x++) {
                string name;
                cout << "What is strand " << x + 1 << "? ";
                getline(cin, name);
                organism tmp(name, dna[x]);
                orgs.push_back(tmp);
            }
            cout << string(20, '-') << endl;

        }
        string outinfo = calculate_orgs(orgs, given_rna);
        output_info(outinfo, quiet, output_to_file, filename);
    }
    return 0;
}

string compare_rna(vector<organism> &orgs) {
    int n = orgs.size();
    
    vector<bool> v(n);
    fill(v.begin(), v.end() - n + 2, true);

    vector<vector<organism>> permutations;

    do { 
        vector<organism> perm;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                perm.push_back(orgs[i]);
             }
        } 
    permutations.push_back(perm);
     } while (prev_permutation(v.begin(), v.end()));

    stringstream output;

    for (int x = 0; x < permutations.size(); x++) {
        int pt_muts = 0;
        output << "Comparison between " << 
          permutations[x][0].get_name() << " and " <<
          permutations[x][1].get_name() << ":" << endl;
        vector<string> rna1 = permutations[x][0].get_rna_sep();
        vector<string> rna2 = permutations[x][1].get_rna_sep();
        for (int i = 0; i < rna1.size(); i++) {
            if (rna1[i] != rna2[i]) {
                output << "The strands differ at codon " <<
                  i + 1 << ": " << rna1[i] << 
                  " vs " << rna2[i] << endl;
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

string calculate_orgs(vector<organism> &orgs, bool r) {
    vector<organism> new_orgs;
    stringstream output;
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

void output_info(string &output, bool q, bool oput, string fn) {
    if (oput) {
        ofstream outfile;
        outfile.open(fn);
        outfile << output;
        outfile.close();
    }
    if (q) {
        return;
    }
    cout << output;

}

void default_process(bool t, bool q, bool r, bool otf, string fn) {
    int times;
    int x=0;
    
    if (t) 
        times = 10000;
    else
        times = 1;
    
    do {
        vector<organism> orgs;
        organism afri("African", "CTACGTTCATCTGGTCAGAACTGGTTA");
        organism ande("Andes", "TCACCTACGTCGATCTGGTCAGGACTT");
        organism anta("Antarctic", "CACCTACGTTGATCCGGTCAGGACTGGTTA");
        
        orgs.push_back(afri);
        orgs.push_back(ande);
        orgs.push_back(anta);

        string outinfo = calculate_orgs(orgs, r);
        output_info(outinfo, q, otf, fn);
        x++;
    } while (x < times);
}

void showhelpinfo(char *s) {
  cout<<"Usage:   "<< s <<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-d  use default values, useful in testing"<<endl;
  cout<<"         "<<"-q  quiet"<<endl;
  cout<<"         "<<"-o  output file"<<endl;
  cout<<"         "<<"-r  input is rna"<<endl;
  cout<<"         "<<"-t  enable test mode"<<endl;
  cout<<"         "<<"-v  show version infomation"<<endl;
  cout<<"example: "<< s <<" -o output.txt <DNA1> <DNA2> ..."<<endl;
}


