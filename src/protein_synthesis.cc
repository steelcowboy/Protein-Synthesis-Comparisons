#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<algorithm>

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
        int start;
        vector<string> rna_sep;
    
    public:
        organism(string n, string d) {
            name = n;
            dna = d;
        }

        string get_name() {
            return this->name;
        }

        string get_rna() {
            return this->rna;
        }

        vector<string> get_rna_sep() {
            return this->rna_sep;
        }

        void translate() {
            int x = this->dna.length();
            for(int i=0; i<x; i++) {
                this->rna += rna_matrix[this->dna[i]];
            }
        }

        void sequence_protein() {
            int i;
            for (i = 0; i < this->rna.length()-3; i++) {
                if (this->rna.substr(i, 3) == "AUG") {
                    this->start = i;
                    cout << "Found a start sequence at " << i << endl;
                    break;
                }
            }
            this->rna = this->rna.substr(this->start);
            int len = this->rna.length();
            int codon_count = len/3;
            this->rna = this->rna.substr(0, codon_count*3);
            rna_sep.push_back("AUG");
            for (i = 1; i < codon_count; i++) {
                string codon = this->rna.substr(3*i, 3);
                if (codon == "UAA" || codon == "UAG" || codon == "UGA") {
                    rna_sep.push_back(codon);
                    cout << "Found end sequence " << codon << "\n\n"; 
                    break;
                }
                else { 
                    rna_sep.push_back(codon);
                }
            }
            for (int x = 0; x < rna_sep.size(); x++) {
                cout << rna_sep[x] << " ";
            }
            cout << "\n\n";
        }
};

void compare_rna(vector<organism> &orgs) {
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

    for (int x = 0; x < permutations.size(); x++) {
        int pt_muts = 0;
        cout << "Comparison between " << permutations[x][0].get_name() << " and " << permutations[x][1].get_name() << ":" << endl;
        vector<string> rna1 = permutations[x][0].get_rna_sep();
        vector<string> rna2 = permutations[x][1].get_rna_sep();
        for (int i = 0; i < rna1.size(); i++) {
            if (rna1[i] != rna2[i]) {
                cout << "The strands differ at codon " << i + 1 << ": " << rna1[i] << " vs " << rna2[i] << endl;
                for (int j =0; j < 3; j++) {
                    if (rna1[i][j] != rna2[i][j]) {
                        pt_muts++;
                    }
                }
            }
        }
        cout << "\tNumber of point mutations: " << pt_muts << "\n\n";
    }
}

int main(int argc, char *argv[]) {
    vector<organism> orgs;
    int x;
    vector<string> dna;
    vector<string> names;
    
    if (string(argv[1]) == "-ni") {  
        organism afri("African", "CTACGTTCATCTGGTCAGAACTGGTTA");
        organism ande("Andes", "TCACCTACGTCGATCTGGTCAGGACTT");
        organism anta("Antarctic", "CACCTACGTTGATCCGGTCAGGACTGGTTA");
        
        orgs.push_back(afri);
        orgs.push_back(ande);
        orgs.push_back(anta);
    }
    else {
        for (x=1; x < argc; x++) {
            dna.push_back(string(argv[x]));
        }
        for (x=1; x < argc; x++) {
            string name;
            cout << "What is strand " << x << "? ";
            getline(cin, name);
            names.push_back(name);
        }
        cout << string(20, '-') << endl;
        for (x=0; x < argc -1; x++) {
            organism tmp(names[x], dna[x]);
            orgs.push_back(tmp);
        }
    }
    for (x=0; x < orgs.size(); x++) {
        orgs[x].translate();
        cout << "\nRNA sequence " << x + 1 << ": " << orgs[x].get_rna() << "\n\n";
        orgs[x].sequence_protein();
        cout << string(20, '-') << endl;
    }
    compare_rna(orgs);
    return 0;
}
    





