#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int sol_size;
vector<string> names, solNames;
set<string> solSet;
map<string, int> namesMap;

vector<int> p, w;


const string CODE_BOLD ="\e[1m";
const string CODE_OK ="\e[32m";
const string CODE_WRONG ="\e[31m";
const string CODE_BETTER = "\e[33m";
const string CODE_NEWOPT = "\e[36m";
const string CODE_END = "\e[0m";


bool check_file(const string& name) {
    if (FILE* file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}


void readSolution(string out) {
    ifstream outfile(out);
    string line;
    getline(outfile, line);
    istringstream iss(line);
    iss >> sol_size;

    iss.clear();
    iss.str(string());

    getline(outfile, line);
    iss.str(line);

    string cur;

    for (int i = 0; i < sol_size; i++) {
        iss >> cur;
        solNames.push_back(cur);
        solSet.insert(cur);
    }

    cout << sol_size << "\n";
    for (int i = 0; i < sol_size; i++) {
        cout << solNames[i] << " ";
    }
    cout << "\n-------------------------\n\n";
    outfile.close();
}

int getNr(string name) {
    if (namesMap.count(name))
        return namesMap[name];
    names.push_back(name);
    p.push_back(-1);
    w.push_back(-1);
    int nr = names.size() - 1;
    namesMap[name] = nr;
    return nr;
}

int ffind(int nr) {
    return (p[nr] < 0) ? nr : p[nr] = ffind(p[nr]);
}

void funion(int nr1, int nr2) {
    nr1 = ffind(nr1);
    nr2 = ffind(nr2);

    if (nr1 == nr2) return;

    if (w[nr1] > w[nr2])
        p[nr2] = nr1;
    else
        p[nr1] = nr2;
    if (w[nr1] == w[nr2])
        w[nr2]++;
}

void saveOutToOpt(string out, string opt) {
    ifstream outfile(out);
    ofstream optfile(opt);
    string line;
    while(getline(outfile, line))
        optfile << line << endl;
}

int main(int argc, char** argv) {
    // Args --->>> Input file Output File

    if (argc < 3) {
        cerr << CODE_WRONG + "Not enough arguments. \n";
        cerr << "Usage: ./Verifier InFile OutFile [Opt]\n" + CODE_END;
        return 1;
    }

    string in = argv[1];

    ifstream infile(in);

    readSolution(argv[2]);

    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        string a, b;
        if (!(iss >> a >> b)) {
            break;
        }  // error
        if (solSet.find(a) != solSet.end() || solSet.find(b) != solSet.end()) {
        	// At least one vertex is in solution
            continue;
        }
        int nr_a = getNr(a), nr_b = getNr(b);

        int ccNra = ffind(nr_a), ccNrb = ffind(nr_b);

        if (ccNra == ccNrb && ccNra != -1) {
            cout << CODE_BOLD + CODE_WRONG + "ERROR: Graph is not acyclic " + CODE_END + "\n";
            infile.close();
            return 2;
        }

        funion(ccNra, ccNrb);
    }
    infile.close();


    int optSize = -1;
    if (argc > 3) {
        string opt = argv[3];
        string line;

        if (!check_file(opt)){
            cout << CODE_BOLD + CODE_NEWOPT << "Creating new opt for " << argv[1] << CODE_END << "\n";
            saveOutToOpt(argv[2], opt);
        }
        else {
            ifstream optfile(opt);
            getline(optfile, line);
            istringstream iss(line);

            iss >> optSize;
            if (optSize < sol_size)
                cout << CODE_BOLD + CODE_WRONG + "[WORSE] " + CODE_END 
                    + "Solution provided is worse than opt = " << optSize << "\n"; 
            else if (optSize > sol_size) {
                cout << CODE_BOLD + CODE_BETTER + "[BETTER] " + CODE_END 
                    + " Solution provided is better than opt = " << optSize << "\n"; 
                saveOutToOpt(argv[2], opt);
            }
            else 
                cout << CODE_BOLD + CODE_OK + "[EQUAL] " + CODE_END;
            optfile.close();
        }
    }
    else
       cout << "OPT not provided\n";

    cout << "Solution for " << in << " is correct, size = " << sol_size << "\n";

    return 0;
}