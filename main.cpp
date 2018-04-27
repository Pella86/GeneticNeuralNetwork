#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h> //getopt

#include "GeneticAlgorithm.h"
#include "TextParseHelper.h"
#include "StringtoDouble.h"
#include "TestHelper.h"

using namespace std;

constexpr const char help_message[] = "\
- Usage -\nThe program runs a genetic algorithm to refine a neural network.\n\
The genetic algorithm parameters have to be specified into the file:\n\
    GA_config.txt file\n\
The input and output must be a (nxm) matrix of double numbers named:\n\
    GA_inputs.txt and GA_outputs.txt.\nIf the files are given in a folder the \
program runs by giving the option:\n\
    -d path/to/folder/.\n\
The program can run also by giving the files separately with the options:\
    -c path/to/config/file.txt\n\
    -i path/to/inputs/file.txt\n\
    -o path/to/outputs/file.tx\n\
The program runs a self generated test and saves self generated input and\
outputs if the program is run with option:\n\
    -t\n\
Any other option will call for this help (also available with the option -h).";

vector<vector<double>> read_iofile_txt(string filename);
vector<vector<double>> read_iofile_bin(string filename);
bool is_text_file(string filename);

int main(int argc, char** argv)
{
    char c;
    bool optd = false;
    string darg;
    bool opti = false;
    string iarg;
    bool opto = false;
    string oarg;
    bool optc = false;
    string carg;
    bool opth = false;
    bool optt = false;

    while( (c = getopt(argc, argv, "htd:i:o:c:") ) != -1 ){
        switch(c){
            case 'd':
                optd = true;
                darg = optarg;
            break;
            case 'i':
                opti = true;
                iarg = optarg;
            break;
            case 'o':
                opto = true;
                oarg = optarg;
            break;
            case 'c':
                optc = true;
                carg = optarg;
            break;
            case 't':
                optt = true;
            break;
            case 'h':
            default:
                opth = true;
            break;

        }
    }

    if( optd || (optc && opti && opto) || optt){
        // continue
    }
    else{
        opth = true;
    }

    vector<vector<double>> inputs;
    vector<vector<double>> outputs;

    GeneticAlgorithm ga;

    if(optd){
        cout << "Directory option" << endl;
        cout << darg << endl;

        string darg_s(darg);

        string filename_inputs = darg + "GA_inputs.txt";
        ifstream file_inputs(filename_inputs);
        if(file_inputs.good()){
           inputs = read_iofile_txt(filename_inputs);
        }
        else{
            filename_inputs = darg + "GA_inputs.nni";
            inputs = read_iofile_bin(filename_inputs);
        }

        string filename_outputs = darg + "GA_outputs.txt";
        ifstream file_outputs(filename_outputs);
        if(file_outputs.good()){
           inputs = read_iofile_txt(filename_outputs);
        }
        else{
            filename_outputs = darg + "GA_outputs.nni";
            outputs = read_iofile_bin(filename_outputs);
        }

        ga.read_from_file(darg_s + "GA_config.txt");
    }

    if(opti && opto && optc){
        cout << "File specification option" << endl;
        cout << "inputs file: " << iarg << endl;
        cout << "outputs file: " << oarg << endl;
        cout << "config file: " << carg << endl;

        bool read_txt;

        string iarg_s(iarg);
        read_txt = is_text_file(iarg_s);
        inputs = (read_txt)? read_iofile_txt(iarg_s) : read_iofile_bin(iarg_s);

        string oarg_s(oarg);
        read_txt = is_text_file(oarg_s);
        outputs =(read_txt)? read_iofile_txt(oarg_s) : read_iofile_bin(oarg_s);

        ga.read_from_file((string)carg);
    }

   if(ga.ga_initialized && (optd || (optc && opti && opto))){

        if(inputs.size() > 0 && outputs.size() > 0 && inputs.size() == outputs.size()){
            ga.run(inputs, outputs);
        }
        else{
            cout << "Inputs or outputs badly formatted." << endl;
        }

    }
    else if(!ga.ga_initialized && (optd || (optc && opti && opto))){
        cout << "Genetic Algorithm not initialized correctly" << endl;
    }

    if(optt){
        cout << "Test run" << endl;

        cout << "Testing node" << endl;
        /*Node testing functions*/
        test_node();

        cout << "Testing algorithm" << endl;
        /*Network testing functions*/
        test_network();

        cout << "Test genetic algorithm" << endl;
        /*Genetic algorithm tests*/
        test_genetic_algorithm();

        cout << "Saving examples files" << endl;
        create_example_files();
    }

    if(opth){

        cout << help_message << endl;
    }

    return 0;
}

bool is_text_file(string filename){
    vector<string> filenameext;

    filenameext = str_split(filename, '.');

    if(filenameext.size() >= 2){
        if(filenameext[1] == "txt"){
            return true;
        }
    }
    return false;
}

vector<vector<double>> read_iofile_txt(string filename){
    ifstream file(filename, ios::in);

    vector<vector<double>> inputs;

    string line;

    if(file.is_open()){

        while(getline(file, line)){

            // skip empty lines and comments
            if(line.size() > 1 && line[0] == '#'){
                continue;
            }

            vector<double> input_line;

            vector<string> chunks = str_split(line, ' ');

            for(auto s : chunks){
                if(s.size() > 0){
                    input_line.push_back(stod(s));
                }
            }

            if(input_line.size() > 0){
                inputs.push_back(input_line);
            }
        }

        file.close();
    }

    return inputs;
}

vector<vector<double>> read_iofile_bin(string filename){

    ifstream file(filename, ios::in | ios::binary);
    streampos co = 512;

    vector<vector<double>> io_pairs;

    if(file.is_open()){
        file.seekg(co);

        size_t n_pairs;
        file.read(reinterpret_cast<char*> (&n_pairs), sizeof(size_t));
        co += sizeof(size_t);

        size_t n_io;
        file.read(reinterpret_cast<char*> (&n_io), sizeof(size_t));
        co += sizeof(size_t);



        for(size_t i = 0; i < n_pairs; ++i){
            vector<double> line;
            for(size_t j = 0; j < n_io; ++j){
                double data;
                file.read(reinterpret_cast<char*> (&data), sizeof(double));
                line.push_back(data);

            }
            io_pairs.push_back(line);
        }


        file.close();
    }

    return io_pairs;

}
