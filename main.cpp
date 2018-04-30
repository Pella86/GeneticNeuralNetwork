#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
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
The program runs a self generated test and saves self generated input, \
outputs and a config file if the program is run with option:\n\
    -t\n\
Any other option will call for this help (also available with the option -h).";


// helper functions to read the files
io_matrix read_iofile_txt(string filename);
io_matrix read_iofile_bin(string filename);
bool is_text_file(string filename);

int main(int argc, char** argv)
{
    // available options + options with arguments
    string const options = "htdioc";
    string const options_with_args = "dioc";

    map<char, bool> opt_found;
    map<char, string> opt_argument;

    // initialize the maps
    for(string::const_iterator c = options.begin(); c != options.end(); ++c){
        opt_found[*c] = false;
        if(options_with_args.find(*c) != string::npos){
            opt_argument[*c] = "";
        }
    }

    // getopt cycle
    char c;
    while((c = getopt(argc, argv, "htd:i:o:c:")) != -1){
        opt_found[c] = true;
        if(options_with_args.find(c) != string::npos){
            opt_argument[c] = optarg;
        }
    }

    // check if the user gave an option, if is not the case, show help
    if( opt_found['d'] || (opt_found['c'] && opt_found['i'] && opt_found['o']) || opt_found['t']){
        // continue
    }
    else{
        opt_found['h'] = true;
    }

    // initialize genetic algorithm elements
    io_matrix inputs;
    io_matrix outputs;

    GeneticAlgorithm ga;

    if(opt_found['d']){
        cout << "Directory option" << endl;
        cout << opt_argument['d'] << endl;

        string initial_folder(opt_argument['d']);

        // search for the inputs in text format or binary format
        string filename_inputs = initial_folder + "GA_inputs.txt";
        ifstream file_inputs(filename_inputs);
        if(file_inputs.good()){
           inputs = read_iofile_txt(filename_inputs);
        }
        else{
            filename_inputs = initial_folder + "GA_inputs.nni";
            inputs = read_iofile_bin(filename_inputs);
        }

        // search for the outputs either in text or binary format
        string filename_outputs = initial_folder + "GA_outputs.txt";
        ifstream file_outputs(filename_outputs);
        if(file_outputs.good()){
            outputs = read_iofile_txt(filename_outputs);
        }
        else{
            filename_outputs = initial_folder + "GA_outputs.nni";
            outputs = read_iofile_bin(filename_outputs);
        }

        // parse the configuration file
        ga.read_from_file(initial_folder + "GA_config.txt");
    }

    bool is_file_inputs = opt_found['i'] && opt_found['o'] && opt_found['c'];

    if(is_file_inputs){
        cout << "File specification option" << endl;
        cout << "inputs file: " << opt_argument['i'] << endl;
        cout << "outputs file: " << opt_argument['o'] << endl;
        cout << "config file: " << opt_argument['c'] << endl;

        bool read_txt;

        // read the inputs file, either text (.txt) or binary (anything else
        string iarg_s(opt_argument['i']);
        read_txt = is_text_file(iarg_s);
        inputs = (read_txt)? read_iofile_txt(iarg_s) : read_iofile_bin(iarg_s);

        string oarg_s(opt_argument['o']);
        read_txt = is_text_file(oarg_s);
        outputs =(read_txt)? read_iofile_txt(oarg_s) : read_iofile_bin(oarg_s);

        // read the configuration
        ga.read_from_file((string)opt_argument['c']);
    }


    // if ga initialized and only one of the two options (dir init xor file init)
    if(ga.ga_initialized && (opt_found['d'] != is_file_inputs)){

        if(inputs.size() > 0 && outputs.size() > 0 && inputs.size() == outputs.size()){
            ga.run(inputs, outputs);
        }
        else{
            cout << "Inputs or outputs badly formatted." << endl;
        }

    }
    else if(!ga.ga_initialized && (opt_found['d'] != is_file_inputs)){
        cout << "Genetic Algorithm not initialized correctly" << endl;
    }
    else if(opt_found['t']){
        // continue
    }
    else{
        opt_found['h'] = true;
    }

    // test option chosen
    if(opt_found['t']){
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

    // help option triggered or chosen
    if(opt_found['h']){
        cout << help_message << endl;
    }

    return 0;
}

// controls the extention in a path and returns true if is txt
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

// reading the text file
io_matrix read_iofile_txt(string filename){
    ifstream file(filename, ios::in);

    // File format:
    // #: comment
    // empty lines ignored
    // then space separated values

    io_matrix inputs;

    string line;

    if(file.is_open()){

        while(getline(file, line)){

            // skip empty lines and comments
            if(line.size() > 1 && line[0] == '#'){
                continue;
            }

            vector<double> input_line;

            // split the line at spaces
            vector<string> chunks = str_split(line, ' ');

            // keep only the lines that have significant numbers
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

// read the custom bianry format
io_matrix read_iofile_bin(string filename){
    cout << "Reading binary file..." << endl;
    ifstream file(filename, ios::in | ios::binary);

    // File format
    // 512 bytes header
    // size_t m -> row number (io pair number)
    // size_t n -> column number (length of in/output vector)

    streampos co = 512; // skip header

    io_matrix io_pairs;

    if(file.is_open()){
        file.seekg(co);

        size_t n_pairs;
        file.read(reinterpret_cast<char*> (&n_pairs), sizeof(size_t));
        co += sizeof(size_t);
        file.seekg(co);

        size_t n_io;
        file.read(reinterpret_cast<char*> (&n_io), sizeof(size_t));
        co += sizeof(size_t);
        file.seekg(co);

        cout << "Io matrix sizes: " << n_pairs << "x" << n_io << endl;

        // read matrix values
        for(size_t i = 0; i < n_pairs; ++i){
            vector<double> line;
            for(size_t j = 0; j < n_io; ++j){
                double data;
                file.read(reinterpret_cast<char*> (&data), sizeof(double));
                line.push_back(data);
                co += sizeof(double);
                file.seekg(co);
            }
            io_pairs.push_back(line);
        }
        file.close();
    }
    return io_pairs;
}
