#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <unistd.h>


#include "GeneticAlgorithm.h"
#include "TextParseHelper.h"
#include "StringtoDouble.h"

using namespace std;

template<class T>
string vec2str(vector<T> v){
    stringstream ss;
    ss << "(";
    for(auto i : v){ ss << i << ", ";}
    ss << ")";
    return ss.str();
}

void test_node(){
    cout << "Test node..." << endl;

    // generate random node
    random_device rd;
    default_random_engine rng(rd());

    // create a node with 10 weights and the
    // normal distr param mu as 0 and std as 3
    Node node(10, rng, 0, 3);

    cout << node << endl;

    // save the node
    node.save_bin("./test_rnd_node.ndf");

    // read the saved random node
    Node read_node("./test_rnd_node.ndf");

    cout << read_node << endl;
}

void test_network(){
    cout << "Testing network..." << endl;
    random_device rd;
    default_random_engine rng(rd());

    Network nn({2, 3, 4, 5}, rng, 0, 4);

    cout << nn << endl;

    nn.save_to_file("./test_net.nnf");

    Network nn_read("./test_net.nnf");

    cout << nn_read << endl;

}

void test_genetic_algorithm(){
    int input_n = 1;
    int output_n = 4;

    io_nn_type inputs;
    io_nn_type outputs;
    int n_io_pairs = 8;
    for(int i = 0; i < n_io_pairs; ++i){
        /*input definition*/
        vector<double> input;
         // put the number as input
        for(int j = 0; j < input_n; ++ j){
            input.push_back(i);
        }
        inputs.push_back(input);


        /*output definition*/
        vector<double> output;
        int mask = 0b00000001;

        for(int j = 0; j < output_n; ++j){
            if( (i & mask) > 0){ // scan i with the mask to show the sig bit
                output.push_back(1);
            }
            else{
                output.push_back(0);
            }
            mask <<= 1;
        }
        //reverse output so it follows a big endian layout
        reverse(output.begin(), output.end());

        outputs.push_back(output);

        // clear the input/output place holder (is anyhow copied in the corr vec
        input.clear();
        output.clear();
    }

    cout << "Final vectors" << endl;

    for(size_t i = 0; i < inputs.size(); ++i){
        cout << vec2str<double>(inputs[i]) << endl;
        cout << vec2str<double>(outputs[i]) << endl;
    }

    vector<int> n_layers = {input_n, 32, 16, 8, output_n};
    string output_folder = "./test_run/";
    GeneticAlgorithm ga(n_layers, output_folder);

    ga.run(inputs, outputs);
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

int main(int argc, char** argv)
{

    // se � una cartella scan for the files -d

    int c;

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

    if(optc || optd || opti || opto || optt){
        // continue
    }
    else{
        opth = true;
    }

    vector<vector<double>> inputs;
    vector<vector<double>> outputs;

    GeneticAlgorithm ga;


    if(optd){
        cout << "d option: " << darg << endl;

        string darg_s(darg);

        // look in the folder for config file, input, output

        inputs = read_iofile_txt(darg_s + "GA_inputs.txt");
        outputs = read_iofile_txt(darg_s + "GA_outputs.txt");
        ga.read_from_file(darg_s + "GA_config.txt");


    }

    if(opti && opto && optc){

        inputs = read_iofile_txt((string)iarg);
        outputs = read_iofile_txt((string)oarg);
        ga.read_from_file((string)carg);
    }

   if(ga.ga_initialized){

        if(inputs.size() > 0 && outputs.size() > 0 && inputs.size() == outputs.size()){
            ga.run(inputs, outputs);
        }
        else{
            cout << "Inputs or outputs badly formatted." << endl;
        }

    }
    else{
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

    }

    if(opth){

        cout << "- Usage -" << endl;
    }


    // se � -i inputs -o outputs -c config
    // -h help
    // -t test




//    bool test = false;
//
//
//
//    cout << "test input read" << endl;
//
//    vector<vector<double>> inputs = read_iofile_txt("./data/GA_inputs.txt");
//
//    cout << "retrived inputs" << endl;
//    for(auto v : inputs){
//        cout << vec2str(v) << endl;
//    }
//
//    vector<vector<double>> outputs = read_iofile_txt("./data/GA_outputs.txt");
//
//    cout << "retrived outputs" << endl;
//    for(auto v : outputs){
//        cout << vec2str(v) << endl;
//    }
//
//    cout << "test genetic algorithm" << endl;
//
//    GeneticAlgorithm ga("./data/GA_config.txt");
//
//    ga.run(inputs, outputs);
//
//    if(test){
//        /*Node testing functions*/
//        test_node();
//
//        /*Network testing functions*/
//        test_network();
//
//        /*Genetic algorithm tests*/
//        test_genetic_algorithm();
//    }

    return 0;
}
