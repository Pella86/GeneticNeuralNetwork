#include "TestHelper.h"
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include "Node.h"
#include "Network.h"
#include "GeneticAlgorithm.h"

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

    io_matrix inputs;
    io_matrix outputs;
    int n_io_pairs = 8;
    for(int i = 0; i < n_io_pairs; ++i){
        /*input definition
            as input there is a number
        */

        vector<double> input;
         // put the number as input
        for(int j = 0; j < input_n; ++ j){
            input.push_back(i);
        }
        inputs.push_back(input);

        /*output definition
            as output there's the binary representation of that number
        */
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

void create_example_files(){

    // config file
    string config_example = "\
# genetic algorithm configuration\n\npop_len = 20\nretain_n = 5\n\
retain_chance = 0.1\nmut_chance = 0.5\ncross_chance = 0.5\n\
network_layers = 4 32 16 2\nrounds = 10\noutput_folder = ./results/\n";

    ofstream config_file("./GA_config.txt", ios::out);
    if(config_file.is_open()){
        config_file << config_example << endl;
        config_file.close();
    }

    // inputs file
    string inputs_example = "\
# input table\n\n\n0 0 0 1\n0 0 1 0\n0 1 0 0\n1 0 0 0\n";

    ofstream inputs_file("./GA_inputs.txt", ios::out);
    if(inputs_file.is_open()){
        inputs_file << inputs_example << endl;
        inputs_file.close();
    }

    // output file
    string outputs_example = "\
# output table\n\n\n0 0\n0 1\n1 0\n1 1\n";

    ofstream outputs_file("./GA_outputs.txt", ios::out);
    if(outputs_file.is_open()){
        outputs_file << outputs_example << endl;
        outputs_file.close();
    }

    // create input binary file
    ofstream bin_inputs_file("./GA_bin_inputs.nni");

    size_t n_pairs = 4;
    size_t n_inputs = 4;

    streampos co = 512;
    if(bin_inputs_file.is_open()){
        bin_inputs_file.seekp(co);

        //first is n
        bin_inputs_file.write(reinterpret_cast<char*> (&n_pairs), sizeof(size_t));
        co += sizeof(size_t);

        bin_inputs_file.seekp(co);

        //second is m
        bin_inputs_file.write(reinterpret_cast<char*> (&n_inputs), sizeof(size_t));
        co += sizeof(size_t);
        bin_inputs_file.seekp(co);


        for(size_t i = 0; i < n_pairs; ++i){

            for(size_t j = 0; j < n_inputs; ++j){
                double data = (i == j)? 1 : 0;
                bin_inputs_file.write(reinterpret_cast<char*> (&data), sizeof(double));
                co += sizeof(double);
                bin_inputs_file.seekp(co);
            }
        }
        bin_inputs_file.close();
    }

    // create binary output file
    // mxn matrix
    ofstream bin_outputs_file("./GA_bin_outputs.nni");

    size_t n_outputs = 2;

    co = 512; // 512 header
    if(bin_outputs_file.is_open()){
        bin_outputs_file.seekp(co);

        //first integer is n
        bin_outputs_file.write(reinterpret_cast<char*> (&n_pairs), sizeof(size_t));
        co += sizeof(size_t);
        bin_outputs_file.seekp(co);

        //second is m
        bin_outputs_file.write(reinterpret_cast<char*> (&n_outputs), sizeof(size_t));
        co += sizeof(size_t);
        bin_outputs_file.seekp(co);

        for(size_t i = 0; i < n_pairs; ++i){
            int mask = 0b00000001;
            for(size_t j = 0; j < n_outputs; ++j){
                double data = ( (i&mask) > 0)? 1 : 0;
                bin_outputs_file.write(reinterpret_cast<char*> (&data), sizeof(double));
                co += sizeof(double);
                bin_outputs_file.seekp(co);
            }
        }
        bin_outputs_file.close();
    }
}
