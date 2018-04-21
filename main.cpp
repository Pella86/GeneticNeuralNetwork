#include <iostream>
#include <vector>
#include <string>
#include <sstream>


//#include "Node.h"
//#include "Network.h"
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
    string folder_name = "C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks ref/data/ga_cpp_run";
    GeneticAlgorithm ga(n_layers, folder_name);

    ga.run(inputs, outputs);

}


int main()
{
    bool test = false;


    GeneticAlgorithm ga("./data/GA_config.txt");

    if(test){
        /*Node testing functions*/
        test_node();

        /*Network testing functions*/
        test_network();

        /*Genetic algorithm tests*/
        test_genetic_algorithm();
    }
    return 0;
}
