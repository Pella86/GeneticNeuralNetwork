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
    // generate random node
    Node node = Node();

}


int main()
{
/*Node testing functions*/

//    cout << "loadain node from file" << endl;
//
//    Node n = Node("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/my_node");
//
//    cout << n << endl;
//
//    cout << "saving node to file" << endl;
//
//    n.save_bin("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/my_c_node");
//
//    cout << "Random node" << endl;
//
//    random_device rd;
//    default_random_engine rng(rd());
//
//    Node random_node = Node(10, rng);
//    cout << random_node << endl;
//    random_node.save_bin("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/my_c_node");

/*Network testing functions*/

//    Network nn = Network("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/nn_test");
//
//    cout << nn << endl;
//
//    nn.save_to_file("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/nn_test_CPP_write");
//
//    Network nn_test = Network("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/nn_test_CPP_write");
//
//    cout << nn_test << endl;
//
//    vector<double> output = nn_test.calculate({1,2});
//    cout << "output" << endl;
//    for(auto d : output) {
//            cout << d << endl;
//    }
//
//    random_device rd;
//    default_random_engine rng(rd());
//
//    vector<int> nodes_per_layer({2, 3, 3, 1});
//    Network rand_net(nodes_per_layer, rng);
//
//    cout << rand_net << endl;
//
//    output = rand_net.calculate({1,2});
//    cout << "output" << endl;
//    for(auto d : output) {
//            cout << d << endl;
//    }

/*Genetic algorithm tests*/

//    vector<int> n_layers = {4, 50, 50, 50, 50, 1};
//    GeneticAlgorithm ga(n_layers);
//
//    io_nn_type inputs;
//    io_nn_type outputs;
//    int number_of_inputs = 4;
//    for(int i = 0; i < number_of_inputs; ++i){
//        vector<double> res;
//        for(int j = 0; j < 4; ++j){
//            if(j  == i % 4){
//                res.push_back(1.0);
//            }
//            else{
//                res.push_back(0.0);
//            }
//        }
//        inputs.push_back(res);
//        vector<double> output;
//        double perc = i % 4;
//        output.push_back( perc / 4);
//        outputs.push_back(output);
//        res.clear();
//    }
//
//    ga.run(inputs, outputs);

/*Genetic algorithm byte converter test*/

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

    return 0;
}
