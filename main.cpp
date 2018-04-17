#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <random>
#include <limits>
#include <cassert>
#include <ctime>
#include <map>
#include <algorithm>

using namespace std;

typedef vector<vector<double>> io_nn_type;



double sigmoid(double z){
    if(z < 0){
        double s = 1 - 1 / ( 1 + exp(z));
        return s;
    }
    else{
        double s = 1 / ( 1 + exp(-z));
        return s;
    }
}

double dot(vector<double> l1, vector<double> l2){
    double d = 0;

    for(unsigned int i = 0; i < l1.size(); ++i){
        d += l1[i] * l2[i];
    }
    return d;
}

class Node{
public:
    Node(std::vector<double> inw, double inb){
        w = inw;
        b = inb;
    }

    Node(){
        cout << "Node: empty initialization" << endl;
    }

    Node(string filename){
        from_file(filename);
    }

    double z(vector<double> input){
        return dot(input, w) + b;
    }

    double output(vector<double> input){
        return sigmoid(z(input));
    }

    string to_string(){
        stringstream ss;

        ss << "|w: (" << w.size() << "| ";

        uint32_t i = 0;

        for(; i < w.size() - 1; ++i){
            ss << w[i] << ", ";
        }
        ss << w[i];

        ss << ") |b: " << b;

        return ss.str();
    }

    void from_file(string filename){
        ifstream file(filename, ios::in|ios::binary);
        if(file.is_open()){
            streampos co = 0;
            read_byte_chunk(file, co);


        }
        file.close();
    }

    void read_byte_chunk(ifstream& file, streampos& co){
        streampos init_co = co;

        file.seekg(co);
        file.read(reinterpret_cast<char*>(&b), sizeof(double));
        cout << "Node: read b: " << b << endl;

        co += sizeof(double);
        file.seekg(co);

        int n_weights = 0;
        file.read(reinterpret_cast<char*>(&n_weights), sizeof(int));
        cout << "Node: read n of weights:" << n_weights << endl;

        co += sizeof(int);

        for(int i = 0; i < n_weights; ++i){
            file.seekg(co);

            double weight;
            file.read(reinterpret_cast<char*>(&weight), sizeof(double));
            w.push_back(weight);
            co += sizeof(double);
        }

        cout << "total bytes read: " << co - init_co << endl;
    }

    void save_bin(string filename){
        ofstream file(filename, ios::out|ios::binary);
        if(file.is_open()){
            streampos co = 0;
            write_byte_chunk(file, co);
            file.close();

        }

    }

    void write_byte_chunk(ofstream& file, streampos& co){
        streampos init_co = co;

        file.seekp(co);

        file.write(reinterpret_cast<char*>(&b),sizeof(b));
        co += sizeof(b);

        int n_weights = (int)w.size();
        file.write(reinterpret_cast<char*>(&n_weights), sizeof(n_weights));
        co += sizeof(n_weights);

        cout << "number of weights to write:" << n_weights << endl;
        int i = 0;
        for(; i < n_weights; ++i){
            file.seekp(co);
            file.write(reinterpret_cast<char*>(&w[i]), sizeof(w[i]));
            co += sizeof(w[i]);

            cout << "co: " << co << endl;
        }

        cout << "total bytes written: " << co - init_co << endl;

    }

    std::vector<double> w;
    double b;

};

class Network{
public:

    Network(){
        cout << "Network: empty initialization" << endl;
    }

    Network(string filename){
        cout << "Network: file inizialization" << endl;

        streampos co = 0; // cumulative offset

        ifstream file(filename, ios::in|ios::binary);
        if(file.is_open()){
            // read the first int -> lengths of layers
            int n_layers;
            file.read(reinterpret_cast<char*>(&n_layers), sizeof(int));

            cout << "Network: has n layers: " << n_layers << endl;

            co += sizeof(int);

            // read how many nodes per layers are there
            vector<int> nodes_layer;
            for(int i = 0; i < n_layers; ++i){
                file.seekg(co);
                int n_nodes;
                file.read(reinterpret_cast<char*>(&n_nodes), sizeof(int));
                nodes_layer.push_back(n_nodes);
                co += sizeof(int);
                cout << n_nodes << ", ";
            }
            cout << endl;

            // for each layer
            for(int i = 0; i < n_layers; i += 1){

                vector<Node> nv;
                layers.push_back(nv);
                // for each node
                for(int inode = 0; inode < nodes_layer[i]; ++inode){
                    Node n;
                    n.read_byte_chunk(file, co);
                    cout << "node read:" << endl;
                    cout << n.to_string() << endl;
                    layers[i].push_back(n);
                }
            }
        }
        cout << "total bytes read: " << co << endl;
        file.close();
    }

    Network(vector<int> node_layers, default_random_engine& rng){

        cout << "Network: vector initialization" << endl;



        normal_distribution<double> distribution(0.0, 1.0);




        for(size_t ilayer = 1; ilayer < node_layers.size(); ++ilayer){
            vector<Node> nv;
            layers.push_back(nv);

            int n_nodes = node_layers[ilayer];

            for(int inode = 0; inode < n_nodes; ++inode){

                double b = distribution(rng);
                vector<double> weights;

                for(int iweight = 0; iweight < node_layers[ilayer - 1]; ++iweight){
                    double w = distribution(rng);
                    weights.push_back(w);
                }

                Node n = Node(weights, b);
                layers[ilayer - 1].push_back(n);
            }
        }

    }

    string to_string(){
        stringstream ss;

        ss << "- Network -" << endl;

        ss << "n layers:" << layers.size() << endl;
        ss << "nodes per layer" << endl;
        for(auto i : layers) ss << i.size() << " ";
        ss << endl;

        for(auto l : layers){
            ss << "------" << endl;
            for(auto n : l){
                ss << n.to_string() << endl;
            }
        }

        return ss.str();
    }


    void save_to_file(string filename){
        ofstream file(filename, ios::out | ios::binary);
        if(file.is_open()){
            streampos co = 0;
            int n_layers = layers.size();
            cout << "n_layers:" << n_layers << endl;
            file.write(reinterpret_cast<char*> (&n_layers), sizeof(n_layers));
            co += sizeof(n_layers);

            for(vector<Node> i : layers) {
                int n_nodes = i.size();

                cout << n_nodes << endl;

                file.seekp(co);
                file.write(reinterpret_cast<char*> (&n_nodes), sizeof(n_nodes));
                co += sizeof(n_nodes);

            }

            for(auto l : layers){
                for(auto n : l){
                    cout << "Writing node" << endl;
                    cout << n.to_string() << endl;
                    n.write_byte_chunk(file, co);
                }
            }

            cout << "total bytes written: " << co <<endl;

        }
        file.close();
    }

    vector<double> calculate(vector<double> input){


        vector<double> dyn_output;
        vector<double> dyn_input = input;
        for (auto l : layers){
            dyn_output.clear();
            for(auto n : l){
                dyn_output.push_back(n.output(dyn_input));
            }
            dyn_input.clear();
            dyn_input = dyn_output;
        }

        return dyn_output;

    }

    vector<vector<Node>> layers;
};



class GeneticAlgorithm{
public:
    GeneticAlgorithm(vector<int> n_layers){
        pop_len = 20;
        retain_n = 10;
        retain_chance = 0.1;
        mut_chance = 0.5;
        cross_chance = 0.5;
        network_layers = n_layers;
        rounds = 1000;
    }

    void run(io_nn_type inputs, io_nn_type exp_outputs){
        cout << "Running genetic algorithm..." << endl;


        cout << "Input/expected output pairs for nn test" << endl;
        assert(inputs.size() == exp_outputs.size());

        // print to console inputs vs expected outputs
        for(size_t i = 0; i < inputs.size(); ++i){
            cout << "input: " << endl;
            for(auto j : inputs[i]) cout << j;
            cout << endl;

            cout << "expected output" << endl;
            for(auto j : exp_outputs[i]) cout << j;
            cout << endl;

        }

        // create a population of random initialized neural networks
        cout << "Creating population..." << endl;

        vector<Network> population;
        for(int i = 0; i < pop_len; ++i){
            Network n = Network(network_layers, rng);
            population.push_back(n);
        }

        init_pop = population;

        // score the population
        cout << "Scoring initial population..." << endl;

        vector<pair<double, Network>> scored_pop;
        for(auto nn : population){
            double nn_score = score(nn, inputs, exp_outputs);
            pair<double, Network> p(nn_score, nn);
            scored_pop.push_back(p);
        }

        // sort based on score
        cout << "Sorting based on score..." << endl;

        // custom score sorter to access the pair variables
        struct{
            bool operator()(pair<double, Network> a, pair<double, Network> b){
                return a.first < b.first;
            }
        } CustomScoreSorting;
        sort(scored_pop.begin(), scored_pop.end(), CustomScoreSorting);

        // run the simulation
        cout << "Running the generations..." << endl;

        for(int i = 0; i < rounds; ++i){
            cout << "--------- generation " << i << "-------------" << endl;
            vector<pair<double, Network>> next_gen;

            // add the first retain_n neural networks to next generation
            assert(retain_n < scored_pop.size());

            for(size_t i = 0; i < retain_n; ++i){
                next_gen.push_back(scored_pop[i]);
            }

            // pick random neural networks from bottm
            uniform_real_distribution<double> rand_num(0.0, 1.0);

            for(size_t i = retain_n; i < scored_pop.size(); ++i){
                if(rand_num(rng) < retain_chance){
                    next_gen.push_back(scored_pop[i]);
                }
            }

            // mutate the population
            int mut_count = 0;
            for(auto pnn : next_gen){
                if(rand_num(rng) < mut_chance){
                    mutate(pnn);
                    mut_count += 1;
                }
            }
            cout << "mutated " << mut_count << " individuals" << endl;

            // cross population

            // calculate number of children to come up to the pop length
            size_t n_children = pop_len - next_gen.size();
            cout << "generating " << n_children << " children" << endl;

            // prepare vector to contain the children
            vector<pair<double, Network>> children;

            // inizialize the random parent number distribution
            uniform_int_distribution<size_t> rand_child(0, next_gen.size() - 1);

            while(children.size() < n_children){
                // chose a male and a female
                size_t male_n = rand_child(rng);
                size_t female_n = rand_child(rng);

                // skip if male and female are the same
                if(male_n != female_n){
                    pair<double, Network> male = next_gen[male_n];
                    pair<double, Network> female = next_gen[female_n];

                    pair<double, Network> child = female;
                    // reset score
                    child.first = -1;

                    // chose which node to give, the child is default female
                    // but by random choice inherits from male
                    for(size_t lit = 0; lit < child.second.layers.size(); ++lit){
                        for(size_t nit = 0; nit < child.second.layers[lit].size(); ++nit){
                            if(rand_num(rng) < cross_chance){
                                child.second.layers[lit][nit] = male.second.layers[lit][nit];
                            }
                        }
                    }
                    children.push_back(child);
                }

            }
            next_gen.insert(next_gen.end(), children.begin(), children.end());

            // score and sort
            cout << "scoring new generation" << endl;

            for(size_t ipnn = 0; ipnn < next_gen.size(); ++ipnn){
                if(next_gen[ipnn].first < 0){
                    next_gen[ipnn].first = score(next_gen[ipnn].second, inputs, exp_outputs);
                }
            }
            cout << "sorting new generation" << endl;

            sort(next_gen.begin(), next_gen.end(), CustomScoreSorting);

            scored_pop.clear();
            scored_pop = next_gen;

        }

        // proof of concept
        vector<Network> result_pop;
        for(auto pnn : scored_pop){
            result_pop.push_back(pnn.second);
        }

        for(size_t i = 0; i < inputs.size(); ++i){
            cout << "---- input: ------" << endl;
            for(auto j : inputs[i]) cout << j << " ";
            cout << endl;

            cout << "expected output" << endl;
            for(auto j : exp_outputs[i]) cout << j;
            cout << endl;

            vector<double> output;
            for(size_t j = 0; j < 1/*scored_pop.size()*/; ++j){
                output = init_pop[j].calculate(inputs[i]);
                cout << "init pop" << endl;
                cout << "score: " << score(init_pop[j], inputs, exp_outputs) << endl;
                cout << "output: ";
                for(auto o : output){
                    cout << o << " ";
                }
                cout << endl;
                output.clear();

                cout << "----" << endl;

                output = result_pop[j].calculate(inputs[i]);
                cout << "refined pop" << endl;
                cout << "score: " << score(result_pop[j], inputs, exp_outputs) << endl;
                cout << "output: ";
                for(auto o : output) {
                    cout << o << " ";
                }
                cout << endl;
                output.clear();


            }

        }

    }

    void mutate(pair<double, Network>& pnn){
        // reset score
        pnn.first = -1;

        // proceed with the mutation
        Network& nn = pnn.second;

        // pick a random layer
        uniform_int_distribution<size_t> rand_layer(0, nn.layers.size() - 1);
        size_t layer_n = rand_layer(rng);

        // pick a random node in that layer
        vector<Node> nodes = nn.layers[layer_n];
        uniform_int_distribution<size_t> rand_node(0, nodes.size() - 1);
        size_t node_n = rand_node(rng);

        // pick the randomly chosed node
        Node& n = nn.layers[layer_n][node_n];

        // create the normal distribution for the weights and biases
        normal_distribution<double> rand_n(0.0, 1.0);

        // assign bias
        n.b = rand_n(rng);

        // assign weight
        size_t nweights = n.w.size();
        for(size_t iw = 0; iw < nweights; ++iw){
            n.w[iw] = rand_n(rng);
        }

    }

    double score(Network nn, io_nn_type inputs, io_nn_type exp_results){
        // scores a given neural network nn using inputs and comparing them to
        // expected results with the least squared difference method

        double d = 0;
        for(size_t i = 0; i < inputs.size(); ++i){
            // calculate the output from the elements of inputs
            vector<double> network_output = nn.calculate(inputs[i]);

            // assign to a place holder the right expected result
            vector<double> exp_res = exp_results[i];

            assert(network_output.size() == exp_res.size());

            // calculate the differences
            for(size_t j = 0; j < network_output.size(); ++j){
                d += pow(exp_res[j] - network_output[j], 2);
            }
        }
        return d;
    }

    int rounds;
    int pop_len;
    size_t retain_n;
    double retain_chance;
    double mut_chance;
    double cross_chance;
    vector<Network> init_pop;
    vector<int> network_layers;

private:
    random_device rd{};
    default_random_engine rng{rd()};
};


int main()
{
//    float memblock;
//
//    ifstream file("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/test_bin_file_2", ios::in|ios::binary);
//
//    if(file.is_open()){
//        cout << "file open" << endl;
//        for(unsigned int i = 0; i < 10; ++i){
//            file.read(reinterpret_cast<char*>(&memblock), sizeof(float));
//            cout << memblock << endl;
//            cout << i*sizeof(float) << endl;
//            streampos size = i*sizeof(float);
//            file.seekg(size);
//
//        }
//
//
//        file.close();
//    }

//    Node n = Node("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/my_node");
//
//    cout << "bias" << endl;
//    cout << n.b << endl;
//
//    cout << "weights" << endl;
//    for(size_t i = 0; i < n.w.size(); ++i){
//
//        cout << n.w[i] << endl;
//    }
//
//    n.save_bin("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/my_c_node");

//    Network nn = Network("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/nn_test");
//
//    cout << nn.to_string() << endl;
//
//    nn.save_to_file("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/nn_test_CPP_write");
//
//    Network nn_test = Network("C:/Users/Mauro/Desktop/Vita Online/Programming/Neural Networks/nn_test_CPP_write");
//
//    cout << nn_test.to_string() << endl;
//
//    vector<double> output = nn_test.calculate({1,2});
//    cout << "output" << endl;
//    for(auto d : output) {
//            cout << d << endl;
//    }
//
//    vector<int> nodes_per_layer({2, 3, 3, 1});
//    Network rand_net(nodes_per_layer);
//
//    cout << rand_net.to_string() << endl;



    vector<int> n_layers = {4, 50, 50, 1};
    GeneticAlgorithm ga(n_layers);

    io_nn_type inputs;
    io_nn_type outputs;
    int number_of_inputs = 4;
    for(int i = 0; i < number_of_inputs; ++i){
        vector<double> res;
        for(int j = 0; j < 4; ++j){
            if(j  == i % 4){
                res.push_back(1.0);
            }
            else{
                res.push_back(0.0);
            }
        }
        inputs.push_back(res);
        vector<double> output;
        double perc = i % 4;
        output.push_back( perc / 4);
        outputs.push_back(output);
        res.clear();
    }

    ga.run(inputs, outputs);





    cout << "Hello world!" << endl;
    return 0;
}
