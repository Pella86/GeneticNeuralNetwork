#include "GeneticAlgorithm.h"
#include "StringtoDouble.h"
#include "TextParseHelper.h"

#include <sstream>
#include <cmath>

using namespace std;

class ScoredNetwork{
public;
    Network nn;
    double score;
};


GeneticAlgorithm::~GeneticAlgorithm()
{
    //dtor
}

GeneticAlgorithm::GeneticAlgorithm(vector<int> n_layers, string folder_name){
    pop_len = 20;
    retain_n = 5;
    retain_chance = 0.1;
    mut_chance = 0.5;
    cross_chance = 0.5;
    network_layers = n_layers;
    rounds = 10;
    output_folder = folder_name;
}


GeneticAlgorithm::GeneticAlgorithm(string config_file){
    cout << "reading config file" << endl;

    ifstream file(config_file, ios::in);

    map<string, string> data; // contains the couple identifier value

    if(file.is_open()){
        string line;

        // each line has a parameter
        while(getline(file, line)){

            // skip empty lines and comments
            if(line.size() > 1 && line[0] != '#'){
                vector<string> sstring = str_split(line, '=');

                string identifier = strip(sstring[0]);
                string value = strip(sstring[1]);

                cout << "[" << identifier << "]" << "[" << value << "]" << endl;

                data[identifier] = value;
            }
        }

        file.close();
    }

    this->pop_len = mystoi<int>(data["pop_len"]);
    this->retain_n = mystoi<size_t>(data["retain_n"]);
    this->rounds = mystoi<int>(data["rounds"]);
    this->mut_chance = stod(data["mut_chance"]);
    this->cross_chance = stod(data["cross_chance"]);
    this->retain_chance = stod(data["retain_chance"]);

    vector<string> nodes_layer_str = str_split(data["network_layers"], ' ');

    for(auto s : nodes_layer_str){
        this->network_layers.push_back(mystoi<int>(s));
        cout << s << endl;
    }
    this->output_folder = data["output_folder"];
}

void GeneticAlgorithm::run(io_nn_type inputs, io_nn_type exp_outputs){
    cout << "Running genetic algorithm..." << endl;
    double first_best;


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

    first_best = scored_pop.begin()->first;

    cout << "Saving networks..." << endl;

    int name_counter = 0;
    for(auto pnn : scored_pop){
        string nn_name = output_folder + "/init_nn_" + to_string(name_counter);
        pnn.second.save_to_file(nn_name);
        ++name_counter;
    }

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


        for(size_t i = 0; i < retain_n; ++i){
            next_gen.push_back(mutate(next_gen[i]));
        }
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
        cout << "Algorithm performance, first best " << first_best << endl;
        cout << "First: " << next_gen.begin()->first << " Last: " << next_gen.back().first << endl;

    }

    // proof of concept
    vector<Network> result_pop;
    for(auto pnn : scored_pop){
        result_pop.push_back(pnn.second);
    }

    cout << "Saving resulting networks" << endl;
    name_counter = 0;
    for(auto nn : result_pop){
        string nn_name = output_folder + "/final_nn_" + to_string(name_counter);
        nn.save_to_file(nn_name);
        ++name_counter;
    }

    cout << "score: " << score(init_pop[0], inputs, exp_outputs) << endl;
    cout << "score: " << score(result_pop[0], inputs, exp_outputs) << endl;

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

            cout << "output: ";
            for(auto o : output){
                cout << o << " ";
            }
            cout << endl;
            output.clear();

            cout << "----" << endl;

            output = result_pop[j].calculate(inputs[i]);
            cout << "refined pop" << endl;

            cout << "output: ";
            for(auto o : output) {
                cout << o << " ";
            }
            cout << endl;
            output.clear();


        }

    }

}

pair<double, Network> GeneticAlgorithm::mutate(pair<double, Network> pnn){

    pair<double, Network> res_nn;
    // reset score
    res_nn.first = -1;

    // proceed with the mutation
    Network nn = pnn.second;

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
    normal_distribution<double> rand_n(0.0, 2.5);

    // assign bias
    n.b = rand_n(rng);

    // assign weight
    size_t nweights = n.w.size();
    for(size_t iw = 0; iw < nweights; ++iw){
        n.w[iw] = rand_n(rng);
    }
    res_nn.second = nn;
    return res_nn;

}

double GeneticAlgorithm::score(Network nn, io_nn_type inputs, io_nn_type exp_results){
    // scores a given neural network nn using inputs and comparing them to
    // expected results with the least squared difference method

    double d = 0;
    for(size_t i = 0; i < inputs.size(); ++i){
        // calculate the output from the elements of inputs
        vector<double> network_output = nn.calculate(inputs[i]);

        // assign to a place holder the right expected result
        vector<double> exp_res = exp_results[i];



        // calculate the differences
        assert(network_output.size() == exp_res.size());
        for(size_t j = 0; j < network_output.size(); ++j){
            d += pow(exp_res[j] - network_output[j], 2);
        }
    }
    return d;
}


