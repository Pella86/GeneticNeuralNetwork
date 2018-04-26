#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include<vector>
#include<string>
#include<random>
#include<map>
#include<cassert>
#include<algorithm>

#include "Network.h"

typedef std::vector<std::vector<double>> io_nn_type;

class ScoredNetwork{
public:
    ScoredNetwork(double iscore, Network inn);
    ScoredNetwork();

    double score;
    Network nn;
};


class GeneticAlgorithm
{
    public:

        GeneticAlgorithm();
        GeneticAlgorithm(std::vector<int> n_layers, std::string folder_name);
        GeneticAlgorithm(std::string config_filename);

        virtual ~GeneticAlgorithm();

        void read_from_file(std::string filename);

        void run(io_nn_type inputs, io_nn_type exp_outputs);

        void mutate(ScoredNetwork& pnn);

        double score(Network nn, io_nn_type inputs, io_nn_type exp_results);

        void save_networks(std::vector<ScoredNetwork> population, std::string network_name);

        bool ga_initialized;

    protected:

        int rounds;
        int pop_len;
        size_t retain_n;
        double retain_chance;
        double mut_chance;
        double cross_chance;
        std::vector<Network> init_pop;
        std::vector<int> network_layers;
        std::string output_folder;


    private:
        std::random_device rd{};
        std::default_random_engine rng{rd()};
};

#endif // GENETICALGORITHM_H
