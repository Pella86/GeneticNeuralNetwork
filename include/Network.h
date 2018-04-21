#ifndef NETWORK_H
#define NETWORK_H

#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <random>

#include "Node.h"

class Network
{
    public:
        Network();
        // open the network from file
        Network(std::string filename);

        // initialize network by layers ex: {2, 3, 4, 5}
        Network(std::vector<int> node_layers, std::default_random_engine& rng, double mu = 0.0, double sigma = 1.0);

        virtual ~Network();

        // save network to file
        void save_to_file(std::string filename);

        // calculate network output
        std::vector<double> calculate(std::vector<double> input);

        friend std::ostream& operator<<(std::ostream& os,  const Network& nn);

        // members
        std::vector<std::vector<Node>> layers;
    protected:

    private:

};

#endif // NETWORK_H
