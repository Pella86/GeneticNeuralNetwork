#ifndef NODE_H
#define NODE_H

#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

class Node
{
    public:
        Node();
        Node(std::vector<double> inw, double inb);
        Node(std::string filename);
        Node(size_t n_weights, std::default_random_engine& rng, double mu = 0.0, double sigma = 1.0);

        virtual ~Node();

        double z(std::vector<double> input);
        double output(std::vector<double> input);
        void from_file(std::string filename);
        void read_byte_chunk(std::ifstream& file, std::streampos& co);
        void save_bin(std::string filename);
        void write_byte_chunk(std::ofstream& file, std::streampos& co);
        friend std::ostream& operator<< (std::ostream& os, const Node& node);

        std::vector<double> w;
        double b;

    protected:

    private:

};

#endif // NODE_H
