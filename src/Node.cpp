#include "Node.h"

#include <cmath>

using namespace std;

/*******************************************************************************
 Sigmoid function
   converts the input in a 0 to 1 result
   the separation between negative z is to avoid exp will have too large numbers
   todo: insert a math exception to handle overflows
********************************************************************************/
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

/*******************************************************************************
 Dot product fuction
   between two double vectors
*******************************************************************************/


double dot(vector<double> l1, vector<double> l2){
    double d = 0;

    for(unsigned int i = 0; i < l1.size(); ++i){
        d += l1[i] * l2[i];
    }
    return d;
}

/*******************************************************************************
Class Node
    The class manages the node
    can be initialized by giving the weights and bias or
    default initialization or
    random initialization or
    or read directly from a file.
*******************************************************************************/

Node::Node()
{
    //ctor
}

Node::~Node()
{
    //dtor
}

Node::Node(std::vector<double> inw, double inb){
    w = inw;
    b = inb;
}

Node::Node(string filename){
    from_file(filename);
}

/*******************************************************************************
 Node(n_weights, rng, mu, sigma)
   randomly initialize the node
   it initializes with n_weights
   needs a random number generator
   by default is a normal distribution with mean mu = 0.0 and sigma = 1.0
*******************************************************************************/
Node::Node(size_t n_weights, default_random_engine& rng, double mu, double sigma){
    normal_distribution<double> rand_n(mu, sigma);

    b = rand_n(rng);

    for(size_t iweight = 0; iweight < n_weights; ++iweight){
        double weight = rand_n(rng);
        w.push_back(weight);
    }
}

// calculates the z parameter
double Node::z(vector<double> input){
    return dot(input, w) + b;
}

// transform the z parameter in a 0..1 output
double Node::output(vector<double> input){
    return sigmoid(z(input));
}

// read node from filename
void Node::from_file(string filename){
    ifstream file(filename, ios::in|ios::binary);
    if(file.is_open()){
        streampos co = 0;
        // the core of the reading is delegated to this member function
        read_byte_chunk(file, co);
        file.close();
    }
}

// read the chunk corresponding to a node from a file
void Node::read_byte_chunk(ifstream& file, streampos& co){
    w.clear();

    // co is the cumulative offset, the cumulative offset will be moved
    // at the end of the byte chunk that correspond to the node read
    // streampos init_co = co;

    // set cumulative offset to the beginning of the node chunk
    file.seekg(co);

    // read the bias
    file.read(reinterpret_cast<char*>(&b), sizeof(double));
    co += sizeof(double);
    file.seekg(co);

    // read how many weights
    int n_weights = 0;
    file.read(reinterpret_cast<char*>(&n_weights), sizeof(int));
    co += sizeof(int);

    // read as many doubles from the file as there are nodes
    for(int i = 0; i < n_weights; ++i){
        file.seekg(co);

        double weight;
        file.read(reinterpret_cast<char*>(&weight), sizeof(double));
        w.push_back(weight);
        co += sizeof(double);
    }
}

// saves the file in a binary format
void Node::save_bin(string filename){
    ofstream file(filename, ios::out|ios::binary);
    if(file.is_open()){
        streampos co = 0;
        write_byte_chunk(file, co);
        file.close();
    }
}

// function that writes in the file the byte version of this node
void Node::write_byte_chunk(ofstream& file, streampos& co){
    //set the cumulative offset to were to start writing
    file.seekp(co);

    // write the bias
    file.write(reinterpret_cast<char*>(&b),sizeof(b));
    co += sizeof(b);
    //cout << "Node: written bias " << endl;

    // write the number of nodes
    int n_weights = (int)w.size();
    file.write(reinterpret_cast<char*>(&n_weights), sizeof(n_weights));
    co += sizeof(n_weights);

    // writes the weights
    int i = 0;
    for(; i < n_weights; ++i){
        file.seekp(co);
        file.write(reinterpret_cast<char*>(&w[i]), sizeof(w[i]));
        co += sizeof(w[i]);
    }
}

std::ostream& operator<< (std::ostream& os, const Node& node){
    os << "w: (" << node.w.size() << ") |(";

    uint32_t i = 0;
    for(; i < node.w.size() - 1; ++i){
        os << node.w[i] << ", ";
    }
    os << node.w[i];

    os << ") |b: " << node.b;

    return os;

}

