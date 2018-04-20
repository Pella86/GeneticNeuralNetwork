#include "Network.h"

using namespace std;

Network::Network()
{
    //ctor
}

Network::~Network()
{
    //dtor
}


Network::Network(string filename){
    streampos co = 0; // cumulative offset

    ifstream file(filename, ios::in|ios::binary);
    if(file.is_open()){
        cout << "Reading network..." << endl;
        // read the first int -> lengths of layers
        int n_layers;
        file.read(reinterpret_cast<char*>(&n_layers), sizeof(int));
        co += sizeof(int);

        // read how many nodes per layers are there and stores them in
        // nodes layer
        vector<int> nodes_layer;
        for(int i = 0; i < n_layers; ++i){
            file.seekg(co);
            int n_nodes;
            file.read(reinterpret_cast<char*>(&n_nodes), sizeof(int));
            nodes_layer.push_back(n_nodes);
            co += sizeof(int);
        }
        cout << endl;

        // for each layer
        for(int i = 0; i < n_layers; i += 1){

            // create a placeholder for the node vector
            vector<Node> nv;
            layers.push_back(nv);

            // for each node
            for(int inode = 0; inode < nodes_layer[i]; ++inode){
                Node n;
                n.read_byte_chunk(file, co);
                layers[i].push_back(n);
            }
        }
    }
    cout << "total bytes read: " << co << endl;
    file.close();
}

Network::Network(vector<int> node_layers, default_random_engine& rng){
    // create the layers
    for(size_t ilayer = 1; ilayer < node_layers.size(); ++ilayer){
        vector<Node> nv;
        layers.push_back(nv);

        // append each layers n nodes
        int n_nodes = node_layers[ilayer];
        for(int inode = 0; inode < n_nodes; ++inode){
            int n_nodes = node_layers[ilayer - 1];
            Node n(n_nodes, rng, 0, 4);
            layers[ilayer - 1].push_back(n);
        }
    }
}


void Network::save_to_file(string filename){
    ofstream file(filename, ios::out | ios::binary);
    if(file.is_open()){
        streampos co = 0; // cumulative offset set to the beginning of file

        // write number of layers
        int n_layers = layers.size();
        cout << "written number layers:" << n_layers << endl;
        file.write(reinterpret_cast<char*> (&n_layers), sizeof(n_layers));
        co += sizeof(n_layers);

        // write layers per node
        cout << "writing layer per node..." << endl;
        for(vector<Node> layer : layers) {
            int n_nodes = layer.size();
            file.seekp(co);
            file.write(reinterpret_cast<char*> (&n_nodes), sizeof(n_nodes));
            co += sizeof(n_nodes);
        }

        //write the nodes
        int node_counter = 0;
        for(auto l : layers){
            for(auto n : l){
                n.write_byte_chunk(file, co);
                ++node_counter;
            }
        }
        cout << "written " << node_counter << " nodes" << endl;
        cout << "total bytes written: " << co <<endl;
    }
    file.close();
}


vector<double> Network::calculate(vector<double> input){
    //the input is copied in a dynamic input array that will be reused in
    //the calculation cycle
    vector<double> dyn_output;
    vector<double> dyn_input = input;
    for (auto l : layers){
        // reset the output vector
        dyn_output.clear();

        // populate output vector with results from layers
        for(auto n : l){
            dyn_output.push_back(n.output(dyn_input));
        }

        // copy the dynamic input so that it can be used again
        dyn_input.clear();
        dyn_input = dyn_output;
    }

    return dyn_output;
}

std::ostream& operator<<(std::ostream& os, const Network& nn){
    os << "--- Network ---" << endl;
    os << "n layers:" << nn.layers.size() << endl;
    os << "nodes per layer: ";
    for(auto i : nn.layers) os << i.size() << " ";
    os << endl;

    for(auto l : nn.layers){
        os << "------" << endl;
        for(auto n : l){
            os << n << endl;
        }
    }
    return os;

}
