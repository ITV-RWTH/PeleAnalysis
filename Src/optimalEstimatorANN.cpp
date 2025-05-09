#include <optimalEstimatorANN.H>

Net::Net(int input_size, amrex::Vector<int> neurons, int output_size) {
  n_layers = neurons.size();
  layers.resize(n_layers+1,nullptr);
  if (n_layers == 0) {
    layers[0] = register_module("fc",torch::nn::Linear(input_size,output_size));
    layers[0]->to(torch::kDouble);
  } else {
    layers[0] = register_module("fc0", torch::nn::Linear(input_size,neurons[0]));
    layers[0]->to(torch::kDouble);
    layers[n_layers] = register_module("fc"+std::to_string(n_layers),torch::nn::Linear(neurons[n_layers-1],output_size));
    layers[n_layers]->to(torch::kDouble);
    for (int n = 1; n < n_layers; n++) {
      layers[n] = register_module("fc"+std::to_string(n),torch::nn::Linear(neurons[n-1],neurons[n]));
      layers[n]->to(torch::kDouble);
    }
  }
  initializeWeights();
}

// Implement the forward pass
torch::Tensor Net::forward(torch::Tensor x) {
  //following Berger et al. (2018), linear activation function, followed by tansig function
  /*
  x = layers[0]->forward(x);
  for (int n = 1; n <= n_layers; n++) {
    x = tansig(layers[n]->forward(x));
  }
  */
  
  for (int n = 0; n < n_layers; n++) {
    x = tansig(layers[n]->forward(x));
  }
  x = layers[n_layers]->forward(x);
  
  return x;
}

torch::Tensor Net::tansig(torch::Tensor x) {
    return 2.0 / (1.0 + torch::exp(-2.0 * x)) - 1.0;
}

  
  // Function to initialize weights
void Net::initializeWeights() {
  // Calculate the scaling factor
  for (int n = 0; n <= n_layers; n++) {
    torch::nn::init::xavier_uniform_(layers[n]->weight);
    /*
    auto weight = layers[n]->weight;
    auto bias = layers[n]->bias;
    int output_neurons = weight.size(0);
    amrex::Print() << n << std::endl;
    int input_neurons = weight.size(1);
    //Xavier initialisation
    amrex::Real beta = std::sqrt(6.0)/std::sqrt(output_neurons+input_neurons);
    //amrex::Real beta = 0.7*std::pow(output_neurons,1.0/input_neurons);//std::sqrt(output_neurons)/std::sqrt(input_neurons);
    amrex::Print() << beta << std::endl;
    torch::nn::init::normal_(weight,0,beta);
    torch::nn::init::zeros_(bias);
    */
  }
}
    
