#include <optimalEstimatorANN.H> //contains and libraries for prescribed architecture

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

torch::Tensor compute_regularisation(const torch::nn::Module& model) {
  torch::Tensor l2_reg = torch::tensor(0.0);
  for (const auto param : model.parameters()) {
    l2_reg += torch::sum(torch::pow(param, 2));
  }
  return l2_reg;
}
/*
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
split_data(torch::Tensor input, torch::Tensor target, amrex::Real train_ratio, amrex::Real val_ratio) {
    auto dataset_size = input.size(0);
    auto train_size = static_cast<int>(dataset_size * train_ratio);
    auto val_size = dataset_size-train_size;

    std::vector<int> indices(dataset_size);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g);

    auto train_indices = torch::tensor(std::vector<int>(indices.begin(), indices.begin() + train_size), torch::kLong);
    auto val_indices = torch::tensor(std::vector<int>(indices.begin() + train_size, indices.end()), torch::kLong);
    
    auto train_input = input.index_select(0, train_indices);
    auto train_target = target.index_select(0, train_indices);
    auto val_input = input.index_select(0, val_indices);
    auto val_target = target.index_select(0, val_indices);

    return std::make_tuple(train_input, train_target, val_input, val_target);
}
*/

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    // Open plotfile header and create an amrData object pointing into it
    std::string plotFileName; pp.get("infile",plotFileName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Set up input field data names, and destination components to load data upon read.
    int nFeatures = pp.countval("features");
    int nTargets = pp.countval("targets");
    int nCompIn = nFeatures + nTargets;
    int nEpochs = 1000; pp.query("nEpochs",nEpochs);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    Vector<std::string> features(nFeatures);
    pp.getarr("features",features);
    Vector<std::string> targets(nTargets);
    pp.getarr("targets",targets);
    for (int n = 0; n < nFeatures; n++) {
      inNames[n] = features[n];
      destFillComps[n] = n;
    }
    for (int n = 0; n < nTargets; n++) {
      inNames[n+nFeatures] = targets[n];
      destFillComps[n+nFeatures] = n+nFeatures;
    }
    
    // Loop over AMR levels in the plotfile, read the data and do work
    int finestLevel = amrData.FinestLevel();
    int minLevel = 0;
    pp.query("finestLevel",finestLevel);
    pp.query("minLevel",minLevel);
    
    std::string path="optimal_estimator"; pp.query("model_path",path);    
    path += ".pt";
    int Nlev = finestLevel + 1;
    const int nGrow = 0;
    int n_layers = pp.countval("neurons");
    Vector<int> neurons(n_layers); pp.getarr("neurons",neurons);
    torch::manual_seed(42); //ensure same seeding used each time

    auto dtype0 = torch::kDouble;

    auto model = std::make_shared<Net>(nFeatures,neurons,nTargets);
    //Net model(nFeatures,n_neurons,nTargets);
    
#ifdef AMREX_USE_CUDA
    torch::Device device0(torch::kCUDA);
    model->to(device0);
    amrex::Print() << "Copying model to GPU." << std::endl;
    // set tensor options
    auto tensoropt = torch::TensorOptions().dtype(dtype0).device(device0);
#else
    auto tensoropt = torch::TensorOptions().dtype(dtype0);//.device(device0);
#endif
    
    Vector<Real> f_max(nFeatures,-1e100);
    Vector<Real> f_min(nFeatures,1e100);
    Vector<Real> t_max(nTargets,-1e100);
    Vector<Real> t_min(nTargets,1e100);
    int batch_size = -1; pp.query("batch_size",batch_size);
    Vector<MultiFab> indata(Nlev);
    
    for (int lev = minLevel; lev < Nlev; lev++) {
      //Get the array of boxes for this level
      BoxArray ba = amrData.boxArray(lev);
      if (batch_size > 0) {
	//IntVect boxshape = {AMREX_D_DECL(box_size,box_size,box_size)};
	//ba.minmaxSize(boxshape,boxshape); //everything comes in in equally sized boxes (box_size^AMREX_SPACEDIM)
	//batch_size += ba.size();
	ba.maxSize(batch_size);
      }
      //Distribution mapping i.e. how are boxes distributed across processors
      const DistributionMapping dm(ba);
      //Vector<torch::Tensor> lev_batches(amrData.boxArray(lev).size());
      
      indata[lev].define(ba,dm,nCompIn,nGrow);
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata[lev],lev,inNames,destFillComps); //magic IO call
      Print() << "Data has been read for level " << lev << std::endl;
      
      for (int nv = 0; nv < nFeatures; nv++) {
	f_max[nv] = std::max(f_max[nv],indata[lev].max(nv));
	f_min[nv] = std::min(f_min[nv],indata[lev].min(nv));
      }
      for (int nv = 0; nv < nTargets; nv++) {
	t_max[nv] = std::max(t_max[nv],indata[lev].max(nv+nFeatures));
	t_min[nv] = std::min(t_min[nv],indata[lev].min(nv+nFeatures));
      }
      if (lev == Nlev - 1) {
	for (int nv = 0; nv < nFeatures; nv++) {
	  Print() << "f_max[" << nv << "] = " << f_max[nv] << std::endl;
	  Print() << "f_min[" << nv << "] = " << f_min[nv] << std::endl;
	}
	for (int nv = 0; nv < nTargets; nv++) {
	  Print() << "t_max[" << nv << "] = " << t_max[nv] << std::endl;
	  Print() << "t_min[" << nv << "] = " << t_min[nv] << std::endl;
	}
      }
    }
    Vector<amrex::Gpu::ManagedVector<Real>> training_batches, target_batches;
    Vector<int> ncell_batch;
    
    for (int lev = minLevel; lev < Nlev; lev++) {
      for (MFIter mfi(indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)  {
	const Box& bx = mfi.tilebox();

	Array4<Real> const& training = indata[lev].array(mfi,0);
	Array4<Real> const& target = indata[lev].array(mfi,nFeatures);
	
	const IntVect bx_lo = bx.smallEnd();
	const IntVect nbox = bx.size();
#if AMREX_SPACEDIM == 2
	int ncell = nbox[0] * nbox[1];
#else
	int ncell = nbox[0] * nbox[1] * nbox[2];
#endif
	
	// create a temporary array to store MultiFab data
	// this is needed to create a tensor from a contiguous block of memory
	amrex::Gpu::ManagedVector<Real> trainingVec(ncell*nFeatures);
	amrex::Gpu::ManagedVector<Real> targetVec(ncell*nTargets);
       
	Real* AMREX_RESTRICT trainingPtr = trainingVec.dataPtr();
	Real* AMREX_RESTRICT targetPtr = targetVec.dataPtr();
	AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
				{
				  int ii = i - bx_lo[0];
				  int jj = j - bx_lo[1];
				  int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
				  int kk = k - bx_lo[2];
				  index += kk*nbox[0]*nbox[1];
#endif
				  
				  for (int n = 0; n < nFeatures; n++) {
				    trainingPtr[index*nFeatures + n] = -1.0 + 2.0*(training(i, j, k, n)-f_min[n])/(f_max[n]-f_min[n]);
				  }
				  for (int n = 0; n < nTargets; n++) {
				    targetPtr[index*nTargets + n] = -1.0 + 2.0*(target(i,j,k,n)-t_min[n])/(t_max[n]-t_min[n]);
				  }
				});
	

	ncell_batch.push_back(ncell);
	training_batches.push_back(trainingVec);
	target_batches.push_back(targetVec);
	
      }
    }
    
    if (training_batches.empty()) {
      amrex::Abort("Some processes have no data!");
    }
    
    int num_batches = training_batches.size();
    /*
    //pick a random number (<num_batches)
    std::random_device rd;         // Obtain a random number from hardware
    std::mt19937 gen(rd());        // Seed the generator
    
    // Create a uniform integer distribution
    std::uniform_int_distribution<> dis(0, num_batches-1);
    
    // Generate a random number within the range
    int r_batches = dis(gen);
    */
    Vector<int> indices(num_batches);
    for (int i = 0; i < num_batches; i++) {
      indices[i] = i;
    }
    // Initialize a random number generator
    std::random_device rd;
    std::mt19937 g(rd());

    // Shuffle the indices
    std::shuffle(indices.begin(), indices.end(), g);


    std::cout << "Number of batches on process " << ParallelDescriptor::MyProc() << " = " << num_batches << std::endl;

    Real split = 0.7;
    int num_train = (int)(split*num_batches);
    int num_val = num_batches-num_train;
    Vector<amrex::Gpu::ManagedVector<Real>> features_t(num_train), target_t(num_train);
    Vector<amrex::Gpu::ManagedVector<Real>> features_e(num_val), target_e(num_val);
    Vector<int> ncell_t(num_train);
    Vector<int> ncell_e(num_val);
    
    
    for (int i = 0; i < num_train; i++) {
      features_t[i] = training_batches[indices[i]];
      target_t[i] = target_batches[indices[i]];
      ncell_t[i] = ncell_batch[indices[i]];
    }
    for (int i = num_train; i < num_batches; i++) {
      features_e[i-num_train] = training_batches[indices[i]];
      target_e[i-num_train] = target_batches[indices[i]];
      ncell_e[i-num_train] = ncell_batch[indices[i]];
    }
    
    Real learning_rate = 1e-3; pp.query("learning_rate",learning_rate);
    
    //at this point we have a bunch of distributed tensors, and a model
    // Create an optimizer
    //auto optimizer = torch::optim::SGD(model->parameters(), learning_rate);
    auto optimizer = torch::optim::Adam(model->parameters(), learning_rate);
    //auto optimizer = torch::optim::LBFGS(model->parameters(), learning_rate);

    // Training loop
    torch::Tensor loss;
    Real alpha = 0.01; pp.query("alpha",alpha);
    Real beta = 1.0-alpha;
    Real prev_epoch_training_loss, prev_epoch_validation_loss;

    
    for (int epoch = 0; epoch < nEpochs; ++epoch) {   
      Real epoch_training_loss = 0.0;
      Real epoch_validation_loss = 0.0;      
      //int counter = 0;
      for (int nb = 0; nb < num_train; nb++) {
	torch::Tensor train_input = torch::from_blob(features_t[nb].dataPtr(),{ncell_t[nb],nFeatures},tensoropt);
	torch::Tensor train_target = torch::from_blob(target_t[nb].dataPtr(),{ncell_t[nb],nTargets},tensoropt);

	//torch::Tensor val_input = torch::from_blob(features_e[nb].dataPtr(),{ncell_e[nb],nFeatures},tensoropt);
	//torch::Tensor val_target = torch::from_blob(target_e[nb].dataPtr(),{ncell_e[nb],nTargets},tensoropt);
	
	//torch::Tensor training_batch = torch::from_blob(training_batches[nb].dataPtr(),{ncell_batch[nb],nFeatures},tensoropt);
	//torch::Tensor target_batch = torch::from_blob(target_batches[nb].dataPtr(),{ncell_batch[nb],nTargets},tensoropt);
	
	//auto [train_input, train_target, val_input, val_target] = split_data(training_batch, target_batch, 0.7, 0.3);
	//if (torch::std(train_target).item<Real>() < 0.001) {
	//  counter++;
	//  continue;
	//}
	//Print() << torch::std(train_target) << std::endl;
							 
	//auto [train_input,train_target] = removeConstantRows(train_input_all,train_target_all);
	
	model->train();

	auto output = model->forward(train_input);
	auto mse_loss = torch::mse_loss(output, train_target);
	auto reg_loss = compute_regularisation(*model);
	
	// Forward pass
	//if (!loss.defined()) {
	loss = beta*mse_loss + alpha * reg_loss;
	//}
	
	// Compute loss
	//auto loss = custom_loss(output, train_target, *model);
	// Backward pass and optimization step
	optimizer.zero_grad();
	/*auto trHinv= compute_sum_inv_d2loss_dw2(loss,*model);
	int n_params = model->parameters().size();
	auto gamma = n_params*(1-trHinv/(2.0*(reg_loss + trHinv)));
	auto alpha = gamma/(2.0*reg_loss);
	auto beta = (nTargets-gamma)/(2.0*mse_loss);
	loss = beta * mse_loss + alpha * reg_loss;*/
	loss.backward();
	optimizer.step();
		
	//if (nb == num_batches - 1) { //just compute on the last batch
	epoch_training_loss += loss.item<Real>();	  
	// Validation phase
	//}
      }
      for (int nb = 0; nb < num_val; nb++) {
	torch::Tensor val_input = torch::from_blob(features_e[nb].dataPtr(),{ncell_e[nb],nFeatures},tensoropt);
	torch::Tensor val_target = torch::from_blob(target_e[nb].dataPtr(),{ncell_e[nb],nTargets},tensoropt);

	model->eval();
	auto val_output = model->forward(val_input);
	auto val_loss = torch::mse_loss(val_output, val_target);
	epoch_validation_loss += val_loss.item<Real>();

      }
      //Print() << num_batches << std::endl;
      //Print() << counter << std::endl;
      epoch_training_loss /= num_train; // split
      epoch_validation_loss /= num_val; // 1- split
      /*
      if (epoch > 0) {
	//we're starting to overfit, bias the regularisation more
	if (prev_epoch_training_loss > epoch_training_loss && prev_epoch_validation_loss < epoch_validation_loss) {
	  Print() << "Overfitting, emphasising regularisation... " << std::endl;
	  beta *= 0.9;
	  alpha = 1.0-beta;
	} else {
	  Print() << "Not overfitting, emphasising fitting... " << std::endl;
	  alpha *= 0.9;
	  beta = 1.0-beta;
	}
      }
      */
      //alpha -= alpha/(Real)nEpochs;
      //beta += beta/(Real)nEpochs;
      prev_epoch_training_loss = epoch_training_loss;
      prev_epoch_validation_loss = epoch_validation_loss;
      
      for (auto& param : model->parameters()) {
	param = param.contiguous();
	ParallelDescriptor::ReduceRealSum((param.grad()).data_ptr<Real>(),param.numel());
      }
      ParallelDescriptor::ReduceRealSum(epoch_training_loss);
      ParallelDescriptor::ReduceRealSum(epoch_validation_loss);
      Print() << "Epoch [" << epoch+1 << "/" << nEpochs << "], Training Loss: " << epoch_training_loss << ", Validation Loss: " << epoch_validation_loss << std::endl;
     
      //std::cout << "Process: " << ParallelDescriptor::MyProc() << ", Epoch [" << epoch << "/" << nEpochs << "], Training Loss: " << epoch_training_loss << ", Validation Loss: " << epoch_validation_loss << std::endl;
      
    }
    
    
    
    if (ParallelDescriptor::IOProcessor()) {
      //save the model
      torch::save(model,path);
      std::string minmax_path="minmax"; pp.query("minmax_path",minmax_path);
      minmax_path += ".bin";
      std::ofstream file(minmax_path,std::ios::binary);      
      file.write((char*)f_min.dataPtr(),sizeof(Real)*nFeatures);
      file.write((char*)f_max.dataPtr(),sizeof(Real)*nFeatures);
      file.write((char*)t_min.dataPtr(),sizeof(Real)*nTargets);
      file.write((char*)t_max.dataPtr(),sizeof(Real)*nTargets);
      file.close();
    }
    
    
  }
  Finalize();
  return 0;
}
