#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

//-----------------------------------------------------------------
// DIM = 3
//-----------------------------------------------------------------

#if AMREX_SPACEDIM==3
// Computes the integral along dir
void integrate1d(int dir, int dir1, int dir2, Vector<Vector<Vector<Real>>>& outdata, Vector<Real>& x, Vector<Real>& y,  PlotFileData& pf, Vector<MultiFab*> indata, Vector<MultiFab*> volFracData, int nVars, int finestLevel, int cComp, Real cMin, Real cMax, int avg) {

  Box probDomain = pf.probDomain(finestLevel);
  int const ldir1 = probDomain.length(dir1);
  int const ldir2 = probDomain.length(dir2);                                          
  int refRatio = 1;

  int out_size = (nVars + 1) * ldir1 * ldir2;
  Gpu::DeviceVector<Real> gpu_outdata(out_size, 0.0);
  Real* d_out = gpu_outdata.data();

  // Loop from highest level to lowest
  for (int lev = finestLevel; lev >= 0; lev--) {
    Real dzLev = pf.cellSize(lev)[dir];      
    if (lev < finestLevel) refRatio *= pf.refRatio(lev);  // Akkumulates the refRatio -> Highest Level has ratio of 1 
    Print() << "Current refRatio "<< refRatio << std::endl;
    Print() << "Integrating level "<< lev << std::endl;
    
    for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {       
      const Box& bx = mfi.tilebox();                             
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      Array4<Real> const& volFracBox = (*volFracData[lev]).array(mfi);

      
      amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

	  if (inbox(i,j,k,nVars) > 1e-8 && (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
	    // d[0] = i;
	    // d[1] = j;
	    // d[2] = k;
      Real volFrac = volFracBox(i,j,k,0);

      // GPU-safe
      for (int rx = 0; rx < refRatio; ++rx)
                    {
                        for (int ry = 0; ry < refRatio; ++ry)
                        {
                            Array<int, 3> idx = {i, j, k};
                            int idx1 = refRatio * idx[dir1] + rx;
                            int idx2 = refRatio * idx[dir2] + ry;

                            if (idx1 < ldir1 && idx2 < ldir2)
                            {
                                int flat = idx1 + idx2 * ldir1;

                                Gpu::Atomic::Add(&d_out[0 + flat], dzLev * volFrac);  // Adding entry dzLev*volFrac to d_out for outdata[0][i][j]
                                for (int n = 1; n <= nVars; ++n)
                                {
                                    Gpu::Atomic::Add(&d_out[n * ldir1 * ldir2 + flat],
                                        dzLev * volFrac * inbox(i,j,k,n - 1));  // Integral for every Var
                                }
                            }
                        }
                      }

	  //   for (int rx = 0; rx < refRatio; rx++) {
	  //     for (int ry = 0; ry < refRatio; ry++) {
		// outdata[0][refRatio*d[dir1]+rx][refRatio*d[dir2]+ry] += dzLev*volFrac; 
		// for (int n = 1; n < nVars+1; n++) {
		//   outdata[n][refRatio*d[dir1]+rx][refRatio*d[dir2]+ry] += dzLev*volFrac*inbox(i,j,k,n-1);
		// } // for nVars-loop
	  //     } //for ry-loop
	  //   } // for rx-loop
	  } // if
  }); //PARALLEL_FOR
	
	} // MFIter
	
    } // lev-loop
  //}


  // Host copy & reduction
  Vector<Real> host_outdata(out_size);
  Gpu::copy(Gpu::deviceToHost, gpu_outdata.begin(), gpu_outdata.end(), host_outdata.begin());

  for (int n = 0; n <= nVars; ++n) {
      ParallelDescriptor::ReduceRealSum(&host_outdata[n * ldir1 * ldir2], ldir1 * ldir2);
  }

  // Repack into nested output
  for (int n = 0; n <= nVars; ++n) {
      for (int j = 0; j < ldir2; ++j) {
          for (int i = 0; i < ldir1; ++i) {
              int flat = i + j * ldir1;
              outdata[n][i][j] = host_outdata[n * ldir1 * ldir2 + flat];
          }
      }
  }

  // for (int i = 0; i < ldir1; i++) {
  //   for (int n = 0; n<nVars+1; n++) {
  //     ParallelDescriptor::ReduceRealSum(outdata[n][i].data(),ldir2);
  //   }
  // }
  if (avg) {
    for (int n = 1; n<nVars+1; n++) {
      for (int i = 0; i < ldir1; i++) {
	for (int j = 0; j < ldir2; j++) {
	  if (outdata[0][i][j] > 0.0) outdata[n][i][j] /= outdata[0][i][j];
	}
      }
    }
  }
  Array< Real, AMREX_SPACEDIM > plo = pf.probLo();
  Array< Real, AMREX_SPACEDIM > phi = pf.probHi();

  // Real dxFine = amrData.DxLevel()[finestLevel][dir1];
  // Real dyFine = amrData.DxLevel()[finestLevel][dir2];
  Real dxFine = pf.cellSize(finestLevel)[dir1];
  Real dyFine = pf.cellSize(finestLevel)[dir2];
  for (int i = 0; i < ldir1; i++) {
    x[i] = plo[dir1] + (i+0.5)*dxFine;
  }
  for (int i = 0; i < ldir2; i++) {
    y[i] = plo[dir2] + (i+0.5)*dyFine;
  }  
  return;
} // integrate1d

// Computes the integral for each position in dir 
void integrate2d(int dir, int dir1, int dir2, Vector<Vector<Real>>& outdata, Vector<Real>& x, PlotFileData& pf, Vector<MultiFab*> indata, Vector<MultiFab*> volFracData, int nVars, int finestLevel, int cComp, Real cMin, Real cMax, int avg) {
  Box probDomain = pf.probDomain(finestLevel);
  int ldir = probDomain.length(dir);
  int refRatio = 1;

  int out_size = (nVars + 1) * ldir;
  Gpu::DeviceVector<Real> gpu_outdata(out_size, 0.0);
  Real* d_out = gpu_outdata.data();

  for (int lev = finestLevel; lev >= 0; lev--) {
    Real dxLev = pf.cellSize(lev)[dir1];
    Real dyLev = pf.cellSize(lev)[dir2];
    Real areaLev = dxLev*dyLev;
    if (lev < finestLevel) refRatio *= pf.refRatio(lev);
    Print() << "Integrating level "<< lev << std::endl;

    for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      Array4<Real> const& volFracBox = (*volFracData[lev]).array(mfi);

      amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (inbox(i,j,k,nVars) > 1e-8 &&  (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {

          Real volFrac = volFracBox(i,j,k,0);

          // GPU-safe
          for (int rx = 0; rx < refRatio; ++rx)
                        {
                            
                          Array<int, 3> idx = {i, j, k};
                          int idx1 = refRatio * idx[dir] + rx;

                          if (idx1 < ldir)
                          {
                              int flat = idx1;

                              Gpu::Atomic::Add(&d_out[0 + flat], areaLev * volFrac);  // Adding entry areaLev*volFrac to d_out for outdata[0][i]
                              for (int n = 1; n <= nVars; ++n)
                              {
                                  Gpu::Atomic::Add(&d_out[n * ldir + flat],
                                      areaLev * volFrac * inbox(i,j,k,n - 1));  // Integral for every Var
                              }
                          }
                      
                        }


      //     d[0] = i;
      //     d[1] = j;
      //     d[2] = k;
      //     Real volFrac = volFracBox(i,j,k,0);
      //     for (int rx = 0; rx < refRatio; rx++) {
      //       outdata[0][refRatio*d[dir]+rx] += areaLev*volFrac;     
      //       for (int n = 1; n < nVars+1; n++) {
      //   outdata[n][refRatio*d[dir]+rx] += areaLev*volFrac*inbox(i,j,k,n-1);
	    //   }
	    // }
	  }
	}); // ParallelFor
    } //MFIter
  } //lev-loop

  // Host copy & reduction
  Vector<Real> host_outdata(out_size);
  Gpu::copy(Gpu::deviceToHost, gpu_outdata.begin(), gpu_outdata.end(), host_outdata.begin());

  for (int n = 0; n <= nVars; ++n) {
      ParallelDescriptor::ReduceRealSum(&host_outdata[n * ldir], ldir);
  }

  // Repack into nested output
  for (int n = 0; n <= nVars; ++n) {
      for (int i = 0; i < ldir; ++i) {
          int flat = i;
          outdata[n][i] = host_outdata[n * ldir + flat];
      }
      
  }


  // for (int n = 0; n < nVars+1; n++) {
  //   ParallelDescriptor::ReduceRealSum(outdata[n].data(),ldir);
  // }
  Real dzFine = pf.cellSize(finestLevel)[dir];
  if (avg) {
    for (int n = 1; n<nVars+1; n++) {
      for (int i = 0; i < ldir; i++) {
	if (outdata[0][i] > 0.0) outdata[n][i] /=  outdata[0][i]; // Value / overall Area
      }
    }
  }
  Array< Real, AMREX_SPACEDIM > plo = pf.probLo();
  Array< Real, AMREX_SPACEDIM > phi = pf.probHi();
  for (int i = 0; i < ldir; i++) {
    x[i] = plo[dir] + (i+0.5)*dzFine;
  }
  return;
} // integrate2d
  
void integrate3d(Vector<Real>& outdata, PlotFileData& pf, Vector<MultiFab*> indata, Vector<MultiFab*> volFracData, int nVars, int finestLevel, int cComp, Real cMin, Real cMax, int avg) {

  int out_size = (nVars + 1);
  Gpu::DeviceVector<Real> gpu_outdata(out_size, 0.0);
  Real* d_out = gpu_outdata.data();

  for (int lev = 0; lev <= finestLevel; lev++) {
    Real dxLev = pf.cellSize(lev)[0];
    Real dyLev = pf.cellSize(lev)[1];
    Real dzLev = pf.cellSize(lev)[2];
    Real volLev = dxLev*dyLev*dzLev;
    Print() << "Integrating level "<< lev << std::endl;

    for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      Array4<Real> const& volFracBox = (*volFracData[lev]).array(mfi);
      amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

	  if (inbox(i,j,k,nVars) > 1e-8 &&  (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
      Real volFrac = volFracBox(i,j,k,0);

      // GPU-safe
      Gpu::Atomic::Add(&d_out[0], volLev * volFrac);  // Adding entry volLev*volFrac to d_out for outdata[0]
      for (int n = 1; n <= nVars; ++n)
      {
          Gpu::Atomic::Add(&d_out[n],
              volLev * volFrac * inbox(i,j,k,n - 1));  // Integral for every Var
      }
	    // outdata[0] += volLev*volFrac;
	    // for (int n = 1; n < nVars+1; n++) {
	    //   outdata[n] += volLev*volFrac*inbox(i,j,k,n-1); 
	    // }
	  }
	}); //ParallelFor
    } // MFIter
  } // lev-loop


  // Host copy & reduction
  Vector<Real> host_outdata(out_size);
  Gpu::copy(Gpu::deviceToHost, gpu_outdata.begin(), gpu_outdata.end(), host_outdata.begin());

  ParallelDescriptor::ReduceRealSum(host_outdata.data(), nVars + 1);

  // Repack into nested output
  for (int n = 0; n <= nVars; ++n) {
    outdata[n] = host_outdata[n];
  }
  // ParallelDescriptor::ReduceRealSum(outdata.data(),nVars+1);

  if (avg) {
    for (int n = 1; n<nVars+1; n++) {
      if(outdata[0]> 0.0) outdata[n] /= outdata[0]; // Value / overall Volume
    }
  }
  return;
} // integrate3d


//-----------------------------------------------------------------
// DIM = 2
//-----------------------------------------------------------------

#elif AMREX_SPACEDIM==2 

// 1D Integral along dirInt 
void integrate1d(int dirInt, int dir, Vector<Vector<Real>>& outdata, Vector<Real>& x, PlotFileData& pf, Vector<MultiFab*> indata, Vector<MultiFab*> volFracData, int nVars, int finestLevel, int cComp, Real cMin, Real cMax, int avg) {
  Box probDomain = pf.probDomain(finestLevel);
  int ldir = probDomain.length(dir);
  int refRatio = 1;

  int out_size = (nVars + 1) * ldir;
  Gpu::DeviceVector<Real> gpu_outdata(out_size, 0.0);
  Real* d_out = gpu_outdata.data();

  for (int lev = finestLevel; lev >= 0; lev--) {
    Real dxLev = pf.cellSize(lev)[dirInt];
    if (lev < finestLevel) refRatio *= pf.refRatio(lev);
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      Array4<Real> const& volFracBox = (*volFracData[lev]).array(mfi);

      amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	  if (inbox(i,j,k,nVars) > 1e-8 &&  (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
      Real volFrac = volFracBox(i,j,k,0);


      // GPU-safe
      for (int rx = 0; rx < refRatio; ++rx)
                    {
                        
                      Array<int, 3> idx = {i, j, k};
                      int idx1 = refRatio * idx[dir] + rx;

                      if (idx1 < ldir)
                      {
                          int flat = idx1;

                          Gpu::Atomic::Add(&d_out[0 + flat], dxLev * volFrac);  // Adding entry dxLev*volFrac to d_out
                          for (int n = 1; n <= nVars; ++n)
                          {
                              Gpu::Atomic::Add(&d_out[n * ldir + flat],
                                  dxLev * volFrac * inbox(i,j,k,n - 1));  // Integral for every Var
                          }
                      }
                  
                    }


	  //   for (int rx = 0; rx < refRatio; rx++) {
	  //     outdata[0][refRatio*d[dir]+rx] += dxLev*volFrac;
	  //     for (int n = 1; n < nVars+1; n++) {
		// outdata[n][refRatio*d[dir]+rx] += dxLev*volFrac*inbox(i,j,k,n-1);
	  //     }
	  //   }
	  }
	}); // ParallelFor
    } // MFIter
  }


  // Host copy & reduction
  Vector<Real> host_outdata(out_size);
  Gpu::copy(Gpu::deviceToHost, gpu_outdata.begin(), gpu_outdata.end(), host_outdata.begin());

  // for (int n = 0; n <= nVars; ++n) {
  //     ParallelDescriptor::ReduceRealSum(&host_outdata[n * ldir], ldir);
  // }

  ParallelDescriptor::ReduceRealSum(host_outdata.data(), out_size);

  // Repack into nested output
  for (int n = 0; n <= nVars; ++n) {
      for (int i = 0; i < ldir; ++i) {
          int flat = i;
          outdata[n][i] = host_outdata[n * ldir + flat];
      }
      
  }

  // for (int n = 0; n < nVars+1; n++) {
  //   ParallelDescriptor::ReduceRealSum(outdata[n].data(),ldir);
  // }
  Real dyFine = pf.cellSize(finestLevel)[dir];
  if (avg) {
    for (int n = 1; n<nVars+1; n++) {
      for (int i = 0; i < ldir; i++) {
	if (outdata[0][i] > 0.0) outdata[n][i] /=  outdata[0][i];
      }
    }
  }
  Array< Real, AMREX_SPACEDIM > plo= pf.probLo();
  Array< Real, AMREX_SPACEDIM > phi= pf.probHi();
  for (int i = 0; i < ldir; i++) {
    x[i] = plo[dir] + (i+0.5)*dyFine;
  }
  return;
} // integrate1d
  

// 2D Integral
void integrate2d(Vector<Real>& outdata, PlotFileData& pf, Vector<MultiFab*> indata, Vector<MultiFab*> volFracData, int nVars, int finestLevel, int cComp, Real cMin, Real cMax, int avg) {

  int out_size = (nVars + 1);
  Gpu::DeviceVector<Real> gpu_outdata(out_size, 0.0);
  Real* d_out = gpu_outdata.data();

  for (int lev = 0; lev <= finestLevel; lev++) {
    Real dxLev = pf.cellSize(lev)[0];
    Real dyLev = pf.cellSize(lev)[1];
    Real areaLev = dxLev*dyLev;
    Print() << "Integrating level "<< lev << std::endl;
    for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& inbox  = (*indata[lev]).array(mfi);
      Array4<Real> const& volFracBox = (*volFracData[lev]).array(mfi);
      amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {



	  if (inbox(i,j,k,nVars) > 1e-8 &&  (cComp < 0 || (inbox(i,j,k,cComp) >= cMin && inbox(i,j,k,cComp) < cMax))) {
      Real volFrac = volFracBox(i,j,k,0);

      // GPU-safe
      Gpu::Atomic::Add(&d_out[0], areaLev * volFrac);  // Adding entry areaLev*volFrac to d_out for outdata[0]
      for (int n = 1; n <= nVars; ++n)
      {
          Gpu::Atomic::Add(&d_out[n],
              areaLev * volFrac * inbox(i,j,k,n - 1));  // Integral for every Var
      }

	    // outdata[0] += areaLev*volFrac;
	    // for (int n = 1; n < nVars+1; n++) {
	    //   outdata[n] += areaLev*volFrac*inbox(i,j,k,n-1);
	    // }
	  }
	});
    }
  }

  // Host copy & reduction
  Vector<Real> host_outdata(out_size);
  Gpu::copy(Gpu::deviceToHost, gpu_outdata.begin(), gpu_outdata.end(), host_outdata.begin());

  ParallelDescriptor::ReduceRealSum(host_outdata.data(), nVars + 1);

  // Repack into nested output
  for (int n = 0; n <= nVars; ++n) {
    outdata[n] = host_outdata[n];
  }

  // ParallelDescriptor::ReduceRealSum(outdata.data(),nVars+1);
  if (avg) {
    for (int n = 1; n<nVars+1; n++) {
      if(outdata[0]> 0.0) outdata[n] /= outdata[0];
    }
  }
  return;
} // integrate2d
#endif


//-----------------------------------------------------------------
// Write data functions
//-----------------------------------------------------------------

void writeDat1D(Vector<Real> vect, std::string filename, std::string folder, int dim) {
  std::string fullPath = folder + "/" + filename;
  FILE *file = fopen(fullPath.c_str(),"w");
  for (int i = 0; i < dim; i++) {
    fprintf(file,"%e ",vect[i]);
  }
  fclose(file);
  return;
}

void writeDat2D(Vector<Vector<Real>> vect, std::string filename, std::string folder, int dim1, int dim2) {

  std::string fullPath = folder + "/" + filename;

  FILE *file = fopen(fullPath.c_str(),"w");
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      fprintf(file,"%e ",vect[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return;
}

void writePPM(Vector<Vector<Real>> vect, std::string filename, int dim1, int dim2, int goPastMax, Real vMin, Real vMax) {
  unsigned char *buff=(unsigned char*)malloc(3*dim1*dim2*sizeof(char));
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      int bc   = ((dim1-i-1)*dim2+j)*3;
      Real val = vect[i][j];
      Real colour = fmax(0.,fmin(1.5,(val-vMin)/(vMax-vMin)));
      if (colour<0.125) {
	buff[bc]   = 0;
	buff[bc+1] = 0;
	buff[bc+2] = (int)((colour+0.125)*1020.);
      } else if (colour<0.375)  {
	buff[bc]   = 0;
	buff[bc+1] = (int)((colour-0.125)*1020.);
	buff[bc+2] = 255;
      } else if (colour<0.625)  {
	buff[bc]   = (int)((colour-0.375)*1020.);
	buff[bc+1] = 255;
	buff[bc+2] = (int)((0.625-colour)*1020.);
      } else if (colour<0.875)  {
	buff[bc]   = 255;
	buff[bc+1] = (int)((0.875-colour)*1020.);
	buff[bc+2] = 0;
      } else if (colour<1.000)  {
	buff[bc]   = (int)((1.125-colour)*1020.);
	buff[bc+1] = 0;
	buff[bc+2] = 0;
      } else if (goPastMax==1) {
	if (colour<1.125)  {
	  buff[bc]   = (int)((colour-0.875)*1020.);
	  buff[bc+1] = 0;
	  buff[bc+2] = (int)((colour-1.000)*1020.);
	} else if (colour<1.250) {
	  buff[bc]   = 255;
	  buff[bc+1] = 0;
	  buff[bc+2] = (int)((colour-1.000)*1020.);
	} else if (colour<1.500)  {
	  buff[bc]   = 255;
	  buff[bc+1] = (int)((colour-1.250)*1020.);
	  buff[bc+2] = 255;
	} else { // default if above 1.5 with goPastMax==1
	  buff[bc]   = 255;
	  buff[bc+1] = 255;
	  buff[bc+2] = 255;
	}
      } else { // default if above 1 with goPastMax==0
	buff[bc]   = 128;
	buff[bc+1] = 0;
	buff[bc+2] = 0;
      }
    }
  }
  FILE *file = fopen(filename.c_str(),"w");
  fprintf(file,"P6\n%i %i\n255\n",dim2,dim1);
  fwrite(buff,dim1*dim2*3,sizeof(unsigned char),file);
  fclose(file);
  return;
}

void findMinMax(Vector<Vector<Real>> vect, int dim1, int dim2, Real &min, Real &max) {
  min = vect[0][0];
  max = vect[0][0];
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      if (vect[i][j] < min) min = vect[i][j];
      if (vect[i][j] > max) max = vect[i][j];
    }
  }
  return;
}





//-----------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------

int main(int argc, char *argv[])
{
  amrex::Initialize(argc, argv);
  {
  ParmParse pp;
  
  // Init of input infile
  std::string infile;
  pp.get("infile",infile);
  Print() << "infile = " << infile << std::endl; 

  // Init of path for output
  std::string path;
  pp.get("path",path);
  Print() << "output path = " << path << std::endl; 

  PlotFileData pf(infile);


  // read in the variable names to use
  int const nVars= pp.countval("vars");    
  Vector<std::string> vars(nVars);
  pp.getarr("vars", vars);
  Vector<int> destFillComps(nVars);
  Print() << "nVars= " << nVars << std::endl;
  for (int n = 0; n<nVars; n++) {
    destFillComps[n] = n;
    Print() << "vars[" << n << "]= " << vars[n] << std::endl;
  }
  

  // Init of integralDimension
  int integralDimension;
  pp.get("integralDimension",integralDimension);
  AMREX_ALWAYS_ASSERT(integralDimension<=AMREX_SPACEDIM);

  // Init of finest Level 
  int finestLevel = pf.finestLevel();
  pp.query("finestLevel", finestLevel);
  int Nlev = finestLevel + 1;

  // Init for optional conditons
  std::string cVar;
  Real cMin,cMax;
  int cComp=-1;
  pp.query("cVar",cVar);
  pp.query("cMin",cMin);
  pp.query("cMax",cMax);
  // Sets cComp to number of choosen variable, also checks if cVar is in listed vars
  if (!cVar.empty()) {             
    for (int n = 0; n<nVars; n++) {
      if (vars[n] == cVar) {
	cComp = n;
	break;
      }
    }
    if (cComp < 0) {
      Abort("cVar not in list of vars!");
    }
  }

  
  int avg = 0;
  pp.query("avg",avg);


  // Init of integral directions + output-file string
  int dir, dir1, dir2;
  std::string format="dat";
  
  Print() << "integralDimension = " << integralDimension << std::endl;
#if AMREX_SPACEDIM==3
  switch(integralDimension) {
  case 1:
    {
      pp.get("dir",dir);
      dir1 = (dir+1)%3;     
      dir2 = (dir+2)%3;
      pp.query("format",format);
      AMREX_ALWAYS_ASSERT(format=="ppm" || format=="dat");
      break;
    }
  case 2:
    {
      pp.get("dir1",dir1);   
      pp.get("dir2",dir2);
      dir = 3-dir1-dir2;    
      break;
    }
    //case 3 doesn't care about directions
  }
#elif AMREX_SPACEDIM==2
  if (integralDimension == 1) {
    pp.get("dir",dir);
    dir1 = (dir+1)%2;
  }
#endif

// Extract timestep number from string
  int start = infile.find("plt");
  std::string timestep = infile.substr(start, 8);
  Print() << "Timestep = " << timestep << std::endl;

  std::string outfile= timestep +"_integral";
  if(integralDimension < AMREX_SPACEDIM) {
    outfile += "_dir"+std::to_string(dir);
  }
  if(!cVar.empty()) {
    outfile+="_c"+cVar+"_"+std::to_string(cMin)+"_"+std::to_string(cMax);
  }
  if(avg) {
    outfile+= "_avg";
  }
    

//-----------------------------------------------------------------
// IO
//-----------------------------------------------------------------


  Vector<MultiFab*> indata(Nlev);   
  Vector<MultiFab*> volFracData(Nlev);
  for (int lev = 0; lev < Nlev; lev++) {
    BoxArray probBoxArray = pf.boxArray(lev);
    DistributionMapping dmap = pf.DistributionMap(lev);
    indata[lev] = new MultiFab(probBoxArray,dmap,nVars+1,0);
    Print() << "Loading data on level " << lev << std::endl;
    for (int n=0; n< nVars; ++n) {
        const MultiFab& src = pf.get(lev, vars[n]);
        MultiFab::Copy(*indata[lev], src, 0, n, 1, 0);
      }
    Print() << "Data loaded" << std::endl;
    indata[lev]->setVal(1.0,nVars,1);  

  #if AMREX_USE_EB
    volFracData[lev] = new MultiFab(probBoxArray, dmap, 1, 0);
    Print() << "Loading volFrac data on level " << lev << std::endl;
    MultiFab::Copy(*volFracData[lev], pf.get(lev, "volFrac"), 0, 0, 1, 0);
      
    Print() << "volFrac data loaded" << std::endl;
  #else 
    volFracData[lev] = new MultiFab(probBoxArray, dmap, 1, 0);
    volFracData[lev]->setVal(1.0);
  #endif
  }


//-----------------------------------------------------------------
// Cell Intersections
//-----------------------------------------------------------------


// Finding the intersections of cells
// Forced to run on CPU 
#ifndef AMREX_USE_GPU
  Print() << "Determining intersects..." << std::endl;
  for (int lev = 0; lev < finestLevel; lev++) {
    BoxArray baf = (*indata[lev+1]).boxArray();       
    baf.coarsen(pf.refRatio(lev));	            
    for (MFIter mfi(*indata[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      FArrayBox& myFab = (*indata[lev])[mfi];
      int idx = mfi.index();                                 
      std::vector< std::pair<int,Box> > isects = baf.intersections((*indata[lev]).boxArray()[idx]);
      for (int ii = 0; ii < isects.size(); ii++) {
	myFab.setVal(0.0,isects[ii].second,nVars,1);  
      }
    }
  }
  Print() << "Intersects determined" << std::endl;
#endif


//-----------------------------------------------------------------
// Real work
//-----------------------------------------------------------------

#if AMREX_SPACEDIM==3  
  switch(integralDimension) {
  case 1: //1D integral, results in 2D data (assuming 3D plotfile), output either dat or ppm
    {
      Box probDomain = pf.probDomain(finestLevel);
      int ldir1 = probDomain.length(dir1);
      int ldir2 = probDomain.length(dir2);
      Vector<Real> x(ldir1);
      Vector<Real> y(ldir2);
      Vector<Real> tmp1(ldir2,0.0);       
      Vector<Vector<Real>> tmp2(ldir1,tmp1);  
      Vector<Vector<Vector<Real>>> outdata(nVars+1,tmp2); 
      //do 1d integration
      integrate1d(dir,dir1,dir2,outdata,x,y,pf,indata,volFracData,nVars,finestLevel,cComp,cMin,cMax,avg); 
      Print() << "Integration completed" << std::endl;
      //output data in desired format
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(x,outfile+"_x.dat",path,ldir1);
	  writeDat1D(y,outfile+"_y.dat",path,ldir2);
	  writeDat2D(outdata[0],outfile+"_length.dat",path,ldir1,ldir2);
	  for (int n = 1; n < nVars+1; n++) {
	    writeDat2D(outdata[n],outfile+"_"+vars[n-1]+".dat",path,ldir1,ldir2);
	  }
	} else if (format == "ppm") {
	  int goPastMax = 1;
	  pp.query("goPastMax",goPastMax);
	  Vector<Real> vMin(nVars+1);
	  Vector<Real> vMax(nVars+1);
	  findMinMax(outdata[0],ldir1,ldir2,vMin[0],vMax[0]);
	  for (int n=1; n<nVars+1; n++) {
	    char argName[12];
	    sprintf(argName,"useminmax%i",n);
	    int nMinMax = pp.countval(argName);
	    if (nMinMax > 0) {
	      Print() << "Reading min/max from command line" << std::endl;
	      if (nMinMax != 2) {
		Abort("Need to specify 2 values for useMinMax");
	      } else {
		pp.get(argName, vMin[n], 0);
		pp.get(argName, vMax[n], 1);
	      }
	    } else {
	      Print() << "Using file values for min/max" << std::endl;
	      findMinMax(outdata[n],ldir1,ldir2,vMin[n],vMax[n]);
	    }
	  }
	  
	  writePPM(outdata[0],outfile+"_length.ppm",ldir1,ldir2,goPastMax,vMin[0],vMax[0]);
	  for (int n = 1; n < nVars+1; n++) { 
	    writePPM(outdata[n],outfile+"_"+vars[n-1]+".ppm",ldir1,ldir2,goPastMax,vMin[n],vMax[n]);
	  }
	} //can add more formats here if we want - add to assert above
      }
      break;
    }
  case 2: //2D integral
    {
      Box probDomain = pf.probDomain(finestLevel);
      int ldir = probDomain.length(dir);
      Vector<Real> x(ldir);
      Vector<Real> tmp(ldir,0.0);
      Vector<Vector<Real>> outdata(nVars+1,tmp);
      integrate2d(dir,dir1,dir2,outdata,x,pf,indata,volFracData,nVars,finestLevel,cComp,cMin,cMax,avg);
      Print() << "Integration completed" << std::endl;
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(x,outfile+"_x.dat",path,ldir);
    std::string stringAttachment = "_allVars.dat";
    pp.query("stringAttachment",stringAttachment);
	  writeDat2D(outdata,outfile + stringAttachment + ".dat",path,nVars+1,ldir);
	} //can add more formats here if we want
      }
      break;
    }
  case 3: //3D integral
    {
      format="dat"; //probably add an option for binary output
      Vector<Real> outdata(nVars+1,0.0);
      integrate3d(outdata,pf,indata,volFracData,nVars,finestLevel,cComp,cMin,cMax,avg);
      Print() << "Integration completed" << std::endl;
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
    std::string stringAttachment = "_allVars.dat";
    pp.query("stringAttachment",stringAttachment);
	  writeDat1D(outdata,outfile+stringAttachment+".dat",path,nVars+1);
	} //can add more formats here
      }
      break;
    }
  }

// For 2D Domain 
#elif AMREX_SPACEDIM==2
  switch(integralDimension) {
  case 1: // 1D Integral 
    {
      Box probDomain = pf.probDomain(finestLevel);
      int ldir = probDomain.length(dir1);
      Vector<Real> x(ldir);
      Vector<Real> tmp(ldir,0.0);
      Vector<Vector<Real>> outdata(nVars+1,tmp);
      integrate1d(dir,dir1,outdata,x,pf,indata,volFracData,nVars,finestLevel,cComp,cMin,cMax,avg);
      Print() << "Integration completed" << std::endl;
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(x,outfile+"_x.dat",path,ldir);
	  writeDat2D(outdata,outfile+"_allVars.dat",path,nVars+1,ldir);
	} //can add more formats here if we want
      }
      break;
    }
  case 2: // 2D Integral 
    {
      format="dat"; //probably add an option for binary output
      Vector<Real> outdata(nVars+1,0.0);
      integrate2d(outdata,pf,indata,volFracData,nVars,finestLevel,cComp,cMin,cMax,avg);
      Print() << "Integration completed" << std::endl;
      Print() << "Writing data as "+format << std::endl;
      if (ParallelDescriptor::IOProcessor()) {
	if (format == "dat") {
	  writeDat1D(outdata,outfile+"_allVars.dat",path,nVars+1);
	} //can add more formats here
      }
      break;
    }
  }
#endif
  }

  Finalize();
  return 0;  
}
















