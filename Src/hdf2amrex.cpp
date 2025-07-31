#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include "H5Cpp.h"
#include <hdf2amrex.H>

using namespace amrex;
using namespace H5;

void
read_hdf5(std::string infile, Vector<std::string> const& vars, MultiFab& data, Geometry& geom,const int nGrow)
{
  ParmParse pp;
  std::string datagroupname,coordgroupname;
  datagroupname = "fields"; pp.query("datagroup",datagroupname);
  coordgroupname = "coords"; pp.query("coordgroup",coordgroupname);

  int nvars = vars.size();
  Print() << "opening file " << infile << std::endl;
  H5File file(infile.c_str(), H5F_ACC_RDONLY);
  Group datagroup = file.openGroup(datagroupname);
  Group coords = file.openGroup(coordgroupname);
  DataSet xdat = coords.openDataSet("xm");
  DataSet ydat = coords.openDataSet("ym");
  DataSet zdat = coords.openDataSet("zm");
  DataSpace space = xdat.getSpace();
  hsize_t xdim[1],ydim[1],zdim[1];
  int ndims = space.getSimpleExtentDims(xdim,NULL);
  space = ydat.getSpace();
  ndims = space.getSimpleExtentDims(ydim,NULL);
  space = zdat.getSpace();
  ndims = space.getSimpleExtentDims(zdim,NULL);
  Vector<int> ncell = {static_cast<int>(xdim[0]),static_cast<int>(ydim[0]),static_cast<int>(zdim[0])};
  Vector<int> per = {0,0,0}; pp.queryarr("per",per);
  //grid group contains x, y, z coords
  Vector<Real> x,y;
  x.resize(ncell[0]);
  y.resize(ncell[1]);
#if AMREX_SPACEDIM == 3
  Vector<Real> z;
  z.resize(ncell[2]);
#endif
  
  DataType xtype = xdat.getDataType();
  DataType ytype = ydat.getDataType();
 
  xdat.read(&x[0],xtype);
  ydat.read(&y[0],ytype);
#if AMREX_SPACEDIM == 3
  DataType ztype = zdat.getDataType();
  zdat.read(&z[0],ztype);
#endif
  
  Print() << "x = " << x[0] << " " << x[ncell[0]-1] << std::endl;
  Print() << "y = " << y[0] << " " << y[ncell[1]-1] << std::endl;
  Print() << "z = " << z[0] << " " << z[ncell[2]-1] << std::endl;
  //fill problo and probhi with top and bottom x,y
  Vector<Real> probLo, probHi;
  probLo.resize(AMREX_SPACEDIM);
  probHi.resize(AMREX_SPACEDIM);
  probLo[0] = x[0];
  probLo[1] = y[0];
  probHi[0] = x[ncell[0]-1];
  probHi[1] = y[ncell[1]-1];
#if AMREX_SPACEDIM == 3
  //fortran to c++ switches
  probLo[2] = z[0];
  probHi[2] = z[ncell[2]-1];
#endif
  
  //create box with nx,ny,nz
  IntVect dlo = IntVect::Zero;
  IntVect dhi(ncell);
  dhi[0] -= 1;
  dhi[1] -= 1;
  dhi[2] -= 1;
  Box bx(dlo,dhi);
  Print() << "We have made a box which looks like: " << bx << std::endl;
  Print() << "Box lo: " << bx.smallEnd() << std::endl;
  Print() << "Box hi: " << bx.bigEnd() << std::endl;
  //create box array with the box
  BoxArray ba(bx);
  //will divide up box when parallel
  //int nprocs = ParallelDescriptor::NProcs();
  //int totalcells = AMREX_D_TERM(ncell[0],* ncell[1],* ncell[2]);
  //Real approxsize = std::pow((Real)totalcells/(Real)nprocs,1.0/(Real)AMREX_SPACEDIM);
  //Real approxpower = std::log2(approxsize);
  //int power = (int)(approxpower);
  //int size = std::pow(2,power);

  int size = 32;
  pp.query("max_grid_size",size);
  Print() << "Max box size shrank to " << size << std::endl;
  ba.maxSize(size);
  
  //create realbox
  RealBox rb(&probLo[0],&probHi[0]);
  
  Print() << "We have a RealBox: " << rb << std::endl;
  //coordinate system (0 => cartesian)
  int coord = 0;    
  //create geometry object, needed to write the file
  geom = Geometry(bx,&rb, coord, &(per[0]));
  
  //create appropriate mf
  data.define(ba, DistributionMapping(ba),nvars,nGrow);
  Print() << "About to start processing with " << ParallelDescriptor::NProcs() << " processors" << std::endl;
  for (int n = 0; n < nvars; n++) {
    Print() << "Reading field " << vars[n] << std::endl;
    DataSet dataset = datagroup.openDataSet(vars[n]);
    DataSpace dataspace = dataset.getSpace();
    DataType datatype = dataset.getDataType();
    for (MFIter mfi(data); mfi.isValid(); ++mfi)
      {
	const Box& myBox = mfi.validbox();
	//FArrayBox& myFab = data[mfi];
	const int *lo = myBox.loVect();
	const int *hi = myBox.hiVect();
	hsize_t count[3] = {1,static_cast<hsize_t>(hi[1]-lo[1]+1),static_cast<hsize_t>(hi[0]-lo[0]+1)};	
	Array4<Real> const& arr = data.array(mfi);
#if AMREX_SPACEDIM == 3
	//fortran to c++ switches
	count[0]  = hi[2]-lo[2]+1;
#endif
	Vector<Real> tmpbox(count[0]*count[1]*count[2]);
	DataSpace memspace = H5Screate_simple(3, count, nullptr);
	hsize_t offset[3] = {0,static_cast<hsize_t>(lo[1]),static_cast<hsize_t>(lo[0])};
#if AMREX_SPACEDIM == 3
	offset[0] = lo[2];
#endif
	
	dataspace.selectHyperslab(H5S_SELECT_SET,count,offset);
	dataset.read(tmpbox.data(),datatype,memspace,dataspace);
	//std::cout << "On processor: " << ParallelDescriptor::MyProc() << ", count = "<< count[0] << " " << count[1] << " " << count[2] << ", offset = " << offset[0] << " " << offset[1] << " " << offset[2] << " " << std::endl; 

	amrex::ParallelFor(myBox,[arr,&tmpbox,n,count,lo] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	  int ii = i - lo[0];
	  int jj = j - lo[1];
	  int kk = k - lo[2];
	  int idx = kk * count[1] * count[2] + jj * count[2] + ii;
	  arr(i,j,k,n) = tmpbox[idx];//[k][j][i];
	});
      }
  }
}
