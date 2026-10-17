#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<name> [options] \n\tOptions:\n";
    std::cerr << "\t     outfile=<name>\n";
    exit(1);
}

static std::string parseTitle(std::istream& is);
static std::vector<std::string> parseVarNames(std::istream& is);
static std::string rootName(const std::string in);

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool verbose = false; pp.query("verbose",verbose);
    std::string infile;
    pp.get("infile",infile);
    std::string outfile = rootName(infile)+"_surfAvg.dat";
    pp.query("outfile",outfile);
    std::string outfile_x = rootName(infile)+"_surfAvg_x.dat";
    int dir;
    pp.get("dir",dir);
    int nbins = 100;
    pp.query("nbins",nbins);
    std::ifstream ifs;
    std::istream* is = (infile=="-" ? (std::istream*)(&std::cin) : (std::istream*)(&ifs) );
    if (infile!="-")
      {
	if (verbose)
	  std::cerr << "Opening " << infile << std::endl;
	ifs.open(infile.c_str(),std::ios::in|std::ios::binary);
	if (ifs.fail())
	  std::cerr << "Unable to open file : " << infile << std::endl;
      }
    else
      {
	if (verbose)
	  std::cerr << "Reading from stream" << std::endl;
      }
    const std::string title = parseTitle(*is);
    const std::vector<std::string> names = parseVarNames(*is);
    const int nComp = names.size();
    Vector<Real> tmp(nComp-AMREX_SPACEDIM+1,0.0);
    Vector<Vector<Real>> dataArr(nbins,tmp);
    Vector<Real> dirArr(nbins);
    Real xmin, xmax;
    pp.get("xmin",xmin);
    pp.get("xmax",xmax);
    for (int ibin =0; ibin < nbins;ibin++) {
      dirArr[ibin] = xmin + ibin*(xmax-xmin)/(Real)nbins;
    }
    int nElts;
    int MYLEN;
    (*is) >> nElts;
    (*is) >> MYLEN;
    
    AMREX_ALWAYS_ASSERT(MYLEN == AMREX_SPACEDIM);
    
    FArrayBox nodeFab;
    nodeFab.readFrom((*is));
    Real* nodeData = nodeFab.dataPtr();

    Vector<int> connData(nElts*MYLEN,0);
    (*is).read((char*)connData.dataPtr(),sizeof(int)*connData.size());
    constexpr Real dimWeight = 1.0/static_cast<Real>(AMREX_SPACEDIM);
    
    for (int i=0; i<nElts; ++i)
      {
	//get element offset
	const int offsetElt = i*MYLEN;
	//offset for each node on that element
	Array<int,AMREX_SPACEDIM> offsetData;
	//coordinate x[node][dim] 
	Array<Array<Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> x;
	//which location bin is the node in?
	Array<int,AMREX_SPACEDIM> idx;	
	for (int inode = 0; inode < AMREX_SPACEDIM; inode++) {
	  offsetData[inode] = connData[offsetElt+inode]-1;
	  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
	    x[inode][idim] = nodeData[offsetData[inode]*nComp+idim];
	  }
	  idx[inode] = (int)(nbins * (x[inode][dir]-xmin)/(xmax-xmin));
	}
	//weight = length or area depending on dimension
#if AMREX_SPACEDIM == 2
	const Real weight = std::hypot(x[1][0]-x[0][0],x[1][1]-x[0][1]);
#else
	Array<Real,AMREX_SPACEDIM> edge;
	Real s = 0.0;
	for (int iedge = 0; iedge < AMREX_SPACEDIM; iedge++) {
	  edge[iedge] = std::hypot(x[(iedge+1)%AMREX_SPACEDIM][0]-x[iedge][0],x[(iedge+1)%AMREX_SPACEDIM][1]-x[iedge][1],x[(iedge+1)%AMREX_SPACEDIM][2]-x[iedge][2]);
	  s += 0.5*edge[iedge];
	}
	//Heron's formula for triangle area
	const Real weight = std::sqrt(s*(s-edge[0])*(s-edge[1])*(s-edge[2]));
#endif
	for (int inode = 0; inode < AMREX_SPACEDIM; inode++) {
	  dataArr[idx[inode]][0] += dimWeight * weight;
	}
	for (int icomp = 0; icomp < nComp-AMREX_SPACEDIM;icomp++) {
	  Real locAvg = 0.0;
	  for (int inode = 0; inode < AMREX_SPACEDIM; inode++) {
	    locAvg += dimWeight * nodeData[offsetData[inode]*nComp + AMREX_SPACEDIM + icomp];
	  }
	  for (int inode = 0; inode < AMREX_SPACEDIM; inode++) {	    	  
	    dataArr[idx[inode]][icomp+1] += dimWeight * weight * locAvg;
	  }
	}
      }
    //Write ASCII file
    std::ofstream os(outfile.c_str(),std::ios::out);

    for (int ibin = 0; ibin < nbins; ibin++) {
      //x bins 
      os << dirArr[ibin] << " ";
      //area bins
      os << dataArr[ibin][0] << " "; 
      for (int n = 1; n < nComp-AMREX_SPACEDIM+1; n++) {
	//variable bins (divide by area for avg)
	if (dataArr[ibin][0] > 0) {
	  dataArr[ibin][n] /= dataArr[ibin][0];
	}
	os << dataArr[ibin][n] << " ";
      }
      os << "\n";
    }
    os.close();

}
amrex::Finalize();
return 0;
}

static
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return amrex::Tokenize(line,std::string(" "));
}

static std::string parseTitle(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return line;
}

static std::string rootName(const std::string inStr)
{
#ifdef WIN32
    const std::string dirSep("\\");
#else
    const std::string dirSep("/");
#endif
    std::vector<std::string> res = amrex::Tokenize(inStr,dirSep);
    res = amrex::Tokenize(res[res.size()-1],std::string("."));
    std::string out = res[0];
    for (int i=1; i<res.size()-1; ++i)
        out = out + std::string(".") + res[i];
    return out;
    //return (res.size()>1 ? res[res.size()-2] : res[res.size()-1]);
}

