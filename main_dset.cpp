#include "dsetxi_netsites.hpp"
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <omp.h>
using std::cout;
int main(int argc, char **argv) {
    srand(static_cast<unsigned>(getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    if (argc < 2){
      std::cerr<<"Requires three args: L trials num_threads";
      exit(1);
    }
    
    unsigned L = atoi(argv[1]);
    int trials=atoi(argv[2]);
    omp_set_num_threads(atoi(argv[3]));
    
    double to_p=0.995;
    double to_k=4;
    double kavg=4;
    double lambda= 0.1;
    int measure_every = L/10;
    std::ofstream pxi_outfile,ppinf_outfile;
    std::vector<double> lambda_vec{  0.02, 0.1, 0.2, 0.5, 1, 2, 5};
    bool firstlambda=true;
    std::stringstream pxifname,ppinffname;
    int timestamp = time(0)%100000;
    pxifname <<"pxi_L"<<L<<"_multiplelambda."<<timestamp<<".json";
    ppinffname <<"ppinf_L"<<L<<"_multiplelambda."<<timestamp<<".json";
    pxi_outfile.open(pxifname.str());
    ppinf_outfile.open(ppinffname.str());
    pxi_outfile <<"{";
    ppinf_outfile <<"{";
    
for (auto& lambda : lambda_vec){    
    if (!firstlambda){
        pxi_outfile <<",\n";
        ppinf_outfile <<",\n";
    }
    firstlambda=false;
    pxi_outfile << "\""<<lambda<<"\" : [";
    ppinf_outfile << "\""<<lambda<<"\" : [";    
    bool firsttrial=true;
#pragma omp parallel shared(firsttrial)
    {
#pragma omp for
	for(int i=0; i<trials; i++){
	    try{
	      DSetXiNetSites dsX(L,kavg,lambda);
	      dsX.incrementalComponents(to_p,measure_every);

	    #pragma omp critical 
		{
		  if(false){
		    std::ofstream cmap_file,cjson_file;
		    cmap_file.open("/tmp/clustermap.json");
		    cjson_file.open("/tmp/clusterdata.json");
		    dsX.print_cluster_json(cjson_file, 1);
		    dsX.print_cluster_map(cmap_file);
		    cmap_file.close();
		    cjson_file.close();
		    std::cerr << "Found a bad one!\n";  
		  }
		  
	    #pragma omp flush(firsttrial)

		if(!firsttrial){
		    pxi_outfile <<",\n";
		    ppinf_outfile <<",\n";
		}
		firsttrial=false;
		dsX.print_k_xi(pxi_outfile);
		dsX.print_Pinf(ppinf_outfile);
		}
	    } catch (exception& e){
	      std::cerr << "With L = "<<L<<", lambda = "<< lambda << ", trial "<<i;
	      std::cerr << "caught:\n " <<e.what() <<std::endl;
	    }
	}
    }
    pxi_outfile << "]";

    ppinf_outfile << "]";
    
}
    pxi_outfile << "}";
    pxi_outfile.close();
    
    ppinf_outfile << "}";
    ppinf_outfile.close();
    /*
    std::ofstream cluster_outfile;
    cluster_outfile.open("/tmp/clusters.json");
    cluster_outfile.close();
    dsX.print_cluster_json(cluster_outfile,1);
    */
    
}

