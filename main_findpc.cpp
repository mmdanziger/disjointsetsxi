#include "dsetxi_netsites.hpp"
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <omp.h>
#include <iomanip>
using std::cout;
int main(int argc, char **argv) {
    srand(static_cast<unsigned>(getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    omp_set_num_threads(8);
    unsigned L = 2000;
    int trials=10;
    double kavg = 3;
    std::ofstream pinf_file;
    std::vector<double> lambda_vec{ 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5};
    bool firstlambda=true;
    pinf_file.open("pc_L2000_k3_multiplelambda.json");
    pinf_file <<"{";
    
for (auto& lambda : lambda_vec){    
    if (!firstlambda){
        pinf_file <<",\n";
    }
    firstlambda=false;
    pinf_file << "\""<<lambda<<"\" : [";
    bool firsttrial=true;
#pragma omp parallel shared(firsttrial)
    {
#pragma omp for
    for(int i=0; i<trials; i++){
        DSetXiNetSites dsX(L,kavg,lambda);
        double pc = dsX.find_pc_only();
     #pragma omp critical 
        {
     #pragma omp flush(firsttrial) 

        if(!firsttrial){
            pinf_file <<",\n";
        }
        firsttrial=false;
        pinf_file << std::setprecision(10) << pc;
        pinf_file.flush();
        }
      }
    }
    pinf_file << "]";    
}
    pinf_file << "}";
    pinf_file.close();

}

