#include <iostream>
#include "DSetXiNetSites8.hpp"
#include <omp.h>
int main(int argc, char **argv) {
    long L =1000;
    long node_res = 100;
    double until_p = 0.5;
    int tests = 1000;
    int steps = (L*L)*until_p/node_res;
    std::vector<std::vector<double>> xi_output(steps);
    for(int i=0; i<steps; i++){
      xi_output[i].resize(tests);
    }
#pragma omp parallel for
    for(long test=0; test<tests; test++){
      DSetXiNetSites8 D(L,time(0)*omp_get_thread_num()+test);
      D.incrementalComponents(until_p,node_res);
      //D.writeDataToFile(test, "/tmp/xi_results");
      int i=0;
      for(auto &pxi : D.get_xi_pair_vec()){
	xi_output[i][test] = pxi.second;
	i++;
      }
      
    }
    int i=0;
    std:ofstream of("/tmp/newout.txt");
    for(auto &xivec : xi_output){
      of << static_cast<double>(++i*node_res)/static_cast<double>(L*L) << "\t";
      for(auto &xival : xivec){
	of << xival <<"\t";
      }
      of << std::endl;
      
    }
	
    return 0;
}
