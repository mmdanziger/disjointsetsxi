#include <iostream>
#include "corrlength.hpp"
#include <boost/program_options.hpp>
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <sys/types.h>
#include <unistd.h>

namespace po = boost::program_options;
int main(int argc, char **argv) {
    srand(static_cast<unsigned>(getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    unsigned L = 200;
    double lambda = 10;
    double kavg=4;
    NetworkType network_type=Lambda;
    std::string config_file;
    try {

        po::options_description desc("Supported options",1024);
        desc.add_options()
        ("help", "produce help message")
        ("L,L",po::value<unsigned>(&L)->default_value(100), "L")
        ("lambda,l",po::value<double>(&lambda)->default_value(-1), "lambda")
        ("kavg,k",po::value<double>(&kavg)->default_value(4),"kavg (initial)")
        ("lattice","check holes on lattice")
        ("random","generate ER network")
        ("config",po::value<string>(&config_file)->default_value("interdep.cfg"),"path to config file");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        ifstream ifs(config_file.c_str());
        if (ifs){
            cout<<"Config file read ("<<config_file<<").\n";
            po::store(po::parse_config_file(ifs,desc),vm);
            po::notify(vm);
        }
            if ((vm.count("help"))){
                cout << desc << "\n";
                return 1;
            }
            if (vm.count("random") || lambda < 0)
                network_type = ER;
            if (vm.count("lattice"))
                network_type = InverseLattice;
            
           }catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    return 1;
    } 
        
    XiMeasurer xm(L,lambda,kavg,network_type);
    double p_step=0.01;
    vector<pair<double,double>> xi_results;
    double p = 0.6;
    while(p<0.7){
     xi_results.emplace_back(std::make_pair(p,xm.measure_xi(p)));
     p+=p_step;
    }
    std::sort(xi_results.begin(), xi_results.end(), [](pair<double,double> a, pair<double,double> b){return a.first<b.first;});
    bool first=true;
    std::cout<<"[";
    for(auto pr : xi_results){
        if (!first)
            std::cout<<',';
        std::cout << "\n["<<pr.first<<","<<pr.second<<"]";
    }
    std::cout<<"]";
    return 0;
}
