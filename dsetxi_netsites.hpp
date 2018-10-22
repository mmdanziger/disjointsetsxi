#include "dsetxi.hpp"
#include <chrono>
#include <omp.h>
#include <exception>
#include <boost/random/exponential_distribution.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "../../jsonutils.hpp"
#define verbosity 1

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef std::pair<int,std::pair<int,int> > DistancePair;
int randsign() { return rand() > (RAND_MAX/2 +1)? 1 : -1;}
int randint(int N) { return rand()%N; }

using std::endl; 
using std::exception;
using std::vector;

class DSetXiNetSites : public DSetXi {
private:
    vector<int> node_order;
    vector<bool> added_yet;
    Graph G;
    
    //Lambda net stuff
    boost::random::exponential_distribution<double> exp;
    double kavg,lambda;
    boost::mt19937 gen;
    std::vector<DistancePair> pencils;
    
public:
    DSetXiNetSites(int Linput, double kavg_input, double lambda_input);  
    void calculate_distances();
    void make_connectivity_links();
    int draw_target(int sNode);
    void link_targets(int s);
    void add_link(int s, int t);
    void incrementalComponents(double to_p, int node_res);
    double find_pc_only();
    vector<vector<int>> get_adjacency_list();

};

DSetXiNetSites::DSetXiNetSites(int Linput, double kavg_input, double lambda_input): 
DSetXi(Linput),kavg(kavg_input),lambda(lambda_input),
gen((static_cast<unsigned>((omp_get_thread_num() +1)* getpid()) * std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()))
{
        node_order.resize(N);
        added_yet.resize(N);
	disconnected.resize(N);
	avoid_disconnected = true;
        for (int i =0; i<N; i++){
            node_order[i] = i;
            added_yet[i] = false;
        }
        periodic_bc=false;
        std::shuffle(node_order.begin(), node_order.end(), gen);
        stop_at_giant=false;
        G=Graph(N);
        exp.param(lambda);
        calculate_distances();
        make_connectivity_links();

}


void DSetXiNetSites::calculate_distances()
{
    int uniquePairs = static_cast<int>(L*(L+1)/2) -1; //-1 to disallow (0,0) -> force at least (0,1)
    pencils.resize( uniquePairs );
    //cout << "Resized pencils" <<endl;
    int m=0;
    for (int i = 1; i<L;++i){//start at ONE: no dist=0 links
        for (int j=0; j<=i;++j){
            try{
                if (m >= uniquePairs){
                    cout <<"I'm sorry Dave. I can't let "<<m<<" do that to ("<<i<<","<<j<<")..."<<endl;
                    return;
                }
           pencils[m++] = std::make_pair(i*i+j*j, std::make_pair(i,j));
            }catch(exception &e){
                cout << "Pencil construction failed" <<endl;
                cout <<e.what() <<endl;
            }
        }
    }
    std::sort(pencils.begin(), pencils.end(), [](const DistancePair& dpA, const DistancePair& dpB){ return dpA.first < dpB.first;});

}

int DSetXiNetSites::draw_target(int sNode)
{
  double rexp,rexp2;
  do{
      rexp = exp(gen);
      rexp2 = rexp*rexp;}while(rexp > L/2 );
  auto it = std::upper_bound(pencils.cbegin(), pencils.cend(), rexp2, [](const int A, const DistancePair B){return A<B.first;});
  if(it!=pencils.cbegin() && ((*it).first + (*(it-1)).first > 2*rexp2 )  ){
    it--;
    if(it!=pencils.cbegin() && (*(it-1)).first == (*(it)).first){
      if (randsign() > 0)
	it--;
    }
  } else if ((it+1) != pencils.cend() && (*(it+1)).first == (*(it)).first){
      if (randsign() > 0)
	it++;
  }
  int dx=randsign()*(*it).second.first;
  int dy=randsign()*(*it).second.second;
  if (randsign()>0){//switch i and j with 50% chance because only one pair appears in the pencils list
    std::swap(dx,dy);
    }
  int tNode_x = sNode%L+dx;
  int tNode_y = sNode/L+dy;
  if(periodic_bc){
      tNode_x = (tNode_x +L)%L;
      tNode_y = (tNode_y +L)%L;
  } else {
      if (tNode_y >= L || tNode_y <0)
          tNode_y-=(2*dy);
      if (tNode_x >= L || tNode_x <0)
          tNode_x-=(2*dx);
      if (tNode_y >= L || tNode_y <0 || tNode_x >= L || tNode_x <0) //if you still can't find something, return failure
          return -1;
  } 
  return tNode_y*L+tNode_x;
}

void DSetXiNetSites::make_connectivity_links(){
    int give_up_count = N/2;
    for (int i =0; i<kavg*N/2; ++i){
        bool added=false;
        int attempt_count=0;
        while(!added){
        int sNode = randint(N);
        int tNode = draw_target(sNode);
          if (tNode > 0 && !boost::edge(sNode,tNode,G).second){
            boost::add_edge(sNode,tNode,G);//should i make sure that sNode<tNode? don't think it matters
            added=true;
          }
          if (++attempt_count > give_up_count){
              if(verbosity>0)
                  cout<<"Unable to find a neighbor for node "<< sNode <<std::endl;
            break;
          }
            
        }
    }
	  std::vector<int> component_vec(boost::num_vertices(G));
	  int n_comps = boost::connected_components(G,&component_vec[0]);
	  std::vector<int> compsizes(n_comps);
	  std::for_each(component_vec.cbegin(), component_vec.cend(),[&compsizes](const int &comp_id){compsizes[comp_id]++;});
// 	  jsonArray(compsizes,std::cout);
	  int giant_id = std::distance(compsizes.cbegin(), std::max_element(compsizes.cbegin(),compsizes.cend()));
	  sumM2 = compsizes[giant_id];
	  std::cout << "Disconnected : " << N - sumM2 << std::endl;
	  for(uint i=0; i<static_cast<uint>(N); i++)
	    disconnected[i] = component_vec[i] == giant_id ? false : true;
	  
// 	  std::ofstream components_file("/tmp/components.json");
// 	  jsonArray(component_vec,components_file);
// 	  components_file.close();
}



void DSetXiNetSites::add_link(int s, int t)
{
    if(added_yet[s] && added_yet[t] && !(found_giant &&stop_at_giant))
        merge_components(s,t);

}


void DSetXiNetSites::link_targets(int s)
{
    auto neighbors = boost::adjacent_vertices(s,G);
    for(auto it=neighbors.first; it!=neighbors.second; ++it){
        add_link(s,*it);
    }
}

vector<vector<int>> DSetXiNetSites::get_adjacency_list(){
	  vector<vector<int>> adjacency_list(N);
	  auto epair=boost::edges(G);
	  for(auto ei = epair.first; ei!=epair.second; ++ei)
	  {
	    int sNode = boost::source(*ei,G);
	    int tNode = boost::target(*ei,G);
	    adjacency_list[sNode].push_back(tNode);
	    adjacency_list[tNode].push_back(sNode);
	  }
	return adjacency_list;
}


double DSetXiNetSites::find_pc_only()
{
    int s,node_count=0;
    while(node_count<N){
        s = node_order[node_count];
        added_yet[s] = true;
        link_targets(s);
        if(found_giant)
            break;
        node_count++;
    }
    return static_cast<double>(node_count)/static_cast<double>(N);
}


void DSetXiNetSites::incrementalComponents(double to_p, int node_res)
{
    long to_node = to_p * N;
    int s,node_count;
    node_count=0;
    while(node_count < to_node){
        s = node_order[node_count];
        added_yet[s]=true;
	if (disconnected[s]){ // nodes that are not in the GCC at p = 1 are not included at all
	 node_count++;
	  continue;
	}
        link_targets(s);
        if( (found_giant && stop_at_giant) || node_res < 2 || node_count%node_res==0 || (node_count+1 >= to_node)){
            double current_p = static_cast<double>(node_count)/static_cast<double>(N);
            double currentXi = get_xi();
            xi_pair_vec.emplace_back(make_pair(current_p,currentXi)); 
            pinf_pair_vec.emplace_back(make_pair(current_p,get_Pinf()));
	    
            if (currentXi > L ){
                std::cerr<<"Bad Xi encountered, check output.\n";
                break;
            }
            if (found_giant && stop_at_giant){
                break;
            }
	    if(current_p > 0.85 && currentXi > 20*(1/lambda)){
		    std::ofstream cmap_file,cjson_file;
		    cmap_file.open("/tmp/clustermap.json");
		    cjson_file.open("/tmp/clusterdata.json");
		    print_cluster_json(cjson_file, 1);
		    print_cluster_map(cmap_file);
		    cmap_file.close();
		    cjson_file.close();
		    std::cerr << "Found a bad one!\n";  
		    throw logic_error("Xi is too large");
	    }
        }
       
    node_count++;
        
    }
}
