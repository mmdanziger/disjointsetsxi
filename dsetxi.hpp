#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <exception>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <bitset>
#include <cassert>

using std::pair;
using std::make_pair;
using std::vector;
using std::logic_error;
using std::map;
using std::cout;
using std::cerr;
using std::setprecision;
long double eps = 1e-7;

struct node{
 int x_rel,y_rel;
 int next_id_in_cluster=-1;
 int containing_cluster_key;
 node(){};
 node(int index, int L):x_rel(index%L),y_rel(index/L),containing_cluster_key(index){}
 //node(int index, int L):x_rel(0),y_rel(0),containing_cluster_key(index){}
};

struct cluster{
  double xcm=0,ycm=0;
  long m=0;
  long double inertia=0;
  int first_node_id = -1;
  int last_node_id = -1;
  std::bitset<4> touching_edges;//ordered: top right bottom left
  cluster(){}
  cluster(int x, int y, int index):xcm(x),ycm(y),m(1),first_node_id(index),last_node_id(index){}
  void initialize_edges(int L) {
      touching_edges[0]=ycm==0;
      touching_edges[1]=xcm==L-1;
      touching_edges[2]=ycm==L-1;
      touching_edges[3]=xcm==0;
  }
};



class DSetXi {
protected:
    const int L,N;
    map<int,cluster> clusters;
    vector<node> nodes;
    vector<pair<double,double>> xi_pair_vec;
    vector<pair<double,double>> pinf_pair_vec;
    bool periodic_bc=true;
    bool found_giant=false;
    int giant_id;
    long double sumMI;
    long double sumM2;
    float zero_fraction;
    bool stop_at_giant;
    std::vector<bool> disconnected;
    bool avoid_disconnected=false;
public:
    DSetXi(int Linput);
    void merge_components(int s, int t);
    pair<int,int> distance_pair(int s, int t);
    void update_cm_inertia(cluster &sCluster, const cluster tCluster);
    double get_xi(int print_output=0);
    double get_Pinf();
    template<typename Stream> void print_cluster_map(Stream& stream);
    template <typename Stream> void print_k_xi(Stream& stream);
    template <typename Stream> void print_Pinf(Stream& stream);
    template <typename Stream> void print_cluster_info(Stream& stream, int min_size);
    template <typename Stream> void print_cluster_json(Stream& stream, int min_size, bool giant_coords=false);
    template <typename Stream> void print_fractal_dimension(Stream& stream, int min_size);
};

DSetXi::DSetXi(int Linput)
    :L(Linput),N(Linput*Linput){
        sumMI=0;
        sumM2=N;
        nodes.resize(N);
        for (int i =0; i<N; i++){
            nodes[i] = node(i,L);
            clusters[i] = cluster(i%L,i/L,i);
            clusters[i].initialize_edges(L);
        }
    zero_fraction = L>=1000 ? 10.0 / L : 0.05;
    }
 

void DSetXi::merge_components(int s, int t)
{
#ifdef CHECKXI
  double lastxi = get_xi();
#endif
  
    auto s_it = clusters.find( nodes[s].containing_cluster_key );
    auto t_it = clusters.find( nodes[t].containing_cluster_key );

        if (t_it->second.m > s_it->second.m){ //we want the algorithm to treat s as the bigger set to which t is merged
            std::swap(s,t);
            std::swap(s_it,t_it);
        }
        cluster& sCluster = s_it->second;
        cluster& tCluster = t_it->second;
        int sKey = s_it->first;
        auto sx = nodes[s].x_rel;
        auto sy = nodes[s].y_rel;
        auto tx = nodes[t].x_rel;
        auto ty = nodes[t].y_rel;
        auto dl = distance_pair(s,t);
        int Cx = sx - tx + dl.first;
        int Cy = sy - ty + dl.second;
    if (s_it != t_it){
        int next_node_id = tCluster.first_node_id;
        while (next_node_id != -1){
            node& child = nodes[next_node_id];
            child.x_rel += Cx;
            child.y_rel += Cy;
            child.containing_cluster_key = sKey;
            next_node_id = child.next_id_in_cluster;
        }
        tCluster.xcm += Cx;
        tCluster.ycm += Cy;
        nodes[sCluster.last_node_id].next_id_in_cluster = tCluster.first_node_id; 
        sCluster.last_node_id = tCluster.last_node_id;
	
	if(found_giant && s_it->first == giant_id){ //if s is the giant, just remove the t-contributions
	  sumM2 -= (tCluster.m*tCluster.m);
	  sumMI -= (tCluster.m*tCluster.inertia);
	  update_cm_inertia(sCluster, tCluster);  ///<--- the value of sCluster.inertia is unreliable because cm is undefined
	  
	} else {
	
	  sumM2 -= (tCluster.m*tCluster.m + sCluster.m*sCluster.m);
	  sumMI -= (tCluster.m*tCluster.inertia + sCluster.m*sCluster.inertia);
	  update_cm_inertia(sCluster, tCluster);
          sumM2 += sCluster.m*sCluster.m;
          sumMI += sCluster.m*sCluster.inertia;
        
	}
        //cout << sCluster.touching_edges << " + " <<tCluster.touching_edges << " = ";
        sCluster.touching_edges |= tCluster.touching_edges;
        //cout <<sCluster.touching_edges <<std::endl;
        clusters.erase(t_it);
        
    } else {
        if((periodic_bc && (Cx!=0 || Cy !=0)) || (!periodic_bc && (sCluster.touching_edges.all() && sCluster.m > zero_fraction*N))) {
//             std::cout<<"Giant @ "<< s_it->first<<"!\n";
//             print_cluster_info(std::cout,3);
//             print_cluster_map(std::cout);
            if(!found_giant){
                giant_id=s_it->first;
                found_giant=true;
                sumM2 -= sCluster.m * sCluster.m;
                sumMI -= sCluster.m * sCluster.inertia;
                //cout << "Found giant @ " << giant_id << "\n";
            }
            else if (giant_id != s_it->first){
                std::cerr << "After "<< giant_id << " Found second spanning component at "<<s_it->first << "... Not expected in lattice.\n";
                std::cerr << "Discovered by link : "<<s<<" -- "<<t<<std::endl;
                //print_cluster_info(std::cout, 3);
                std::ofstream cmap_file;
                cmap_file.open("/tmp/clustermap.json");
                print_cluster_map(cmap_file);
                cmap_file.close();
                throw logic_error("Two giant components detected");
            }
        }

        }
#ifdef CHECKXI
 double xi = get_xi();
 if (xi > lastxi && found_giant){
   
                std::cerr << "Xi increased above pc.\n";
                std::cerr << "Discovered by link : "<<s<<" -- "<<t<<std::endl;
                print_cluster_info(std::cout, 3);
                std::ofstream cmap_file;
                cmap_file.open("/tmp/clustermap.json");
                print_cluster_map(cmap_file);
                cmap_file.close();
 }
#endif
        
    }
void DSetXi::update_cm_inertia(cluster &sCluster, const cluster tCluster)
{
    int M = sCluster.m + tCluster.m;
    double xcm = (sCluster.m * sCluster.xcm + tCluster.m * tCluster.xcm)/M;
    double ycm = (sCluster.m * sCluster.ycm + tCluster.m * tCluster.ycm)/M;
    double I = sCluster.inertia + tCluster.inertia + 
        sCluster.m * (xcm - sCluster.xcm)*(xcm - sCluster.xcm) + 
        sCluster.m * (ycm - sCluster.ycm)*(ycm - sCluster.ycm) +
        tCluster.m * (xcm - tCluster.xcm)*(xcm - tCluster.xcm) + 
        tCluster.m * (ycm - tCluster.ycm)*(ycm - tCluster.ycm);
    sCluster.m = M;
    sCluster.xcm = xcm;
    sCluster.ycm = ycm;
    sCluster.inertia = I;

}

    
pair< int, int > DSetXi::distance_pair(int s, int t)
{
    int dx =  t%L - s%L; //had this switched for reasons i can't remember, this is correct c-ordering
    int dy =  t/L - s/L;
    if(periodic_bc){
        if (dx > L/2)
            dx=dx-L;//L-dx?
        else if (dx < -L/2)
            dx=L+dx;
        if (dy > L/2)
            dy=dy-L;//L- dy?
        else if (dy < -L/2)
            dy = L +dy;   
    }
    return make_pair(dx,dy);
}

double DSetXi::get_Pinf(){
    if(found_giant){
       auto giant = clusters.find(giant_id)->second;
       return static_cast<long double>(giant.m) / static_cast<long double>(N);
    }
    return 0;
}

double DSetXi::get_xi(int print_output)
{
    long double correctSumM2=0,correctSumMI=0;
    /*if(found_giant){
       cluster giant = clusters.find(giant_id)->second;
       correctSumM2 = sumM2 - ((giant.m)*(giant.m));
       correctSumMI = sumMI - ((giant.m)*(giant.inertia));
       //cout << "Above p_c.  Returning " << "std::sqrt("<<nogiantmi<<"/"<<nogiantm2<<") "<< " instead of "<< "std::sqrt("<<correctSumMI<<"/"<<sumM2<<") "<<"\n";

    } else {*/
      correctSumM2 = sumM2;
      correctSumMI = sumMI;
    //}
    if (correctSumM2<1 || correctSumMI <0)
        return 0;
    if (print_output == 0){
      double xi = std::sqrt(correctSumMI/correctSumM2);
      if (xi<0 || correctSumM2 < 0 || correctSumMI < 0){
        std::cerr<<"Xi < 0: sqrt("<<correctSumMI << " / " << correctSumM2 << ")\n";
	std::ofstream clusterjsonfile("/tmp/clusterdata.json");
	print_cluster_json(clusterjsonfile,1);
	clusterjsonfile.close();
        throw logic_error("Negative Xi");
      } else if (isnan(xi)){
        std::cerr<<"Xi is nan: sqrt("<<correctSumMI << " / " << correctSumM2<< ")\n";
        throw logic_error("Xi is nan");        
      }
     return xi; 
    }
    long double m2sum=0;
    long double mIsum=0;
    if(avoid_disconnected && disconnected.size() != static_cast<std::size_t>(N)){
      std::cerr << "Expected disconnected vec size: " << N << " Actual size: "<<disconnected.size() <<std::endl;
      throw std::logic_error("Trying to access disconnected when not initialized.");
    }
    long contributing_clusters = 0;
    for(auto it = clusters.cbegin(); it!=clusters.cend(); it++){
	   
           if( (found_giant && it->first == giant_id) || (avoid_disconnected && disconnected[it->second.first_node_id]))
               continue;
            const cluster& this_cluster = it->second;
            mIsum += (this_cluster.m * this_cluster.inertia);
            m2sum += (this_cluster.m * this_cluster.m);
	   contributing_clusters++;
    }
    long double xi2 = mIsum / static_cast<long double>(m2sum);
    double xi =  std::sqrt(xi2);
    if (xi<0 || m2sum < 0 || mIsum < 0){
        std::cerr<<"Xi < 0: sqrt("<<mIsum << " / " << m2sum << ")\n";
        throw logic_error("Negative Xi");
    } else if (isnan(xi)){
        std::cerr<<"Xi is nan: sqrt("<<mIsum << " / " << m2sum << ")\n";
        throw logic_error("Xi is nan");        
    }
    if (fabs(mIsum - correctSumMI)/mIsum > eps || fabs(m2sum - correctSumM2)/m2sum > eps){
        cerr<<"Discrepancy between tracked Xi and summed Xi:\n";
        cerr<<"Tracked Xi is : sqrt("<<correctSumMI << " / " << correctSumM2 << ")\n";
        cerr<<"Summed  Xi is : sqrt("<<mIsum << " / " << m2sum << ")\n";
    }
    if (print_output==1)
        cout << "Xi ("<< contributing_clusters <<")= sqrt("<<mIsum << " / " << m2sum << ") = "<< xi <<"\n";
    return xi;
}

template <typename Stream>
void DSetXi::print_cluster_map(Stream& stream)
{
    stream<<"[";
    for(int i=0; i<L;i++){
        if(i!=0)
            stream<<",\n";
        stream<<"[";
        for(int j=0; j<L; j++){
            if(j!=0)
                stream<<",";
            stream <<nodes[i*L + j%L].containing_cluster_key;
        }
        stream<<"]";
    }
    stream<<"]\n";
}



template <typename Stream>
void DSetXi::print_k_xi(Stream& stream)
{
stream << "[";
bool firstVal=true;
for(auto dpair : xi_pair_vec){
    if(!firstVal)
        stream << ",\n";
    firstVal=false;
    stream << "[" << setprecision(10) << dpair.first << "," << setprecision(10) << dpair.second << "]";
}
stream <<"]\n";
}

template <typename Stream>
void DSetXi::print_Pinf(Stream& stream)
{
stream << "[";
bool firstVal=true;
for(auto dpair : pinf_pair_vec){
    if(!firstVal)
        stream << ",\n";
    firstVal=false;
    stream << "[" << setprecision(10) << dpair.first << "," << setprecision(10) << dpair.second << "]";
}
stream <<"]\n";
}

template <typename Stream>
void DSetXi::print_cluster_info(Stream& stream, int min_size)
{
    for(auto& clst_kv : clusters){
        cluster this_cluster = clst_kv.second;
        if (this_cluster.m<min_size)
            continue;
        stream << "Cluster with root: "<< this_cluster.first_node_id << " : \n";
        stream <<
        "\tm : " << this_cluster.m << 
        " x_cm : " << this_cluster.xcm <<
        " y_cm : " << this_cluster.ycm <<
        " I : " << this_cluster.inertia << "\n";
        int next_node_id = this_cluster.first_node_id;
        stream<<"abs:\t";
        while(next_node_id!=-1){
            node this_node = nodes[next_node_id];
            stream << "(" << next_node_id%L << "," << next_node_id/L<<")\t";
            next_node_id = this_node.next_id_in_cluster;
        }
        stream<<"\nrel:\t";
        next_node_id = this_cluster.first_node_id;
        while(next_node_id!=-1){
            node this_node = nodes[next_node_id];
            stream << "(" << this_node.x_rel << "," << this_node.y_rel<<")\t";
            next_node_id = this_node.next_id_in_cluster;
        }
        stream<<"\n";
    }
}

template <typename Stream>
void DSetXi::print_cluster_json(Stream& stream, int min_size, bool giant_coords)
{
    stream << "[";
    bool firstval0=true;
    for(auto& clst_kv : clusters){
        cluster this_cluster = clst_kv.second;
        if (this_cluster.m<min_size)
            continue;
        if (!firstval0)
            stream << ",\n";
        firstval0=false;
        std::string giant_string = (found_giant && giant_id == clst_kv.first)? "true" : "false";
        stream << "{\"root\": " << this_cluster.first_node_id ;
        stream <<
        ", \"m\" : " << this_cluster.m << 
        ", \"x_cm\" : " << this_cluster.xcm <<
        ", \"y_cm\" : " << this_cluster.ycm <<
        ", \"I\" : " << this_cluster.inertia <<
        ", \"giant\" : "<< giant_string;
	if(avoid_disconnected){
	  if (disconnected[this_cluster.first_node_id])
	    stream << ", \"disconnected\" : true";
	  else
	    stream << ", \"disconnected\" : false";
	}
	if(!( found_giant && giant_id == clst_kv.first && !giant_coords)){
	  int next_node_id = this_cluster.first_node_id;
	  stream<<", \"abs coords\": [";
	  bool firstval1=true;
	  while(next_node_id!=-1){
	      if (!firstval1)
		  stream << ",";
	      firstval1=false;
	      node this_node = nodes[next_node_id];
	      stream << "[" << next_node_id%L << "," << next_node_id/L<<"]\t";
	      next_node_id = this_node.next_id_in_cluster;
	  }
	  stream<<"],\n \"rel coords\" : [";
	  next_node_id = this_cluster.first_node_id;
	  firstval1=true;
	  while(next_node_id!=-1){
	      if (!firstval1)
		  stream << ",";
	      firstval1=false;            
	      node this_node = nodes[next_node_id];
	      stream << "[" << this_node.x_rel << "," << this_node.y_rel<<"]\t";
	      next_node_id = this_node.next_id_in_cluster;
	  }
	  stream <<"]";
	}
	stream <<"}";
    }
    stream<<"]\n";
}

template <typename Stream>
void DSetXi::print_fractal_dimension(Stream& stream, int min_size)
{
    stream << "[";
    bool firstval0=true;
    for(auto& clst_kv : clusters){
        cluster this_cluster = clst_kv.second;
        if (this_cluster.m<min_size)
            continue;
	if(avoid_disconnected && disconnected[this_cluster.first_node_id])
	    continue;
        if (!firstval0)
            stream << ",\n";
        firstval0=false;
	stream << "{ \"I\" : "<< this_cluster.inertia <<", \"m\" : "<<this_cluster.m<<"}";
    stream<<"]\n";
}