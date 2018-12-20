

#ifndef SRC_DSETXI_HPP_
#define SRC_DSETXI_HPP_


#include "my_structs.hpp"
#include <math.h>
using namespace std;
long double eps = 1e-2;
#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

class DSetXi {
protected:
    int L,N;
    //map<int,cluster> clusters;
    vector<node> nodes;
    vector<cluster> clusters;
    vector<int> alive;
    vector<pair<double,double>> xi_pair_vec;
    vector<pair<double,double>> xi_with_giant_pair_vec;
    vector<pair<double,double>> pinf_pair_vec;
    vector<pair<double,double>> largest_pair_vec;
    vector<pair<double,double>> second_pair_vec;
    vector<pair<double,double>> mean_mass_vec;


    vector<pair<double,double>> m2_pair_vec;
    bool periodic_bc=true;
    bool found_giant=false;
    int giant_id;
    long double sumMI;
    //long double sumMI;
    long double sumM2;
    float zero_fraction;
    //for test
    long double sumSN;
    double largest_size=1, second_size=1;
    int largest_id = 0 ,second_id = 1;
    bool merge_flag = false;
 
    bool ERflag = true;
    
    int clusters_number=0;
    int node_number=0;

public:
    DSetXi(int Linput);
    bool merge_components(int s, int t);
    void distance_pair(int s, int t, int &Cx, int &Cy);
    void update_cm_inertia(cluster &sCluster, const cluster tCluster);
    double get_xi(int print_output);
    double get_xi_with_giant();
    double get_Pinf();
    void recalculate_largest_clusters(int l_id);
    virtual void calc_radius(int l_id);
    template<typename Stream> void print_cluster_map(Stream& stream);
    template <typename Stream> void print_k_xi(Stream& stream);
    template <typename Stream> void print_Pinf(Stream& stream);
    template <typename Stream> void print_cluster_info(Stream& stream, int min_size);
    template <typename Stream> void print_cluster_json(Stream& stream, int min_size,vector<int>* connected=0);
    void check_data();
};
/*#ifdef 64BIT
using real=double;
#else
using real=float;
#endif*/
//using real=std::setprecision(17)double;

DSetXi::DSetXi(int Linput)
    :L(Linput),N(Linput*Linput),nodes(N),clusters(N),alive(N){
        sumMI=0;
        sumM2=0;
        sumSN=0;
        for (int i =0; i<N; i++){
            nodes[i] = node(i,L);
            clusters[i] = cluster(i%L,i/L,i);
            //clusters.push_back(cluster(i%L,i/L,i));
            alive[i] = 0;
            clusters[i].initialize_edges(L);
        }
    //zero_fraction = L>=1000 ? 10.0 / L : 0.05;
    }
void DSetXi::calc_radius(int l_id){
	cerr<<"execute function in the base class!!!!!, should be redefined in the deriven class!!!"<<endl;
};
bool DSetXi::merge_components(int s, int t)
{
    /*auto s_it = clusters.find( nodes[s].containing_cluster_key );
    auto t_it = clusters.find( nodes[t].containing_cluster_key );*/
	cluster sCluster = clusters[nodes[s].containing_cluster_key];
	cluster tCluster = clusters[nodes[t].containing_cluster_key];

       /* if (t_it->second.m > s_it->second.m){ //we want the algorithm to treat s as the bigger set to which t is merged
            std::swap(s,t);
            std::swap(s_it,t_it);
        }*/
	      if (tCluster.m > sCluster.m){ //we want the algorithm to treat s as the bigger set to which t is merged
	    	std::swap(s,t);
	      }
	    node s_node=nodes[s] ,t_node=nodes[t];
        /*cluster& sCluster = s_it->second;
        cluster& tCluster = t_it->second;*/
        //int sKey = s_it->first;
	    int sKey = s_node.containing_cluster_key;
	    int tKey = t_node.containing_cluster_key;
	    auto sCluster_it = clusters.begin() + sKey;
	    auto tCluster_it = clusters.begin() + tKey;
        auto sx = s_node.x_rel;
        auto sy = s_node.y_rel;
        auto tx = t_node.x_rel;
        auto ty = t_node.y_rel;
        int Cx = sx - tx;
        int Cy = sy - ty;
        distance_pair(s,t,Cx,Cy);
    //if (s_it != t_it){
      if (sKey != tKey){
    	  sumSN = sumSN - sCluster.m*sCluster.m - tCluster.m*tCluster.m +
    	  			(tCluster.m + sCluster.m)*(tCluster.m + sCluster.m);
    	  clusters_number--;
        //int next_node_id = tCluster.first_node_id;
    	  int next_node_id = tCluster_it->first_node_id;
        while (next_node_id != -1){
            node& child = nodes[next_node_id];
            child.x_rel += Cx;
            child.y_rel += Cy;
            child.containing_cluster_key = sKey;
            next_node_id = child.next_id_in_cluster;
        }

        tCluster_it->xcm += Cx;
        tCluster_it->ycm += Cy;
        nodes[sCluster_it->last_node_id].next_id_in_cluster = tCluster_it->first_node_id;
        sCluster_it->last_node_id = tCluster_it->last_node_id;
        double long sm2 = sCluster_it->m*sCluster_it->m;
        double long  smi = sCluster_it->m*sCluster_it->inertia;
       /* if(found_giant && s_it->first == giant_id){ //if s is the giant, don't subtract it
            sm2=0;
            smi=0;
        }*/
        if(found_giant && sKey == giant_id){ //if s is the giant, don't subtract it
            sm2=0;
            smi=0;
        }


        sumM2 -= (tCluster_it->m*tCluster_it->m +sm2);
        sumMI -= (tCluster_it->m*tCluster_it->inertia + smi);
        update_cm_inertia(*sCluster_it, *tCluster_it);


        /*if (! (found_giant && s_it->first == giant_id)){ //if s is the giant, don't add it to the tally
            sumM2 += sCluster.m*sCluster.m;
            sumMI += sCluster.m*sCluster.inertia;
        }*/
        if (! (found_giant && sKey == giant_id)){ //if s is the giant, don't add it to the tally
            sumM2 += sCluster_it->m*sCluster_it->m;
            sumMI += sCluster_it->m*sCluster_it->inertia;
        }
        //cout << sCluster.touching_edges << " + " <<tCluster.touching_edges << " = ";
        sCluster_it->touching_edges |= tCluster_it->touching_edges;
        //cout <<sClusperiodic_bcter.touching_edges <<std::endl;
        //clusters.erase(t_it);
        alive[tKey] = 0;
        
        if(sCluster_it->touching_edges.all() && sKey==largest_id && tKey==second_id){        
        		merge_flag = true;
         }
        
        if(sKey==largest_id && tKey==second_id)
        {
        	largest_size = sCluster_it->m;
        	recalculate_largest_clusters(sKey);
        }
        else
        {
        	if(sKey==largest_id)
        	{
        		largest_size=sCluster_it->m;
        	}
        	else if(sKey==second_id)
        	{
        		if(sCluster_it->m>largest_size)
        		{
        			second_size = largest_size;
        			second_id = largest_id;
        			largest_id = sKey;
        			largest_size = sCluster_it->m;
        		}
        		else
        		{
        			second_size = sCluster_it->m;
        		}
        	}
        	else
        	{
        		if(sCluster_it->m>largest_size)
        		{
        			second_size = largest_size;
        			second_id = largest_id;
        			largest_id = sKey;
        			largest_size = sCluster_it->m;
        		}
        		else if(sCluster_it->m>second_size)
        		{
        			second_size = sCluster_it->m;
        			second_id = sKey;
        		}
        	}
        }
  
       
           
        
    } else {
    	//if((periodic_bc && (Cx!=0 || Cy !=0)) || (!periodic_bc && (sCluster.touching_edges.all() && sCluster.m > zero_fraction*N))) {
    	if(( periodic_bc && (Cx!=0 || Cy !=0) ) || (!periodic_bc  && sCluster_it->touching_edges.all() && merge_flag)){
//             std::cout<<"Giant @ "<< s_it->first<<"!\n";
//             print_cluster_info(std::cout,3);
//             print_cluster_map(std::cout);
            if(!found_giant){
                //giant_id=s_it->first;
            	clusters_number--;
            	giant_id=sKey;
                found_giant=true;
                sumM2 -= sCluster_it->m * sCluster_it->m;
                sumMI -= sCluster_it->m * sCluster_it->inertia;
                //cout << "Found giant @ " << giant_id << "\n";
            }
            //else if (giant_id != s_it->first){
              //  std::cerr << "After "<< giant_id << " Found second spanning component at "<<s_it->first << "... Not expected in lattice.\n";
              else if (giant_id != sKey && !ERflag){
            	  //cout<<"two giants"<<endl;
            	  return false;
            	  std::cerr << "After "<< giant_id << " Found second spanning component at "<<sKey << "... Not expected in lattice.\n";
                std::cerr << "Discovered by link : "<<s<<" -- "<<t<<std::endl;
              //  print_cluster_info(std::cout, 3);
                std::ofstream cmap_file;
                cmap_file.open("/tmp/clustermap.json");
            //    print_cluster_map(cmap_file);
                cmap_file.close();
                throw logic_error("Two giant components detected");
            }
        }

        }
      return true;
    }
void DSetXi::recalculate_largest_clusters(int l_id)
{
	double s_size=0;
	int s_id;
	for(int k=0; k<clusters.size() ;k++)
	   {
			if(alive[k] && k!=l_id)
			{
				if(clusters[k].m>s_size)
				{
					s_size=clusters[k].m;
					s_id = k;
				}
			}
	   }
	second_size = s_size;
	second_id = s_id;	
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
inline void DSetXi::distance_pair(int s, int t, int &Cx, int &Cy)
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
    Cx+=dx;
    Cy+=dy;
}
    
double DSetXi::get_Pinf(){
   /* if(found_giant){
       auto giant = clusters.find(giant_id)->second;
       return static_cast<long double>(giant.m) / static_cast<long double>(N);

    }*/
	if(found_giant){
	  auto giant = clusters[giant_id];
	  return static_cast<long double>(giant.m) / static_cast<long double>(N);
	 }
    return 0;
}

double DSetXi::get_xi(int print_output=0)
{
	/*if(found_giant){
	       cluster giant = clusters.find(giant_id)->second;
	       double long nogiantm2 = sumM2 - ((giant.m)*(giant.m));
	       double long nogiantmi = sumMI - ((giant.m)*(giant.inertia));

		       //cout << "Above p_c.  Returning " << "std::sqrt("<<nogiantmi<<"/"<<nogiantm2<<") "<< " instead of "<< "std::sqrt("<<sumMI<<"/"<<sumM2<<") "<<"\n";
	       if (nogiantm2 <1 || nogiantmi<1)
	    	   return 0;
	       return sqrt(nogiantmi/nogiantm2);

	}*/
   
    if (sumM2<1 || sumMI<1)
        return 0;
   return sqrt(sumMI/sumM2);
  
	
    /*//long long m2sum=0;
	long double m2sum=0;
	long double mIsum=0;
    for(auto it = clusters.begin(); it!=clusters.end(); it++){
           if((found_giant && it->first == giant_id))
               continue;
            const cluster& this_cluster = it->second;
            mIsum += (this_cluster.m * this_cluster.inertia);
            m2sum += (this_cluster.m * this_cluster.m);
    }
    long double xi2 = mIsum / (m2sum);
    long double xi =  sqrt(xi2);
    if (xi<0 || m2sum < 0 || mIsum < 0){
        std::cerr<<"Xi < 0: sqrt("<<mIsum << " / " << m2sum << ")\n";
        throw logic_error("Negative Xi");
    } else if (std::isnan(xi)){
    	return 0;
        std::cerr<<"Xi is nan: sqrt("<<mIsum << " / " << m2sum << ")\n";
        throw logic_error("Xi is nan");
    }
    if (fabs(mIsum - sumMI) > eps || fabs(m2sum - sumM2) > eps){
        cerr<<"Discrepancy between tracked Xi and summed Xi:\n";
        cerr<<"Tracked Xi is nan: sqrt("<<sumMI << " / " << sumM2 << ")\n";
        cerr<<"Summed  Xi is nan: sqrt("<<mIsum << " / " << m2sum << ")\n";
    }
    if (print_output==1)
        cout << "Xi = sqrt("<<mIsum << " / " << m2sum << ") = "<< xi <<"\n";
    return xi;*/
}
double DSetXi::get_xi_with_giant()
{
		if(found_giant){
	       cluster giant = clusters[giant_id];
	       double long withgiantm2 = sumM2 + ((giant.m)*(giant.m));
	       double long withgiantmi = sumMI + ((giant.m)*(giant.inertia));


	       //cout << "Above p_c.  Returning " << "std::sqrt("<<nogiantmi<<"/"<<nogiantm2<<") "<< " instead of "<< "std::sqrt("<<sumMI<<"/"<<sumM2<<") "<<"\n";
	       if (withgiantm2 <1 || withgiantmi<1)
	           return 0;
	       return sqrt(withgiantmi/withgiantm2);

	    }
	  
	    if (sumM2<1 || sumMI<1)
	        return 0;
	   return sqrt(sumMI/sumM2);
	
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
    stream << "[" << std::setprecision(10) << dpair.first << "," << setprecision(10) << dpair.second << "]";
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
		for(int i=0;i<N;i++){
    /*for(auto& clst_kv : clusters){
        cluster this_cluster = clst_kv.second;
        if (this_cluster.m<min_size)
            continue;*/
    	cluster this_cluster = clusters[i];
    	/*if (this_cluster.m<min_size || !alive[i])
    	     continue;*/
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
//void DSetXi::print_cluster_json(Stream& stream, int min_size, vector<int>* connected)
void DSetXi::print_cluster_json(Stream& stream, int min_size, vector<int>* connected)
{
    stream << "[";
    bool firstval0=true;
    for(int i=0; i<N; i++){
    /*for(auto& clst_kv : clusters){
        cluster this_cluster = clst_kv.second;
        if (this_cluster.m<min_size || !connected->at(this_cluster.first_node_id))
            continue;*/
    	cluster this_cluster = clusters[i];
    /*if(connected!=0)
    	if (this_cluster.m<min_size || !connected->at(this_cluster.first_node_id) ||!alive[i])
	            continue;
    else
    	if (this_cluster.m<min_size ||!alive[i])
    		            continue;*/
        if (!firstval0)
            stream << ",\n";
        firstval0=false;
        //std::string giant_string = (found_giant && giant_id == clst_kv.first)? "true" : "false";
        std::string giant_string = (found_giant && giant_id == i)? "true" : "false";
        stream << "{\"root\": " << this_cluster.first_node_id ;
        stream <<
        ", \"m\" : " << this_cluster.m <<
        ", \"x_cm\" : " << this_cluster.xcm <<
        ", \"y_cm\" : " << this_cluster.ycm <<
        ", \"I\" : " << this_cluster.inertia <<
        ", \"giant\" : "<< giant_string;
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
        stream <<"]}";
    }
    stream<<"]\n";
}
void DSetXi::check_data(){
	double xi, yi, xi_2, yi_2,m;
	double R, I, R2, m2 =0, im = 0;
	for(int i=0;i<clusters.size();i++){
		xi = 0;
		yi = 0;
		xi_2 = 0;
		yi_2 = 0;
		m = clusters[i].m;
		I = clusters[i].inertia;
  	    int next_node_id = clusters[i].first_node_id;
        while (next_node_id != -1){
            node& child = nodes[next_node_id];
            xi += child.x_rel;
            yi += child.y_rel;
            xi_2 += child.x_rel*child.x_rel;
            yi_2 += child.y_rel*child.y_rel;
            next_node_id = child.next_id_in_cluster;
        }
        xi/=m;
        yi/=m;
        yi_2/=m;
        xi_2/=m;


        R = sqrt((double)(xi_2 + yi_2 - xi*xi - yi*yi));
        R2 = sqrt((double)I/m);

        if(alive[i] and i !=giant_id){
        	m2 +=m*m;
        	im +=I*m;
        }


        if(abs(R  - R2) > 0.00001 and alive[i]){
        	cout<<"Error!"<<endl;
        	cout<<"R1 = "<<R<<"   R2 = "<<R2<<endl;
        	cout<<"m = "<<clusters[i].m<<endl;
        	cout<<"xi = "<<xi<<endl;
        	cout<<"yi = "<<yi<<endl;
        	cout<<"xi_2 = "<<xi_2<<endl;
        	cout<<"yi_2 = "<<yi_2<<endl;

      	    int next_node_id = clusters[i].first_node_id;
            while (next_node_id != -1){
                node& child = nodes[next_node_id];
                cout<<"node"<<next_node_id<<"  ("<<child.x_rel<<","<<child.y_rel<<")"<<endl;
                next_node_id = child.next_id_in_cluster;
            }

        }

	}

	if(abs(m2 - sumM2) > 0.0001 or abs(im - sumMI) > 0.0001){
		cout<<"Error!"<<endl;
		cout<<"m2="<<m2<<"     sumM2="<<sumM2<<endl;
		cout<<"im="<<im<<"     sumMI="<<sumMI<<endl;
	}

}


#endif
