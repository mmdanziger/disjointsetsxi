#include "dsetxi.hpp"


class DSetXiSites : public DSetXi {
private:
    vector<int> node_order;
    vector<bool> added_yet;
public:
    DSetXiSites(int Linput);  
    void link_targets(int s);
    void add_link(int s, int t);
    void incrementalComponents(double to_p, int node_res);
};

DSetXiSites::DSetXiSites(int Linput) : DSetXi(Linput)
    {

        node_order.resize(N);
        added_yet.resize(N);
        for (int i =0; i<N; i++){
            node_order[i] = i;
            added_yet[i] = false;
        }
	periodic_bc = true;
	stop_at_giant = false;
        std::random_shuffle(node_order.begin(), node_order.end());
    }

void DSetXiSites::incrementalComponents(double to_p, int node_res)
{
    long to_node = to_p * N;
    int s,node_count;
    node_count=0;
    while(node_count < to_node){
        s = node_order[node_count]; 
        added_yet[s]=true;
        link_targets(s);
        if( (found_giant && stop_at_giant)|| node_res < 2 || node_count%node_res==0){
            double current_p = static_cast<double>(node_count)/static_cast<double>(N);
            double currentXi = get_xi();
	    xi_pair_vec.emplace_back(make_pair(current_p,currentXi)); 
	    pinf_pair_vec.emplace_back(make_pair(current_p,get_Pinf()));

            if (currentXi > L ){
                std::cerr<<"Bad Xi encountered, check output.\n";
                break;
            }
            if (found_giant && stop_at_giant){
                //std::cout<<"Found giant at p = "<<current_p<<"\n";
                break;
            }
            
        }
    node_count++;
        
    }
}


void DSetXiSites::add_link(int s, int t)
{
    if(added_yet[s] && added_yet[t] && !(found_giant && stop_at_giant))
        merge_components(s,t);
}

    
void DSetXiSites::link_targets(int s)
{//add links for all 8 neighbors of s
   
    if ( s == 0){
            add_link(0,1);
            add_link(0,L-1);
            add_link(0,L);
            add_link(0,N-L);
            add_link(0,L+1);
            add_link(0,2*L-1);
            add_link(0,N-L+1);
            add_link(0,N-1);
    } else if ( s == L-1) {
            add_link(L-1,0);
            add_link(L-1,L-2);
            add_link(L-1,2*L-1);            
            add_link(L-1,N-1);
            add_link(L-1,L);
            add_link(L-1,2*L-2);
            add_link(L-1,N-L);
            add_link(L-1,N-2);            
    } else if ( s ==  N-L) {
            add_link(N-L,N-L+1);
            add_link(N-L,N-1);
            add_link(N-L,0);
            add_link(N-L,N-2*L);
            add_link(N-L,1);
            add_link(N-L,L-1);
            add_link(N-L,N-2*L+1);
            add_link(N-L,N-L-1);
    } else if ( s ==  N-1){
            add_link(N-1,N-L);
            add_link(N-1,N-2);            
            add_link(N-1,L-1);
            add_link(N-1,N-L-1);            
            add_link(N-1,0);
            add_link(N-1,L-2);            
            add_link(N-1,N-2*L);
            add_link(N-1,N-L-2);            
    } else {
            add_link(s,s+1);
            add_link(s,s-1);
            add_link(s,(s+L)%N);
            add_link(s,(s+N-L)%N);
            add_link(s,(s+L-1)%N);
            add_link(s,(s+N-L+1)%N);
            if(s%L == L-1)
                add_link(s,(s-2*L+1)%N);
            else
                add_link(s,(s+L+1)%N);
            if(s%L == 0)
                add_link(s,(s+2*L-1)%N);
            else
                add_link(s,(s+N-L-1)%N);
       
    }

}
