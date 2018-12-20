#include "dsetxi.hpp"

class DSetXiBonds : public DSetXi {
protected:
    vector<pair<int,int>> edge_list;

public:
    DSetXiBonds(int L);
    void make_NNN_lattice();
    void add_link_to_edge_list(int s, int t);
    void incrementalComponents(double toK, int link_res);
};

DSetXiBonds::DSetXiBonds(int L) :DSetXi(L)
    {
        sumMI=0;
        sumM2=N;
        nodes.resize(N);
        for (int i =0; i<N; i++){
            nodes[i] = node(i,L);
            clusters[i] = cluster(i%L,i/L,i);
        }
        make_NNN_lattice();
        std::random_shuffle(edge_list.begin(), edge_list.end());
    }
    
void DSetXiBonds::make_NNN_lattice()
{
    for(int i =0; i<N; i++){
       // n.neighbor ---- 
       if (i%L == L-1){
           add_link_to_edge_list(i,i-L+1);
       } else {
           add_link_to_edge_list(i,i+1);
       }
       // n.neighbor   |
       //              |
       add_link_to_edge_list(i,(i+L)%N);//the modulo corrects the bottom row automatically
       // n.n.neighbor /
       //             /
       if (i/L == 0){//first row, need to go up
           if(i == L-1){
               add_link_to_edge_list(L-1,N-L);//only one point like that
           }else{
               add_link_to_edge_list(i,N-L+i+1);//in last row (N-L) i+1 steps from the edge
           }
       } else{
           add_link_to_edge_list(i,i-L+1);
       }
       /* n.n. neighbor \
                         \         */
       add_link_to_edge_list(i,(i+1)%N); 
    }
}

void DSetXiBonds::add_link_to_edge_list(int s, int t){
    edge_list.push_back(make_pair(s,t));
}

void DSetXiBonds::incrementalComponents(double toK, int link_res)
{
    long to_link = toK*N/2;
    int s,t,link_count;
    link_count=0;
    while(link_count < to_link){
        s = edge_list[link_count].first;
        t = edge_list[link_count].second;
        merge_components(s,t);
        link_count++;
        if(found_giant || link_res <2 || link_count%link_res==0){
            double currentK = static_cast<double>(link_count)*2/static_cast<double>(N);
            double currentXi = get_xi();
            xi_pair_vec.emplace_back(make_pair(currentK, currentXi));
            if (currentXi > L || (currentXi < 1e-5) ){
                std::cerr<<"Bad Xi encountered, check output.\n";
                break;
            }
        }
        if(found_giant)
            break;
    }
}

