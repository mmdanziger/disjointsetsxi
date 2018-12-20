#include "DSetXi.hpp"


class DSetXiNetSites8 : public DSetXi {
private:
    vector<int> node_order;
    vector<int> added_yet;
    bool stop_at_giant = false;


public:
    DSetXiNetSites8(int Linput, int seed);
    void link_targets(int s);
    void add_link(int s, int t);
    void incrementalComponents(double to_p, int node_res);
    void writeDataToFile(int TestNum, string where_to);
    void runSergeyCode(int node_count, double & sergey_MI, double & sergey_M2);
    const vector<pair<double,double>>& get_xi_pair_vec() const{return xi_pair_vec;};

};
DSetXiNetSites8::DSetXiNetSites8(int Linput, int seed) : DSetXi(Linput)
    {

        node_order.resize(N);
        added_yet.resize(N);
        for (int i =0; i<N; i++){
            node_order[i] = i;
            added_yet[i] = 0;
        }
        mt19937 gen;
        gen.seed(seed);
        std::shuffle(node_order.begin(), node_order.end(), gen);
        //std::random_shuffle(node_order.begin(), node_order.end());
    }
void DSetXiNetSites8::incrementalComponents(double to_p, int node_res)
{
    long to_node = to_p * N;
    //long to_node =  0.4073 * N;
    int s,node_count;
    node_count=0;
    while(node_count < to_node){
        s = node_order[node_count];
        added_yet[s]=1;
        sumM2++;
        alive[s] = 1;
        link_targets(s);
        double current_p = static_cast<double>(node_count)/static_cast<double>(N);
        if(node_res < 2 || node_count%node_res==0){
		double sergey_MI,sergey_M2;
		//double currentXi = get_xi();
        	double currentXi = get_xi(current_p);
		/*runSergeyCode(node_count,sergey_MI,sergey_M2);
		cout<<"p= "<<current_p<<" MI= "<<sumMI<< " (" << sergey_MI <<")\t"
		  <<"M2= "<<sumM2<< " (" << sergey_M2 <<")\n";
		*/
	    
	      //xi_pair_vec.push_back(make_pair(current_p,currentXi));
	      xi_pair_vec.push_back(make_pair(current_p,currentXi));
	      m2_pair_vec.push_back(make_pair(current_p,sumM2));
	      largest_pair_vec.push_back(make_pair(current_p,clusters[largest_id].m));
             //   pinf_pair_vec.push_back(make_pair(current_p,get_Pinf()));
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
/*    if(node_count % L == 0){
    	cout<<"node_count="<<node_count<<endl;
    	check_data();
    }*/

    }
}
void DSetXiNetSites8::add_link(int s, int t)
{
    if(added_yet[s]==1 && added_yet[t]==1 && !(found_giant && stop_at_giant))
        merge_components(s,t);


}
void DSetXiNetSites8::link_targets(int s)
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
void DSetXiNetSites8::writeDataToFile(int TestNum, string where_to)
{
	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = where_to + "dataxi";
	ss<<TestNum;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<xi_pair_vec.size();i++)
		dataInfo<<xi_pair_vec[i].first<<"	"<<xi_pair_vec[i].second<<"	"<<m2_pair_vec[i].second<<"	"<<largest_pair_vec[i].second<<endl;
	dataInfo.close(); // p, mi, m2
}

void DSetXiNetSites8::runSergeyCode(int node_count, double & sergey_MI, double & sergey_M2)
{
  int *burn,*burn1; 
  int *disc,*x1_shell,*y1_shell,*x_shell,*y_shell,*dummy;
  int nburn,nburn1; 
  int top;
  int cln; //cluster number
  double cls,clnmax,clsmax,R2,R2max;
  double s1x,s1y,s2x,s2y;
  double scn,sR2,sn2;
  int i,j,k;
  int x,y,x0,y0,x1; 
  long L2;//table's side length and area
  unsigned int rn; //rundom number
  double res;
  char Test_num[3];
  

  //strcpy(Test_num,argv[3]);
  //char fname[100]="/tmp/sergey_test"; //file name 
  //strcat(fname,Test_num);
  //strcat(fname, ".txt");

  
  //ff=fopen(fname,"w");

  L2=L*L; //number of cells in the lattice
  
  disc=(int *)malloc(L2*sizeof(int));
  burn=(int *)malloc(L2*sizeof(int));
  burn1=(int *)malloc(L2*sizeof(int));
  x1_shell=(int *)malloc(L2*sizeof(int));
  y1_shell=(int *)malloc(L2*sizeof(int));
  x_shell=(int *)malloc(L2*sizeof(int));
  y_shell=(int *)malloc(L2*sizeof(int));
  sR2=0; // sum of R^2
  sn2=0; // sum of cluster size squared
  scn=0; // sum center of nodes?
  
  for(k=0;k<L2;k++)
    disc[k]=L2;
  for(k=0;k<=node_count;k++) 
	disc[node_order[k]]=-1;
  
      top=0;//disc[top] always contain node we care about, after we finish handle a node, 'top' goes for the next node.
      cln=0;// cluster number
      clsmax=0;//all the __max paramter just keep the max value of __, the algorithem don't actually use them.
      clnmax=0;
      while(1){
	while(disc[top]>=0)//exit condition: have we went through all the nodes we needed?
	  {
	    top++;
	    if(top==L2)goto finish;
	  }
	disc[top]=cln; 
        burn[0]=top; 
	x_shell[0]=0;
	y_shell[0]=0;
	nburn=1; //run variable of burn
	cls=0; // cluster size
	s1x=0; //source (x0,y0) and squared Coordinates
	s1y=0;
	s2x=0;
	s2y=0;
	while(nburn)
	  {
	    nburn1=0;//run variable of burn1
	    while(nburn)
	      {
		nburn--;
		k=burn[nburn]; //source node(x0,y0) , target node(x,y)
		x0=x_shell[nburn];
		y0=y_shell[nburn];
		cls++;
		s1x+=x0;
		s2x+=x0*x0;
		s1y+=y0;
		s2y+=y0*y0;
		//defines the cell coordiantes from the address
		x1=k%L;
		for(i=0;i<8;i++)//distributes 4 grains to the four neibors
		  {
		    x=x0;
		    y=y0;
		  switch(i) //each neibor must be treated individually since
		      //it may be in the vicinity of different table's edges
		      {//4 nearest Neighbours
		      case 0: //Eastern neigbor
			{
			  j=k+1;
			  x++;
			  if(x1+1==L)
			    j-=L;
			  break;
			}
		      case 1: //Western neibor
			{
			  j=k-1;
			  x--;
			  if(x1==0)
			    j+=L;//j is out of the table
			  break;
			}
		      case 2://Northern neibor
			{
			  j=k+L;
			  y++;
			  if(j>=L2)
			    j-=L2;
			  break;
			}
		      case 3: //Southern neighbor
			{
			  j=k-L;
			  y--;
			if(j<0)
			  j+=L2;
			break;
			}//4 secend nearest Neighbours
		      case 4: //Eastern neigbor (E+N)
			{
			  x++;
                          y++;
			  j=k+1;
			  if(x1+1==L)
			    j-=L;
			  j+=L;
			  if(j>=L2)
			    j-=L2;
			  break;
			}
		      case 5: //Western neibor ((W+N)
			{
			  x--;
			  y++;
			  j=k-1;
			  if(x1==0)
			    j+=L;//j is out of the table
			  j+=L;
			  if(j>=L2)
			    j-=L2;
			  break;
			}
		      case 6://Northern neibor (i dont think the comment is correct, it should be E+S)
			{
			  x++;
			  y--;
			  j=k+1;
			  if(x1+1==L)
			    j-=L;
			  j-=L;
			  if(j<0)
			    j+=L2;
			  break;
			}
		      case 7: //Southern neighbor (S+W)
			{
			  x--;
			  y--;
			  j=k-1;
			  if(x1==0)
			    j+=L;
			  j-=L;
			  if(j<0)
			    j+=L2;	
			break;
			}

		      }
		    if(disc[j]<0)// disc[j]==-1 it means that the node at j is node we care about and not connected to cluster.method
		      {
			disc[j]=cln;//connect it to the cluster
			burn1[nburn1]=j;
			x1_shell[nburn1]=x;
			y1_shell[nburn1]=y;
			nburn1++;
		      }
		  }
	      }
	    nburn=nburn1;
	    dummy=x1_shell;x1_shell=x_shell;x_shell=dummy;//lots of swap
	    dummy=y1_shell;y1_shell=y_shell;y_shell=dummy;
	    dummy=burn1;burn1=burn;burn=dummy;
	  }
	R2=cls*(s2y+s2x)-s1x*s1x-s1y*s1y; 
	//printf("%lf %lf\n",R2,cls);
	sR2+=R2;
	sn2+=cls*cls;
	if(cls>clsmax)
	  {
	    clsmax=cls;
	    clnmax=cln;
	    R2max=R2;
	  }
	cln++;
      }
    finish:
      //      sR2-=R2max;
      //      sn2-=clsmax*clsmax;
      //      scn+=cln-1;
      scn+=cln;
      sergey_MI = sR2;
      sergey_M2 = sn2;
    }

  //if(sn2 == 0 || scn == 0)
  //	fprintf(ff,"%lf %lf %lf %lf %lf\n",p,0.0,0.0,scn,clsmax); 
  //else
  //  	fprintf(ff,"%lf %lf %lf %lf %lf\n",p,sR2/sn2,sn2/scn,scn,clsmax); 
  
