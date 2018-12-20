

#ifndef SRC_MYSTRUCTS_HPP_
#define SRC_MYSTRUCTS_HPP_

#include <vector>
#include <bitset>
#include <random>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <omp.h>

using namespace std;
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
  //double xcm=0,ycm=0;
  long m=0;
  long double inertia=0;
  //long double inertia=0;
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
struct PitStruct{
	int i;
	int j;
	int length; //i^2 + j^2
	PitStruct(){};
	PitStruct(int x, int y):i(x) ,j(y) ,length(x*x + y*y){};
};

struct MyCompareStruct{
	bool operator()(PitStruct p, double l) const
		{
			return p.length<l;
		}
};


#endif
