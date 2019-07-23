#ifndef _DATABASE_HPP_
#define _DATABASE_HPP_

#include "nuclear.hpp"

#define PARTICLE_NUM 350//登録できる粒子数の上限

using namespace std;

class database
{
private:
  nucl data[PARTICLE_NUM];
  int part_num;
  int A;
  int Z;
  int N;
  string name;
  double delta;
  double AMU;
  void inputdata();
  void setdata(int n,int z,string na,double del);
  int det_part(int N,int Z);//determine the nucl. from N,Z
  int det_part(int A,string ss);//determin the nucl. from name,A
  
public:
  database();
  double get_delta(int N,int Z);//get delta of nucl. from N,Z
  double get_delta(int A,string ss);//get delta of nucl. from name,A
  double get_mass(int N, int Z);
  double get_mass(int A, string ss);
  string get_name(int N,int Z);//get name of nucl.
  int read_stdin();//if faild return 1,success return 0
  int read_comd(char*);//if faild return 1,success return 0
  int get_N();
  int get_N(int A,string ss);
  int get_Z();
  int get_Z(int A,string ss);
  int get_A();
  string get_name();
};

#endif
