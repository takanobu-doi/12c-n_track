#ifndef _NUCLEAR_HPP_
#define _NUCLEAR_HPP_

#include <string>

#define amu 931.494013

using namespace std;

class nucl
{
protected:
  int N;
  int Z;
  string name;//nucl. name
  double delta;//nucl. mass(MeV)
  double E;//nucl. energy
  double Ex;//excited energy(MeV)

public:
  nucl();
  void setdata(int n,int z,string na,double del);
  void setdata(int n,int z,string na,double del,double e);
  void set_N(int n);
  void set_Z(int z);
  void set_name(string na);
  void set_delta(double del);
  void set_E(double K);
  void set_Ex(double ex);
  int get_N();
  int get_Z();
  int get_A();
  string get_name();
  double get_delta();
  double get_mass();
  double get_E();
  double get_K();
  double get_Ex();
};

#endif
