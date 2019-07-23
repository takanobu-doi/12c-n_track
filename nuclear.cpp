#include "nuclear.hpp"
#include <iostream>
#include <string>

using namespace std;

nucl::nucl()
{
  N = 0;
  Z = 0;
  name = "";
  delta = 0;
  E = 0;
  Ex = 0;
}

void nucl::setdata(int n,int z,string na,double del)
{
  set_N(n);
  set_Z(z);
  set_name(na);
  set_delta(del);
  set_E(0);
}

void nucl::setdata(int n,int z,string na,double del,double K)
{
  set_N(n);
  set_Z(z);
  set_name(na);
  set_delta(del);
  set_E(K);
}

void nucl::set_N(int n)
{
  N = n;
}

void nucl::set_Z(int z)
{
  Z = z;
}

void nucl::set_name(string na)
{
  name = na;
}

void nucl::set_delta(double del)
{
  delta = del;
}

void nucl::set_E(double K)
{
  E = K+get_mass();
}

void nucl::set_Ex(double ex)
{
  Ex = ex;
}

int nucl::get_N()
{
  return N;
}

int nucl::get_Z()
{
  return Z;
}

int nucl::get_A()
{
  return N+Z;
}

string nucl::get_name()
{
  return name;
}

double nucl::get_delta()
{
  return delta;
}

double nucl::get_mass()
{
  return (N+Z)*amu+delta+get_Ex();
}

double nucl::get_E()
{
  return E;
}

double nucl::get_K()
{
  return E-get_mass();
}

double nucl::get_Ex()
{
  return Ex;
}
