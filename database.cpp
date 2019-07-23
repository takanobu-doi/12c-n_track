#include <iostream>
#include <string>
#include <cstdlib>
#include <typeinfo>
#include <fstream>
#include <algorithm>
#include "database.hpp"
#include "nuclear.hpp"

using namespace std;

database::database()
{
  inputdata();
  cout  << part_num << " particles are registed." << endl;
  AMU = 931.495; // [MeV]
}

void database::setdata(int n,int z,string na,double del)
{
  if(part_num >= PARTICLE_NUM){
    cout << "Dataset is too many." << endl;
    cout << "You should more space for data." << endl;
    exit(EXIT_FAILURE);
  }
  data[part_num].setdata(n,z,na,del);
  part_num++;
}

int database::det_part(int N,int Z)
{
  for(int i=0;i<part_num;i++){
    if(data[i].get_N()==N&&data[i].get_Z()==Z){
      return i;
    }
  }
  cout << "Such nuclear does not registed in this database." << endl;
  exit(EXIT_FAILURE);
}

int database::det_part(int A,string ss)
{
  for(int i=0;i<part_num;i++){
    if(data[i].get_N()+data[i].get_Z()==A
       &&data[i].get_name()==ss){
      return i;
    }
  }
  cout <<"Such nuclear does not registed in this database." << endl;
  exit(EXIT_FAILURE); 
}

int database::read_stdin()
{
  cout << ":";
  if(getline(cin,name)){
    
    A = (int)strtol(name.c_str(),NULL,10);
    
    if(A==0){
      if(name=="n"){
	A = 1;
      }else if(name=="p"){
	A = 1;
      }else if(name=="d"){
	A = 2;
      }else if(name=="t"){
	A = 3;
      }else if(name=="r"){
	A = 0;
	}else{
	cout << "Illigal input format." << endl;
	return 1;
      }
    }
    else{
      for(double keta=1;keta<=A;keta=keta*10){
	name.erase(name.begin());
      }
      transform(name.begin(),name.begin()+1,name.begin(),::toupper);
    }
  }else{
    cout << "Illigal forlmat." << endl;
    return 1;
  }
  return 0;
}

int database::read_comd(char *n)
{
  A = (int)strtol(n,NULL,10);
  name = n;
  
  if(A==0){
    if(name=="n"){
      A = 1;
    }else if(name=="p"){
      A = 1;
    }else if(name=="d"){
      A = 2;
    }else if(name=="t"){
      A = 3;
    }else if(name=="r"){
      A = 0;
    }else{
      cout << "Illigal input format." << endl;
      return 1;
    }
  }
  else{
    for(double keta=1;keta<=A;keta=keta*10){
      name.erase(name.begin());
    }
    transform(name.begin(),name.begin()+1,name.begin(),::toupper);
  }
  return 0;
}

double database::get_delta(int N,int Z)
{
  return data[det_part(N,Z)].get_delta();
}

double database::get_delta(int A,string ss)
{
  return data[det_part(A,ss)].get_delta();
}

double database::get_mass(int N, int Z)
{
  return (N+Z)*AMU+data[det_part(N,Z)].get_delta();
}

double database::get_mass(int A, string ss)
{
  return A*AMU+data[det_part(A, ss)].get_delta();
}

string database::get_name(int N,int Z){
  return data[det_part(N,Z)].get_name();
}

int database::get_N()
{
  return N;
}

int database::get_N(int A,string ss)
{
  return data[det_part(A,ss)].get_N();
}

int database::get_Z()
{
  return Z;
}

int database::get_Z(int A,string ss)
{
  return data[det_part(A,ss)].get_Z();
}

int database::get_A()
{
  return A;
}

string database::get_name()
{
  return name;
}
