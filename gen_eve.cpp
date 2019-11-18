#include "gen_eve.hpp"
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <random>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include "database.hpp"

gen_eve::gen_eve(std::string BEAM_NAME, std::string TARGET_NAME, std::vector<std::vector<std::string>> PARTICLE_NAME,
		 double T_beam/*[MeV]*/, std::vector<std::vector<double>> PARTICLE_EX/*[MeV]*/ )
{
  // get paticle imformation
  database dataset;
  int A;
  std::string particle_name = BEAM_NAME;
  if(particle_name.size()>1){
    A = stoi(particle_name);
    particle_name.erase(particle_name.begin(), particle_name.begin()+int(log10(A))+1);
  }else{
    if(particle_name=="p"){
      A = 1;
    }else if(particle_name=="d"){
      A = 2;
    }else if(particle_name=="t"){
      A = 3;
    }else if(particle_name=="n"){
      A = 1;
    }else{
      A = 0;
    }
  }
  mass_beam = dataset.get_mass(A, particle_name.c_str())/1000.;

  particle_name = TARGET_NAME;
  if(particle_name.size()>1){
    A = stoi(particle_name);
    particle_name.erase(particle_name.begin(), particle_name.begin()+int(log10(A))+1);
  }else{
    if(particle_name=="p"){
      A = 1;
    }else if(particle_name=="d"){
      A = 2;
    }else if(particle_name=="t"){
      A = 3;
    }else if(particle_name=="n"){
      A = 1;
    }else{
      A = 0;
    }
  }
  mass_target = dataset.get_mass(A, particle_name.c_str())/1000.;
  
  for(auto it1=PARTICLE_NAME.begin();it1!=PARTICLE_NAME.end();++it1){
    std::vector<double> mass_temp;
    for(auto it=(*it1).begin();it!=(*it1).end();++it){
      particle_name = *it;
      if(particle_name.size()>1){
	A = stoi(particle_name);
	particle_name.erase(particle_name.begin(), particle_name.begin()+int(log10(A))+1);
      }else{
	if(particle_name=="p"){
	  A = 1;
	}else if(particle_name=="d"){
	  A = 2;
	}else if(particle_name=="t"){
	  A = 3;
	}else if(particle_name=="n"){
	  A = 1;
	}else{
	  A = 0;
	}
      }
      mass_temp.push_back(dataset.get_mass(A, particle_name.c_str())/1000.);
    }
    mass.push_back(mass_temp);
  }
  E_beam = T_beam/1000.+mass_beam; // total energy of incident particle [GeV]
  P_beam = TMath::Sqrt(E_beam*E_beam-mass_beam*mass_beam); // P [GeV/c]
  beam = TLorentzVector(0., 0., -P_beam, E_beam); // assume z-axis direction
  target = TLorentzVector(0., 0., 0., mass_target);
  W = beam+target;

  Ex = PARTICLE_EX;

  event = new TGenPhaseSpace();
  rndm = new TRandom3();
  rndm->SetSeed(seed_gen());
}

gen_eve::~gen_eve()
{
  delete event;
  delete rndm;
}

double  gen_eve::Generate()
{
  double weight;
  double uniform_rndm;
  double weight_max = 10;

  particles.clear();

  TLorentzVector w = W;
  for(unsigned int jj=0;jj!=mass.size();++jj){
    std::vector<double> Mass;
    for(unsigned int ii=0;ii!=mass[jj].size();++ii){
      double m = mass[jj][ii];
      double ex = Ex[jj][ii];
      m = m+ex/1000.;
      Mass.push_back(m);
    }
    event->SetDecay(w, Mass.size(), Mass.data());
    do{
      weight = event->Generate();
      uniform_rndm = rndm->Uniform(0., weight_max);
//    uniform_rndm = rndm->Uniform(0., event->GetWtMax());
    }while(uniform_rndm > weight);

    if(particles.size()){
      particles.pop_back();
    }
    for(unsigned int ii=0;ii<Mass.size();ii++){
      particles.push_back(*event->GetDecay(ii));
    }
    w = *(particles.end()-1);
  }
  return 7.65; // this number has no mean.
}

double gen_eve::GetBeamMass()
{
  return mass_beam;
}

double gen_eve::GetParticleMass(int i, int j)
{
  return mass[i][j];
}

TLorentzVector gen_eve::GetBeamVector()
{
  return beam;
}

TLorentzVector gen_eve::GetParticleVector(int i)
{
  return particles[i];
}

unsigned int gen_eve::GetParticleNumber()
{
  return particles.size();
}
