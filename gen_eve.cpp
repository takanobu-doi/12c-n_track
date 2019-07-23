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

gen_eve::gen_eve(std::string BEAM_NAME, std::string TARGET_NAME, std::vector<std::string> PARTICLE_NAME,
		 double T_beam/*[MeV]*/)
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
  
  for(auto it=PARTICLE_NAME.begin();it!=PARTICLE_NAME.end();++it){
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
    mass.push_back(dataset.get_mass(A, particle_name.c_str())/1000.);
  }
  
  E_beam = T_beam/1000.+mass_beam; // total energy of incident particle [GeV]
  P_beam = TMath::Sqrt(E_beam*E_beam-mass_beam*mass_beam); // P [GeV/c]
  beam = TLorentzVector(0., 0., -P_beam, E_beam); // assume z-axis direction
  target = TLorentzVector(0., 0., 0., mass_target);
  W = beam+target;

  event = new TGenPhaseSpace();
  rndm = new TRandom3();
  rndm->SetSeed(seed_gen());
}

gen_eve::~gen_eve()
{
  delete event;
  delete rndm;
}

void gen_eve::Generate()
{
  double weight;
  double uniform_rndm;
  double weight_max = 10;

  particles.clear();

  event->SetDecay(W, mass.size(), mass.data());
  do{
    weight = event->Generate();
    uniform_rndm = rndm->Uniform(0., weight_max);
//    uniform_rndm = rndm->Uniform(0., event->GetWtMax());
  }while(uniform_rndm > weight);

  for(unsigned int ii=0;ii<mass.size();ii++){
    particles.push_back(*event->GetDecay(ii));
  }

  return;
}

double gen_eve::GetBeamMass()
{
  return mass_beam;
}

double gen_eve::GetParticleMass(int i)
{
  return mass[i];
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
