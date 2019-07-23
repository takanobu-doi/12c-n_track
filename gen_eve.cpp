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
		 std::vector<std::vector<double>> PARTICLE_EX, double T_beam/*[MeV]*/)
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

  int ii = 0;
  for(auto it1=PARTICLE_NAME.begin();it1!=PARTICLE_NAME.end();++it1){
    std::vector<double> mass_temp;
    int jj = 0;
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
      mass_temp.push_back(dataset.get_mass(A, particle_name.c_str())/1000.+PARTICLE_EX[ii][jj]/1000.);
      jj++;
    }
    mass.push_back(mass_temp);
    ii++;
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

  event->SetDecay(W, mass[0].size(), mass[0].data());
  do{
    weight = event->Generate();
    uniform_rndm = rndm->Uniform(0., weight_max);
//    uniform_rndm = rndm->Uniform(0., event->GetWtMax());
  }while(uniform_rndm > weight);

  for(unsigned int ii=0;ii<mass[0].size();ii++){
    particles.push_back(*event->GetDecay(ii));
  }

  for(unsigned int ii=1;ii<mass.size();ii++){
    event->SetDecay(*(particles.end()-1), mass[ii].size(), mass[ii].data());
    particles.pop_back();
    do{
      weight = event->Generate();
      uniform_rndm = rndm->Uniform(0., weight_max);
    }while(uniform_rndm > weight);

    for(unsigned int jj=0;jj<mass[ii].size();jj++){
      particles.push_back(*event->GetDecay(jj));
    }
  }

  return;
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
