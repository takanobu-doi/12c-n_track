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
//  beam = TLorentzVector(0., 0., P_beam, E_beam); // assume z-axis direction
//  beam = TLorentzVector(0.,
//			-P_beam*TMath::Sin(30*TMath::DegToRad()),
//			P_beam*TMath::Cos(30*TMath::DegToRad()),
//			E_beam);
  beam = TLorentzVector(0.,
			P_beam*TMath::Sin(30*TMath::DegToRad()),
			P_beam*TMath::Cos(30*TMath::DegToRad()),
			E_beam);
  target = TLorentzVector(0., 0., 0., mass_target);
  W = beam+target;

  Ex = PARTICLE_EX;

//  event = new TGenPhaseSpace();
  rndm = new TRandom3();
  rndm->SetSeed(seed_gen());
}

gen_eve::~gen_eve()
{
//  delete event;
  delete rndm;
}

double  gen_eve::Generate()
{
  double weight;
  double uniform_rndm;
  double weight_max = 100;

  particles.clear();

  TLorentzVector w = W;
  w.Boost(0, 0, -W.Beta());
  double Ea, Pa;
  Ea = (w.E()*w.E()+(mass[0][0]+Ex[0][0]/1000)*(mass[0][0]+Ex[0][0]/1000)
	-(mass[0][1]+Ex[0][1]/1000)*(mass[0][1]+Ex[0][1]/1000))/(2*w.E());
  Pa = TMath::Sqrt(Ea*Ea-(mass[0][0]+Ex[0][0]/1000)*(mass[0][0]+Ex[0][0]/1000));
  TLorentzVector Eject(0, 0, Pa, Ea);
  Ea = (w.E()*w.E()-(mass[0][0]+Ex[0][0]/1000)*(mass[0][0]+Ex[0][0]/1000)
	+(mass[0][1]+Ex[0][1]/1000)*(mass[0][1]+Ex[0][1]/1000))/(2*w.E());
  Pa = TMath::Sqrt(Ea*Ea-(mass[0][1]+Ex[0][1]/1000)*(mass[0][1]+Ex[0][1]/1000));
  TLorentzVector Recoil(0, 0, -Pa, Ea);

  double Theta = TMath::ACos(rndm->Uniform(-1, 1));
  double Phi = rndm->Uniform(0, 2*TMath::Pi());
  Eject.RotateX(Theta);
  Eject.RotateZ(Phi);
  Eject.Boost(0, 0, W.Beta()); 
  particles.push_back(Eject); 
  Recoil.RotateX(Theta);
  Recoil.RotateZ(Phi);
  Recoil.Boost(0, 0, W.Beta());
  if(mass.size()==1){
    particles.push_back(Recoil);
  }else{
    w = Recoil;
  }
  
  for(unsigned int jj=1;jj!=mass.size();++jj){
    std::vector<double> Mass;
    TGenPhaseSpace event;
    for(unsigned int ii=0;ii!=mass[jj].size();++ii){
      double m = mass[jj][ii];
      double ex = Ex[jj][ii];
      m = m+ex/1000.;
      Mass.push_back(m);
    }
    event.SetDecay(w, Mass.size(), Mass.data());
    do{
      weight = event.Generate();
      uniform_rndm = rndm->Uniform(0., weight_max);
    }while(uniform_rndm > weight);
    
//    if(particles.size()){
//      particles.pop_back();
//    }
//    for(unsigned int ii=0;ii<Mass.size();ii++){
//      particles.push_back(*event->GetDecay(ii));
//    }
    for(unsigned int ii=0;ii!=Mass.size()-1;++ii){
//      TLorentzVector vec = *event.GetDecay(ii);
//      particles.push_back(vec);
      particles.push_back(*event.GetDecay(ii));
    }
    if(jj!=mass.size()-1){
      w = *(event.GetDecay(Mass.size()-1));
    }else{
      particles.push_back(*event.GetDecay(Mass.size()-1));
    }
  }
  return (Recoil.M()-mass[0][1])*1000; // this number has no mean.
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
