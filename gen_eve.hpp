#ifndef _GEN_EVE_HPP_
#define _GEN_EVE_HPP_

#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <random>
#include <string>
#include <vector>

class gen_eve
{
public:
  gen_eve(std::string, std::string, std::vector<std::vector<std::string>>, double, std::vector<std::vector<double>>);
  ~gen_eve();
  double Generate();
  double GetBeamMass();
  double GetParticleMass(int, int);
  TLorentzVector GetBeamVector();
  TLorentzVector GetParticleVector(int);
  unsigned int GetParticleNumber();
private:
  double E_beam;
  double P_beam;
  double mass_beam;
  double mass_target;
  std::vector<std::vector<double>> mass;
  std::vector<std::vector<double>> Ex;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector W;
  std::vector<TLorentzVector> particles;

  TGenPhaseSpace *event;

  std::random_device seed_gen;
  TRandom3 *rndm;
};
#endif
