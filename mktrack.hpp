#ifndef _MKTRACK_HPP_
#define _MKTRACK_HPP_

#include "gen_eve.hpp"

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include <TGraph.h>
#include <TSpline.h>
#include <TRandom3.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "ComponentAnalyticField.hh"
#include "GeometrySimple.hh"
#include "SolidBox.hh"
#include "Sensor.hh"
#include "AvalancheMC.hh"
#include "TrackSrim.hh"
using namespace Garfield;

class mktrack
{
public:
  mktrack(int, int, int);
  ~mktrack();
  void Initialize(int, int);
  int SetGasFile();
  int SetSrimFile();
  int SetWaveFile();
  int SetRangeFile();
  int SetDriftFile();
  void SetParameters(int, int, int);
//  std::vector<std::string> GetParticleName();
  int ModTrack();
  double GetFlush(int, int, int);
  int GetTOT(int, int, int);
  int DefineDetector();
  int Generate(int&, double&);
  int GenTrack(TrackSrim*, TLorentzVector, double[3], double[3][2], int);
  void ClearBuffer();
  void AddRawWave(double [], double, int);
  void ShowSrim();
  int ShowIdealValues(std::ostream&, int);
  int ShowTeacherValues(std::ostream&, int);
  void ElectronDrift(double, double, double, double, double&, double&, double&, double&);
  double GetDriftv();
private:
  MediumMagboltz *gas = nullptr;
  GeometrySimple *geo = nullptr;
  SolidBox *box = nullptr;
  ComponentAnalyticField *comp = nullptr;
  Sensor *sensor = nullptr;
  AvalancheMC *drift = nullptr;
  TrackSrim *srim_beam = nullptr;
  std::vector<TrackSrim*> srim_particle;
  gen_eve *event = nullptr;
  TGraph *wave_temp = nullptr;
  TSpline5 *wave = nullptr;
  std::vector<TGraph*> EnetoRange_temp;
  std::vector<TSpline5*> EnetoRange;
  TRandom3 *rndm = nullptr;
  // general parameters
  std::string beam_name;
  std::vector<std::string> target_name;
  std::vector<std::vector<std::vector<std::string>>> particle_name;
  std::vector<std::vector<std::vector<double>>> particle_ex;
  std::vector<std::vector<bool>> particle_flag;
  std::string srim_name;
  std::string dirname;
  int event_id;
  int pressure;
  double beam_energy;
  double Ex_min;
  double Ex_max;
  double BEAM_RADIUS;
  double BEAM_X_CENTER;
  double BEAM_Y_CENTER;
  double VTX_X_MEAN;
  double VTX_X_SIGMA;
  double VTX_Y_MEAN;
  double VTX_Y_SIGMA;
  double VTX_Y_START;
  double VTX_Y_STOP;
  double VTX_Z_START;
  double VTX_Z_STOP;
  // gas parameters
  double W_Val;
  double Fano_Factor;
  double Mass_Gas;
  double Charge_Gas;
  double density;
  int Cluster_Size;
  // detector parameters
  double center[3];
  double half[3];
  double y_plate;
  double v_plate;
  double y_grid;
  double v_grid;
  double E_FIELD;
  double vtx[3];
  double start_point[3];
  double area[3][2];
  double beam_area[3][2];
  
  // TPC data buffer
  std::vector<std::vector<std::vector<double>>> flush;
  // raw wave

  double gain;
  int ie_step;
  double cmTomm;
  double mmTocm;
  double threshold[2];

  std::vector<std::vector<std::vector<double>>> point;
  std::vector<double> range;

  //*** dirft paraeters ***//
  double driftv;
  double diff_long;
  double diff_tra;

  // triger modify
  int buff;
};

#endif
