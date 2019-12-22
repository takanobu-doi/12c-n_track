#include "mktrack.hpp"
#include "gen_eve.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <random>

#include <TSpline.h>
#include <TGraph.h>
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

//#define ONE_ELE 1

void mktrack::SetParameters(int Event_id, int Pressure, int Gas)
{
  // general parameters
  beam_name = "4He";
  // define event-id
  target_name = {"n"};
  particle_name = {{{"n", "4He"}}
  }; // last particle each bracket is not used for tracking, but all particle of the last bracket is used.
  particle_ex = {{{0, 0}}
  };
  particle_flag = {{false, false}
  };
  dirname = "table/";
  event_id = Event_id;
  pressure = Pressure;
//  beam_energy = 14; // [MeV]
  beam_energy = 0.5; // [MeV]
//  beam_energy = 4.2; // [MeV]
  Ex_min = 0.; // [MeV]
  Ex_max = 20.; // [MeV]
  BEAM_RADIUS = 10;         // mm
  BEAM_X_CENTER = 102.4/2.; // mm
  BEAM_Y_CENTER = 140.*0.2;  // mm
  VTX_X_MEAN = 102.4/2.;    // mm
  VTX_X_SIGMA = 0.1;	    // mm
  VTX_Y_MEAN = 140./2.;	    // mm
  VTX_Y_SIGMA = 0.1;	    // mm
  VTX_Y_START = 140.*1/8.;  // mm
  VTX_Y_STOP = 140.*7/8;    // mm
  VTX_Z_START = 102.4*0/8.; // mm
  VTX_Z_STOP = 102.4*9/8;   // mm       

  // detector parameters
  center[0] = 10.24/2;
  center[1] = 14./2;
  center[2] = 10.24/2;
  half[0] = 15.;
  half[1] = 14./2;
  half[2] = 14./2;
  y_plate = 14.;
  y_grid = 0.;
  v_grid = -1250.;
  area[0][0] = 0.;
  area[1][0] = 0.;
  area[2][0] = 0.;
  area[0][1] = 102.4;
  area[1][1] = 140.;
  area[2][1] = 102.4;
  beam_area[0][0] = 0.;
  beam_area[1][0] = 0.;
  beam_area[2][0] = 0.;
  beam_area[0][1] = 102.4;
  beam_area[1][1] = 140.;

  // gas parameters
  switch(Gas){
  case 0: // detection gas is CH4
    srim_name = "CH4_";
    W_Val = 30.0; // for CH4
    Mass_Gas = 16.;   // for CH4 
    Charge_Gas = 10.; // for CH4
    density = 0.00065819*pressure/1000.; // for CH4
    gain = 360.;
    if(Pressure==50){
      v_plate = -70+v_grid; // for CH4 50hPa
    }else if(Pressure==100){
      v_plate = -120.+v_grid; // for CH4 100hPa
    }else{
      exit(0);
    }
    break;
  case 1: // detection gas is CH4(3)+H2(7)
    srim_name = "CH4_H2_";
    W_Val = 30.0*0.3+13.6*0.7; // for CH4 3 H2 4
    Mass_Gas = 16.*0.3+2.*0.7;   // for CH4 3 H2 7
    Charge_Gas = 10.*0.3+2.*0.7; // for CH4 3 H2 7
    density = 0.000025535*pressure/100.; // for CH4 3 H2 7
    gain = 355;
    if(Pressure==100){
      v_plate = -605.+v_grid; // for CH4 3 H2 7 100hPa
    }else{
      exit(0);
    }    
    break;
  case 2: // detection gas is CH4(4)+He(6)
    srim_name = "CH4_He_";
    W_Val = 30.0*0.4+41.3*0.6; // for CH4 4 He 6
    Mass_Gas = 16.*0.4+4*0.6;   // for CH4 4 He 6
    Charge_Gas = 10.*0.4+2*0.6; // for CH4 4 He 6
    density = 0.000036181*pressure/100.; // for CH4 4 He 6
    gain = 320;
    if(Pressure==100){
      v_plate = -265.+v_grid; // for CH4 4 He 6 100hPa
    }else{
      exit(0);
    }
    break;
  case 3: // detection gas is iC4H10
    srim_name = "iC4H10_";
    W_Val = 26;
    Mass_Gas = 58;
    Charge_Gas = 34;
    density = 0.000035769*pressure/15.;
    gain = 0;
    if(Pressure==15){
      v_plate = -90+v_grid;
    }else{
      exit(0);
    }
    break;
  case 4: // detection gas is iC4H10(1)+H2(9)
    srim_name = "iC4H10_1_H2_9_";
    W_Val = 26*0.1+13.6*0.9;
    Mass_Gas = 58*0.1+2*0.9;
    Charge_Gas = 34*0.1+2*0.9;
    density = 0.00003129*pressure/100.;
    gain = 270;
    if(Pressure==100){
      v_plate = -955+v_grid;
    }else{
      exit(0);
    }
    break;
  case 5: // detection gas is iC4H10(1)+He(9)
    srim_name = "iC4H10_1_He_9_";
    W_Val = 26*0.1+41.3*0.9;
    Mass_Gas = 58*0.1+4*0.9;
    Charge_Gas = 34*0.1+2*0.9;
    density = 0.000038627*pressure/100.;
    gain = 390;
    if(Pressure==100){
      v_plate = -460+v_grid;
    }else{
      exit(0);
    }
    break;
  case 6:
    srim_name = "iC4H10_3_H2_7_";
    W_Val = 26*0.3+13.6*0.7;
    Mass_Gas = 58*0.3+2*0.7;
    Charge_Gas = 34*0.3+2*0.7;
    density = 0.000038664*Pressure/50;
    gain = 0;
    if(Pressure==50){
      v_plate = -330+v_grid;
    }else{
      exit(0);
    }
    break;
  case 7:
    srim_name = "iC4H10_2_H2_8_";
    W_Val = 26*0.2+13.6*0.8;
    Mass_Gas = 58*0.2+2*0.8;
    Charge_Gas = 34*0.2+2*0.8;
    density = 0.000027155*Pressure/50;
    gain = 0;
    if(Pressure==50){
      v_plate = -370+v_grid;
    }else{
      exit(0);
    }
    break;
  case 8:
    srim_name = "iC4H10_3_He_7_";
    W_Val = 26*0.3+41.3*0.7;
    Mass_Gas = 58*0.3+4*0.7;
    Charge_Gas = 34*0.3+2*0.7;
    density = 0.000041517*Pressure/50;
    gain = 500;
    if(Pressure==50){
      v_plate = -210+v_grid;
    }else{
      exit(0);
    }
    break;
  default:
    exit(0);
    break;
  }
  Fano_Factor = 1.0;
  Cluster_Size = 1; // default is 30
  Beam_Cluster_Size = 1;
  Particle_Cluster_Size = 20;
  E_FIELD = (v_grid-v_plate)/(y_plate-y_grid);

//  gain = 100.;
//  gain = 300.; // default is 1000.
  ie_step = 1; // default is 100

  cmTomm = 10.;
  mmTocm = 0.1;
  threshold[0] = 0.1; // anode
  threshold[1] = 0.1;  //cathode

  buff = 50;
  return;
}

mktrack::mktrack(int Event_id, int Pressure, int Gas)
{
  SetParameters(Event_id, Pressure, Gas);
  Initialize(Event_id, Pressure);
}

mktrack::~mktrack()
{
  if(gas!=nullptr){
    delete gas;
  }
  if(geo!=nullptr){
    delete geo;
  }
  if(box!=nullptr){
    delete box;
  }
  if(comp!=nullptr){
    delete comp;
  }
  if(sensor!=nullptr){
    delete sensor;
  }
  if(drift!=nullptr){
    delete drift;
  }
  if(srim_beam!=nullptr){
    delete srim_beam;
  }
  srim_particle.clear();
  if(event!=nullptr){
    delete event;
  }
  if(wave_temp!=nullptr){
    delete wave_temp;
  }
  if(wave!=nullptr){
    delete wave;
  }
  EnetoRange_temp.clear();
  EnetoRange.clear();
  if(rndm!=nullptr){
    delete rndm;
  }
}

void mktrack::Initialize(int Event_id, int Pressure){
  if(SetGasFile()==0){
    std::cerr << "Cannot open gasfile" << std::endl;
    exit(0);
  }
  if(DefineDetector()==0){
    std::cerr << "Cannot define detector" << std::endl;
    exit(0);
  }
  if(SetSrimFile()==0){
    std::cerr << "Cannot open srimfile" << std::endl;
    exit(0);
  }
  event = new gen_eve(beam_name, *(target_name.begin()+event_id), *(particle_name.begin()+event_id), beam_energy, *(particle_ex.begin()+event_id));
  if(SetWaveFile()==0){
    std::cerr << "Cannot open wavefile" << std::endl;
    exit(0);
  }
  if(SetRangeFile()==0){
    std::cerr << "Cannot open rangefile" << std::endl;
    exit(0);
  }
  if(SetDriftFile()==0){
    std::cerr << "Cannot set driftfile" << std::endl;
    exit(0);
  }
  
  std::random_device seed_gen;
  
  rndm = new TRandom3();
  rndm->SetSeed(seed_gen());

  std::vector<double> ch(256);
  std::vector<std::vector<double>> clk;
  for(int ii=0;ii<1024;++ii){
    clk.push_back(ch);
  }
  for(int ii=0;ii<2;++ii){
    flush.push_back(clk);
  }
}

int mktrack::SetGasFile()
{
  std::string magfname = dirname+srim_name+std::to_string(pressure)+".gas";
  gas = new MediumMagboltz();
  std::cout << "Loading gasfile: " << magfname << std::endl;
  if(gas->LoadGasFile(magfname)==0){
    return 0;
  }
  return 1;
}

int mktrack::SetSrimFile()
{
  std::string srimfname;
  if(beam_name=="n"){
    srim_beam = nullptr;
  }else{
    srimfname = dirname+beam_name+"_"+srim_name+std::to_string(pressure)+".srim";
  
    srim_beam = new TrackSrim();
    if(sensor == nullptr){
      std::cerr << "Sensor is not defined yet." << std::endl;
      return 0;
    }
    srim_beam->SetSensor(sensor);
//  std::cout << "Loading srimfile: " << srimfname << std::endl;
    if(srim_beam->ReadFile(srimfname)==0){
      return 0;
    }
    srim_beam->SetWorkFunction(W_Val);
    srim_beam->SetFanoFactor(Fano_Factor);
    srim_beam->SetModel(4);
    srim_beam->SetAtomicMassNumbers(Mass_Gas,
				    Charge_Gas);
    srim_beam->SetDensity(density);
    srim_beam->SetTargetClusterSize(Cluster_Size);
    srim_beam->DisableTransverseStraggling();
    srim_beam->DisableLongitudinalStraggling();
//    srim_beam->EnableTransverseStraggling();
//    srim_beam->EnableLongitudinalStraggling();
  }
  
  for(auto it1=(*(particle_name.begin()+event_id)).begin();it1!=(*(particle_name.begin()+event_id)).end();++it1){
    for(auto it=(*it1).begin();it!=(*it1).end();++it){
      if(it==(*it1).begin()&&srim_particle.size()){
	srim_particle.pop_back();
      }
      if(*it=="n"){
	srim_particle.push_back(nullptr);
	continue;
      }
      srimfname = dirname+(*it)+"_"+srim_name+std::to_string(pressure)+".srim";
      srim_particle.push_back(new TrackSrim());
      (*(srim_particle.end()-1))->SetSensor(sensor);
//    std::cout << "Loading srimfile: " << srimfname << std::endl;
      if((*(srim_particle.end()-1))->ReadFile(srimfname)==0){
	return 0;
      }
      (*(srim_particle.end()-1))->SetWorkFunction(W_Val);
      (*(srim_particle.end()-1))->SetFanoFactor(Fano_Factor);
      (*(srim_particle.end()-1))->SetModel(4);
      (*(srim_particle.end()-1))->SetAtomicMassNumbers(Mass_Gas,
						       Charge_Gas);
      (*(srim_particle.end()-1))->SetDensity(density);
      (*(srim_particle.end()-1))->SetTargetClusterSize(Cluster_Size);
//      (*(srim_particle.end()-1))->DisableTransverseStraggling();
//      (*(srim_particle.end()-1))->DisableLongitudinalStraggling();
      (*(srim_particle.end()-1))->EnableTransverseStraggling();
      (*(srim_particle.end()-1))->EnableLongitudinalStraggling();
    }
  }

  return 1;
}

int mktrack::SetWaveFile()
{
  std::string wavefilename = "table/wave_temp.dat";
//  std::cout << "Loading wavefile: " << wavefilename << std::endl;
  std::ifstream ifile(wavefilename);
  if(ifile.fail()){
    std::cerr << "There is not " << wavefilename << std::endl;;
    return 0;
  }
  std::vector<double> t;
  std::vector<double> x;
  std::string str;
  std::stringstream stream;  
  double t_temp;
  double x_temp;
  while(getline(ifile, str)){
    stream.clear();
    stream << str;
    stream >> t_temp >> x_temp;
    t.push_back(t_temp);
    x.push_back(x_temp);
  }
  wave_temp = new TGraph(t.size(), t.data(), x.data());
  wave_temp->SetBit(1);
  wave = new TSpline5("wave", wave_temp);
  
  return 1;
}

int mktrack::SetRangeFile()
{
  std::string rangefname;
  for(auto it1=(*(particle_name.begin()+event_id)).begin();it1!=(*(particle_name.begin()+event_id)).end();++it1){
    if(EnetoRange.size()){
      EnetoRange.pop_back();
    }
    for(auto it=(*it1).begin();it!=(*it1).end();++it){
      if((*it)=="n"){
	EnetoRange.push_back(nullptr);
	continue;
      }
      rangefname = dirname+(*it)+"_"+srim_name+std::to_string(pressure)+"_ene_to_range.dat";
      std::ifstream ifile(rangefname);
      if(ifile.fail()){
	std::cerr << "There is not " << rangefname << std::endl;
	return 0;
      }
      std::vector<double> e;
      std::vector<double> r;
      std::string str;
      double e_temp;
      double r_temp;
      while(getline(ifile, str)){
	std::stringstream stream;
	stream << str;
	stream >> e_temp >> r_temp;
	e.push_back(e_temp);
	r.push_back(r_temp);
      }
      EnetoRange.push_back(new TSpline5("EnetoRange", e.data(), r.data(), e.size()));
    }
  }
  return 1;
}

int mktrack::SetDriftFile()
{
  std::vector<double> efields, bfields, angles;
  gas->GetFieldGrid(efields, bfields, angles);
  double velocity[efields.size()];
  double diff_t[efields.size()];
  double diff_l[efields.size()];
  double town_a[efields.size()];
  double attach[efields.size()];

  unsigned int i=0;
  for(i=0;i<efields.size();++i){
    gas->GetElectronVelocityE(i, 0, 0, velocity[i]); // cm/ns
    gas->GetElectronTransverseDiffusion(i, 0, 0, diff_t[i]); // sqrt(cm)
    gas->GetElectronLongitudinalDiffusion(i, 0, 0, diff_l[i]); // sqrt(cm)
    gas->GetElectronTownsend(i, 0, 0, town_a[i]);
    gas->GetElectronAttachment(i, 0, 0, attach[i]);

    if(town_a[i]<0){
      town_a[i] = 0;
    }
    if(attach[i]<0){
      attach[i] = 0;
    }
  }

  if(i==0){
    return 0;
  }

  TGraph gr_driftv(efields.size(), efields.data(), velocity);
  TGraph gr_diff_tra(efields.size(), efields.data(), diff_t);
  TGraph gr_diff_long(efields.size(), efields.data(), diff_l);

  driftv = gr_driftv.Eval(E_FIELD);
  diff_tra = gr_diff_tra.Eval(E_FIELD);
  diff_long = gr_diff_long.Eval(E_FIELD);

  return 1;
}

int mktrack::DefineDetector()
{
  geo = new GeometrySimple();
  box = new SolidBox(center[0], center[1], center[2], half[0], half[1], half[2]);
  if(gas==nullptr){
    std::cerr << "MediumMagboltz is not defined yet." << std::endl;
    return 0;
  }
  geo->AddSolid(box, gas);
  comp = new ComponentAnalyticField();
  comp->AddPlaneY(y_plate, v_plate, "plate");
  comp->AddPlaneY(y_grid, v_grid, "grid");
  comp->SetGeometry(geo);
  sensor = new Sensor();
  sensor->AddComponent(comp);
  sensor->SetArea(center[0]-half[0], center[1]-half[1], center[2]-half[2],
		  center[0]+half[0], center[1]+half[1], center[2]+half[2]);
  drift = new AvalancheMC();
  drift->SetSensor(sensor);

  return 1;
}

void mktrack::ShowSrim()
{
  std::cout << "**********************************************************" << std::endl;
  srim_beam->Print();
  std::cout << "**********************************************************" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  for(auto it=srim_particle.begin();it!=srim_particle.end();++it){
    std::cout << "**********************************************************" << std::endl;
    (*it)->Print();
    std::cout << "**********************************************************" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }

  return;
}

//std::vector<std::string> mktrack::GetParticleName()
//{
//  return *(particle_name.begin()+event_id);
//}

int mktrack::Generate(int &status, double &ex)
{
  ClearBuffer();

  status = 1;
  ex = event->Generate();

  vtx[0] = BEAM_X_CENTER;
  vtx[1] = BEAM_Y_CENTER;
  vtx[2] = VTX_Z_STOP;
  start_point[0] = vtx[0];
  start_point[1] = vtx[1];
  start_point[2] = 0.;
  beam_area[2][1] = vtx[2];
  
//  if(event->GetParticleVector(1).Theta()>92.*TMath::DegToRad()){
//    return 0;
//  }

  // judge if particle is stopped inside
  unsigned int N_particle = event->GetParticleNumber();
  for(unsigned int i_particle=1;i_particle<N_particle;++i_particle){
    if(EnetoRange[i_particle]==nullptr){
      continue;
    }
    double r = EnetoRange[i_particle]->Eval((event->GetParticleVector(i_particle).E()-
					     event->GetParticleVector(i_particle).M())*1000.);
    double dr = TMath::Sqrt(event->GetParticleVector(i_particle).Px()*event->GetParticleVector(i_particle).Px()+
			    event->GetParticleVector(i_particle).Py()*event->GetParticleVector(i_particle).Py()+
			    event->GetParticleVector(i_particle).Pz()*event->GetParticleVector(i_particle).Pz());
    double dx[3] = {event->GetParticleVector(i_particle).Px()/dr,
		    event->GetParticleVector(i_particle).Py()/dr,
		    event->GetParticleVector(i_particle).Pz()/dr};

    if(particle_flag[event_id][i_particle] &&
       ((r*dx[0]+vtx[0]<area[0][0] || r*dx[0]+vtx[0]>area[0][1] ||
	 r*dx[1]+vtx[1]<area[1][0] || r*dx[1]+vtx[1]>area[1][1] ||
	 r*dx[2]+vtx[2]<area[2][0] || r*dx[2]+vtx[2]>area[2][1])||
	r<0)
       ){
      return 0;
    }
  }

  status = 2;
  // beam track
  if(srim_beam!=nullptr){
    GenTrack(srim_beam, event->GetBeamVector(), start_point, beam_area, -1);
  }
  point.clear();
  range.clear();

  // other particles track
  for(unsigned int i_particle=0;i_particle<N_particle;++i_particle){
    if(srim_particle[i_particle]==nullptr){
      continue;
    }
    if(GenTrack(srim_particle[i_particle], event->GetParticleVector(i_particle), vtx, area, i_particle)==0){
      if(i_particle==0){
	point.clear();
	range.clear();
      }else if(particle_flag[event_id][i_particle]){
	return 0;
      }else{
	continue;
      }
    }
  } // end of for(int i_part...

//  if(point.size()!=event->GetParticleNumber()-1 || range.size()!=event->GetParticleNumber()-1){
//    return 0;
//  }

  // modify triger timing
  int shifted = ModTrack();

//  // modify labeled point "y-shifted"
//  for(unsigned int i=0;i<point.size();i++){
//    point[i][0][1] = point[i][0][1]-shifted;
//    point[i][1][1] = point[i][1][1]-shifted;
//    point[i][2][1] = point[i][2][1]-shifted;
//    point[i][3][1] = point[i][3][1]-shifted;
//  }

  return 1;
}

int mktrack::GenTrack(TrackSrim *srim, TLorentzVector particle_vec, double VTX[3], double sensed_area[3][2], int i_particle)
{
  srim->SetKineticEnergy((particle_vec.E()-particle_vec.M())*1.e9); // [GeV] -> [eV]

  double t_0 = 0.;

  std::vector<std::vector<double>> temp_point;
  
  if(!srim->NewTrack(VTX[0]*mmTocm, VTX[1]*mmTocm, VTX[2]*mmTocm, t_0,
		     particle_vec.Px(), particle_vec.Py(), particle_vec.Pz())){
    return 0;
  }

  int n_cluster = 0;
  int tot_ne = 0;
  double cluster_pos[4];
  double ele_end_pos[4];
  int ne;
  double ec;
  double ekin;
  int first_flag = 1;
  double drift_time = 0;
  std::vector<double> sca_a, sca_c, end_a, end_c;
  
  while(srim->GetCluster(cluster_pos[1], cluster_pos[2], cluster_pos[3], cluster_pos[0],
			 ne, ec, ekin)){
    if(cluster_pos[1]<sensed_area[0][0]*mmTocm || cluster_pos[1]>sensed_area[0][1]*mmTocm ||
       cluster_pos[2]<sensed_area[1][0]*mmTocm || cluster_pos[2]>sensed_area[1][1]*mmTocm ||
       cluster_pos[3]<sensed_area[2][0]*mmTocm || cluster_pos[3]>sensed_area[2][1]*mmTocm){
      return 0;
    }
    n_cluster++;
    tot_ne+=ne;
    
//    int drift_status;
    int add_ele;
    int ie = 0;


    while((ne-ie)!=0){
      if(ie+ie_step>ne){
	add_ele = ne-ie;
      }else{
	add_ele = ie_step;
      }
      
//      if(drift->AvalancheElectron(cluster_pos[1], cluster_pos[2], cluster_pos[3], cluster_pos[0])){
//	int ne_sub = drift->GetNumberOfElectronEndpoints();
      int ne_sub = 1;
	for(int ie_sub=0;ie_sub<ne_sub;++ie_sub){
//	  drift->GetElectronEndpoint(ie_sub,
//				     cluster_pos[1], cluster_pos[2], cluster_pos[3], cluster_pos[0],
//				     ele_end_pos[1], ele_end_pos[2], ele_end_pos[3], ele_end_pos[0],
//				     drift_status);
	  ElectronDrift(cluster_pos[1], cluster_pos[2], cluster_pos[3], cluster_pos[0],
			ele_end_pos[1], ele_end_pos[2], ele_end_pos[3], ele_end_pos[0]);
//	  if(ele_end_pos[2]<0 && drift_status==-5){
	    drift_time = (ele_end_pos[0]-cluster_pos[0]);
	    AddRawWave(ele_end_pos, drift_time, add_ele);
	    if(first_flag){
	      sca_a = std::vector<double>{(cluster_pos[3]*cmTomm/0.4), (drift_time/10.)};
	      sca_c = std::vector<double>{(cluster_pos[1]*cmTomm/0.4), (drift_time/10.)};
	      first_flag = 0;
	    } // end of if(first...
	    end_a = std::vector<double>{(cluster_pos[3]*cmTomm/0.4), (drift_time/10.)};
	    end_c = std::vector<double>{(cluster_pos[1]*cmTomm/0.4), (drift_time/10.)};
//	  } // end of if(ele_e...
	} // end of for(int ie_sub...
//      } // end of if(drift->Ava...

      ie+=add_ele;
    } // end of while(ne...

  } // end of while(srim...

  if(i_particle == -1){
    return 0;
  }else if(!particle_flag[event_id][i_particle]){
    return 0;
  }else if(first_flag == 1){
    return 0;
  }
  temp_point.push_back(sca_a);
  temp_point.push_back(sca_c);
  temp_point.push_back(end_a);
  temp_point.push_back(end_c);
  point.push_back(temp_point);
  range.push_back(TMath::Sqrt((VTX[0]-cluster_pos[1]*cmTomm)*(VTX[0]-cluster_pos[1]*cmTomm)+
			      (VTX[1]-cluster_pos[2]*cmTomm)*(VTX[1]-cluster_pos[2]*cmTomm)+
			      (VTX[2]-cluster_pos[3]*cmTomm)*(VTX[2]-cluster_pos[3]*cmTomm)));
  return 1;
}

int mktrack::ModTrack()
{
  int triger = -1;
  for(int jj=0;jj<1024;++jj){
    for(int kk=0;kk<256;++kk){
      if(flush[0][jj][kk]>threshold[0]){
	triger = jj;
	break;
      }
    }
    if(triger!=-1){
      break;
    }
  }

  if(triger>buff){
    for(int ii=0;ii<2;ii++){
      for(int kk=0;kk<256;kk++){
	for(int jj=0;jj<1024-triger+buff;jj++){
	  flush[ii][jj][kk] = flush[ii][jj+triger-buff][kk];
	}
	for(int jj=1024-triger+buff;jj<1024;jj++){
	  flush[ii][jj][kk] = 0;
	}
      }
    }
  }else{
    for(int ii=0;ii<2;ii++){
      for(int kk=0;kk<256;kk++){
	for(int jj=1023;jj>=buff-triger;jj--){
	  flush[ii][jj][kk] = flush[ii][jj-buff+triger][kk];
	}
	for(int jj=buff-triger-1;jj>=0;jj--){
	  flush[ii][jj][kk] = 0;
	}
      }
    }
  }

  return triger-buff;
}

double mktrack::GetFlush(int ii, int jj, int kk)
{
  return flush[ii][jj][kk];
}

int mktrack::GetTOT(int ii, int jj, int kk)
{
  if(flush[ii][jj][kk]>threshold[ii]){
    return 1;
  }else{
    return 0;
  }
}

void mktrack::ClearBuffer()
{
  for(int ii=0;ii<2;++ii){
    for(int jj=0;jj<1024;++jj){
      for(int kk=0;kk<256;++kk){
	flush[ii][jj][kk] = 0;
      }
    }
  }
}

void mktrack::AddRawWave(double ele_end_pos[], double drift_time, int ne)
{
  int strp[2];
  double ns;
  double pulse_height;
  double temp_height;
  double min_height;
  double max_height;
  int break_flag;
  int max_flag;
  
//  strp[0] = (int)(ele_end_pos[3]*cmTomm/0.4) + (int)rndm->Uniform(-1.2,1.2);
//  strp[1] = (int)(ele_end_pos[1]*cmTomm/0.4) + (int)rndm->Uniform(-1.2,1.2);
  strp[0] = (int)(ele_end_pos[3]*cmTomm/0.4);
  strp[1] = (int)(ele_end_pos[1]*cmTomm/0.4);
  
  min_height = wave->Eval(0)*ne*gain;
  temp_height = wave->Eval(-drift_time)*ne*gain;
//  temp_height = 0;
  for(int ac=0;ac<2;++ac){
    break_flag = 0;
    max_flag = 0;
    if(strp[ac]>=0 && strp[ac]<256){
      for(int clk=0;clk<1024;++clk){
	ns = clk*10.0;
	pulse_height = wave->Eval(ns-drift_time)*ne*gain;
	
	if(max_flag==0 && pulse_height<temp_height && pulse_height>min_height*2){
	  max_flag = 1;
	  max_height = temp_height;
	}
	
	if(pulse_height<0) pulse_height = 0;
	if(max_flag==1 && pulse_height<max_height*0.05){
	  pulse_height = 0;
	  break_flag = 1;
	}
//	flush[ac][clk][strp[ac]] += pulse_height;
	*(flush[ac][clk].begin()+strp[ac]) += pulse_height;
	temp_height = pulse_height;
	if(break_flag==1){
	  break;
	}
      } // end of for(int clk...
    } // end of if(strp..
  } // end of for(int ac...
}

int mktrack::ShowIdealValues(std::ostream& os, int exist)
{
  if(exist==0){
    os << "# ";
    for(unsigned int i=0;i!=event->GetParticleNumber();++i){
      os << "e" << i+2 << "[MeV] m" << i+2 << "[MeV] phi" << i+2 << "[rad] theta" << i+2 << "[rad] ";
    }
    os << std::endl;
  }
  for(unsigned int i=0;i<event->GetParticleNumber();++i){
    os << event->GetParticleVector(i).E()*1000 << " "
       << event->GetParticleVector(i).M()*1000 << " "
       << event->GetParticleVector(i).Phi() << " "
       << event->GetParticleVector(i).Theta() << " ";
  }
  os << std::endl;
  
  return 1;
}

int mktrack::ShowTeacherValues(std::ostream& os, int exist)
{
  if(exist==0){
    os << "# ";
//    for(unsigned int i=1;i<event->GetParticleNumber();i++){
//      os << "range" << i+2 << "[mm] phi" << i+2 << "[rad] theta" << i+2 << "[rad] ";
//      os << "sca_a_x sca_a_y sca_c_x sca_c_y end_a_x end_a_y end_c_x end_c_y [pixel]";
//    }
//    os << std::endl;
    for(unsigned int i=0;i<range.size();++i){
      os << "range" << i+3 << "[mm] phi" << i+3 << "[rad] theta" << i+3 << "[rad] ";
      os << "sca_a_x sca_a_y sca_c_x sca_c_y end_a_x end_a_y end_c_x end_c_y [pixel]";
    }
    os << std::endl;
  }
//  for(unsigned int i=0;i<event->GetParticleNumber()-1;i++){
  for(unsigned int i=0;i<range.size();++i){
    os << range[i] << " " << std::flush;
    os << event->GetParticleVector(i+1).Phi() << " " << std::flush;
    os << event->GetParticleVector(i+1).Theta() << " " << std::flush;
    os << point[i][0][0] << " " << std::flush;
    os << point[i][0][1] << " " << std::flush;
    os << point[i][1][0] << " " << std::flush;
    os << point[i][1][1] << " " << std::flush;
    os << point[i][2][0] << " " << std::flush;
    os << point[i][2][1] << " " << std::flush;
    os << point[i][3][0] << " " << std::flush;
    os << point[i][3][1] << " " << std::flush;
  }
  os << std::endl;
  
  return 1;
}

void mktrack::ElectronDrift(double cluster_pos1, double cluster_pos2, double cluster_pos3, double cluster_pos0,
			double &ele_end_pos1, double &ele_end_pos2, double &ele_end_pos3, double &ele_end_pos0)
{
  double drift_len = cluster_pos2-0;
  double sigma_tra = diff_tra*TMath::Sqrt(drift_len);
  double sigma_long = diff_long*TMath::Sqrt(drift_len);
//  double sigma_tra = 0.08; // [cm]
//  double sigma_long = 0.08; // [cm]

  ele_end_pos0 = (drift_len+rndm->Gaus(0, sigma_long))/(driftv); // ns
  if(ele_end_pos0<0){
    ele_end_pos0 = 0;
  }
  ele_end_pos1 = cluster_pos1+rndm->Gaus(0, sigma_tra); // cm
  ele_end_pos2 = 0;
  ele_end_pos3 = cluster_pos3+rndm->Gaus(0, sigma_tra); // cm
  
  return;
}

double mktrack::GetDriftv()
{
  return driftv;
}
