#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <thread>
#include <mutex>
#include <cmath>

#include "mktrack.hpp"

int check_cmd(int argc, char *argv[], std::vector<std::string> reaction, unsigned int &event_num, std::string &fname, int &scatter, int &pressure);
void show_progress(int &sum_num, unsigned int event_num, int &status, double &Ex);

int main(int argc, char *argv[]){

  unsigned int event_num;
  std::string fname;
  int scatter;
  int pressure;

  std::vector<std::string> reaction = {"12C(n,n)12C",
				       "12C(n,n)3a (7.65)",
				       "12C(n,n)3a (7.65)",
				       "12C(n,n)3a (9.64)",
				       "12C(n,n)3a (9.64)",
				       "p(n,n)p"};

  if(check_cmd(argc, argv, reaction, event_num, fname, scatter, pressure) == 0){
    return 1;
  }
  
  std::string datdir = "data/";

  std::ofstream flush_out(datdir+fname+"_flush.dat");
  std::ofstream tot_out(datdir+fname+"_tot.dat");
  std::ofstream idealvalue_out(datdir+fname+"_idealvalue.dat");
  std::ofstream teachervalue_out(datdir+fname+"_teachervalue.dat");
  std::mutex mtx;
  int exist_ideal = 0;
  int exist_teacher = 0;

//  size_t mp = std::thread::hardware_concurrency();
//  if(mp==0){
//    mp = 1;
//  }else{
//    mp = (size_t)(mp-2);
//    if(mp>event_num){
//      mp = event_num;
//    }
//  }
  size_t mp = 1;

  std::vector<std::thread> ths(mp);

//  std::cout << ths.size() << std::endl;
  
  int step_num = event_num/ths.size();
  int sum_num = 0;
  int tot_num;
  int event_num_temp = event_num;
  int status = 0;
  double Ex = 0;
  double drift_v = -1;
//  for(int ii=0;ii<ths.size();ii++){
  for(auto itr=ths.begin();itr!=ths.end();++itr){
//    std::cerr << ii << std::endl;
    if(event_num_temp-step_num*ths.size()>0){
      tot_num = step_num+1;
      event_num_temp--;
    }else{
      tot_num = step_num;
    }

    // generate events by multi threads
    *itr = std::thread([&mtx, &flush_out, &tot_out, &idealvalue_out, &teachervalue_out, &exist_ideal, &exist_teacher, &sum_num, &status, &Ex, &drift_v](int scatter, int pressure, int tot_num){
//	std::string fname_2 = "temp_tot_"+std::to_string(ii++)+".dat";
//	std::ofstream tot_out_2(datdir+fname_2);
	mtx.lock();
	mktrack *MAIKo = new mktrack(scatter, pressure);
	drift_v = MAIKo->GetDriftv();
	mtx.unlock();
	for(int num=0;num<tot_num;){
	  // generate event & get picutres
	  while(MAIKo->Generate(status, Ex)==0){}
	  
	  mtx.lock();
	  status = 3;
//	  std::cout << num << std::endl;
	  exist_ideal = MAIKo->ShowIdealValues(idealvalue_out, exist_ideal);
	  exist_teacher = MAIKo->ShowTeacherValues(teachervalue_out, exist_teacher);
	  
	  // output event
	  for(int ii=0;ii<2;ii++){
	    for(int jj=0;jj<1024;jj++){
	      for(int kk=0;kk<256;kk++){
		tot_out << MAIKo->GetTOT(ii, jj, kk) << " ";
//		tot_out_2 << MAIKo.GetTOT(ii, jj, kk) << " ";
		flush_out << MAIKo->GetFlush(ii, jj, kk) << " ";
	      }
	    }
	  }
	
	  tot_out << std::endl;
//	  tot_out_2 << std::endl;
	  flush_out << std::endl;
	  num++;
	  sum_num++;
	  mtx.unlock();
	  status = 0;
//	  tot_out_2.close();
	}
	delete MAIKo;
      }, scatter, pressure, tot_num);
  }

  while(drift_v<0){
    std::this_thread::sleep_for(std::chrono::seconds(1));
  }
  std::ofstream param_out(datdir+fname+"_param.dat");
  for(int ii=0;ii<argc;ii++){
    param_out << argv[ii] << " ";
  }
  param_out << std::endl;
  param_out << "Reaction              : " << reaction[scatter] << std::endl;
  param_out << "Number of events      : " << event_num << std::endl;
  param_out << "Pressure (hPa)        : " << pressure << std::endl;
  param_out << "Drift velocity (cm/ns): " << drift_v << std::endl;
  param_out.close();
  
  std::cout << "\e[?25l" << std::flush; // disappear cursol
  show_progress(sum_num, event_num, status, Ex);
      
  for(unsigned int ii=0;ii<ths.size();ii++){
    ths[ii].join();
  }
  
  flush_out.close();
  tot_out.close();

  std::cout << "Generated " << sum_num << " events." << std::endl;
  std::cout << "\e[?25h" << std::flush; // appear cursol

  return 0;
}

int check_cmd(int argc, char *argv[], std::vector<std::string> reaction, unsigned int &event_num, std::string &fname, int &scatter, int &pressure)
{
  if(argc<2){
    std::cout << "-t, --test: run by default values" << std::endl;
    std::cout << "-n, --number: number to generate" << std::endl;
    std::cout << "-f, --file: filename" << std::endl;
    std::cout << "-p, --pressure: pressure" << std::endl;
    std::cout << "-s, --scatter: select pattern" << std::endl;
    std::cout << "-l, --list: show scattering pattern" << std::endl;
    return 0;
  }else{
    std::map<std::string, int> list = {{"-t", 1},
				       {"--test", 1},
				       {"-n", 2},
				       {"--number", 2},
				       {"-f", 3},
				       {"--file", 3},
				       {"-p", 4},
				       {"--pressure", 4},
				       {"-s", 5},
				       {"--scatter", 5},
				       {"-l", 6},
				       {"--list", 6}
    };
    int i = 1;
    event_num = 1;
    fname = "temp";
    scatter = 0;
    pressure = 500;
    int reaction_num=0;
    while(i<argc){
      std::string opt = argv[i++];
      switch(list[opt]){
      case 0:
	std::cout << "-t, --test: run by default values" << std::endl;
	std::cout << "-n, --number: number to generate" << std::endl;
	std::cout << "-f, --file: filename" << std::endl;
	std::cout << "-p, --pressure: pressure" << std::endl;
	std::cout << "-s, --scatter: select pattern" << std::endl;
	std::cout << "-l, --list: show scattering pattern" << std::endl;
	return 0;
      case 1:
	break;
      case 2:
	if(argv[i][0] == '-'){
	  break;
	}else{
	  event_num = atoi(argv[i++]);
	}
	break;
      case 3:
	if(argv[i][0] == '-'){
	  break;
	}else{
	  fname = argv[i++];
	}
	break;
      case 4:
	if(argv[i][0] == '-'){
	  break;
	}else{
	  pressure = atoi(argv[i++]);
	}
	break;
      case 5:
	if(argv[i][0] == '-'){
	  break;
	}else{
	  scatter = atoi(argv[i++]);
	}
	break;
      case 6:
	for(auto itr=reaction.begin();itr!=reaction.end();++itr){
	  std::cout << reaction_num++ << ": " << *itr << std::endl;
	}
	return 0;
      default:
	std::cout << "-t, --test: run by default values" << std::endl;
	std::cout << "-n, --number: number to generate" << std::endl;
	std::cout << "-f, --file: filename" << std::endl;
	std::cout << "-p, --pressure: pressure" << std::endl;
	std::cout << "-s, --scatter: select pattern" << std::endl;
	std::cout << "-l, --list: show scattering pattern" << std::endl;
	return 0;
      }
    }
  }

  return 1;
}

void show_progress(int &sum_num, unsigned int event_num, int &status, double &Ex)
{
  double percent;
  unsigned int max_chr = 50;
  int step = 0;
  int num_size = log10(event_num)+1;
  std::string progress;
  time_t start = time(NULL);
  time_t pass;

  std::cout << std::endl << std::endl;
  
  while(percent!=100){
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    std::cout << "Status: ";
    switch(status){
    case 0:
      std::cout << "ready to generate track";
      break;
    case 1:
      std::cout << "generating event       ";
      break;
    case 2:
      std::cout << "generating track       ";
      break;
    case 3:
      std::cout << "outputing to file      ";
      break;
    default:
      std::cout << "unknown                ";
      break;
    }

    std::cout << "(Ex = " << std::setfill(' ') << std::setw(6) << std::right << Ex;
    std::cout << std::setfill(' ') << std::setw(10) << std::left << " MeV)" << std::endl;

    progress = "";
    percent = (sum_num*100)/(double)event_num;
    for(unsigned int ii=0;ii<(sum_num*max_chr)/event_num;ii++){
      progress += "#";
    }
    if(progress.size()<max_chr){
      switch(step){
      case 0:
	step = 1;
	progress += "-";
	break;
      case 1:
	step = 2;
	progress += "\\";
	break;
      case 2:
	step = 3;
	progress += "|";
	break;
      case 3:
	step = 0;
	progress += "/";
	break;
      default:
	step = 0;
	progress += "-";
	break;
      }
    }

    pass = time(NULL)-start;
    
    std::cout << "[" << "\e[1;33;49m"
	      << std::setfill(' ') << std::setw(max_chr) << std::left << progress << "\e[0m" << "]"
	      << std::setfill(' ') << std::setw(5) << std::right << std::fixed << std::setprecision(1) << percent << "%"
	      <<"(" << std::setfill(' ') << std::setw(num_size) << std::right << sum_num
	      << "/" << event_num << ")"
	      << " " << pass/(60*60) << ":"
	      << std::setfill('0') << std::setw(2) << std::right << (pass/60)%60 << ":"
	      << std::setfill('0') << std::setw(2) << std::right << pass%60 <<  std::flush;
    if(progress.size() == max_chr && percent == 100){
      std::cout << std::endl;
      break;
    }
    std::cout << std::endl;
    std::cout << "\r" << "\e[2A" << std::flush;
  }
  return;
}
