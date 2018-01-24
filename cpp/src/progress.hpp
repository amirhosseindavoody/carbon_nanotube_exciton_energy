#ifndef _progress_h_
#define _progress_h_

#include <iostream>
#include <ctime>
#include <string>
#include <iomanip>

// helper class to monitor progress of the loops
class progress_bar
{
private:
  const int _barWidth=100;
  int _pos=0;
  float _progress=0;
  std::string _title="";
  int _i_max=0;
  int _i=0;

  bool _start_chrono = true;
  std::time_t _start_time;
  std::time_t _current_time;


public:
  // constructor to initialize internal state of the progress bar
  progress_bar(const int& i_max, const std::string& title="")
  {
    _title = title;
    _i_max = i_max-1;
    _i=0;
    std::cout << _title << ":\n";
  };

  std::string estimate_remaining_time()
  {

    std::string remaining_time;
    if (_start_chrono){
      _start_time = std::time(nullptr);
      _current_time = std::time(nullptr);
      _start_chrono = false;
    }
    else {
      _current_time = std::time(nullptr);
    }

    int sec = int(std::difftime(_current_time,_start_time)*(1-_progress)/(_progress+0.001));

    int hour = sec/3600;
	  sec = sec%3600;
	  int min = sec/60;
	  sec = sec%60;
	  std::stringstream ss;
    ss << "remaining time " << std::setw(2) << std::setfill('0') << hour << ":" 
                            << std::setw(2) << std::setfill('0') << min << ":" 
                            << std::setw(2) << std::setfill('0') << sec;
    return ss.str();
  };



  // stepping function to keep track of progress internally
  void step()
  {
    _progress = float(_i)/float(_i_max);
    _pos = _progress*_barWidth;
    std::cout << "[";
    for (int j=0; j<_barWidth; ++j) {
      if (j < _pos) std::cout << "=";
      else if (j == _pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int((_progress) * 100.0) << "% " << estimate_remaining_time() << "\r";
    if (_i == _i_max)
    {
      std::cout << "\n";
    }
    std::cout.flush();
    _i++;
  };

  // stepping function to pass progress level explicitly and do the rest internally
  void step(const int& i)
  {
    _progress = float(i)/float(_i_max);
    _pos = _progress*_barWidth;
    std::cout << "[";
    for (int j=0; j<_barWidth; ++j) {
      if (j < _pos) std::cout << "=";
      else if (j == _pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int((_progress) * 100.0) << "% " << estimate_remaining_time() << "\r";
    if (i == _i_max)
    {
      std::cout << "\n";
    }
    std::cout.flush();
  };

  // constructor without any internal state initialization
  progress_bar() {};

  // stepping function that does not need initialization
  void step(const int& i, const int& i_max, const std::string& title)
  {
    _progress = float(i)/float(i_max-1);
    _pos = _progress*_barWidth;
    std::cout << title << ": [";
    for (int j=0; j<_barWidth; ++j) {
      if (j < _pos) std::cout << "=";
      else if (j == _pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int((_progress) * 100.0) << "% " << estimate_remaining_time() << "\r";
    std::cout.flush();
    if (i == i_max-1)
    {
      std::cout << std::endl;
    }
  };
};

#endif //_progress_h_