#ifndef _progress_h_
#define _progress_h_

#include <iostream>

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

public:
  // constructor to initialize internal state of the progress bar
  progress_bar(const int& i_max, const std::string& title="")
  {
    _title = title;
    _i_max = i_max-1;
    _i=0;
    std::cout << _title << ":\n";
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
    std::cout << "] " << int((_progress) * 100.0) << "%\r";
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
    std::cout << "] " << int((_progress) * 100.0) << "%\r";
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
    std::cout << "] " << int((_progress) * 100.0) << "%\r";
    std::cout.flush();
    if (i == i_max-1)
    {
      std::cout << std::endl;
    }
  };
};

#endif //_progress_h_