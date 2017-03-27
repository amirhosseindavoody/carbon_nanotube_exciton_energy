#ifndef functions_h
#define functions_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include "dirent.h"
#include <list>
#include "CNT.h"
#include <vector>
#include <memory>
#include "tableElem.h"
#include <locale>
#include <math.h>
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include <regex>
#include <omp.h>
#include <stdint.h>

#include "exciton.h"


//method declarations
void updateExcitonList(int numExcitonsAtCont, vector<exciton> &excitons, vector<int> &currCount, vector<shared_ptr<segment>> inContact);
void writeStateToFile(ofstream &file, vector<int> &currCount, double time);
void markCurrentExcitonPosition(vector<CNT> &cnt_list, exciton &curr_exciton, vector<int> &currCount, vector<double> &regionBdr);
bool hasMovedToOutContact(exciton &curr_exciton, vector<double> &regionBdr, vector<CNT> &cnt_list);
void injectExciton(exciton &curr_exciton, vector<shared_ptr<segment>> &inContact);
void assignNextState(vector<CNT> &cnt_list, exciton &curr_exciton, double gamma, vector<double> &regionBdr);
void add_self_scattering(vector<CNT> &cnt_list, double maxGam);
double getRand(bool excludeZero);
double make_rate_table(vector<shared_ptr<segment>> &seg_list, segment &seg, double maxDist);
double convert_units(string unit, double val);
void init_random_number_generator();


#endif // functions_h