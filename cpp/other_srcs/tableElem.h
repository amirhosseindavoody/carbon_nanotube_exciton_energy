#pragma once

#include <Eigen>

#include "CNT.h"
#include "segment.h"

using namespace std;
using namespace Eigen;

class tableElem
{ 

	double r; //distance between segment and segment with ID: tubeidx, segidx
	double theta; //Angle in radians between two segments
	double gamma; // this is the constant that is used in calculating the conventional Forster rates
	double rate; //Transistion rate to be used in simulation
	int tubeidx; //tube interacting with from current tube
	int segidx; //segment number on the tube of tubeidx

private:
	void set_rate();

public:
	tableElem(double distance=1, double angle=0, double forster_const=0, int tube_idx=0, int seg_idx=0);
	double getRate();
	double getr();
	double getTheta();
	double getGamma();
	int getTubeidx();
	int getSegidx();
};

