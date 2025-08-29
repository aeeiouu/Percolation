#pragma once

#ifndef LFILE_H
#define LFILE_H

#ifndef FSTREAM
#define FSTREAM
#include <fstream>
#endif

#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef STRING
#define STRING
#include <string>
#endif

class LIBRARY {
private:
	class PAGE* page;
	int PAGECOUNT;
	LIBRARY(const LIBRARY& );
	LIBRARY();
public:
	LIBRARY(const char* path);
	~LIBRARY();
	int GetPagecount();
	PAGE* GetPage();
};

class PAGE
{
private:
	int TIMESTEP;
	int NUMBER;
	double BOUNDARY[6];
	int* ID;
	int* TYPE;
	double* POSITION[3];
	double* VELOCITY[3];
	PAGE(const PAGE& );
public:
	PAGE();
	~PAGE();
	void SetTimeStep(int time);
	int GetTimeStep();
	void SetNumber(int number);
	int GetNumber();
	void SetBoundary(double value, int index);
	double GetBoundary(int index);
	void SetMemory();
	int* GetID();
	int* GetType();
	double* GetPosition(int index);
	double* GetVelocity(int index);
	double* GetBoundary();
};
#endif