#pragma once
#ifndef RADIAL_H
#define RADIAL_H

#define PI 3.14159265359
constexpr auto NA = 6.0221409e+23;

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif 

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef ALGORITHM
#define ALGORITHM
#include <algorithm>
#endif

#ifndef VECTOR
#define VECTOR
#include <vector>
#endif

#ifndef MUTEX
#define MUTEX
#include <mutex>
#endif

#ifndef THREAD
#define TRHEAD
#include <thread>
#endif

#ifndef CHRONO
#define CHRONO
#include <chrono>
#endif

#ifndef QUEUE
#define QUEUE
#include <queue>
#endif

#ifndef CONDITION_VARIABLE
#define CONDITION_VARIABLE
#include <condition_variable>
#endif

#ifndef CSTDINT
#define CSTDINT
#include <cstdint>
#endif

#include "LFile.h"
using namespace std;

class RADIAL
{
private:
	LIBRARY* PLIB;
	vector<int>**  HB;
	int Num_Atom;
	RADIAL(const RADIAL& radial);

public:
	RADIAL()
	{
		PLIB = nullptr;
		Num_Atom = 0;
		HB = nullptr;
	}
	RADIAL(LIBRARY* plib,int num_atom);
	void RDF(int type1, int type2, double r, double dr);
	void hydrogenbond(int definition);
	double OSP(int num_atom);
	void SOP();
	vector<int>** GetHydrobond();
	~RADIAL();

	double periodic(double pos1, double pos2, double* BOUNDARY);
	void periodic(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int j, int k, double *delta);
	void checkdelta(int* const neighbor, double * const Ndelta, int j, double dis2delta);
	double dis2(double* delta);
	void crossproduct(double* delta1, double* delta2, double* delta3);
	bool oxygencheck(double* del,int definition);
	bool hydrogencheck(double* del,int definition);
	double cos_angle(double* delta1, double* delta2);
	bool anglecheck(double* vec1, double* vec2,int definition);
	//void HBlifetime();
	void HBlifetime2();

};

struct inf_chain
{
	double Coord[3];
	bool able = true;
};

class PERCOLATION : public RADIAL
{
private:
	PERCOLATION();
	LIBRARY* PLIB;
	void set_step(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int count, int*& order_number, bool*& memo, int Molnum,int& memo_count);
	bool CheckPercolation(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int* order_number,int* number_of_cluster, int MAX_order_number,int Molnum);
	bool CheckPercolation_step(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int> index);
	bool Check_remained_path(int** M_Path, int start, int number);
	bool Check_remained_path(vector<int>** M_N_Path, int start, int number);
	bool Check_remained_path(bool** memo_M_N_Path, int start, int number);
	int Check_remained_path(vector<int>** M_N_Path, int start, int number, vector<int>* D_path, bool* memo_Path, vector<int>* node);
	bool Belong_To_Node(vector<int>* node, vector<int>* D_path, int chain_number);
	bool Check_loop(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int>* index, vector<int>* memo);
	bool Check_Matrix(int** M_path, int number);
	bool Check_Matrix(vector<int>** M_N_path, int number);
	void Trim_Matrix(int** M_path, int* number_path, int number);
	bool Determination(int** M_path, int* number_path, double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int number, vector<int> index);
	void Calculate_Vector(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int> index, vector<int>* memo, double* Coord);
	bool Check_Inner(double* memo, double* Coord, bool reverse);
	bool Check_vector(double* memo, double* Coord, bool reverse);
	void Chain(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int>* index, vector<int>** M_N_Path, bool** memo_M_N_Path, int number, vector<int>* D_path, int N_chain, inf_chain* I_D_Path, bool* flag,queue<int>* master_job, int address);
	void set(int start, int end);
	void Push_job_into_Master(queue<int>* master_job, queue<int>* job, int index);

	int synchro = 0;
	int num_Percolation = 0;
	int Set_Num_threads = 0;
	int Num_threads = 0;
	mutex M, job_m;
	std::condition_variable cv, cv2, cv_wait;
	int pass_count = 0;
	std::chrono::milliseconds time = std::chrono::milliseconds(0);
	std::chrono::milliseconds time2 = std::chrono::milliseconds(0);
	bool Check_path = true;


public:
	PERCOLATION(LIBRARY* plib);
	bool hydrogenbondcheck(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int i, int j);
	int Get_Percolation(int num_threads);
};
#endif