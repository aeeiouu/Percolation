#include "LFile.h"

const int MAX = 256;

using namespace std;

LIBRARY::LIBRARY(const char* path)
{
	ifstream FILE;
	FILE.open(path);

	int count = 0;
	char buffer[MAX];
	char* ptr;
	char* context = NULL;

	while (!FILE.eof())
	{
		FILE.getline(buffer, MAX);
		if (!strcmp("ITEM: TIMESTEP\0", buffer))
		{
			count++;
		}
	}

	FILE.close();

	page = new PAGE[count];
	PAGECOUNT = count;

	FILE.open(path);

	count = -1;

	while (!FILE.eof())
	{
		FILE.getline(buffer, MAX);
		if (!strcmp(buffer, "ITEM: TIMESTEP\0"))
		{
			FILE.getline(buffer, MAX);
			count++;
			page[count].SetTimeStep(atoi(buffer));
		}
		FILE.getline(buffer, MAX);
		if (!strcmp(buffer, "ITEM: NUMBER OF ATOMS\0"))
		{
			FILE.getline(buffer, MAX);
			page[count].SetNumber(atoi(buffer));
			page[count].SetMemory();
		}
		FILE.getline(buffer, MAX);
		if (!strcmp(buffer, "ITEM: BOX BOUNDS pp pp pp\0"))
		{
			for (int i = 0; i < 3; i++)
			{
				FILE.getline(buffer, MAX);
				ptr = strtok_s(buffer, " ", &context);
				page[count].SetBoundary(atof(ptr), 2 * i);
				ptr = strtok_s(NULL, " ", &context);
				page[count].SetBoundary(atof(ptr), 2 * i + 1);
			}
		}
		FILE.getline(buffer, MAX);
		if (!strcmp(buffer, "ITEM: ATOMS id type x y z vx vy vz \0"))
		{
			for (int j = 0; j < page[count].GetNumber(); j++)
			{
				FILE.getline(buffer, MAX);
				ptr = strtok_s(buffer, " ", &context);
				page[count].GetID()[j]= atoi(ptr);
				ptr = strtok_s(NULL, " ", &context);
				page[count].GetType()[j]= atoi(ptr);

				for (int m = 0; m < 3; m++)
				{
					ptr = strtok_s(NULL, " ", &context);
					page[count].GetPosition(m)[j]=atof(ptr);
				}
				for (int m = 0; m < 3; m++)
				{
					ptr = strtok_s(NULL, " ", &context);
					page[count].GetVelocity(m)[j]= atof(ptr);
				}
			}
		}
	}
	FILE.close();
}

LIBRARY::~LIBRARY()
{
	delete[] page;
}

int LIBRARY::GetPagecount()
{
	if (this == nullptr)
	{
		return 0;
	}
	return PAGECOUNT;
}

PAGE* LIBRARY::GetPage()
{
	return page;
}

//PAGE
PAGE::PAGE() :TIMESTEP(0), NUMBER(0),ID(nullptr),TYPE(nullptr)
{
	for (int k = 0; k < 3; k++)
	{
		POSITION[k] = nullptr;
		VELOCITY[k] = nullptr;
	}
	memset(BOUNDARY, 0, sizeof(double) * 6);
}

void PAGE::SetTimeStep(int time)
{
	TIMESTEP = time;
}

int PAGE::GetTimeStep()
{
	return TIMESTEP;
}

void PAGE::SetNumber(int number)
{
	NUMBER = number;
}

int PAGE::GetNumber()
{
	return NUMBER;
}

void PAGE::SetBoundary(double value, int index)
{
	BOUNDARY[index] = value;
}

double PAGE::GetBoundary(int index)
{
	return BOUNDARY[index];
}

void PAGE::SetMemory()
{
	if (NUMBER == 0)
	{
		cout << "You should set Number first" << endl;
		return;
	}
	else
	{
		ID = new int[NUMBER];
		TYPE = new int[NUMBER];
		for (int k = 0; k < 3; k++)
		{
			POSITION[k] = new double[NUMBER];
			VELOCITY[k] = new double[NUMBER];
		}
	}
}

int* PAGE::GetID()
{
	return ID;
}

int* PAGE::GetType()
{
	return TYPE;
}

double* PAGE::GetPosition(int index)
{
	return POSITION[index];
}

double* PAGE::GetVelocity(int index)
{
	return VELOCITY[index];
}

double* PAGE::GetBoundary()
{
	return BOUNDARY;
}

PAGE::~PAGE()
{
	delete[] ID;
	delete[] TYPE;
	for (int k = 0; k < 3; k++)
	{
		delete[] POSITION[k];
		delete[] VELOCITY[k];
	}
}