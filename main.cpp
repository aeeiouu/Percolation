#include "LFile.h"
#include "radial.h"

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef VECTOR
#define VECTOR
#include <vector>
#endif

using namespace std;

int main()
{
	LIBRARY lib("h2o.lammpstrj");
	PAGE* page = lib.GetPage();


	double percol;
	PERCOLATION percolation(&lib);
	percol = percolation.Get_Percolation(80);
	cout << "the number of percolation: " << percol << endl;

	cout << "PageCount: " << lib.GetPagecount() << endl;
	
	cout << "the ratio of percolation: " << percol / lib.GetPagecount() << endl;
	
	system("pause");
	return 0;
}
