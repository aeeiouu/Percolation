#include "radial.h"

RADIAL::RADIAL(LIBRARY* plib,int num_atom):PLIB(plib),Num_Atom(num_atom)
{
	int Molnum = PLIB->GetPage()->GetNumber()/Num_Atom;
	HB = new vector<int>* [PLIB->GetPagecount()];
	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		HB[t] = new vector<int>[Molnum];
	}
}

RADIAL::~RADIAL()
{
	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		delete[] HB[t];
	}
	delete[] HB;
}

void RADIAL::RDF(int type1, int type2, double r, double dr)
{
	if (dr == 0)
	{
		cout << "Stepsize Error" << endl;
		return;
	}

	FILE *pfile;
	fopen_s(&pfile, "RDF", "w");

	PAGE* pPage = PLIB->GetPage();
	double* Xpos, *Ypos, *Zpos;
	double* BOUNDARY;
	int* pType;

	double delta = 0;
	int step = (int)(r / dr);
	double rad = 0;
	double delV = 0;
	double volume = 0;
	double radial = 0;

	int count_Type1;
	int count_Type2;

	bool step_stepsize = true;
	double* count;

	if ((r - ((double)step)*dr != 0))
	{
		count = new double[step + 1];
		memset(count, 0, sizeof(double)*(step + 1));
		step_stepsize = true;
	}
	else
	{
		count = new double[step];
		memset(count, 0, sizeof(double)*(step));
		step_stepsize = false;
	}

	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		Xpos = pPage[t].GetPosition(0);
		Ypos = pPage[t].GetPosition(1);
		Zpos = pPage[t].GetPosition(2);
		pType = pPage[t].GetType();
		BOUNDARY = pPage[t].GetBoundary();


		count_Type1 = 0;
		count_Type2 = 0;
		double dis = 0;
		volume = pow(BOUNDARY[1] - BOUNDARY[0], 3);

		for (int m = 0; m < pPage[t].GetNumber(); m++)
		{
			if (pType[m] == type1)
			{
				count_Type1++;
			}
			else if (pType[m] == type2)
			{
				count_Type2++;
			}

			for (int n = m + 1; n < pPage[t].GetNumber(); n++)
			{
				if ((pType[m] == type1) && (pType[n] == type2))
				{
					delta = pow(periodic(Xpos[m], Xpos[n], BOUNDARY), 2);
					delta += pow(periodic(Ypos[m], Ypos[n], (BOUNDARY + 2)), 2);
					delta += pow(periodic(Zpos[m], Zpos[n], (BOUNDARY + 4)), 2);

					dis = 0;
					for (int index = 0; index < step; index++)
					{
						dis = ((double)index) * dr;
						if ((delta >= dis * dis) && (delta < (dis + dr)*(dis + dr)))
						{
							count[index] = count[index] + 2;
							break;
						}
					}
					if (step_stepsize == true)
					{
						dis = ((double)step) * dr;
						if ((delta >= dis * dis) && (delta < r*r))
						{
							count[step] = count[step] + 2;
						}
					}
				}
			}
		}


		for (int index = 0; index < step; index++)
		{
			count[index] /= count_Type1;

			dis = ((double)index) * dr;
			delV = pow(dis + dr, 3) - pow(dis, 3);
			delV *= 4 * PI / 3;
			count[index] /= delV;
			count[index] /= (count_Type2 / volume);

			fprintf_s(pfile, "%f ", count[index]);
		}

		if (step_stepsize == true)
		{
			count[step] /= count_Type1;

			dis = ((double)step) * dr;
			delV = pow(r, 3) - pow(dis, 3);
			delV *= 4 * PI / 3;
			count[step] /= delV;
			count[step] /= (count_Type2 / volume);

			fprintf_s(pfile, "%f ", count[step]);
		}
		fprintf_s(pfile, ";\n");

		if (step_stepsize == true)
		{
			memset(count, 0, sizeof(double)*(step + 1));
		}
		else
		{
			memset(count, 0, sizeof(double)*(step));
		}
	}

	delete[] count;
	fclose(pfile);
}

double RADIAL::periodic(double pos1, double pos2, double* BOUNDARY)
{
	double size = BOUNDARY[1] - BOUNDARY[0];
	double delta = 0;
	delta = pos1 - pos2;
	if (delta > size*0.5)
		delta = delta - size;
	else if (delta < -size * 0.5)
		delta = delta + size;

	return delta;
}

void RADIAL::periodic(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int j, int k, double *delta)
{
	double size = BOUNDARY[1] - BOUNDARY[0];
	delta[0] = Xpos[k] - Xpos[j];
	if (delta[0] > size*0.5)
		delta[0] = delta[0] - size;
	else if (delta[0] < -size * 0.5)
		delta[0] = delta[0] + size;

	size = BOUNDARY[3] - BOUNDARY[2];
	delta[1] = Ypos[k] - Ypos[j];
	if (delta[1] > size*0.5)
		delta[1] = delta[1] - size;
	else if (delta[1] < -size * 0.5)
		delta[1] = delta[1] + size;

	size = BOUNDARY[5] - BOUNDARY[4];
	delta[2] = Zpos[k] - Zpos[j];
	if (delta[2] > size*0.5)
		delta[2] = delta[2] - size;
	else if (delta[2] < -size * 0.5)
		delta[2] = delta[2] + size;
}

double RADIAL::OSP(int num_atom)
{
	PAGE* page = PLIB->GetPage();

	int molNum = page[0].GetNumber();
	molNum /= num_atom;

	double* Xpos;
	double* Ypos;
	double* Zpos;
	double* boundary;

	double ensemble = 0;
	double total_ensemble = 0;

	double temp;

	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		double delta = 0;
		double Vdelta1[3] = { 0, };
		double Vdelta2[3] = { 0, };

		Xpos = page[t].GetPosition(0);
		Ypos = page[t].GetPosition(1);
		Zpos = page[t].GetPosition(2);
		boundary = page[t].GetBoundary();

		for (int i = 0; i < molNum; i++)
		{
			int NI[4] = { -1,-1,-1,-1 };
			double ND[4] = { 0, };

			for (int j = 0; j < molNum; j++)
			{
				if (i != j)
				{
					delta = 0;
					delta += pow(periodic(Xpos[3 * i], Xpos[3 * j], boundary), 2);
					delta += pow(periodic(Ypos[3 * i], Ypos[3 * j], boundary), 2);
					delta += pow(periodic(Zpos[3 * i], Zpos[3 * j], boundary), 2);
					checkdelta(NI, ND, j, delta);
				}
			}

			for (int m = 0; m < 4; m++)
			{
				for (int n = m + 1; n < 4; n++)
				{
					periodic(Xpos, Ypos, Zpos, boundary, 3 * i, 3 * NI[m], Vdelta1);
					periodic(Xpos, Ypos, Zpos, boundary, 3 * i, 3 * NI[n], Vdelta2);

					temp = cos_angle(Vdelta1, Vdelta2);
					temp = temp + (double)1 / (double)3;
					temp = temp * temp * 3 / 8;
					ensemble = ensemble + temp;
				}
			}

			for (int y = 0; y < 4; y++)
			{
				NI[y] = -1;
				ND[y] = 0;
			}
		}

		ensemble /= molNum;
		total_ensemble += 1 - ensemble;
		ensemble = 0;
	}

	return total_ensemble / PLIB->GetPagecount();
}

void RADIAL::checkdelta(int* neighbor, double * Ndelta, int index, double delta)
{
	static int temp;
	static double tempDelta;
	static int Count = 0;

	if (*(neighbor) == -1 && Count != 4)
	{
		*(neighbor) = index;
		*(Ndelta) = delta;
	}
	else
	{
		if (*(Ndelta) > delta)
		{
			tempDelta = *(Ndelta);
			*(Ndelta) = delta;
			temp = *(neighbor);
			*(neighbor) = index;
			if (Count < 4)
			{
				Count = Count + 1;
				checkdelta(neighbor + 1, Ndelta + 1, temp, tempDelta);
			}
			else
			{
				Count = 0;
			}
		}
		else
		{
			if (Count < 4)
			{
				Count = Count + 1;
				checkdelta(neighbor + 1, Ndelta + 1, index, delta);

			}
			else
			{
				Count = 0;
			}
		}
	}
}

double RADIAL::cos_angle(double* delta1, double* delta2)
{
	double value = 0;
	for (int i = 0; i < 3; i++)
	{
		value += delta1[i] * delta2[i];
	}
	value = value / (dis2(delta1)*dis2(delta2));

	return value;

}

double RADIAL::dis2(double* delta)
{
	double value = 0;
	for (int i = 0; i < 3; i++)
	{
		value += delta[i] * delta[i];
	}

	value = sqrt(value);

	return value;
}

bool RADIAL::oxygencheck(double* del, int definition)
{
	switch (definition)
	{
	case 1:
		if (dis2(del) < 3.6)
			return true;
		else
			return false;
		break;

	case 2:
		if (dis2(del) < 3.5)
			return true;
		else
			return false;
		break;

	case 4:
		if (dis2(del) < 3.5)
			return true;
		else
			return false;
		break;

	case 5:
		if (dis2(del) < 3.5)
			return true;
		else
			return false;
		break;
	}
}

bool RADIAL::hydrogencheck(double* del,int definition)
{
	switch (definition)
	{
	case 1:
		if (dis2(del) < 2.4)
			return true;
		else
			return false;
		break;

	case 3:
		if (dis2(del) < 2.31)
			return true;
		else
			return false;
		break;

	case 4:
		if (dis2(del) < 2.5)
			return true;
		else
			return false;
		break;
	}
}

bool RADIAL::anglecheck(double* vec1, double* vec2, int definition)
{
	double value = (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);

	switch (definition)
	{
	case 1:
		value = value / (dis2(vec1));
		value = value / (dis2(vec2));
		if (value >= sqrt(3) / 2)
			return true;
		else
			return false;
		break;

	case 2:
		value = value / (dis2(vec1));
		value = value / (dis2(vec2));
		if (value >= sqrt(3) / 2)
			return true;
		else
			return false;
		break;

	case 3:
		value = value / (dis2(vec1));
		value = value / (dis2(vec2));

		value = acos(value);
		if (value > PI / 2)
		{
			value = PI - value;
		}
		value = value * 180 / PI;

		value = 7.1 - 0.05*value + 0.00021*value*value;

		double temp;
		temp = dis2(vec1);
		value = value * exp(-(temp / 0.343));

		if (value > 0.0085)
			return true;
		else
			return false;
		break;
	}
}

void RADIAL::crossproduct(double* delta1, double* delta2, double* delta3)
{
	delta3[0] = delta1[1] * delta2[2] - delta1[2] * delta2[1];
	delta3[1] = delta1[2] * delta2[0] - delta1[0] * delta2[2];
	delta3[2] = delta1[0] * delta2[1] - delta1[1] * delta2[0];
}

void RADIAL::hydrogenbond(int definition)
{
	PAGE* page = PLIB->GetPage();
	int Molnum = PLIB->GetPage()->GetNumber();
	Molnum = Molnum / Num_Atom;

	double* Xpos;
	double* Ypos;
	double* Zpos;
	double* BOUNDARY;

	double delta1[3] = { 0, };
	double delta2[3] = { 0, };
	double delta3[3] = { 0, };

	switch (definition)
	{
	case 1:

		for (int t = 0; t < PLIB->GetPagecount(); t++)
		{
			Xpos = page[t].GetPosition(0);
			Ypos = page[t].GetPosition(1);
			Zpos = page[t].GetPosition(2);
			BOUNDARY = page[t].GetBoundary();

			for (int i = 0; i < Molnum; i++)
			{
				for (int j = i + 1; j < Molnum; j++)
				{
					periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * i, delta1);
					if (oxygencheck(delta1,definition))
					{
						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 1, 3 * i, delta2);
						if (hydrogencheck(delta2,definition))
						{
							periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * j + 1, delta2);
							if (anglecheck(delta1, delta2,definition))
							{
								HB[t][i].push_back(j);
								HB[t][j].push_back(i);
							}
						}


						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 2, 3 * i, delta2);
						if (hydrogencheck(delta2,definition))
						{
							periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * j + 2, delta2);
							if (anglecheck(delta1, delta2,definition))
							{
								HB[t][i].push_back(j);
								HB[t][j].push_back(i);
							}
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * j, delta1);
						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 1, 3 * j, delta2);
						if (hydrogencheck(delta2,definition))
						{
							periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * i + 1, delta2);
							if (anglecheck(delta1, delta2,definition))
							{
								HB[t][i].push_back(j);
								HB[t][j].push_back(i);
							}
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 2, 3 * j, delta2);
						if (hydrogencheck(delta2,definition))
						{
							periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * i + 2, delta2);
							if (anglecheck(delta1, delta2,definition))
							{
								HB[t][i].push_back(j);
								HB[t][j].push_back(i);
							}
						}
					}
				}
			}
		}
		break;

	case 2:

		for (int t = 0; t < PLIB->GetPagecount(); t++)
		{
			Xpos = page[t].GetPosition(0);
			Ypos = page[t].GetPosition(1);
			Zpos = page[t].GetPosition(2);
			BOUNDARY = page[t].GetBoundary();

			for (int i = 0; i < Molnum; i++)
			{
				for (int j = i + 1; j < Molnum; j++)
				{
					periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * i, delta1);
					if (oxygencheck(delta1,definition))
					{
						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * j + 1, delta2);
						if (anglecheck(delta1, delta2,definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * j + 2, delta2);
						if (anglecheck(delta1, delta2,definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * j, delta1);
						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * i + 1, delta2);
						if (anglecheck(delta1, delta2,definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * i + 2, delta2);
						if (anglecheck(delta1, delta2,definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}
						
					}
				}
			}
		}
		break;

	case 3:

		for (int t = 0; t < PLIB->GetPagecount(); t++)
		{
			Xpos = page[t].GetPosition(0);
			Ypos = page[t].GetPosition(1);
			Zpos = page[t].GetPosition(2);
			BOUNDARY = page[t].GetBoundary();

			for (int i = 0; i < Molnum; i++)
			{
				periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 1, 3 * i, delta1);
				periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 2, 3 * i, delta2);
				crossproduct(delta1, delta2, delta3);

				for (int j = 0; j < Molnum; j++)
				{
					if (i == j)
					{
						continue;
					}

					periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 1, 3 * i, delta1);
					if (hydrogencheck(delta1, definition))
					{
						if (anglecheck(delta1, delta3, definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}
					}

					periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 2, 3 * i, delta1);
					if (hydrogencheck(delta1, definition))
					{
						if (anglecheck(delta1, delta3, definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}
					}
				}

			}

		}
		break;

	case 4:

		for (int t = 0; t < PLIB->GetPagecount(); t++)
		{
			Xpos = page[t].GetPosition(0);
			Ypos = page[t].GetPosition(1);
			Zpos = page[t].GetPosition(2);
			BOUNDARY = page[t].GetBoundary();

			for (int i = 0; i < Molnum; i++)
			{
				for (int j = i + 1; j < Molnum; j++)
				{
					periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * i, delta1);
					if (oxygencheck(delta1, definition))
					{
						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 1, 3 * i, delta2);
						if (hydrogencheck(delta2, definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 2, 3 * i, delta2);
						if (hydrogencheck(delta2, definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * j, delta1);
						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 1, 3 * j, delta2);
						if (hydrogencheck(delta2, definition))
						{
								HB[t][i].push_back(j);
								HB[t][j].push_back(i);
						}

						periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 2, 3 * j, delta2);
						if (hydrogencheck(delta2, definition))
						{
							HB[t][i].push_back(j);
							HB[t][j].push_back(i);
						}
					}
				}
			}
		}
		break;

	case 5:

		for (int t = 0; t < PLIB->GetPagecount(); t++)
		{
			Xpos = page[t].GetPosition(0);
			Ypos = page[t].GetPosition(1);
			Zpos = page[t].GetPosition(2);
			BOUNDARY = page[t].GetBoundary();

			for (int i = 0; i < Molnum; i++)
			{
				for (int j = i + 1; j < Molnum; j++)
				{
					periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * i, delta1);
					if (oxygencheck(delta1, definition))
					{
						HB[t][i].push_back(j);
						HB[t][j].push_back(i);
					}
				}
			}
		}
		break;

	}
}

vector<int>** RADIAL::GetHydrobond()
{
	return HB;
}

void RADIAL::SOP()
{
	PAGE* page = PLIB->GetPage();
	double* Xpos;
	double* Ypos;
	double* Zpos;
	double* BOUNDARY;

	FILE* pfile;
	fopen_s(&pfile, "Structural_Order_Parameter.txt", "w");

	int Molnum = PLIB->GetPage()->GetNumber() / Num_Atom;
	
	bool Hbond = false;
	double delta_HB = 0;
	double delta = 0;

	size_t size = 0;

	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		Xpos = page[t].GetPosition(0);
		Ypos = page[t].GetPosition(1);
		Zpos = page[t].GetPosition(2);
		BOUNDARY = page[t].GetBoundary();

		for (int i = 0; i < Molnum; i++)
		{
			size = HB[t][i].size();
			if (size == 0)
				continue;

			double temp[3] = { 0, };

			for (int k = 0; k < size; k++)
			{
				periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * HB[t][i][k], temp);
				if (delta_HB == 0)
				{
					delta_HB = dis2(temp);
				}
				else if (dis2(temp) > delta_HB)
				{
					delta_HB = dis2(temp);
				}
			}

			for (int j = 0; j < Molnum; j++)
			{
				for (int k = 0; k < size; k++)
				{
					if (j == HB[t][i][k])
						Hbond = true;

					if (j == i)
						Hbond = true;
				}
				if (Hbond == true)
				{
					Hbond = false;
					continue;
				}

				periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * j, temp);
				if (delta == 0)
				{
					delta = dis2(temp);
				}
				else if (dis2(temp) < delta)
				{
					delta = dis2(temp);
				}
			}

			fprintf(pfile, "%f ", delta - delta_HB);
			delta = 0;
			delta_HB = 0;
		}
	}

	fclose(pfile);
}

void RADIAL::HBlifetime()
{	

		int Molnum = PLIB->GetPage()->GetNumber() / Num_Atom;
		vector<int>* counting;
		counting = new vector<int>[Molnum];
		bool* memo;
		memo = new bool[Molnum];
		memset(memo, 0, sizeof(bool)*Molnum);
		vector<int> temp;
		temp.push_back(-1);
		vector<int>::iterator location;
		int* location_memo = nullptr;

		int length;
		int length_counting;
		int index;
		bool modification = false;

		for (int i = 0; i < Molnum; i++)
		{
			for (int t = 0; t < PLIB->GetPagecount(); t++)
			{
				length = HB[t][i].size();
				length_counting = temp.size();
				location = temp.end() - 1;
				if (length_counting == 1)
				{
					for (int k = 0; k < length; k++)
					{
						index = HB[t][i][k];
						if (index > i)
						{
							temp.push_back(index);
							if (memo[index] == false)
							{
								counting[index].push_back(1);
								memo[index] = true;
							}
						}

					}
				}
				else
				{
					while (*(location) != -1)
					{
						for (int k = 0; k < length; k++)
						{
							index = HB[t][i][k];
							if (index > i && index == *(location))
							{
								modification = true;
								counting[index][counting[index].size() - 1]++;

							}
						}
						if (!modification)
						{
							index = *(location);
							memo[index] = false;
							temp.erase(location);
						}
						location_memo = location._Ptr;
						location_memo--;
						location = temp.end() - 1;
						location._Ptr = location_memo;
						modification = false;
					}

					for (int k = 0; k < length; k++)
					{
						index = HB[t][i][k];
						if (index > i && memo[index] == false)
						{
							memo[index] = true;
							temp.push_back(index);
							counting[index].push_back(1);

						}
					}
				}
			}
			temp.clear();
			temp.shrink_to_fit();
			temp.push_back(-1);
			memset(memo, 0, sizeof(bool)*Molnum);
		}

		FILE* pfile;
		fopen_s(&pfile, "HBlifetime.txt", "w");
		for (int i = 0; i < Molnum; i++)
		{
			length = counting[i].size();
			for (int k = 0; k < length; k++)
			{
				fprintf_s(pfile, "%d ", counting[i][k]);
			}
		}

		fclose(pfile);
		delete[] counting;
		delete[] memo;
}

void RADIAL::HBlifetime2()
{
	double* counting;
	counting = new double[PLIB->GetPagecount()];
	memset(counting, 0, sizeof(double)*(PLIB->GetPagecount()));

	int Molnum = PLIB->GetPage()->GetNumber() / Num_Atom;
	int index;
	int length;
	bool Serial = false;
	vector<int> memo;
	vector<int>::iterator location;

	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		if (memo.size() != 0)
		{
			std::cout << endl;
			std::cout << "error" << endl;
		}
		for (int i = 0; i < Molnum; i++)
		{
			length = HB[t][i].size();
			if (length == 0)
			{
				break;
			}

			for (int k = 0; k < length; k++)
			{
				index = HB[t][i][k];
				if (index > i)
				{
					memo.push_back(HB[t][i][k]);
				}
			}

			while (memo.size() != 0)
			{
				location = memo.end() - 1;
				for (int k = t; k < PLIB->GetPagecount(); k++)
				{
					length = HB[k][i].size();
					Serial = false;
					for (int j = 0; j < length; j++)
					{
						if (*(location) == HB[k][i][j])
						{
							counting[k - t]++;
							Serial = true;
							break;
						}
					}
					if (!Serial || k==(PLIB->GetPagecount()-1))
					{
						memo.erase(location);
						break;
					}
				}
			}
		}
	}

	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		counting[t] /= Molnum;
	}

	FILE *pfile;
	fopen_s(&pfile, "HBlifetime2.txt", "w");
	
	for (int t = 0; t < PLIB->GetPagecount(); t++)
	{
		fprintf_s(pfile, "%f ",counting[t]);
	}
	
	fclose(pfile);
	delete[] counting;
}

PERCOLATION::PERCOLATION(LIBRARY* plib)
{
	PLIB = plib;
}

bool PERCOLATION::hydrogenbondcheck(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int i, int j)
{
	double delta1[3] = { 0, };
	double delta2[3] = { 0, };

	periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * i, delta1);
	if (oxygencheck(delta1,1))
	{
		periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 1, 3 * i, delta2);
		if (hydrogencheck(delta2,1))
		{
			periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * j + 1, delta2);
			if (anglecheck(delta1, delta2,1))
			{
				return true;
			}
		}

		periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j + 2, 3 * i, delta2);
		if (hydrogencheck(delta2,1))
		{
			periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * j, 3 * j + 2, delta2);
			if (anglecheck(delta1, delta2,1))
			{
				return true;
			}
		}

		periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * j, delta1);
		periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 1, 3 * j, delta2);
		if (hydrogencheck(delta2,1))
		{
			periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * i + 1, delta2);
			if (anglecheck(delta1, delta2,1))
			{
				return true;
			}
		}

		periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i + 2, 3 * j, delta2);
		if (hydrogencheck(delta2,1))
		{
			periodic(Xpos, Ypos, Zpos, BOUNDARY, 3 * i, 3 * i + 2, delta2);
			if (anglecheck(delta1, delta2,1))
			{
				return true;
			}
		}
	}

	return false;
}

int PERCOLATION::Get_Percolation(int num_threads)
{
	Set_Num_threads = num_threads;
	set(0, PLIB->GetPagecount());

	return num_Percolation;
}

void PERCOLATION::set(int start, int end)
{
	double* Xpos;
	double* Ypos;
	double* Zpos;
	double* BOUNDARY;
	FILE* pfile;
	PAGE* page = PLIB->GetPage();
	int Molnum = page->GetNumber()/3;

	int* order_number = new int[Molnum];
	bool* memo = new bool[Molnum];

	thread::id ID = this_thread::get_id();

	int count = 0;
	int* number_of_cluster;

	for (int t = start; t < end; t++)
	{
		Xpos = page[t].GetPosition(0);
		Ypos = page[t].GetPosition(1);
		Zpos = page[t].GetPosition(2);
		BOUNDARY = page[t].GetBoundary();
		count = 0;

		memset(memo, 0, sizeof(bool)*Molnum);
		for (int index = 0; index < Molnum; index++)
		{
			order_number[index] = index + 1;
		}

		int memo_count = 0;
		while (memo_count != Molnum)
		{
			if (memo[count] != true)
			{
				memo[count] = true;
				memo_count++;
				for (int k = 0; k < Molnum; k++)
				{
					if (memo[k] != true)
					{
						if ((count != k) && hydrogenbondcheck(Xpos, Ypos, Zpos, BOUNDARY, count, k))
						{
							order_number[k] = order_number[count];
							set_step(Xpos, Ypos, Zpos, BOUNDARY, k, order_number, memo, Molnum,memo_count);
						}
					}
				}
			}
			else
			{
				count++;
			}

		}

		int MAX_order_number = 0;
		for (int index = 0; index < Molnum; index++)
		{
			if (order_number[index] == (index + 1))
			{
				MAX_order_number++;
				order_number[index] = MAX_order_number;
			}
			else
			{
				order_number[index] = order_number[(order_number[index]-1)];
			}
		}
		number_of_cluster = new int[MAX_order_number];
		memset(number_of_cluster, 0, sizeof(int)*MAX_order_number);
		
		for (int index = 0; index < Molnum; index++)
		{
			number_of_cluster[(order_number[index] - 1)]++;
		}

		if (CheckPercolation(Xpos, Ypos, Zpos, BOUNDARY, order_number,number_of_cluster, MAX_order_number,Molnum))
		{
			num_Percolation++;
		}

		delete[] number_of_cluster;

		double progress = (double)(t - start + 1) / (end - start) * 100;
		cout << "Thread " << ID << " has calculated " << progress << endl;
	}

	delete[] memo;
	delete[] order_number;
	return;
}

void PERCOLATION::set_step(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int count, int*& order_number, bool*& memo,int Molnum,int& memo_count)
{
	if (memo[count] != true)
	{
		memo[count] = true;
		memo_count++;
		for (int k = 0; k < Molnum; k++)
		{
			if (memo[k] != true)
			{
				if ((count != k) && hydrogenbondcheck(Xpos, Ypos, Zpos, BOUNDARY, count, k))
				{
					order_number[k] = order_number[count];
					set_step(Xpos, Ypos, Zpos, BOUNDARY, k, order_number, memo, Molnum,memo_count);
				}
			}
		}
	}
}

bool PERCOLATION::CheckPercolation(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int* order_number,int* number_of_cluster, int MAX_order_number,int Molnum)
{
	vector<int> index;
	for (int i = 0; i < MAX_order_number; i++)
	{
		if (number_of_cluster[i] > ((BOUNDARY[1] - BOUNDARY[0]) / 3.6 - 1))
		{
			for (int k = 0; k < Molnum; k++)
			{
				if (order_number[k] == (i + 1))
				{
					index.push_back(k);
				}
			}
			
			
			if (CheckPercolation_step(Xpos, Ypos, Zpos, BOUNDARY, index))
			{
				return true;
			}


			index.clear();
		}
	}

	return false;
}

bool PERCOLATION::CheckPercolation_step(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int> index)
{
	int** M_path;
	size_t number = index.size();
	M_path = new int*[number];

	for (int i = 0; i < number; i++)
	{
		M_path[i] = new int[number];
		memset(M_path[i], 0, sizeof(int)*number);
	}

	int* number_path = new int[number];
	memset(number_path, 0, sizeof(int)*number);

	for (int i = 0; i < number; i++)
	{
		for (int j = i + 1; j < number; j++)
		{
			if (hydrogenbondcheck(Xpos, Ypos, Zpos, BOUNDARY, index[i], index[j]))
			{
				M_path[i][j] = 1;
				M_path[j][i] = 1;
				number_path[i] = number_path[i] + 1;
				number_path[j] = number_path[j] + 1;
			}
		}
	}

	
	if (Determination(M_path, number_path, Xpos, Ypos, Zpos, BOUNDARY, number, index))
	{
		for (int k = 0; k < number; k++)
		{
			delete[] M_path[k];
		}
		delete[] M_path;
		delete[] number_path;
		return true;
	}
	else
	{
		for (int k = 0; k < number; k++)
		{
			delete[] M_path[k];
		}
		delete[] M_path;
		delete[] number_path;
		return false;
	}
}

bool PERCOLATION::Check_remained_path(int** M_Path, int start, int number)
{
	for (int i = 0; i < number; i++)
	{
		if (M_Path[start][i] == 1)
		{
			return true;
		}
	}

	return false;
}

bool PERCOLATION::Check_remained_path(vector<int>** M_N_Path, int start, int number)
{
	for (int i = 0; i < number; i++)
	{
		if (!M_N_Path[start][i].empty())
		{
			return true;
		}
	}

	return false;
}

bool PERCOLATION::Check_remained_path(bool** memo_M_N_Path, int start, int number)
{
	for (int i = 0; i < number; i++)
	{
		if (memo_M_N_Path[start][i])
		{
			return true;
		}
	}

	return false;
}


int PERCOLATION::Check_remained_path(vector<int>** M_N_Path, int start, int number,vector<int>* D_path, bool* memo_Path, vector<int>* node)
{
	int chain_number = -1;
	mutex m;

	for (int i = 0; i < number; i++)
	{
		if (!M_N_Path[start][i].empty())
		{
			for (int k = 0; k < M_N_Path[start][i].size(); k++)
			{
				int temp = M_N_Path[start][i][k];
				if (memo_Path[temp] && !Belong_To_Node(node, D_path, temp))
				{
					chain_number = temp;
					return chain_number;
				}
			}
		}
	}

	return chain_number;
}


bool PERCOLATION::Belong_To_Node(vector<int>* node, vector<int>* D_path, int chain_number)
{
	size_t length = node->size();

	for (int i = 1; i < length-1; i++)
	{
		if (D_path[chain_number].front() == (*node)[i] || D_path[chain_number].back() == (*node)[i])
			return true;
	}
	return false;
}

bool PERCOLATION::Check_Matrix(int** M_path, int number)
{
	for (int i = 0; i < number; i++)
	{
		for (int j = 0; j < number; j++)
		{
			if (M_path[i][j] != 0)
			{
				return true;
			}
		}
	}
	return false;
}

bool PERCOLATION::Check_Matrix(vector<int>** M_N_path, int number)
{
	for (int i = 0; i < number; i++)
	{
		for (int j = 0; j < number; j++)
		{
			if (!M_N_path[i][j].empty())
			{
				return true;
			}
		}
	}
	return false;
}


void PERCOLATION::Trim_Matrix(int** M_path, int* number_path, int number)
{
	vector<int> memo;
	int count = 0;
	while (count != number)
	{
		count = 0;
		for (int i = 0; i < number; i++)
		{
			for (int j = 0; j < number; j++)
			{
				if (M_path[i][j] == 1)
				{
					memo.push_back(j);
				}
			}

			if (memo.size() == 1)
			{
				M_path[i][memo[0]] = 0;
				M_path[memo[0]][i] = 0;
				number_path[i] = number_path[i] - 1;
				number_path[memo[0]] = number_path[memo[0]] - 1;
			}
			else
			{
				count++;
			}
			memo.clear();
		}
	}
}


bool PERCOLATION::Check_loop(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int> *index, vector<int>* memo)
{
	
	size_t temp = memo->size();
	double sum = 0;
	double standard;
	for (int i = 0; i < temp - 1; i++)
	{
		sum = sum + periodic(Xpos[3 * (*index)[(*memo)[i + 1]]], Xpos[3 * (*index)[(*memo)[i]]], BOUNDARY);
	}

	standard = BOUNDARY[1] - BOUNDARY[0];
	if (abs(abs(sum) - standard) <1e-8)
		return true;

	sum = 0;
	for (int i = 0; i < temp - 1; i++)
	{
		sum = sum + periodic(Ypos[3 * (*index)[(*memo)[i + 1]]], Ypos[3 * (*index)[(*memo)[i]]], BOUNDARY);
	}

	standard = BOUNDARY[3] - BOUNDARY[2];
	if (abs(abs(sum) - standard) < 1e-8)
		return true;

	sum = 0;
	for (int i = 0; i < temp - 1; i++)
	{
		sum = sum + periodic(Zpos[3 * (*index)[(*memo)[i + 1]]], Zpos[3 * (*index)[(*memo)[i]]], BOUNDARY);
	}

	standard = BOUNDARY[5] - BOUNDARY[4];
	if (abs(abs(sum) - standard) < 1e-8)
		return true;

	return false;
}

bool PERCOLATION::Determination(int** M_path, int* number_path, double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, int number, vector<int> index)
{
	bool IS_Self = false;
	size_t length;
	int N_chain=0;
	
	vector<int>* D_path;
	
	
	Trim_Matrix(M_path, number_path, number);

	for (int i = 0; i < number; i++)
	{
		if (number_path[i] > 2)
		{
			N_chain += number_path[i];
		}
	}

	if (N_chain == 0)
	{
		N_chain = 1;
	}
	else
	{
		N_chain /= 2;
	}

	D_path = new vector<int>[N_chain];
	//Chain_initialization

	inf_chain* I_D_Path = new inf_chain[N_chain];

	for (int i = 0; i < N_chain; i++)
	{
		memset(I_D_Path[i].Coord, 0, sizeof(double) * 3);
		I_D_Path[i].able = true;
	}
	//Chain_information_initialization


	if (N_chain == 1)
	{
		if (!Check_Matrix(M_path, number))
		{
			return false;
		}

		for (int m = 0; m < N_chain; m++)
		{
			int start;
			for (int i = 0; i < number; i++)
			{
				if (Check_remained_path(M_path, i, number))
				{
					start = i;
					D_path[m].push_back(i);
					break;
				}
			}

			while (true)
			{
				for (int i = 0; i < number; i++)
				{
					if (M_path[start][i] != 0)
					{
						M_path[start][i] = 0;
						M_path[i][start] = 0;
						start = i;
						D_path[m].push_back(i);
						break;
					}
				}

				if (!Check_remained_path(M_path, start, number))
				{
					break;
				}
			}

			if (Check_loop(Xpos, Ypos, Zpos, BOUNDARY, &index, &D_path[0]))
			{
				delete[] D_path;
				delete[] I_D_Path;
				return true;
			}
			else
			{
				delete[] D_path;
				delete[] I_D_Path;
				return false;
			}
		}
	}


	for (int m = 0; m < N_chain; m++)
	{
		int start = -1;
		for (int i = 0; i < number; i++)
		{
			if (number_path[i] > 2 && Check_remained_path(M_path, i, number))
			{
				start = i;
				D_path[m].push_back(i);
				break;
			}
		}

		while (true)
		{
			for (int i = 0; i < number; i++)
			{
				if (M_path[start][i] != 0)
				{
					M_path[start][i] = 0;
					M_path[i][start] = 0;
					start = i;
					D_path[m].push_back(i);
					break;
				}
			}

			if (number_path[start] > 2)
			{
				break;
			}
		}

		Calculate_Vector(Xpos, Ypos, Zpos, BOUNDARY, index, &D_path[m], I_D_Path[m].Coord);

		length = D_path[m].size();
		if (D_path[m][0] == D_path[m][length - 1])
		{
			if (Check_loop(Xpos, Ypos, Zpos, BOUNDARY, &index, &D_path[m]))
			{
				delete[] I_D_Path;
				delete[] D_path;
				return true;
			}
			else
			{
				I_D_Path[m].able = false;
			}
		}
	}

	for (int i = 0; i < N_chain; i++)
	{
		size_t length_x = D_path[i].size();
		for (int j = i + 1; j < N_chain; j++)
		{
			size_t length_y = D_path[j].size();
		
			if ((D_path[i][0] == D_path[j][0] && D_path[i][length_x - 1] == D_path[j][length_y - 1]
				&& Check_vector(I_D_Path[i].Coord, I_D_Path[j].Coord, false) && I_D_Path[j].able == true) ||
				(D_path[i][0] == D_path[j][length_y - 1] && D_path[i][length_x - 1] == D_path[j][0]
					&& Check_vector(I_D_Path[i].Coord, I_D_Path[j].Coord, true) && I_D_Path[j].able == true))
			{
				I_D_Path[j].able = false;
			}
		}
	}

	vector<int>** M_N_Path = new vector<int> * [number];
	for (int i = 0; i < number; i++)
	{
		M_N_Path[i] = new vector<int>[number];
	}

	for (int i = 0; i < N_chain; i++)
	{
		if (I_D_Path[i].able)
		{
			length = D_path[i].size();
			M_N_Path[D_path[i][0]][D_path[i][length - 1]].push_back(i);
			M_N_Path[D_path[i][length - 1]][D_path[i][0]].push_back(i);
		}
	}

	/*for (int i = 0; i < N_chain; i++)
	{
		bool able = true;

		if (!I_D_Path[i].able)
			continue;
		else
		{
			int end = D_path[i].front();
			int count = 0;
			for (int j = 0; j < number; j++)
			{
				if (!M_N_Path[end][j].empty())
				{
					break;
				}
				count = j;
			}

			if (count == number - 1)
			{
				able = false;
			}

			end = D_path[i].back();
			count = 0;
			for (int j = 0; j < number; j++)
			{
				if (!M_N_Path[end][j].empty())
				{
					break;
				}
				count = j;
			}

			if (count == number - 1)
			{
				able = false;
			}

			I_D_Path[i].able = able;


		}
	}*/

	//M_N_Path Finished

	bool flag = false;
	bool** memo_M_N_Path = new bool* [number];
	for (int i = 0; i < number; i++)
	{
		memo_M_N_Path[i] = new bool[number];
		for (int j = 0; j < number; j++)
		{
			memo_M_N_Path[i][j] = false;
		}
	}

	for (int i = 0; i < number; i++)
	{
		for (int j = 0; j < number; j++)
		{
			if (!M_N_Path[i][j].empty())
			{
				memo_M_N_Path[i][j] = true;
			}
		}
	}

	Set_Num_threads = N_chain/2;
	if (Set_Num_threads > 20)
		Set_Num_threads = 20;

	queue<int>* chain_exclusion = new queue<int>[Set_Num_threads];

	vector<thread> workers;
	Num_threads = Set_Num_threads;
	workers.reserve(Set_Num_threads);

	for (int i = 0; i < Set_Num_threads; i++)
	{
		int address = (int)&chain_exclusion[i];
		workers.emplace_back([&, address]() 
			{this->Chain(Xpos,Ypos,Zpos,BOUNDARY,&index,M_N_Path,memo_M_N_Path, number, D_path, N_chain,I_D_Path,&flag,chain_exclusion,address); });
	}

	for (int i = 0; i < Set_Num_threads; i++)
	{
		workers[i].join();
	}

	if (flag)
	{
		for (int i = 0; i < Set_Num_threads; i++)
		{
			queue<int> dumping;
			chain_exclusion[i].swap(dumping);
			while (!dumping.empty())
			{
				dumping.pop();
			}
		}


		for (int i = 0; i < number; i++)
		{
			delete[] M_N_Path[i];
			delete[] memo_M_N_Path[i];
		}
		delete[] chain_exclusion;
		printf("%Delete Queue\n");
		delete[] M_N_Path;
		delete[] D_path;
		delete[] I_D_Path;
		return true;
	}
	else
	{
		for (int i = 0; i < Set_Num_threads; i++)
		{
			queue<int> dumping;
			chain_exclusion[i].swap(dumping);
			while (!dumping.empty())
			{
				dumping.pop();
			}
		}


		for (int i = 0; i < number; i++)
		{
			delete[] M_N_Path[i];
			delete[] memo_M_N_Path[i];
		}
		delete[] chain_exclusion;
		printf("%Delete Queue\n");
		delete[] M_N_Path;
		delete[] D_path;
		delete[] I_D_Path;
		return false;
	}
}

void PERCOLATION::Calculate_Vector(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int> index, vector<int>* memo,double* Coord)
{
	size_t temp = memo->size();
	double sum = 0;
	for (int i = 0; i < temp - 1; i++)
	{
		sum = sum + periodic(Xpos[3 * index[(*memo)[i + 1]]], Xpos[3 * index[(*memo)[i]]], BOUNDARY);
	}

	Coord[0] = sum;

	sum = 0;
	for (int i = 0; i < temp - 1; i++)
	{
		sum = sum + periodic(Ypos[3 * index[(*memo)[i + 1]]], Ypos[3 * index[(*memo)[i]]], BOUNDARY);
	}

	Coord[1] = sum;

	sum = 0;
	for (int i = 0; i < temp - 1; i++)
	{
		sum = sum + periodic(Zpos[3 * index[(*memo)[i + 1]]], Zpos[3 * index[(*memo)[i]]], BOUNDARY);
	}

	Coord[2] = sum;

	return;
}

bool PERCOLATION::Check_Inner(double* memo, double* Coord, bool reverse)
{
	double sum = 0;
	for (int i = 0; i < 3; i++)
	{
		sum += memo[i] * Coord[i];
	}

	if (reverse)
	{
		sum = -sum;
	}

	if (sum > 0)
		return true;
	else
		return false;
}

bool PERCOLATION::Check_vector(double* memo, double* Coord, bool reverse)
{
	if (!reverse)
	{
		if ( abs(memo[0]-Coord[0])<1e-8 && abs(memo[1]-Coord[1])<1e-8 && abs(memo[2]-Coord[2])<1e-8)
			return true;
		else
			return false;
	}
	else
	{
		if ( abs(memo[0]+Coord[0])<1e-8 && abs(memo[1]+Coord[1])<1e-8 && abs(memo[2]+Coord[2])<1e-8)
			return true;
		else
			return false;
	}
}

void PERCOLATION::Chain(double* Xpos, double* Ypos, double* Zpos, double* BOUNDARY, vector<int>* index, vector<int>** M_N_Path, bool** memo_M_N_Path, int number, vector<int>* D_path, int N_chain, inf_chain* I_D_Path, bool* flag, queue<int>* master_job, int address)
{
	queue<int>* job = (queue<int>*)address;
	vector<int> memo, chain_order;
	vector<vector<int>> trash;
	size_t trash_length=trash.size();
	vector<int> temp_vector;
	vector<int> node;
	size_t pre_chain = 0;
	//mutex m_lock;
	bool Initialization = true;

	printf("%p Thread Start, %d : Thread Job Address\n", std::this_thread::get_id(), address);

	vector<int>** M_N_Path_Copy = new vector<int>*[number];
	for (int i = 0; i < number; i++)
	{
		M_N_Path_Copy[i] = new vector<int>[number];
		for (int j = 0; j < number; j++)
		{
			M_N_Path_Copy[i][j] = M_N_Path[i][j];
		}
	}


	bool* memo_Path = new bool[N_chain];
	memset(memo_Path, false, sizeof(bool) * N_chain);
	for (int i = 0; i < N_chain; i++)
	{
		if (I_D_Path[i].able)
		{
			memo_Path[i] = true;
		}
	}

	while (Check_Matrix(M_N_Path_Copy, number))
	{
		if (*flag)
		{
			for (int x = 0; x < number; x++)
			{
				for (int y = 0; y < number; y++)
				{
					vector<int> dumping;
					M_N_Path_Copy[x][y].swap(dumping);
					dumping.clear();
					M_N_Path_Copy[x][y].clear();
				}
			}

			for (int t = 0; t < number; t++)
			{
				delete[] M_N_Path_Copy[t];
			}
			delete[] M_N_Path_Copy;
			printf("%Delete M_N_Path_Copy\n");
			delete[] memo_Path;
			return;
		}

		job_m.lock();
		bool job_empty = job->empty();
		job_m.unlock();

		if (!job_empty)
		{

			printf("%p Thread Start, %d : Job Start\n", std::this_thread::get_id(), address);
			
			bool job_removed = false;

			job_m.lock();
			int exclused_index = job->front();
			job_m.unlock();

			M_N_Path_Copy[D_path[exclused_index].front()][D_path[exclused_index].back()].pop_back();
			M_N_Path_Copy[D_path[exclused_index].back()][D_path[exclused_index].front()].pop_back();

			//printf("Job %p: M_N_Path_Copy[%d][%d] pop_back\n", std::this_thread::get_id(), D_path[exclused_index].front(), D_path[exclused_index].back());
			//printf("Job %p: M_N_Path_Copy[%d][%d] pop_back\n", std::this_thread::get_id(), D_path[exclused_index].back(), D_path[exclused_index].front());


			size_t chain_length = chain_order.size();
			for (int i = 0; i < chain_length; i++)
			{
				if (chain_order[i] == exclused_index)
				{
					size_t length;
					while (chain_order.back() != exclused_index)
					{
						length = memo.size();
						while ((memo)[length - 2] != (memo)[length - 1])
						{
							memo.pop_back();
							length = memo.size();
						}
						memo.pop_back();
						node.pop_back();
						chain_order.pop_back();

						trash_length = trash.size();
						if (pre_chain - chain_order.size() == 2)
						{
							size_t length_t = trash[trash_length - 1].size();
							for (int k = 0; k < length_t; k++)
							{
								int index_t = trash[trash_length - 1][k];
								memo_Path[index_t] = true;
							}
							trash.pop_back();
							pre_chain--;
						}
					}

					length = memo.size();
					while ((memo)[length - 2] != (memo)[length - 1])
					{
						memo.pop_back();
						length = memo.size();
					}
					memo.pop_back();
					node.pop_back();
					chain_order.pop_back();

					trash_length = trash.size();
					if (pre_chain - chain_order.size() == 2)
					{
						size_t length_t = trash[trash_length - 1].size();
						for (int k = 0; k < length_t; k++)
						{
							int index_t = trash[trash_length - 1][k];
							memo_Path[index_t] = true;
						}
						trash.pop_back();
						pre_chain--;
					}

					job_m.lock();
					job->pop();
					job_m.unlock();
					job_removed = true;

					break;
				}
			}

			if (!job_removed)
			{
				job_m.lock();
				job->pop();
				job_m.unlock();
			}
		}


		temp_vector.clear();

		if (memo.size() == 0 && Initialization)
		{
			chain_order.clear();
			memo.clear();
			trash.clear();

			int start;
			size_t length;
			int pre_index;

			for (int i = 0; i < N_chain; i++)
			{
				if (I_D_Path[i].able)
				{
					memo_Path[i] = true;
				}
			}


			unique_lock<mutex> lk(M);
			printf("%p : CV2_IN\n", std::this_thread::get_id());
			cv2.wait(lk, [&] {return Check_path || Num_threads == 1; });
			printf("%p : CV2_OUT\n", std::this_thread::get_id());
			Check_path = false;
			printf("%p : Check_path\n", std::this_thread::get_id());

			start = -1;
			for (int i = 0; i < number; i++)
			{
				if (Check_remained_path(memo_M_N_Path, i, number))
				{
					start = i;
					break;
				}

				if (i == number - 1)
				{
					delete[] memo_Path;
					Num_threads--;
					Check_path = true;
					cv2.notify_one();
					printf("%p : Check_path Return\n", std::this_thread::get_id());

					for (int x = 0; x < number; x++)
					{
						for (int y = 0; y < number; y++)
						{
							vector<int> dumping;
							M_N_Path_Copy[x][y].swap(dumping);
							dumping.clear();
							M_N_Path_Copy[x][y].clear();
						}
					}

					for (int t = 0; t < number; t++)
					{
						delete[] M_N_Path_Copy[t];
					}
					delete[] M_N_Path_Copy;
					printf("%Delete M_N_Path_Copy\n");
					return;
				}
			}

			for (int j = 0; j < number; j++)
			{
				if (memo_M_N_Path[start][j])
				{
					job_m.lock();
					memo_M_N_Path[start][j] = false;
					memo_M_N_Path[j][start] = false;
					job_m.unlock();

					pre_index = M_N_Path_Copy[start][j].back();

					temp_vector.push_back(pre_index);
					trash.push_back(temp_vector);
					chain_order.push_back(pre_index);
					
					length = D_path[pre_index].size();
					for (int k = 0; k < length; k++)
					{
						memo.push_back(D_path[pre_index][k]);
					}
					memo_Path[pre_index] = false;

					break;
				}
			}
			Check_path = true;
			cv2.notify_one();

			node.clear();
			node.push_back(memo.front());
			node.push_back(memo.back());
			pre_chain = 1;

			this_thread::sleep_for(50ms);
		}
		else if (memo.size() == 0 && !Initialization)
		{
			int pre_index = M_N_Path_Copy[node.front()][node.back()].back();

			for (int i = 0; i < N_chain; i++)
			{
				if (I_D_Path[i].able)
				{
					memo_Path[i] = true;
				}
			}

			if (D_path[pre_index].front() == node.front())
			{
				size_t length = D_path[pre_index].size();
				for (int k = 0; k < length; k++)
				{
					memo.push_back(D_path[pre_index][k]);
				}
			}
			else
			{
				size_t length = D_path[pre_index].size();
				for (int k = 0; k < length; k++)
				{
					memo.push_back(D_path[pre_index][length - 1 - k]);
				}
			}
			memo_Path[pre_index] = false;
			
			chain_order.push_back(pre_index);

			trash_length = trash.size();
			trash[trash_length - 1].push_back(pre_index);

		}
		else
		{
			size_t length = memo.size();
			int chain_number = -1;

			int chain_index = memo.back();
			chain_number = Check_remained_path(M_N_Path_Copy, chain_index, number, D_path, memo_Path, &node);

			if (chain_number == -1)
			{
				if (chain_order.size() == 1)
				{
					int index_t = chain_order.front();

					unique_lock<mutex> lk(M);
					printf("%p : CV2_IN\n", std::this_thread::get_id());
					cv2.wait(lk, [&] {return Check_path || Num_threads == 1; });
					Check_path = false;
					printf("%p : CV2_OUT\n", std::this_thread::get_id());
					printf("%p : Remove_path\n", std::this_thread::get_id());

					I_D_Path[index_t].able = false;
					Push_job_into_Master(master_job, job, index_t);
					M_N_Path_Copy[memo.front()][memo.back()].pop_back();
					M_N_Path_Copy[memo.back()][memo.front()].pop_back();
					printf("Remove %p: M_N_Path_Copy[%d][%d] pop_back\n", std::this_thread::get_id(), memo.front(), memo.back());
					printf("Remove %p: M_N_Path_Copy[%d][%d] pop_back\n", std::this_thread::get_id(), memo.back(), memo.front());
					printf("%p : %d , %d chain out\n", memo.front(), memo.back());


					Check_path = true;
					cv2.notify_one();




					if (!M_N_Path_Copy[memo.front()][memo.back()].empty())
					{
						Initialization = false;
					}
					else
					{
						Initialization = true;
					}

					memo.clear();
					chain_order.pop_back();
				}
				else
				{
					//printf("%p Thread Start, %d : Chain Dissembly\n", std::this_thread::get_id(), address);
					while ((memo)[length - 2] != (memo)[length - 1])
					{
						memo.pop_back();
						length = memo.size();
					}
					memo.pop_back();
					node.pop_back();
					chain_order.pop_back();
				}



				trash_length = trash.size();
				if (pre_chain - chain_order.size() == 2)
				{
					size_t length_t = trash[trash_length - 1].size();
					for (int k = 0; k < length_t; k++)
					{
						int index_t = trash[trash_length - 1][k];
						memo_Path[index_t] = true;
					}
					trash.pop_back();
					pre_chain--;
				}

			}
			else
			{
				//printf("%p Thread Start, %d : Chain Assembly\n", std::this_thread::get_id(), address);
				length = D_path[chain_number].size();;
				if (memo.back() == D_path[chain_number].front())
				{
					for (int k = 0; k < length; k++)
					{
						memo.push_back(D_path[chain_number][k]);
					}
				}
				else if (memo.back() == D_path[chain_number].back())
				{
					for (int k = 0; k < length; k++)
					{
						memo.push_back(D_path[chain_number][length - 1 - k]);
					}
				}

				memo_Path[chain_number] = false;
				chain_order.push_back(chain_number);
				node.push_back(memo.back());

				trash_length = trash.size();

				if (pre_chain == chain_order.size())
				{
					trash[trash_length - 1].push_back(chain_number);
				}
				else
				{
					temp_vector.clear();
					temp_vector.push_back(chain_number);
					trash.push_back(temp_vector);
					pre_chain = chain_order.size();
				}
			}

			if (chain_order.size() != 0)
			{
				if (memo.back() == memo.front())
				{
					if (Check_loop(Xpos, Ypos, Zpos, BOUNDARY, index, &memo))
					{
						job_m.lock();
						*flag = true;
						job_m.unlock();

						size_t memo_length = memo.size();
						for (int k = 0; k < memo_length; k++)
						{
							cout << memo[k] << " ";
						}
						cout << endl;

						for (int x = 0; x < number; x++)
						{
							for (int y = 0; y < number; y++)
							{
								vector<int> dumping;
								M_N_Path_Copy[x][y].swap(dumping);
								dumping.clear();
								M_N_Path_Copy[x][y].clear();
							}
						}


						delete[] memo_Path;
						for (int t = 0; t < number; t++)
						{
							delete[] M_N_Path_Copy[t];
						}
						delete[] M_N_Path_Copy;
						printf("%Delete M_N_Path_Copy\n");
						return;
					}
					else
					{
						length = memo.size();
						while (memo[length - 2] != memo[length - 1])
						{
							memo.pop_back();
							length = memo.size();
						}
						memo.pop_back();
						node.pop_back();
						chain_order.pop_back();
					}
				}

			}


			//Debugging

			/*size_t memo_length = memo.size();
			for (int k = 0; k < memo_length; k++)
			{
				cout << memo[k] << " ";
			}
			cout << endl;*/
		}
	}
	
	Num_threads--;
	delete[] memo_Path;
	printf("%p : Matrix Return\n", std::this_thread::get_id());

	for (int x = 0; x < number; x++)
	{
		for (int y = 0; y < number; y++)
		{
			vector<int> dumping;
			M_N_Path_Copy[x][y].swap(dumping);
			dumping.clear();
			M_N_Path_Copy[x][y].clear();
		}
	}

	for (int t = 0; t < number; t++)
	{
		delete[] M_N_Path_Copy[t];
	}
	delete[] M_N_Path_Copy;
	printf("%Delete M_N_Path_Copy\n");
	return;
}

void PERCOLATION::Push_job_into_Master(queue<int>* master_job, queue<int>* job, int index)
{
	for (int i = 0; i < Set_Num_threads; i++)
	{
		if (job == &master_job[i])
		{
			continue;
		}
		job_m.lock();
		master_job[i].push(index);
		job_m.unlock();
	}

	return;
}