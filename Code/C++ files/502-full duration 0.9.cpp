// 502Test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include <conio.h>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include <cstring>
#include <string.h>
#include <math.h>
#include <stdio.h>
using namespace std;
#define x 63
#define y 63
int t = 832, landmissingcount = 0;//x-number of rows, y-number of columns, t-no of weeks. Right now, initialized to small data set.
int degree[63][63];
int node[63][63];
int node1d[3969];
double totalcluster = 0;
int edge1d[3969][3969];
double Sxx[63][63];
double cluster[63][63];
double Sxy;
double totaldegree = 0;
double avgdegree;
double corelation;
double degreedistribution[3969];
int mncount, nodecount = 0;
int tempcount = 0;
double kv, ev;
int mtemp[3969], ntemp[3969], x3, x4, x5, x6;
int edge[63][63][63][63];

double data1[63][63][832];//Declare array to hold data set.
double sum, mean[63][63];//Contains mean of all cells in the given time frame - i.e for 832 weeks.
int dist[x*y];     // The output array.  dist[i] will hold the shortest
				   // distance from src to i

bool sptSet[x*y]; // sptSet[i] will true if vertex i is included in shortest
				  // path tree or shortest distance from src to i is finalized

char filename[] = "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw01y1990+landmask";

void filenameupdate(int);
int dijkstra(int src);
int dijkstralag1(int src);
int minDistance(int dist[], bool sptSet[]);
int minDistancelag1(int dist[], bool sptSet[]);

int main()
{
	char val[4];
	int temp[32], bindata[32];
	for (int week = 0; week < t; week++)
	{

		int bindigit = 0;
		int row = 0, column = 0;
		filenameupdate(week);

		ifstream file(filename, ios::binary | ios::in);

		if (!file.is_open()) {
			cout << "Error while opening the file" << filename << endl;
		}

		for (int noofvalues = 0; noofvalues < x*y; noofvalues++)
		{
			bindigit = 0;
			//For loop for reading 4 bytes together
			for (int i = 0; i < 4; i++)
			{
				file.get(val[i]);
				for (int k = 0; k < 8; k++)
				{

					temp[k] = val[i] % 2;
					if (temp[k] == -1)
						temp[k] = (-1)*temp[k];
					val[i] = val[i] >> 1;

				}
				for (int k = 0; k <= 7; k++)
				{
					bindata[bindigit] = temp[k];
					bindigit++;

				}
			}

			/*cout << '\n';
			for (int e = 31; e >= 0; e--)
			cout << bindata[e];

			cout << '\n';
			*/

			//Converting the read 4 bytes to value
			//Finding Mantissa
			double mantissa = 1;
			for (int m = 1; m <= 23; m++)
			{
				mantissa = mantissa + (bindata[23 - m] * pow(2, -m));
			}
			//cout << "mantissa" << mantissa << '\n';

			//Finding Exponent
			double exponenttemp = 0;
			for (int n = 23, l = 0; n <= 30; n++, l++)
				exponenttemp = exponenttemp + (bindata[n] * pow(2, l));

			//cout << "Exponent decimal value" << exponenttemp<<'\n';
			double exponent = pow(2, (exponenttemp - 127));

			//Finding value for these 4 bytes that were read 
			double value = mantissa*exponent;
			if (bindata[31] == 1)
				value = value*(-1);
			//cout << "value" << value << '/n';
			//cout << "row " << row << " column " << column<<"Week"<<week<<'\n';
			data1[row][column][week] = value;//Storing the value in data array
											 //cout <<data1[row][column][week]<<'\n';
											 //if (week==0 && value != 168 && value != 157)
											 //tempcount++;

			if (column == y - 1 && row == x - 1)
				break;
			else if (column == y - 1 && row < x - 1)
			{
				column = 0;
				row++;
			}
			else if (column < y - 1)
				column++;

		}
		/*	cout << "The data values stored in file " << week << "\n";


		for (row = 0; row < x; row++)
		{
		for (column = 0; column < y; column++)
		cout << data1[row][column][week] << " ";
		cout << '\n';

		}
		cout << "\n\n";  */


		file.close();
	}
	//cin.get();
	/*--------------------------End of file read-----------------------------------------------------------------------------------------------------------*/


	int i, j, k, m, n;

	//cout << "Hello World";
	//cout << "\nTemp count = " << tempcount;
	/*Three dimensional array of dimensions 63*63*832
	63*63 cells for each week. 832 weeks in 16 years. All of this corresponds to the small data set.*/

	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++)
		{
			mean[i][j] = 0;
			Sxx[i][j] = 0;
			degree[i][j] = 0;
			node[i][j] = 0;
			cluster[i][j] = 0;
			degreedistribution[(i*x) + j] = 0;
		}
	//To find the mean of the data
	for (i = 0; i < x; i++)
	{
		for (j = 0; j < y; j++)
		{
			landmissingcount = 0;
			sum = 0;
			for (k = 0; k < t; k++)
			{
				if (data1[i][j][k] == 168.00 || data1[i][j][k] == 157.00)
					landmissingcount++;
				else
					sum = sum + data1[i][j][k];
			}
			//if (i == 0 && j == 0)
			//cout << "Sum of the very first element" << sum;
			if (landmissingcount != t)
				mean[i][j] = sum / (t - landmissingcount);
			//printing
			//cout <<"Sum= "<<sum << "mean= "<<mean[i][j] << ' ';

		}
		//cout << '\n';
	}

	//To find Sxx - Sxx is a 2D array of size x*y 
	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++)
		{
			for (k = 0; k < t; k++)
			{
				if (!(data1[i][j][k] == 168.00 || data1[i][j][k] == 157.00))
					Sxx[i][j] = Sxx[i][j] + (pow((data1[i][j][k]-mean[i][j]), 2.0));
			}

		}
	//
	//To find Sxy, Corelation & Edge 
	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++)
			for (m = 0; m < x; m++)
				for (n = 0; n < y; n++)
				{
					Sxy = 0;
					for (k = 0; k < t; k++)
					{
						if (!(data1[i][j][k] == 168.00 || data1[i][j][k] == 157.00 || data1[m][n][k] == 168.00 || data1[m][n][k] == 157.00))
							Sxy = Sxy + ((data1[i][j][k]-mean[i][j])*(data1[m][n][k])-mean[i][j]); //To find Sxy temporarily for 2 cells
					}
					if (Sxx[i][j] != 0 && Sxx[m][n] != 0)
						corelation = (Sxy / sqrt(Sxx[i][j] * Sxx[m][n])); //To find corelation coefficient temporarily for 2 cells
					else
						corelation = 0;
					//if(i==0&&j==57&&m==2)
					//cout << corelation << " is corelation between ( " << i << " , " << j << " ) and " << " ( " << m << " , " << n << " )\n";

					if ((abs(corelation) >= 0.9) && !((i == m) && (j == n)))
					{
						edge[i][j][m][n] = 1;//Edge between cells
						edge[m][n][i][j] = 1;
						node[i][j] = 1;  
						node[m][n] = 1;
						//if (i == 0)
						//{
						//cout << "Element at position ( " << i << " \, " << j << " ) is a node \n";
						//cout << "Element at position ( " << m << " , " << n << " ) is a node \n";
						//cout << "There is edge between the above nodes \n\n";
						//}
						
					}

				}

	//To find the degree & clustering coefficient of each node in graph

	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++)
		{
			mncount = 0;
			kv = 0;
			ev = 0;


			for (m = 0; m < x; m++)
			{
				for (n = 0; n < y; n++)
					if (edge[i][j][m][n] == 1 && node[i][j]==1&&!((i == m) && (j == n)))
					{

						degree[i][j]++;
						mtemp[mncount] = m;
						ntemp[mncount] = n;
						mncount++;
						//cout << i << j << m << n << " is a edge \n";
					}
			}


			kv = degree[i][j] + 1;  //Number of vertices in N(v) can be given by degree of an element+1 inclusive of the element itself.
			totaldegree += degree[i][j];
			degreedistribution[degree[i][j]]++;
		//if (i == 1)
		//cout << "Degree of (" << i << ", " << j << " ) " << degree[i][j] << " kv = "<<kv;

			for (int x1 = 0; x1 < mncount; x1++)
				for (int x2 = x1 + 1; x2 < mncount; x2++)
				{
					x3 = mtemp[x1];
					x4 = ntemp[x1];
					x5 = mtemp[x2];
					x6 = ntemp[x2];
					if (edge[x3][x4][x5][x6] == 1)
						ev++;

				}

			ev = ev + degree[i][j];

			if (node[i][j] == 1)
			{
				cluster[i][j] = (2 * ev) / (kv*(kv - 1));
				totalcluster += cluster[i][j];
				//cout << "Clustering coeff : "<< cluster[i][j];
			}
			
			

				//if(i==1)
				//cout << "Clustering coefficient of (" << i << ", " << j << " ) = " << cluster[i][j]<<" ev = "<< ev <<" node of (i,j)= "<<node[i][j]<<"\n";


		}

	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
			if (node[i][j] == 1)
				nodecount++;
				
	//Average number of degree per node
	//cout << "Total degree : " << totaldegree;
	avgdegree = totaldegree / nodecount;
	//cout << "\nTotal cluster = " << totalcluster;
	avgdegree = totaldegree / nodecount;
	

	//cout << "Total connected nodes: " << nodecount;

	//cout << "Total connected nodes: " << nodecount;
	cout << "\nAverage degree of the nodes in the graph : " << avgdegree;



	totalcluster = totalcluster / nodecount;
	cout << "\nClustering coefficient of the graph: " << totalcluster;

	//Random graph metrices
  double totalclusterrandom = avgdegree / nodecount;
	cout << "\n Clustering coefficient of a random graph: " << totalclusterrandom;
	//nodecount = 0;
	//Converting edges and nodes to 1d...as a sequence of X * Y elements
	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
		{
			if (node[i][j] == 1)
			{
				node1d[i*x + j] = 1;
				//nodecount++;
			}
			for (int m = 0; m < x; m++)
				for (int n = 0; n < y; n++)
					if (edge[i][j][m][n] == 1)
						edge1d[i*x + j][m*x + n] = 1;
		}
	//cout << "Node count in 1d = " << nodecount;
	
	//Characteristic path length
	int sumofdij = 0;
	for (int i = 0; i < x*y; i++)
	{
	if (node1d[i] == 1)
	{
	int dij = dijkstra(i);
	sumofdij = sumofdij + dij;
	}
	}
	//cout << "Sum of dij : " << sumofdij;
	double charpathlength = sumofdij / ((nodecount)*(nodecount - 1));

	ofstream degreewrite;

	degreewrite.open("degreefor09.txt", ios::app | ios::out);
	for (int i = 0; i < 3969; i++)
	{
		degreewrite << degreedistribution[i] << " ";
		//outfile << "\n";
	}
	degreewrite.close();

	cout << "\nDegree distribution array written on a file...";
	


	cout << "\nCharacteristic path length : " << charpathlength;
	double charpathrandom = log(nodecount) / log(avgdegree);
	cout << "\nCharacteristic path length of random graph: " << charpathrandom;
	//Random Graph
	//*/

	
	
	//cout << "Degree of (0,16)" << degree[0][16];
	_getch();
	return 0;
}

// Djikstra shortest path for finding a shortest path from one node to every other node
int dijkstra(int src)
{

	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < x*y; i++)
		dist[i] = INT_MAX, sptSet[i] = false;

	// Distance of source vertex from itself is always 0
	dist[src] = 0;

	// Find shortest path for all vertices
	for (int count = 0; count < (x*y); count++)
	{

		if (node1d[count] == 1)
		{
			// Pick the minimum distance vertex from the set of vertices not
			// yet processed. u is always equal to src in first iteration.
			int u = minDistance(dist, sptSet);

			// Mark the picked vertex as processed
			sptSet[u] = true;

			// Update dist value of the adjacent vertices of the picked vertex.
			for (int v = 0; v < (x*y); v++)

				// Update dist[v] only if is not in sptSet, there is an edge from 
				// u to v, and total weight of path from src to  v through u is 
				// smaller than current value of dist[v]
				if (!sptSet[v] && edge1d[u][v] && dist[u] != INT_MAX
					&& dist[u] + edge1d[u][v] < dist[v])
					dist[v] = dist[u] + edge1d[u][v];
		}
	}


	//here dij is sum of distance from one vertex (i) to every other vertex
	int dij = 0;
	for (int m = 0; m<(x*y); m++)
	{
		if (dist[m] != INT_MAX)
			dij = dij + dist[m];

	}
	return dij;

	
}

int minDistance(int dist[], bool sptSet[])
{
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < (x*y); v++)
		if (node1d[v] == 1 && sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}


void filenameupdate(int filenumber)
{
	if (filenumber == 0)			strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw01y1990+landmask");

	if (filenumber == 1)			strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw02y1990+landmask");

	if (filenumber == 2)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw03y1990+landmask");

	if (filenumber == 3)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw04y1990+landmask");

	if (filenumber == 4)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw05y1990+landmask");

	if (filenumber == 5)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw06y1990+landmask");

	if (filenumber == 6)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw07y1990+landmask");

	if (filenumber == 7)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw08y1990+landmask");

	if (filenumber == 8)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw09y1990+landmask");

	if (filenumber == 9)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw10y1990+landmask");

	if (filenumber == 10)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw11y1990+landmask");

	if (filenumber == 11)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw12y1990+landmask");

	if (filenumber == 12)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw13y1990+landmask");

	if (filenumber == 13)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw14y1990+landmask");

	if (filenumber == 14)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw15y1990+landmask");

	if (filenumber == 15)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw16y1990+landmask");

	if (filenumber == 16)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw17y1990+landmask");

	if (filenumber == 17)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw18y1990+landmask");

	if (filenumber == 18)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw19y1990+landmask");

	if (filenumber == 19)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw20y1990+landmask");

	if (filenumber == 20)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw21y1990+landmask");	 if (filenumber == 21)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw22y1990+landmask");

	if (filenumber == 22)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw23y1990+landmask");

	if (filenumber == 23)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw24y1990+landmask");

	if (filenumber == 24)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw25y1990+landmask");

	if (filenumber == 25)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw26y1990+landmask");

	if (filenumber == 26)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw27y1990+landmask");

	if (filenumber == 27)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw28y1990+landmask");

	if (filenumber == 28)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw29y1990+landmask");

	if (filenumber == 29)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw30y1990+landmask");

	if (filenumber == 30)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw31y1990+landmask");	 if (filenumber == 31)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw32y1990+landmask");

	if (filenumber == 32)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw33y1990+landmask");

	if (filenumber == 33)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw34y1990+landmask");

	if (filenumber == 34)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw35y1990+landmask");

	if (filenumber == 35)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw36y1990+landmask");

	if (filenumber == 36)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw37y1990+landmask");

	if (filenumber == 37)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw38y1990+landmask");

	if (filenumber == 38)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw39y1990+landmask");

	if (filenumber == 39)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw40y1990+landmask");

	if (filenumber == 40)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw41y1990+landmask");	 if (filenumber == 41)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw42y1990+landmask");

	if (filenumber == 42)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw43y1990+landmask");

	if (filenumber == 43)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw44y1990+landmask");

	if (filenumber == 44)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw45y1990+landmask");

	if (filenumber == 45)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw46y1990+landmask");

	if (filenumber == 46)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw47y1990+landmask");

	if (filenumber == 47)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw48y1990+landmask");

	if (filenumber == 48)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw49y1990+landmask");

	if (filenumber == 49)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw50y1990+landmask");

	if (filenumber == 50)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw51y1990+landmask");	 if (filenumber == 51)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1990\\Beaufort_Sea_diffw52y1990+landmask");

	//----------------------------------------------------------1991----------------------------------------------------------------------------------

	if (filenumber == 52)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw01y1991+landmask");

	if (filenumber == 53)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw02y1991+landmask");

	if (filenumber == 54)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw03y1991+landmask");

	if (filenumber == 55)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw04y1991+landmask");

	if (filenumber == 56)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw05y1991+landmask");

	if (filenumber == 57)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw06y1991+landmask");

	if (filenumber == 58)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw07y1991+landmask");

	if (filenumber == 59)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw08y1991+landmask");

	if (filenumber == 60)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw09y1991+landmask");	 if (filenumber == 61)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw10y1991+landmask");

	if (filenumber == 62)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw11y1991+landmask");

	if (filenumber == 63)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw12y1991+landmask");

	if (filenumber == 64)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw13y1991+landmask");

	if (filenumber == 65)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw14y1991+landmask");

	if (filenumber == 66)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw15y1991+landmask");

	if (filenumber == 67)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw16y1991+landmask");

	if (filenumber == 68)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw17y1991+landmask");

	if (filenumber == 69)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw18y1991+landmask");

	if (filenumber == 70)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw19y1991+landmask");	 if (filenumber == 71)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw20y1991+landmask");

	if (filenumber == 72)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw21y1991+landmask");

	if (filenumber == 73)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw22y1991+landmask");

	if (filenumber == 74)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw23y1991+landmask");

	if (filenumber == 75)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw24y1991+landmask");

	if (filenumber == 76)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw25y1991+landmask");

	if (filenumber == 77)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw26y1991+landmask");

	if (filenumber == 78)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw27y1991+landmask");

	if (filenumber == 79)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw28y1991+landmask");

	if (filenumber == 80)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw29y1991+landmask");	 if (filenumber == 81)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw30y1991+landmask");

	if (filenumber == 82)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw31y1991+landmask");

	if (filenumber == 83)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw32y1991+landmask");

	if (filenumber == 84)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw33y1991+landmask");

	if (filenumber == 85)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw34y1991+landmask");

	if (filenumber == 86)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw35y1991+landmask");

	if (filenumber == 87)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw36y1991+landmask");

	if (filenumber == 88)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw37y1991+landmask");

	if (filenumber == 89)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw38y1991+landmask");

	if (filenumber == 90)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw39y1991+landmask");	 if (filenumber == 91)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw40y1991+landmask");

	if (filenumber == 92)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw41y1991+landmask");

	if (filenumber == 93)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw42y1991+landmask");

	if (filenumber == 94)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw43y1991+landmask");

	if (filenumber == 95)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw44y1991+landmask");

	if (filenumber == 96)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw45y1991+landmask");

	if (filenumber == 97)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw46y1991+landmask");

	if (filenumber == 98)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw47y1991+landmask");

	if (filenumber == 99)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw48y1991+landmask");

	if (filenumber == 100)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw49y1991+landmask");	 if (filenumber == 101)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw50y1991+landmask");



	if (filenumber == 102)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw51y1991+landmask");



	if (filenumber == 103)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1991\\Beaufort_Sea_diffw52y1991+landmask");



	//---------------------------------------------------------1992----------------------------------------------------------------------------------



	if (filenumber == 104)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw01y1992+landmask");

	if (filenumber == 105)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw02y1992+landmask");

	if (filenumber == 106)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw03y1992+landmask");

	if (filenumber == 107)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw04y1992+landmask");

	if (filenumber == 108)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw05y1992+landmask");

	if (filenumber == 109)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw06y1992+landmask");

	if (filenumber == 110)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw07y1992+landmask");

	if (filenumber == 111)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw08y1992+landmask");

	if (filenumber == 112)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw09y1992+landmask");	 if (filenumber == 113)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw10y1992+landmask");

	if (filenumber == 114)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw11y1992+landmask");

	if (filenumber == 115)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw12y1992+landmask");

	if (filenumber == 116)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw13y1992+landmask");

	if (filenumber == 117)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw14y1992+landmask");

	if (filenumber == 118)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw15y1992+landmask");

	if (filenumber == 119)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw16y1992+landmask");

	if (filenumber == 120)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw17y1992+landmask");

	if (filenumber == 121)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw18y1992+landmask");

	if (filenumber == 122)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw19y1992+landmask");	 if (filenumber == 123)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw20y1992+landmask");

	if (filenumber == 124)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw21y1992+landmask");

	if (filenumber == 125)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw22y1992+landmask");

	if (filenumber == 126)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw23y1992+landmask");

	if (filenumber == 127)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw24y1992+landmask");

	if (filenumber == 128)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw25y1992+landmask");

	if (filenumber == 129)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw26y1992+landmask");

	if (filenumber == 130)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw27y1992+landmask");

	if (filenumber == 131)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw28y1992+landmask");

	if (filenumber == 132)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw29y1992+landmask");	 if (filenumber == 133)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw30y1992+landmask");

	if (filenumber == 134)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw31y1992+landmask");

	if (filenumber == 135)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw32y1992+landmask");

	if (filenumber == 136)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw33y1992+landmask");

	if (filenumber == 137)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw34y1992+landmask");

	if (filenumber == 138)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw35y1992+landmask");

	if (filenumber == 139)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw36y1992+landmask");

	if (filenumber == 140)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw37y1992+landmask");

	if (filenumber == 141)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw38y1992+landmask");

	if (filenumber == 142)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw39y1992+landmask");	 if (filenumber == 143)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw40y1992+landmask");

	if (filenumber == 144)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw41y1992+landmask");

	if (filenumber == 145)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw42y1992+landmask");

	if (filenumber == 146)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw43y1992+landmask");

	if (filenumber == 147)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw44y1992+landmask");

	if (filenumber == 148)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw45y1992+landmask");

	if (filenumber == 149)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw46y1992+landmask");

	if (filenumber == 150)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw47y1992+landmask");

	if (filenumber == 151)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw48y1992+landmask");

	if (filenumber == 152)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw49y1992+landmask");	 if (filenumber == 153)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw50y1992+landmask");



	if (filenumber == 154)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw51y1992+landmask");



	if (filenumber == 155)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1992\\Beaufort_Sea_diffw52y1992+landmask");



	//---------------------------------------------------------1993----------------------------------------------------------------------------------


	if (filenumber == 156)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw01y1993+landmask");

	if (filenumber == 157)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw02y1993+landmask");

	if (filenumber == 158)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw03y1993+landmask");

	if (filenumber == 159)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw04y1993+landmask");

	if (filenumber == 160)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw05y1993+landmask");

	if (filenumber == 161)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw06y1993+landmask");

	if (filenumber == 162)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw07y1993+landmask");

	if (filenumber == 163)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw08y1993+landmask");

	if (filenumber == 164)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw09y1993+landmask");	 if (filenumber == 165)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw10y1993+landmask");

	if (filenumber == 166)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw11y1993+landmask");

	if (filenumber == 167)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw12y1993+landmask");

	if (filenumber == 168)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw13y1993+landmask");

	if (filenumber == 169)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw14y1993+landmask");

	if (filenumber == 170)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw15y1993+landmask");

	if (filenumber == 171)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw16y1993+landmask");

	if (filenumber == 172)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw17y1993+landmask");

	if (filenumber == 173)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw18y1993+landmask");

	if (filenumber == 174)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw19y1993+landmask");	 if (filenumber == 175)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw20y1993+landmask");

	if (filenumber == 176)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw21y1993+landmask");

	if (filenumber == 177)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw22y1993+landmask");

	if (filenumber == 178)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw23y1993+landmask");

	if (filenumber == 179)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw24y1993+landmask");

	if (filenumber == 180)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw25y1993+landmask");

	if (filenumber == 181)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw26y1993+landmask");

	if (filenumber == 182)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw27y1993+landmask");

	if (filenumber == 183)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw28y1993+landmask");

	if (filenumber == 184)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw29y1993+landmask");	 if (filenumber == 185)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw30y1993+landmask");

	if (filenumber == 186)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw31y1993+landmask");

	if (filenumber == 187)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw32y1993+landmask");

	if (filenumber == 188)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw33y1993+landmask");

	if (filenumber == 189)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw34y1993+landmask");

	if (filenumber == 190)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw35y1993+landmask");

	if (filenumber == 191)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw36y1993+landmask");

	if (filenumber == 192)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw37y1993+landmask");

	if (filenumber == 193)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw38y1993+landmask");

	if (filenumber == 194)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw39y1993+landmask");	 if (filenumber == 195)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw40y1993+landmask");

	if (filenumber == 196)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw41y1993+landmask");

	if (filenumber == 197)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw42y1993+landmask");

	if (filenumber == 198)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw43y1993+landmask");

	if (filenumber == 199)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw44y1993+landmask");

	if (filenumber == 200)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw45y1993+landmask");

	if (filenumber == 201)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw46y1993+landmask");

	if (filenumber == 202)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw47y1993+landmask");

	if (filenumber == 203)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw48y1993+landmask");

	if (filenumber == 204)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw49y1993+landmask");	 if (filenumber == 205)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw50y1993+landmask");



	if (filenumber == 206)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw51y1993+landmask");



	if (filenumber == 207)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1993\\Beaufort_Sea_diffw52y1993+landmask");





	//--------------------------------------------------------------1994----------------------------------------------------------------------------





	if (filenumber == 208)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw01y1994+landmask");

	if (filenumber == 209)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw02y1994+landmask");

	if (filenumber == 210)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw03y1994+landmask");

	if (filenumber == 211)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw04y1994+landmask");

	if (filenumber == 212)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw05y1994+landmask");

	if (filenumber == 213)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw06y1994+landmask");

	if (filenumber == 214)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw07y1994+landmask");

	if (filenumber == 215)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw08y1994+landmask");

	if (filenumber == 216)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw09y1994+landmask");	 if (filenumber == 217)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw10y1994+landmask");

	if (filenumber == 218)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw11y1994+landmask");

	if (filenumber == 219)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw12y1994+landmask");

	if (filenumber == 220)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw13y1994+landmask");

	if (filenumber == 221)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw14y1994+landmask");

	if (filenumber == 222)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw15y1994+landmask");

	if (filenumber == 223)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw16y1994+landmask");

	if (filenumber == 224)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw17y1994+landmask");

	if (filenumber == 225)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw18y1994+landmask");

	if (filenumber == 226)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw19y1994+landmask");	 if (filenumber == 227)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw20y1994+landmask");

	if (filenumber == 228)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw21y1994+landmask");

	if (filenumber == 229)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw22y1994+landmask");

	if (filenumber == 230)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw23y1994+landmask");

	if (filenumber == 231)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw24y1994+landmask");

	if (filenumber == 232)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw25y1994+landmask");

	if (filenumber == 233)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw26y1994+landmask");

	if (filenumber == 234)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw27y1994+landmask");

	if (filenumber == 235)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw28y1994+landmask");

	if (filenumber == 236)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw29y1994+landmask");	 if (filenumber == 237)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw30y1994+landmask");

	if (filenumber == 238)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw31y1994+landmask");

	if (filenumber == 239)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw32y1994+landmask");

	if (filenumber == 240)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw33y1994+landmask");

	if (filenumber == 241)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw34y1994+landmask");

	if (filenumber == 242)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw35y1994+landmask");

	if (filenumber == 243)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw36y1994+landmask");

	if (filenumber == 244)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw37y1994+landmask");

	if (filenumber == 245)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw38y1994+landmask");

	if (filenumber == 246)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw39y1994+landmask");	 if (filenumber == 247)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw40y1994+landmask");

	if (filenumber == 248)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw41y1994+landmask");

	if (filenumber == 249)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw42y1994+landmask");

	if (filenumber == 250)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw43y1994+landmask");

	if (filenumber == 251)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw44y1994+landmask");

	if (filenumber == 252)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw45y1994+landmask");

	if (filenumber == 253)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw46y1994+landmask");

	if (filenumber == 254)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw47y1994+landmask");

	if (filenumber == 255)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw48y1994+landmask");

	if (filenumber == 256)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw49y1994+landmask");	 if (filenumber == 257)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw50y1994+landmask");



	if (filenumber == 258)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw51y1994+landmask");



	if (filenumber == 259)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1994\\Beaufort_Sea_diffw52y1994+landmask");





	//--------------------------------------------------------------------1995---------------------------------------------------------------------





	if (filenumber == 260)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw01y1995+landmask");

	if (filenumber == 261)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw02y1995+landmask");

	if (filenumber == 262)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw03y1995+landmask");

	if (filenumber == 263)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw04y1995+landmask");

	if (filenumber == 264)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw05y1995+landmask");

	if (filenumber == 265)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw06y1995+landmask");

	if (filenumber == 266)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw07y1995+landmask");

	if (filenumber == 267)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw08y1995+landmask");

	if (filenumber == 268)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw09y1995+landmask");	 if (filenumber == 269)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw10y1995+landmask");

	if (filenumber == 270)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw11y1995+landmask");

	if (filenumber == 271)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw12y1995+landmask");

	if (filenumber == 272)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw13y1995+landmask");

	if (filenumber == 273)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw14y1995+landmask");

	if (filenumber == 274)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw15y1995+landmask");

	if (filenumber == 275)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw16y1995+landmask");

	if (filenumber == 276)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw17y1995+landmask");

	if (filenumber == 277)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw18y1995+landmask");

	if (filenumber == 278)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw19y1995+landmask");	 if (filenumber == 279)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw20y1995+landmask");

	if (filenumber == 280)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw21y1995+landmask");

	if (filenumber == 281)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw22y1995+landmask");

	if (filenumber == 282)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw23y1995+landmask");

	if (filenumber == 283)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw24y1995+landmask");

	if (filenumber == 284)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw25y1995+landmask");

	if (filenumber == 285)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw26y1995+landmask");

	if (filenumber == 286)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw27y1995+landmask");

	if (filenumber == 287)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw28y1995+landmask");

	if (filenumber == 288)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw29y1995+landmask");	 if (filenumber == 289)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw30y1995+landmask");

	if (filenumber == 290)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw31y1995+landmask");

	if (filenumber == 291)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw32y1995+landmask");

	if (filenumber == 292)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw33y1995+landmask");

	if (filenumber == 293)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw34y1995+landmask");

	if (filenumber == 294)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw35y1995+landmask");

	if (filenumber == 295)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw36y1995+landmask");

	if (filenumber == 296)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw37y1995+landmask");

	if (filenumber == 297)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw38y1995+landmask");

	if (filenumber == 298)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw39y1995+landmask");	 if (filenumber == 299)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw40y1995+landmask");

	if (filenumber == 300)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw41y1995+landmask");

	if (filenumber == 301)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw42y1995+landmask");

	if (filenumber == 302)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw43y1995+landmask");

	if (filenumber == 303)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw44y1995+landmask");

	if (filenumber == 304)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw45y1995+landmask");

	if (filenumber == 305)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw46y1995+landmask");

	if (filenumber == 306)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw47y1995+landmask");

	if (filenumber == 307)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw48y1995+landmask");

	if (filenumber == 308)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw49y1995+landmask");	 if (filenumber == 309)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw50y1995+landmask");



	if (filenumber == 310)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw51y1995+landmask");



	if (filenumber == 311)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1995\\Beaufort_Sea_diffw52y1995+landmask");





	//----------------------------------------------------------------1996------------------------------------------------------------------------





	if (filenumber == 312)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw01y1996+landmask");

	if (filenumber == 313)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw02y1996+landmask");

	if (filenumber == 314)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw03y1996+landmask");

	if (filenumber == 315)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw04y1996+landmask");

	if (filenumber == 316)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw05y1996+landmask");

	if (filenumber == 317)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw06y1996+landmask");

	if (filenumber == 318)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw07y1996+landmask");

	if (filenumber == 319)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw08y1996+landmask");

	if (filenumber == 320)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw09y1996+landmask");	 if (filenumber == 321)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw10y1996+landmask");

	if (filenumber == 322)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw11y1996+landmask");

	if (filenumber == 323)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw12y1996+landmask");

	if (filenumber == 324)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw13y1996+landmask");

	if (filenumber == 325)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw14y1996+landmask");

	if (filenumber == 326)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw15y1996+landmask");

	if (filenumber == 327)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw16y1996+landmask");

	if (filenumber == 328)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw17y1996+landmask");

	if (filenumber == 329)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw18y1996+landmask");

	if (filenumber == 330)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw19y1996+landmask");	 if (filenumber == 331)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw20y1996+landmask");

	if (filenumber == 332)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw21y1996+landmask");

	if (filenumber == 333)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw22y1996+landmask");

	if (filenumber == 334)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw23y1996+landmask");

	if (filenumber == 335)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw24y1996+landmask");

	if (filenumber == 336)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw25y1996+landmask");

	if (filenumber == 337)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw26y1996+landmask");

	if (filenumber == 338)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw27y1996+landmask");

	if (filenumber == 339)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw28y1996+landmask");

	if (filenumber == 340)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw29y1996+landmask");	 if (filenumber == 341)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw30y1996+landmask");

	if (filenumber == 342)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw31y1996+landmask");

	if (filenumber == 343)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw32y1996+landmask");

	if (filenumber == 344)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw33y1996+landmask");

	if (filenumber == 345)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw34y1996+landmask");

	if (filenumber == 346)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw35y1996+landmask");

	if (filenumber == 347)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw36y1996+landmask");

	if (filenumber == 348)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw37y1996+landmask");

	if (filenumber == 349)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw38y1996+landmask");

	if (filenumber == 350)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw39y1996+landmask");	 if (filenumber == 351)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw40y1996+landmask");

	if (filenumber == 352)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw41y1996+landmask");

	if (filenumber == 353)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw42y1996+landmask");

	if (filenumber == 354)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw43y1996+landmask");

	if (filenumber == 355)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw44y1996+landmask");

	if (filenumber == 356)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw45y1996+landmask");

	if (filenumber == 357)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw46y1996+landmask");

	if (filenumber == 358)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw47y1996+landmask");

	if (filenumber == 359)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw48y1996+landmask");

	if (filenumber == 360)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw49y1996+landmask");	 if (filenumber == 361)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw50y1996+landmask");



	if (filenumber == 362)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw51y1996+landmask");



	if (filenumber == 363)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1996\\Beaufort_Sea_diffw52y1996+landmask");





	//-------------------------------------------------------------1997-----------------------------------------------------------------------------





	if (filenumber == 364)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw01y1997+landmask");

	if (filenumber == 365)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw02y1997+landmask");

	if (filenumber == 366)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw03y1997+landmask");

	if (filenumber == 367)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw04y1997+landmask");

	if (filenumber == 368)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw05y1997+landmask");

	if (filenumber == 369)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw06y1997+landmask");

	if (filenumber == 370)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw07y1997+landmask");

	if (filenumber == 371)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw08y1997+landmask");

	if (filenumber == 372)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw09y1997+landmask");	 if (filenumber == 373)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw10y1997+landmask");

	if (filenumber == 374)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw11y1997+landmask");

	if (filenumber == 375)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw12y1997+landmask");

	if (filenumber == 376)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw13y1997+landmask");

	if (filenumber == 376)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw14y1997+landmask");

	if (filenumber == 377)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw15y1997+landmask");

	if (filenumber == 378)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw16y1997+landmask");

	if (filenumber == 379)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw17y1997+landmask");

	if (filenumber == 380)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw18y1997+landmask");

	if (filenumber == 381)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw19y1997+landmask");	 if (filenumber == 382)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw20y1997+landmask");

	if (filenumber == 383)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw21y1997+landmask");

	if (filenumber == 384)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw22y1997+landmask");

	if (filenumber == 385)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw23y1997+landmask");

	if (filenumber == 386)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw24y1997+landmask");

	if (filenumber == 387)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw25y1997+landmask");

	if (filenumber == 389)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw26y1997+landmask");

	if (filenumber == 390)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw27y1997+landmask");

	if (filenumber == 391)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw28y1997+landmask");

	if (filenumber == 392)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw29y1997+landmask");	 if (filenumber == 393)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw30y1997+landmask");

	if (filenumber == 394)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw31y1997+landmask");

	if (filenumber == 395)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw32y1997+landmask");

	if (filenumber == 396)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw33y1997+landmask");

	if (filenumber == 397)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw34y1997+landmask");

	if (filenumber == 398)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw35y1997+landmask");

	if (filenumber == 399)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw36y1997+landmask");

	if (filenumber == 400)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw37y1997+landmask");

	if (filenumber == 401)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw38y1997+landmask");

	if (filenumber == 402)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw39y1997+landmask");	 if (filenumber == 403)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw40y1997+landmask");

	if (filenumber == 404)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw41y1997+landmask");

	if (filenumber == 405)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw42y1997+landmask");

	if (filenumber == 406)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw43y1997+landmask");

	if (filenumber == 407)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw44y1997+landmask");

	if (filenumber == 408)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw45y1997+landmask");

	if (filenumber == 409)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw46y1997+landmask");

	if (filenumber == 410)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw47y1997+landmask");

	if (filenumber == 411)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw48y1997+landmask");

	if (filenumber == 412)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw49y1997+landmask");	 if (filenumber == 413)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw50y1997+landmask");



	if (filenumber == 414)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw51y1997+landmask");



	if (filenumber == 415)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1997\\Beaufort_Sea_diffw52y1997+landmask");





	//--------------------------------------------------------------1998-------------------------------------------------------------------------





	if (filenumber == 416)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw01y1998+landmask");

	if (filenumber == 417)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw02y1998+landmask");

	if (filenumber == 418)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw03y1998+landmask");

	if (filenumber == 419)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw04y1998+landmask");

	if (filenumber == 420)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw05y1998+landmask");

	if (filenumber == 421)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw06y1998+landmask");

	if (filenumber == 422)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw07y1998+landmask");

	if (filenumber == 423)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw08y1998+landmask");

	if (filenumber == 424)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw09y1998+landmask");	 if (filenumber == 425)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw10y1998+landmask");

	if (filenumber == 426)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw11y1998+landmask");

	if (filenumber == 427)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw12y1998+landmask");

	if (filenumber == 428)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw13y1998+landmask");

	if (filenumber == 429)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw14y1998+landmask");

	if (filenumber == 430)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw15y1998+landmask");

	if (filenumber == 431)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw16y1998+landmask");

	if (filenumber == 432)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw17y1998+landmask");

	if (filenumber == 433)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw18y1998+landmask");

	if (filenumber == 434)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw19y1998+landmask");	 if (filenumber == 435)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw20y1998+landmask");

	if (filenumber == 436)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw21y1998+landmask");

	if (filenumber == 437)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw22y1998+landmask");

	if (filenumber == 438)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw23y1998+landmask");

	if (filenumber == 439)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw24y1998+landmask");

	if (filenumber == 440)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw25y1998+landmask");

	if (filenumber == 441)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw26y1998+landmask");

	if (filenumber == 442)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw27y1998+landmask");

	if (filenumber == 443)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw28y1998+landmask");

	if (filenumber == 444)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw29y1998+landmask");	 if (filenumber == 445)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw30y1998+landmask");

	if (filenumber == 446)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw31y1998+landmask");

	if (filenumber == 447)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw32y1998+landmask");

	if (filenumber == 448)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw33y1998+landmask");

	if (filenumber == 449)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw34y1998+landmask");

	if (filenumber == 450)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw35y1998+landmask");

	if (filenumber == 451)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw36y1998+landmask");

	if (filenumber == 452)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw37y1998+landmask");

	if (filenumber == 453)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw38y1998+landmask");

	if (filenumber == 454)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw39y1998+landmask");	 if (filenumber == 455)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw40y1998+landmask");

	if (filenumber == 456)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw41y1998+landmask");

	if (filenumber == 457)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw42y1998+landmask");

	if (filenumber == 458)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw43y1998+landmask");

	if (filenumber == 459)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw44y1998+landmask");

	if (filenumber == 460)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw45y1998+landmask");

	if (filenumber == 461)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw46y1998+landmask");

	if (filenumber == 462)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw47y1998+landmask");

	if (filenumber == 463)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw48y1998+landmask");

	if (filenumber == 464)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw49y1998+landmask");	 if (filenumber == 465)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw50y1998+landmask");



	if (filenumber == 466)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw51y1998+landmask");



	if (filenumber == 467)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1998\\Beaufort_Sea_diffw52y1998+landmask");





	//---------------------------------------------------------------1999-----------------------------------------------------------------------





	if (filenumber == 468)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw01y1999+landmask");

	if (filenumber == 469)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw02y1999+landmask");

	if (filenumber == 470)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw03y1999+landmask");

	if (filenumber == 471)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw04y1999+landmask");

	if (filenumber == 472)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw05y1999+landmask");

	if (filenumber == 473)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw06y1999+landmask");

	if (filenumber == 474)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw07y1999+landmask");

	if (filenumber == 475)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw08y1999+landmask");

	if (filenumber == 476)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw09y1999+landmask");	 if (filenumber == 477)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw10y1999+landmask");

	if (filenumber == 478)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw11y1999+landmask");

	if (filenumber == 479)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw12y1999+landmask");

	if (filenumber == 480)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw13y1999+landmask");

	if (filenumber == 481)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw14y1999+landmask");

	if (filenumber == 482)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw15y1999+landmask");

	if (filenumber == 483)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw16y1999+landmask");

	if (filenumber == 484)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw17y1999+landmask");

	if (filenumber == 485)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw18y1999+landmask");

	if (filenumber == 486)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw19y1999+landmask");	 if (filenumber == 487)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw20y1999+landmask");

	if (filenumber == 488)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw21y1999+landmask");

	if (filenumber == 489)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw22y1999+landmask");

	if (filenumber == 490)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw23y1999+landmask");

	if (filenumber == 491)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw24y1999+landmask");

	if (filenumber == 492)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw25y1999+landmask");

	if (filenumber == 493)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw26y1999+landmask");

	if (filenumber == 494)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw27y1999+landmask");

	if (filenumber == 495)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw28y1999+landmask");

	if (filenumber == 496)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw29y1999+landmask");	 if (filenumber == 497)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw30y1999+landmask");

	if (filenumber == 498)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw31y1999+landmask");

	if (filenumber == 499)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw32y1999+landmask");

	if (filenumber == 500)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw33y1999+landmask");

	if (filenumber == 501)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw34y1999+landmask");

	if (filenumber == 502)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw35y1999+landmask");

	if (filenumber == 503)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw36y1999+landmask");

	if (filenumber == 504)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw37y1999+landmask");

	if (filenumber == 505)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw38y1999+landmask");

	if (filenumber == 506)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw39y1999+landmask");	 if (filenumber == 507)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw40y1999+landmask");

	if (filenumber == 508)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw41y1999+landmask");

	if (filenumber == 509)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw42y1999+landmask");

	if (filenumber == 510)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw43y1999+landmask");

	if (filenumber == 511)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw44y1999+landmask");

	if (filenumber == 512)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw45y1999+landmask");

	if (filenumber == 513)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw46y1999+landmask");

	if (filenumber == 514)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw47y1999+landmask");

	if (filenumber == 515)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw48y1999+landmask");

	if (filenumber == 516)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw49y1999+landmask");	 if (filenumber == 517)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw50y1999+landmask");



	if (filenumber == 518)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw51y1999+landmask");



	if (filenumber == 519)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\1999\\Beaufort_Sea_diffw52y1999+landmask");





	//-------------------------------------------------------2000------------------------------------------------------------------------------





	if (filenumber == 520)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw01y2000+landmask");

	if (filenumber == 521)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw02y2000+landmask");

	if (filenumber == 522)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw03y2000+landmask");

	if (filenumber == 523)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw04y2000+landmask");

	if (filenumber == 524)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw05y2000+landmask");

	if (filenumber == 525)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw06y2000+landmask");

	if (filenumber == 526)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw07y2000+landmask");

	if (filenumber == 527)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw08y2000+landmask");

	if (filenumber == 528)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw09y2000+landmask");	 if (filenumber == 529)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw10y2000+landmask");

	if (filenumber == 530)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw11y2000+landmask");

	if (filenumber == 531)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw12y2000+landmask");

	if (filenumber == 532)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw13y2000+landmask");

	if (filenumber == 533)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw14y2000+landmask");

	if (filenumber == 534)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw15y2000+landmask");

	if (filenumber == 535)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw16y2000+landmask");

	if (filenumber == 536)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw17y2000+landmask");

	if (filenumber == 537)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw18y2000+landmask");

	if (filenumber == 538)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw19y2000+landmask");	 if (filenumber == 539)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw20y2000+landmask");

	if (filenumber == 540)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw21y2000+landmask");

	if (filenumber == 541)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw22y2000+landmask");

	if (filenumber == 542)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw23y2000+landmask");

	if (filenumber == 543)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw24y2000+landmask");

	if (filenumber == 544)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw25y2000+landmask");

	if (filenumber == 545)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw26y2000+landmask");

	if (filenumber == 546)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw27y2000+landmask");

	if (filenumber == 547)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw28y2000+landmask");

	if (filenumber == 548)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw29y2000+landmask");	 if (filenumber == 549)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw30y2000+landmask");

	if (filenumber == 550)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw31y2000+landmask");

	if (filenumber == 551)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw32y2000+landmask");

	if (filenumber == 552)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw33y2000+landmask");

	if (filenumber == 553)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw34y2000+landmask");

	if (filenumber == 554)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw35y2000+landmask");

	if (filenumber == 555)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw36y2000+landmask");

	if (filenumber == 556)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw37y2000+landmask");

	if (filenumber == 557)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw38y2000+landmask");

	if (filenumber == 558)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw39y2000+landmask");	 if (filenumber == 559)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw40y2000+landmask");

	if (filenumber == 560)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw41y2000+landmask");

	if (filenumber == 561)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw42y2000+landmask");

	if (filenumber == 562)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw43y2000+landmask");

	if (filenumber == 563)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw44y2000+landmask");

	if (filenumber == 564)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw45y2000+landmask");

	if (filenumber == 565)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw46y2000+landmask");

	if (filenumber == 566)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw47y2000+landmask");

	if (filenumber == 567)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw48y2000+landmask");

	if (filenumber == 568)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw49y2000+landmask");	 if (filenumber == 569)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw50y2000+landmask");



	if (filenumber == 570)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw51y2000+landmask");



	if (filenumber == 571)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2000\\Beaufort_Sea_diffw52y2000+landmask");



	//------------------------------------------------------2001--------------------------------------------------------------------------------





	if (filenumber == 572)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw01y2001+landmask");

	if (filenumber == 573)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw02y2001+landmask");

	if (filenumber == 574)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw03y2001+landmask");

	if (filenumber == 575)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw04y2001+landmask");

	if (filenumber == 576)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw05y2001+landmask");

	if (filenumber == 577)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw06y2001+landmask");

	if (filenumber == 578)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw07y2001+landmask");

	if (filenumber == 579)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw08y2001+landmask");

	if (filenumber == 580)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw09y2001+landmask");	 if (filenumber == 581)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw10y2001+landmask");

	if (filenumber == 582)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw11y2001+landmask");

	if (filenumber == 583)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw12y2001+landmask");

	if (filenumber == 584)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw13y2001+landmask");

	if (filenumber == 585)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw14y2001+landmask");

	if (filenumber == 586)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw15y2001+landmask");

	if (filenumber == 587)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw16y2001+landmask");

	if (filenumber == 588)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw17y2001+landmask");

	if (filenumber == 589)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw18y2001+landmask");

	if (filenumber == 590)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw19y2001+landmask");	 if (filenumber == 591)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw20y2001+landmask");

	if (filenumber == 592)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw21y2001+landmask");

	if (filenumber == 593)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw22y2001+landmask");

	if (filenumber == 594)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw23y2001+landmask");

	if (filenumber == 595)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw24y2001+landmask");

	if (filenumber == 596)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw25y2001+landmask");

	if (filenumber == 597)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw26y2001+landmask");

	if (filenumber == 598)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw27y2001+landmask");

	if (filenumber == 599)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw28y2001+landmask");

	if (filenumber == 600)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw29y2001+landmask");	 if (filenumber == 601)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw30y2001+landmask");

	if (filenumber == 602)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw31y2001+landmask");

	if (filenumber == 603)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw32y2001+landmask");

	if (filenumber == 604)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw33y2001+landmask");

	if (filenumber == 605)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw34y2001+landmask");

	if (filenumber == 606)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw35y2001+landmask");

	if (filenumber == 607)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw36y2001+landmask");

	if (filenumber == 608)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw37y2001+landmask");

	if (filenumber == 609)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw38y2001+landmask");

	if (filenumber == 610)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw39y2001+landmask");	 if (filenumber == 611)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw40y2001+landmask");

	if (filenumber == 612)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw41y2001+landmask");

	if (filenumber == 613)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw42y2001+landmask");

	if (filenumber == 614)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw43y2001+landmask");

	if (filenumber == 615)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw44y2001+landmask");

	if (filenumber == 616)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw45y2001+landmask");

	if (filenumber == 617)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw46y2001+landmask");

	if (filenumber == 618)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw47y2001+landmask");

	if (filenumber == 619)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw48y2001+landmask");

	if (filenumber == 620)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw49y2001+landmask");	 if (filenumber == 621)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw50y2001+landmask");



	if (filenumber == 622)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw51y2001+landmask");



	if (filenumber == 623)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2001\\Beaufort_Sea_diffw52y2001+landmask");





	//---------------------------------------------------------2002-----------------------------------------------------------------------------





	if (filenumber == 624)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw01y2002+landmask");

	if (filenumber == 625)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw02y2002+landmask");

	if (filenumber == 626)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw03y2002+landmask");

	if (filenumber == 627)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw04y2002+landmask");

	if (filenumber == 628)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw05y2002+landmask");

	if (filenumber == 629)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw06y2002+landmask");

	if (filenumber == 630)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw07y2002+landmask");

	if (filenumber == 631)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw08y2002+landmask");

	if (filenumber == 632)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw09y2002+landmask");	 if (filenumber == 633)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw10y2002+landmask");

	if (filenumber == 634)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw11y2002+landmask");

	if (filenumber == 635)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw12y2002+landmask");

	if (filenumber == 636)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw13y2002+landmask");

	if (filenumber == 637)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw14y2002+landmask");

	if (filenumber == 638)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw15y2002+landmask");

	if (filenumber == 639)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw16y2002+landmask");

	if (filenumber == 640)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw17y2002+landmask");

	if (filenumber == 641)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw18y2002+landmask");

	if (filenumber == 642)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw19y2002+landmask");	 if (filenumber == 643)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw20y2002+landmask");

	if (filenumber == 644)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw21y2002+landmask");

	if (filenumber == 645)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw22y2002+landmask");

	if (filenumber == 646)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw23y2002+landmask");

	if (filenumber == 647)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw24y2002+landmask");

	if (filenumber == 648)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw25y2002+landmask");

	if (filenumber == 649)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw26y2002+landmask");

	if (filenumber == 650)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw27y2002+landmask");

	if (filenumber == 651)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw28y2002+landmask");

	if (filenumber == 652)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw29y2002+landmask");	 if (filenumber == 653)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw30y2002+landmask");

	if (filenumber == 654)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw31y2002+landmask");

	if (filenumber == 655)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw32y2002+landmask");

	if (filenumber == 656)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw33y2002+landmask");

	if (filenumber == 657)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw34y2002+landmask");

	if (filenumber == 658)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw35y2002+landmask");

	if (filenumber == 659)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw36y2002+landmask");

	if (filenumber == 660)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw37y2002+landmask");

	if (filenumber == 661)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw38y2002+landmask");

	if (filenumber == 662)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw39y2002+landmask");	 if (filenumber == 663)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw40y2002+landmask");

	if (filenumber == 664)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw41y2002+landmask");

	if (filenumber == 665)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw42y2002+landmask");

	if (filenumber == 666)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw43y2002+landmask");

	if (filenumber == 667)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw44y2002+landmask");

	if (filenumber == 668)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw45y2002+landmask");

	if (filenumber == 669)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw46y2002+landmask");

	if (filenumber == 670)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw47y2002+landmask");

	if (filenumber == 671)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw48y2002+landmask");

	if (filenumber == 672)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw49y2002+landmask");	 if (filenumber == 673)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw50y2002+landmask");



	if (filenumber == 674)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw51y2002+landmask");



	if (filenumber == 675)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2002\\Beaufort_Sea_diffw52y2002+landmask");





	//----------------------------------------------------------------2003-----------------------------------------------------------------------





	if (filenumber == 676)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw01y2003+landmask");

	if (filenumber == 677)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw02y2003+landmask");

	if (filenumber == 678)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw03y2003+landmask");

	if (filenumber == 679)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw04y2003+landmask");

	if (filenumber == 680)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw05y2003+landmask");

	if (filenumber == 681)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw06y2003+landmask");

	if (filenumber == 682)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw07y2003+landmask");

	if (filenumber == 683)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw08y2003+landmask");

	if (filenumber == 684)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw09y2003+landmask");	 if (filenumber == 685)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw10y2003+landmask");

	if (filenumber == 686)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw11y2003+landmask");

	if (filenumber == 687)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw12y2003+landmask");

	if (filenumber == 688)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw13y2003+landmask");

	if (filenumber == 689)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw14y2003+landmask");

	if (filenumber == 690)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw15y2003+landmask");

	if (filenumber == 691)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw16y2003+landmask");

	if (filenumber == 692)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw17y2003+landmask");

	if (filenumber == 693)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw18y2003+landmask");

	if (filenumber == 694)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw19y2003+landmask");	 if (filenumber == 695)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw20y2003+landmask");

	if (filenumber == 696)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw21y2003+landmask");

	if (filenumber == 697)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw22y2003+landmask");

	if (filenumber == 698)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw23y2003+landmask");

	if (filenumber == 699)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw24y2003+landmask");

	if (filenumber == 700)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw25y2003+landmask");

	if (filenumber == 701)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw26y2003+landmask");

	if (filenumber == 702)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw27y2003+landmask");

	if (filenumber == 703)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw28y2003+landmask");

	if (filenumber == 704)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw29y2003+landmask");	 if (filenumber == 705)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw30y2003+landmask");

	if (filenumber == 706)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw31y2003+landmask");

	if (filenumber == 707)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw32y2003+landmask");

	if (filenumber == 708)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw33y2003+landmask");

	if (filenumber == 709)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw34y2003+landmask");

	if (filenumber == 710)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw35y2003+landmask");

	if (filenumber == 711)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw36y2003+landmask");

	if (filenumber == 712)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw37y2003+landmask");

	if (filenumber == 713)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw38y2003+landmask");

	if (filenumber == 714)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw39y2003+landmask");	 if (filenumber == 715)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw40y2003+landmask");

	if (filenumber == 716)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw41y2003+landmask");

	if (filenumber == 717)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw42y2003+landmask");

	if (filenumber == 718)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw43y2003+landmask");

	if (filenumber == 719)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw44y2003+landmask");

	if (filenumber == 720)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw45y2003+landmask");

	if (filenumber == 721)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw46y2003+landmask");

	if (filenumber == 722)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw47y2003+landmask");

	if (filenumber == 723)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw48y2003+landmask");

	if (filenumber == 724)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw49y2003+landmask");	 if (filenumber == 725)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw50y2003+landmask");



	if (filenumber == 726)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw51y2003+landmask");



	if (filenumber == 727)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2003\\Beaufort_Sea_diffw52y2003+landmask");





	//----------------------------------------------------------2004----------------------------------------------------------------------------





	if (filenumber == 728)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw01y2004+landmask");

	if (filenumber == 729)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw02y2004+landmask");

	if (filenumber == 730)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw03y2004+landmask");

	if (filenumber == 731)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw04y2004+landmask");

	if (filenumber == 732)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw05y2004+landmask");

	if (filenumber == 733)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw06y2004+landmask");

	if (filenumber == 734)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw07y2004+landmask");

	if (filenumber == 735)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw08y2004+landmask");

	if (filenumber == 736)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw09y2004+landmask");	 if (filenumber == 737)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw10y2004+landmask");

	if (filenumber == 738)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw11y2004+landmask");

	if (filenumber == 739)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw12y2004+landmask");

	if (filenumber == 740)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw13y2004+landmask");

	if (filenumber == 741)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw14y2004+landmask");

	if (filenumber == 742)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw15y2004+landmask");

	if (filenumber == 743)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw16y2004+landmask");

	if (filenumber == 744)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw17y2004+landmask");

	if (filenumber == 745)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw18y2004+landmask");

	if (filenumber == 746)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw19y2004+landmask");	 if (filenumber == 747)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw20y2004+landmask");

	if (filenumber == 748)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw21y2004+landmask");

	if (filenumber == 749)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw22y2004+landmask");

	if (filenumber == 750)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw23y2004+landmask");

	if (filenumber == 751)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw24y2004+landmask");

	if (filenumber == 752)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw25y2004+landmask");

	if (filenumber == 753)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw26y2004+landmask");

	if (filenumber == 754)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw27y2004+landmask");

	if (filenumber == 755)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw28y2004+landmask");

	if (filenumber == 756)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw29y2004+landmask");	 if (filenumber == 757)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw30y2004+landmask");

	if (filenumber == 758)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw31y2004+landmask");

	if (filenumber == 759)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw32y2004+landmask");

	if (filenumber == 760)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw33y2004+landmask");

	if (filenumber == 761)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw34y2004+landmask");

	if (filenumber == 762)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw35y2004+landmask");

	if (filenumber == 763)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw36y2004+landmask");

	if (filenumber == 764)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw37y2004+landmask");

	if (filenumber == 765)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw38y2004+landmask");

	if (filenumber == 766)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw39y2004+landmask");	 if (filenumber == 767)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw40y2004+landmask");

	if (filenumber == 768)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw41y2004+landmask");

	if (filenumber == 769)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw42y2004+landmask");

	if (filenumber == 770)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw43y2004+landmask");

	if (filenumber == 771)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw44y2004+landmask");

	if (filenumber == 772)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw45y2004+landmask");

	if (filenumber == 773)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw46y2004+landmask");

	if (filenumber == 774)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw47y2004+landmask");

	if (filenumber == 775)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw48y2004+landmask");

	if (filenumber == 776)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw49y2004+landmask");	 if (filenumber == 777)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw50y2004+landmask");



	if (filenumber == 778)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw51y2004+landmask");



	if (filenumber == 779)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2004\\Beaufort_Sea_diffw52y2004+landmask");





	//-------------------------------------------------------------2005------------------------------------------------------------------------





	if (filenumber == 780)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw01y2005+landmask");

	if (filenumber == 781)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw02y2005+landmask");

	if (filenumber == 782)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw03y2005+landmask");

	if (filenumber == 783)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw04y2005+landmask");

	if (filenumber == 784)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw05y2005+landmask");

	if (filenumber == 785)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw06y2005+landmask");

	if (filenumber == 786)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw07y2005+landmask");

	if (filenumber == 787)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw08y2005+landmask");

	if (filenumber == 788)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw09y2005+landmask");	 if (filenumber == 789)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw10y2005+landmask");

	if (filenumber == 790)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw11y2005+landmask");

	if (filenumber == 791)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw12y2005+landmask");

	if (filenumber == 792)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw13y2005+landmask");

	if (filenumber == 793)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw14y2005+landmask");

	if (filenumber == 794)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw15y2005+landmask");

	if (filenumber == 795)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw16y2005+landmask");

	if (filenumber == 796)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw17y2005+landmask");

	if (filenumber == 797)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw18y2005+landmask");

	if (filenumber == 798)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw19y2005+landmask");	 if (filenumber == 799)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw20y2005+landmask");

	if (filenumber == 800)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw21y2005+landmask");

	if (filenumber == 801)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw22y2005+landmask");

	if (filenumber == 802)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw23y2005+landmask");

	if (filenumber == 803)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw24y2005+landmask");

	if (filenumber == 804)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw25y2005+landmask");

	if (filenumber == 805)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw26y2005+landmask");

	if (filenumber == 806)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw27y2005+landmask");

	if (filenumber == 807)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw28y2005+landmask");

	if (filenumber == 808)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw29y2005+landmask");	 if (filenumber == 809)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw30y2005+landmask");

	if (filenumber == 810)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw31y2005+landmask");

	if (filenumber == 811)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw32y2005+landmask");

	if (filenumber == 812)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw33y2005+landmask");

	if (filenumber == 813)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw34y2005+landmask");

	if (filenumber == 814)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw35y2005+landmask");

	if (filenumber == 815)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw36y2005+landmask");

	if (filenumber == 816)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw37y2005+landmask");

	if (filenumber == 817)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw38y2005+landmask");

	if (filenumber == 818)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw39y2005+landmask");	 if (filenumber == 819)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw40y2005+landmask");

	if (filenumber == 820)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw41y2005+landmask");

	if (filenumber == 821)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw42y2005+landmask");

	if (filenumber == 822)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw43y2005+landmask");

	if (filenumber == 823)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw44y2005+landmask");

	if (filenumber == 824)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw45y2005+landmask");

	if (filenumber == 825)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw46y2005+landmask");

	if (filenumber == 826)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw47y2005+landmask");

	if (filenumber == 827)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw48y2005+landmask");

	if (filenumber == 828)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw49y2005+landmask");	if (filenumber == 829)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw50y2005+landmask");



	if (filenumber == 830)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw51y2005+landmask");



	if (filenumber == 831)		strcpy_s(filename, "C:\\Sivakumar\\ASU\\SEM 1\\CEN 502\\Project\\CS310_project_subregion\\2005\\Beaufort_Sea_diffw52y2005+landmask");



}
