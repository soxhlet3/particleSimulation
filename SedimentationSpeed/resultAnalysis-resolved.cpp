#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

const double residual = 1e-3;

int get_colNum () {
	double val=0.;
	int colNum=0;
	ifstream dataFile;
	dataFile.open ("gnuplot.dat");
	if (dataFile.is_open()) {
		//cout << "File has been opened" << endl;
		while( dataFile.peek() != '\n' && dataFile >> val)
		{
			++colNum;
		}
		//cout << "There are " << colNum << " columns of data" << endl;
		return colNum;
	} else {
		cout << "File could not be opened" << endl;
		return 0;
	}
}

int get_rowNum (int colNum) {
	double val=0.;
	int itemNum=0;
	ifstream dataFile;
	dataFile.open ("gnuplot.dat");
	if (dataFile.is_open()) {
		//cout << "File has been opened" << endl;
		while( dataFile >> val)
		{
			++itemNum;
		}
		//cout << "There are " << itemNum/colNum << " rows of data" << endl;
		return itemNum/colNum;
	} else {
		cout << "File could not be opened" << endl;
		return 0;
	}
}


void printResults (double t_99, double x_99, double aveSpeed, double stdDev) {

	cout ;
}

int main () {
	int rowNum=0, colNum=0, pos99=0, posStat=0;
	colNum = get_colNum();
	rowNum = get_rowNum(colNum);
	double** dataMatrix = new double*[rowNum];
	for (int i; i<rowNum; i++) {
		dataMatrix[i] = new double[colNum];
		for (int j; j<colNum; j++) {
		dataMatrix[i][j] = 0.;
		}
	}
	double t_99=0., x_99=0., aveSpeed=0.;
	
	
	ifstream dataFile;
	dataFile.open ("gnuplot.dat");
	if (dataFile.is_open()) {
		//cout << "File has been opened" << endl;
		int currentRow=0;
		while(dataFile)
		{
			for (int currentCol=0; currentCol < colNum; ++currentCol) {
				dataFile >> dataMatrix[currentRow][currentCol];
			}
			++currentRow;
		}
		//cout << setprecision(12) << dataMatrix[rowNum-1][2] << endl;
	} else {
		cout << "File could not be opened" << endl;
	}
	//find stationary point with residual
	/*for (int currentRow=0; currentRow < rowNum-1; ++currentRow) {
		//cout << "current Row = " << currentRow+1 << endl;
		if (abs((dataMatrix[currentRow][2]-dataMatrix[currentRow+1][2])/dataMatrix[currentRow][2]) < residual) {
			posStat = currentRow + 1;
			//cout << "stationary point found at row " << posStat + 1<< endl;
			break;
		} else {
			//cout << "stationary point not arrived yet" << endl;
			if (currentRow == rowNum-2) {
				posStat = currentRow + 1;
				cout << "stationary state assumed to be arrived at last time step" << endl;
			}
		}
	}
	//cout << "steady state starts from row " << posStat + 1 << endl;
	//calculate average speed after reaching stationary
	double sumSpeed=0.;
	for (int currentRow=posStat; currentRow < rowNum; ++currentRow) {
		sumSpeed += abs(dataMatrix[currentRow][2]);
	}
	aveSpeed = sumSpeed/(rowNum-1-posStat+1);*/
	aveSpeed = abs( dataMatrix[rowNum-1][2]);
	cout << "Sedimentation speed = " << setprecision(12) << aveSpeed << " m/s" << endl;
	
	
	for (int currentRow=0; currentRow < rowNum-1; ++currentRow) {
		if (abs(abs(dataMatrix[currentRow][2])-aveSpeed)/aveSpeed < 0.01) {
			//cout << "99% speed at row " << currentRow+1 << endl;
			//cout << "residual = " << abs(abs(dataMatrix[currentRow][2])-aveSpeed)/aveSpeed << endl;
			cout << "v_99 = " << abs(dataMatrix[currentRow][2]) << " m/s" << endl;
			t_99 = dataMatrix[currentRow][0];
			cout << "t_99 = " << t_99 << " s" << endl;
			x_99 = dataMatrix[0][1] - dataMatrix[currentRow][1];
			cout << "x_99 = " << x_99*1000. << " mm" << endl;
			break;
		} else if (currentRow == rowNum-1) {
			cout << "Temporal resolution of output data is not small enough to determine x_99 and t_99" << endl;
		}
	}
	
	for (int j=0; j<colNum; j++) {
	delete[] dataMatrix[j];
	}
	delete[] dataMatrix;
	return 0;
}
