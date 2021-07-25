#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <bits/stdc++.h>
using namespace std;

const double residual = 1e-6;

// Delete n-th line from given file
void delete_line(const char *file_name, int n)
{
    // open file in read mode or in mode
    ifstream is(file_name);
  
    // open file in write mode or out mode
    ofstream ofs;
    ofs.open("temp.txt", ofstream::out);
  
    // loop getting single characters
    char c;
    int line_no = 1;
    while (is.get(c))
    {
        // if a newline character
        if (c == '\n')
        line_no++;
  
        // file content not to be deleted
        if (line_no != n)
            ofs << c;
    }
  
    // closing output file
    ofs.close();
  
    // closing input file
    is.close();
  
    // remove the original file
    remove(file_name);
  
    // rename the file
    rename("temp.txt", file_name);
}

void blank_line(const char *file_name)
{   
  ifstream fin(file_name);    
  
  ofstream fout;                
  fout.open("temp.txt", ios::out);
  
  string str;
  while(getline(fin,str))
  { 
    while (str.length()==0 ) 
       getline(fin,str);   
  
    fout<<str<<endl;
  }
  fout.close();  
  fin.close();  
  remove(file_name);        
  rename("temp.txt", file_name);
}


int get_colNum () {
	double val=0.;
	int colNum=0;
	ifstream dataFile;
	dataFile.open ("_output");
	stringstream ss;
	string temp;
	if (dataFile.is_open()) {
		//cout << "File has been opened" << endl;
		getline(dataFile,temp);
		ss << temp;
		while(ss >> val) {
			//cout << val << ' ';
			++colNum;
		}
		/*while( dataFile.peek() != '\n' && dataFile >> val)
		{
			cout << val << endl;
			++colNum;
		}*/
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
	dataFile.open ("_output");
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
	delete_line("_output", 1);
	blank_line("_output");
	int rowNum=0, colNum=0, pos99=0, posStat=0;
	colNum = get_colNum();
	rowNum = get_rowNum(colNum);
	//cout << "sizing done" << endl;
	double** dataMatrix = new double*[rowNum];
	for (int i = 0; i < rowNum; i++) {
		dataMatrix[i] = new double[colNum];
		for (int j = 0; j < colNum; j++) {
			dataMatrix[i][j] = 0.;
		}
	}
	//cout << "2d array done" << endl;
	//double dataMatrix[rowNum][colNum]{};
	//cout << "init done" << endl;
	double t_99=0., x_99=0., aveSpeed=0.;
	
	
	ifstream dataFile;
	dataFile.open ("_output");
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
		//cout << setprecision(12) << dataMatrix[rowNum-1][13] << endl;
	} else {
		cout << "File could not be opened" << endl;
	}
	//find stationary point with residual
	/*for (int currentRow=0; currentRow < rowNum-1; ++currentRow) {
		//cout << "current Row = " << currentRow+1 << endl;
		if (abs((dataMatrix[currentRow][13]-dataMatrix[currentRow+1][13])/dataMatrix[currentRow][13]) < residual) {
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
		sumSpeed += abs(dataMatrix[currentRow][13]);
	}
	aveSpeed = sumSpeed/(rowNum-1-posStat+1);*/
	aveSpeed = abs(dataMatrix[rowNum-1][13]);
	cout << "Sedimentation speed = " << setprecision(12) << aveSpeed << " m/s" << endl;
	
	
	for (int currentRow=0; currentRow < rowNum-1; ++currentRow) {
		if (abs(abs(dataMatrix[currentRow][13])-aveSpeed)/aveSpeed < 0.01) {
			//cout << "99% speed at row " << currentRow+1 << endl;
			//cout << "residual = " << abs(abs(dataMatrix[currentRow][13])-aveSpeed)/aveSpeed << endl;
			cout << "v_99 = " << abs(dataMatrix[currentRow][13]) << " m/s" << endl;
			t_99 = dataMatrix[currentRow][1] - dataMatrix[0][1];
			cout << "t_99 = " << t_99 << " s" << endl;
			x_99 = dataMatrix[0][6] - dataMatrix[currentRow][6];
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
