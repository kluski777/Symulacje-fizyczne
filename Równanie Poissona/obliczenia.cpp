#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

void saveToFile(double**,std::string);
void dealloc2d(double**);
double** alloc2d();
void initialize(double**);
void printmatrix(double**);
void changeU(double**);
double changeA(double**);
void saveToFile(int,std::string);
void saveToFile(double*,std::string);
double** calRoPrim(double**);
double** calcDelta(double**);
void changeU(double**,double);
void MinimizeU(double**);
double minDelta(double**);

const int N = 31;
const int iterations = 501;


double rho(int x, int y){
	return exp(-((x-4)*(x-4)+(y*y))/16)-exp(-((x+4)*(x+4)+(y*y))/16);
}

int main(){
  double** u = alloc2d();
	double* a = new double[iterations];
	initialize(u);
	saveToFile(N,"Nka");
	saveToFile(iterations,"Itery");
	for(int i=0;i<iterations;i++){
		a[i] = changeA(u);
		changeU(u);
		if(i == 100)
			saveToFile(u,"Po100Iteracjach");
	}
	saveToFile(u,"Po500Iteracjach");
	saveToFile(a,"ADlaWR1");
	double** roPrim = calRoPrim(u);
	double** delta = calcDelta(roPrim);
	saveToFile(roPrim,"roPrim1");
	saveToFile(delta,"delta1");
	///// pierwsza czêœæ zadania
	initialize(u);
	for(int i=0;i<iterations;i++){
		a[i] = changeA(u);
		changeU(u,1.9);
		if( i == 100)
			saveToFile(u,"Po100IteracjachP2");
	}
	saveToFile(u,"Po500IteracjachP2");
	saveToFile(a,"ADlaWR2");
	roPrim = calRoPrim(u);
	delta = calcDelta(roPrim);
	saveToFile(roPrim,"roPrim2");
	saveToFile(delta,"delta2");
	///// koniec czêœci d³ugiej
	initialize(u);
	for(int i=0;i<iterations;i++){
		a[i] = changeA(u);
		changeU(u); //mo¿naby uproœciæ tak by by³a tylko 1 funkcja
		MinimizeU(u);
		if( i == 100)
			saveToFile(u,"Po100IteracjachP3");
	}
	saveToFile(u,"Po500IteracjachP3");
	saveToFile(a,"ADlaWR3");
	roPrim = calRoPrim(u);
	delta = calcDelta(roPrim);
	saveToFile(roPrim,"roPrim3");
	saveToFile(delta,"delta3");
	dealloc2d(u);
	dealloc2d(delta);
	delete [] a;
	return 0;
}

void dealloc2d(double** u){
	for(int i=0;i<2*N-1;i++)
		delete[] u[i];
	delete[] u;
}

double** calRoPrim(double** u){
	double** toRet = alloc2d();
	for(int i=1;i<2*N-2;i++){
		for(int j=1;j<2*N-2;j++)
			toRet[i][j] = 4*u[i][j]-u[i][j+1]-u[i][j-1]-u[i+1][j]-u[i-1][j];
	}
	return toRet;
}

double** calcDelta(double** ro2){
	double** temp = alloc2d();
	for(int i=1;i<2*N-2;i++){
		for(int j=1;j<2*N-2;j++)
			temp[i][j] = ro2[i][j] - rho(i-N+1,j-N+1); //na brzegach zostawiamy 0
	}
	return temp;
}

void initialize(double** u){
	for(int i=0;i<2*N-1;i++){
		for(int j=0;j<2*N-1;j++){
			if( i >= 15 && i <= 45 && j >= 15 && j <= 45)
				u[i][j] = 1.0;
			else
				u[i][j] = 0.0;
		}
	}
}

void MinimizeU(double** u){
	double* delta = new double[4];
	int index = 0;
	double min = 0.0;
	for(int i=1;i<2*N-2;i++){
		for(int j=1;j<2*N-2;j++){
			min = std::numeric_limits<double>::max();
			double* x = new double[4];
			for(int a=0;a<=3;a++){
				if(a==3 && abs(x[0]+x[2]-2*x[1]) > 1e-12){
					delta[a] = u[i][j] + 0.25 + 0.5 * (x[0]-x[1])/(x[0]-2*x[1]+x[2]);
				}
				else if(a!=3)
					delta[a] = 0.5*a+u[i][j];
				else
					break;
				if(a!=3 && abs(x[0]+x[2]-2*x[1]) > 1e-12)
					x[a] = a*i*((u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])/2-2*+rho(i-(N-1),j-(N-1)));
				if(x[a] < min){
					min = x[a];
					index = a;
				}
			}
			u[i][j] = delta[index];
		}
	}
	delete [] delta;
}

void changeU(double** u,double w){
	for(int i=1;i<2*N-2;i++){
		for(int j=1;j<2*N-2;j++){
			u[i][j] = (1-w)*u[i][j]+w*(0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]+rho(i-N+1,j-N+1)));
		}
	}
}
		
void changeU(double** u){
	for(int i=1;i<2*N-2;i++){
		for(int j=1;j<2*N-2;j++)
			u[i][j] = u[i+1][j]*0.25 + u[i-1][j]*0.25 + u[i][j+1]*0.25 + u[i][j-1]*0.25 + rho(i-(N-1),j-(N-1))*0.25;
	}
}

void printmatrix(double** u){
	for(int i=0;i<2*N-1;i++){
		printf("\n");
		for(int j=0;j<2*N-1;j++)
			printf("%0.2lf ",u[i][j]);
	}
}

double changeA(double** u){
	double sum = 0.0;
	for(int i=1;i<2*N-2;i++){
		for(int j=1;j<2*N-2;j++)
			sum -= u[i][j]*u[i+1][j]/2+u[i][j]*u[i-1][j]/2+u[i][j]*u[i][j+1]/2+u[i][j]*u[i][j-1]/2-2*u[i][j]*u[i][j] + u[i][j]*rho(i-(N-1),j-(N-1));
	}
	return sum;
}

void saveToFile(int cus,std::string title){
	std::ofstream outdata;
	outdata.open(title);
	if(!outdata){
		std::cout<<"Nie uda³o siê zapisaæ danych do pliku "<<title<<std::endl;
		return;
	}
	outdata << N;
	std::cout<<"Dane zapisano w pliku "<<title<<".txt"<<std::endl;
}

void saveToFile(double** u,std::string title){
	std::ofstream outdata;
	outdata.open(title);
	if(!outdata){
		std::cout<<"Nie uda³o siê zapisaæ danych do pliku "<<title<<std::endl;
		return;
	}
	for(int i=0;i<2*N-1;i++){
		for(int j=0;j<2*N-1;j++)
			outdata << u[i][j] << " ";
		outdata << std::endl;
	}
	std::cout<<"Dane zapisano w pliku "<<title<<".txt"<<std::endl;
}

void saveToFile(double* a,std::string title){
	std::ofstream outdata;
	outdata.open(title);
	if(!outdata){
		std::cout<<"Nie uda³o siê zapisaæ danych do pliku "<<title<<std::endl;
		return;
	}
	for(int i=0;i<iterations;i++)
		outdata << a[i] << " ";
	std::cout<<"Dane zapisano w pliku "<<title<<std::endl;
}

double** alloc2d(){
	double** x = new double*[2*N-1];
	for(int i=0; i<2*N-1;i++){
		x[i] = new double[2*N-1];
	}
	return x;
}