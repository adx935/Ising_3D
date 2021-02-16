// Ising Model in three dimensions

#define _CRT_SECURE_NO_DEPRECATE	//Para evitar que salga error por funcion fprintf()
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

//Generador de números aleatorios entre 0 y 1:
inline double std_rand()
{
	return rand() / (RAND_MAX + 1.0);
}

//Definición de variables:
double J = +1;		//Factor de interacción entre partículas
int Lx, Ly, Lz;		//Número de partículas en X, Y y Z
int N;				//Número de partículas o spines
int ***s;			//Spines
double T;			//Temperatura
double H;			//Campo Magnético
int steps = 0;		//Pasos de Metrópolis
double ratio = 0;	//Razón para aceptar o negar cambio de spin
int MCSteps = 0;	//Pasos de Monte Carlo

//Variables para cálculo de promedios a partir de una simulación a T = 4.5
double Tf = 4;		//Temperatura a obtener los promedios
double beta = 0;		//1/T
double betaCero = 0;	//1/T_o
double deltaBeta = 0;	//beta-beta_o
double M2 = 0;			//Magnetización al cuadrado
double E2 = 0;			//Energía al cuadrado

//Funciones:

//Factores de Boltzmann:
double w[25][3];
void computeBoltzmannFactors()
{
	//printf("Factores de Boltzmann: \n \n");
	for (int i = -12; i <= 12; i += 4) {
		w[i + 12][0] = exp(-(i * J + 2 * H) / T);
		w[i + 12][2] = exp(-(i * J - 2 * H) / T);

		//printf("%f \t %f \n", w[i + 12][0], w[i + 12][2]);
	}
	//printf("\n");
}

//Inicializa el arreglo con distribución aleatoria de spines:
void initializeRandom()
{
	//printf("Initial random spin array: \n \n");

	s = new int**[Lx];
	for (int i = 0; i < Lx; i++)
	{
		s[i] = new int*[Ly];
		for (int j = 0; j < Ly; j++)
		{
			s[i][j] = new int[Lz];

			for (int k = 0; k < Lz; k++)
			{
				double a = std_rand();
				if (a < 0.5) s[i][j][k] = +1;
				else s[i][j][k] = -1;
			//	printf("%d \t", s[i][j][k]);
			}
			//printf("\n");
		}
		//printf("\n");
		//printf("\n");
	}
	//printf("\n");

	computeBoltzmannFactors();
	steps = 0;
}

//Inicializa el arreglo con distribución de spines todos +1
void initializeUp() {

	printf("Initial spin array (all up): \n \n");

	s = new int**[Lx];
	for (int i = 0; i < Lx; i++)
	{
		s[i] = new int*[Ly];
		for (int j = 0; j < Ly; j++)
		{
			s[i][j] = new int[Lz];

			for (int k = 0; k < Lz; k++)
			{
				s[i][j][k] = +1;
				printf("%d \t", s[i][j][k]);
			}
			printf("\n");
		}
		printf("\n");
		printf("\n");
	}
	printf("\n");

	computeBoltzmannFactors();
	steps = 0;
}

//Paso de Metropolis:
bool MetropolisStep()
{
	//Se escoge un spin aleatorio
	int i = int(Lx*std_rand());
	int j = int(Ly*std_rand());
	int k = int(Lz*std_rand());

	//Se consiguen sus vecinos usando condiciones de borde periódicas
	int iPrev = i == 0 ? Lx - 1 : i - 1;
	int iNext = i == Lx - 1 ? 0 : i + 1;
	int jPrev = j == 0 ? Ly - 1 : j - 1;
	int jNext = j == Ly - 1 ? 0 : j + 1;
	int kPrev = k == 0 ? Lz - 1 : k - 1;
	int kNext = k == Lz - 1 ? 0 : k + 1;

	//Se consigue la suma de los vecinos
	int sumNeighbors = s[iPrev][j][k] + s[iNext][j][k] + s[i][jPrev][k] + s[i][jNext][k] + s[i][j][kPrev] + s[i][j][kNext];
	int delta_ss = 2 * s[i][j][k] * sumNeighbors;

	//Razón de los factores de Boltzmann
	ratio = w[delta_ss + 12][1 + s[i][j][k]];

	if (std_rand() < ratio) {
		s[i][j][k] = -s[i][j][k];

		//Imprime el arreglo de spines para visualizar el cambio de spin:
		/*for (int i = 0; i < Lx; i++)
		{
		for (int j = 0; j < Ly; j++)
		{
		printf("%d \t", s[i][j]);
		}
		printf("\n");
		}
		printf("\n");*/

		return true;
	}
	else return false;
}

//Paso de Monte Carlo:
double acceptanceRatio;
void oneMonteCarloStepPerSpin()
{
	int accepts = 0;
	for (int i = 0; i < N; i++)
		if (MetropolisStep())
			++accepts;
	acceptanceRatio = accepts / double(N);
	++steps;
}

//Cálculo de la Magnetización:
double magnetizationPerSpin() {
	int sSum = 0;
	for (int i = 0; i < Lx; i++)
		for (int j = 0; j < Ly; j++)
		{
			for (int k = 0; k < Lz; k++)
			{
				sSum += s[i][j][k];
			}
		}

	double M = (double)sSum / N;
	return M;
}

//Cálculo de la Energía:
double energyPerSpin() {
	int sSum = 0, ssSum = 0;
	for (int i = 0; i < Lx; i++)
	{
		for (int j = 0; j < Ly; j++)
		{
			for (int k = 0; k < Lz; k++)
			{
				sSum += s[i][j][k];
				int iNext = i == Lx - 1 ? 0 : i + 1;
				int jNext = j == Ly - 1 ? 0 : j + 1;
				int kNext = k == Lz - 1 ? 0 : k + 1;

				ssSum += s[i][j][k] * (s[iNext][j][k] + s[i][jNext][k] + s[i][j][kNext]);
			}
		}
	}

	double E = (double)-(J*ssSum + H*sSum) / N;
	return E;
}

FILE *point1;	//Para guardar M y E vs T en un archivo.
FILE *point2;	//Para guardar los promedios en un archivo.


int main(int argc, char *argv[]) {
	cout << " Three-dimensional Ising Model - Metropolis simulation\n"
		<< " ---------------------------------------------------\n"
		<< " Enter number of spins L in each direction: ";
	cin >> Lx;
	Ly = Lx;
	Lz = Lx;
	N = Lx*Ly*Lz;
	cout << " Enter temperature T: ";
	cin >> T;
	cout << " Enter magnetic field H: ";
	cin >> H;
	cout << " Enter number of Monte Carlo steps: ";
	MCSteps;
	cin >> MCSteps;

//	point1 = fopen("MyEvsT_3D.csv", "w");		//Abre archivo con M y E vs T.
//	fprintf(point1, "T \t M \t +-M \t E \t +-E \n");	//Impime el archivo con M y E vs T.
//	point2 = fopen("Promedios_Ising3D.csv", "w");		//Abre archivo con los promedios.
//	fprintf(point2, "T_f \t Mprom \t M^2prom \t Eprom \t E^2prom \t Chi \t C \n");	//Impime el archivo con los promedios.

//Código para cálculos con T fija:

		initializeRandom();
		//initializeUp();

		//Pasos para termalizar el sistema:
		int thermSteps = int(0.2 * MCSteps);
		cout << " Performing " << thermSteps
			<< " steps to thermalize the system ..." << flush;
		for (int s = 0; s < thermSteps; s++)
			oneMonteCarloStepPerSpin();

		//Pasos de producción de la simulación:
		cout << " Done\n Performing production steps ... \n" << flush;
		double mAv = 0, m2Av = 0, eAv = 0, e2Av = 0;
		ofstream file("ising3D.data");

//		printf("T_f \t Mprom \t M^2prom \t Eprom \t E^2prom \t Chi \t C \n");

		//Loop para cálculos de promedios a otras T a partir de una simulación a T=4.5:
/*		for (int l = 0; l < 50; l++)
		{
			betaCero = 1 / T;
			beta = 1 / Tf;
			deltaBeta = beta - betaCero;

			double Mprom = 0;		//Magnetización promedio
			double Eprom = 0;		//Energía promedio
			double M2prom = 0;		//Magnetización al cuadrado promedio
			double E2prom = 0;		//Energía al cuadrado promedio
			double ZetaT = 0;		//Sumatoria en MCsteps de la función de partición
			double Chi = 0;			//Susceptibilidad magnética
			double C = 0;			//Calor específico
*/
			for (int s = 0; s < MCSteps; s++)
			{
				oneMonteCarloStepPerSpin();
				double m = magnetizationPerSpin();
				double e = energyPerSpin();
//				double Etot = e*N;
				mAv += m; m2Av += m * m;
				eAv += e; e2Av += e * e;
				file << m << '\t' << e << '\n';

				//Cálculos de promedios a otras T a partir de una simulación a T=4.5:

				/*M2 = pow(m, 2);
				E2 = pow(e, 2);

				Mprom += abs(m)*exp(-1*deltaBeta*Etot);
				Eprom += e*exp(-1*deltaBeta*Etot);
				M2prom += M2*exp(-1*deltaBeta*Etot);
				E2prom += E2*exp(-1*deltaBeta*Etot);
				ZetaT += exp(-1*deltaBeta*Etot);*/
			}

/*			Mprom = Mprom / ZetaT;
			Eprom = Eprom / ZetaT;
			M2prom = M2prom / ZetaT;
			E2prom = E2prom / ZetaT;

			Chi = beta*(M2prom - pow(Mprom, 2));
			C = pow(beta, 2)*(E2prom - pow(Eprom, 2));

			fprintf(point2, "%f \t %f \t %f \t %f \t %f \t %f \t %f \n", Tf, Mprom, M2prom, Eprom, E2prom, Chi, C);		//Impime el archivo con los promedios.
			printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \n", Tf, Mprom, M2prom, Eprom, E2prom, Chi, C);

			Tf += 0.02;
			betaCero = 0;
			beta = 0;
			deltaBeta = 0;
		}*/

//		fclose(point2);				//Cierra el archivo con los promedios.

		file.close();
		mAv /= MCSteps; m2Av /= MCSteps;
		eAv /= MCSteps; e2Av /= MCSteps;
		double incertM = sqrt(m2Av - mAv*mAv);
		double incertE = sqrt(m2Av - mAv*mAv);

		cout << " \n Magnetization and energy per spin written in file "
			<< " \"ising3D.data\"" << endl;
		cout << " <m> = " << mAv << " +/- " << incertM << endl;
		cout << " <e> = " << eAv << " +/- " << incertE << endl;


//		fprintf(point1, "%f \t %f \t %f \t %f \t %f \n", T, mAv, incertM, eAv, incertE);		//Impime el archivo con M y E vs T.
													
	//Loop para hacer los cálculos para varias T:
/*	for (int k = 0; k < 200; k++)
	{
		T += 0.05;
		printf("\n ----- T=%f ----- \n", T);
		initializeRandom();
		//initializeUp();

		//Pasos para termalizar el sistema:
		int thermSteps = int(0.2 * MCSteps);
		cout << " Performing " << thermSteps
			<< " steps to thermalize the system ..." << flush;
		for (int s = 0; s < thermSteps; s++)
			oneMonteCarloStepPerSpin();

		//Pasos de producción de la simulación:
		cout << " Done\n Performing production steps ... \n" << flush;
		double mAv = 0, m2Av = 0, eAv = 0, e2Av = 0;
		ofstream file("ising3D.data");
		for (int s = 0; s < MCSteps; s++)
		{
			oneMonteCarloStepPerSpin();
			double m = magnetizationPerSpin();
			double e = energyPerSpin();
			mAv += m; m2Av += m * m;
			eAv += e; e2Av += e * e;
			file << m << '\t' << e << '\n';
		}

		file.close();
		mAv /= MCSteps; m2Av /= MCSteps;
		mAv = abs(mAv);							//Queda el valor absoluto de la Magnetización promedio
		eAv /= MCSteps; e2Av /= MCSteps;
		double incertM = sqrt(m2Av - mAv*mAv);
		double incertE = sqrt(m2Av - mAv*mAv);

		cout << " \n Magnetization and energy per spin written in file "
			<< " \"ising3D.data\"" << endl;
		cout << " <m> = " << mAv << " +/- " << incertM << endl;
		cout << " <e> = " << eAv << " +/- " << incertE << endl;


		fprintf(point1, "%f \t %f \t %f \t %f \t %f \n", T, mAv, incertM, eAv, incertE);			//Impime el archivo con M y E vs T.
		
	}*/

//	fclose(point1);				//Cierra el archivo con M y E vs T.

	printf("\n DONE! Click 'Enter' to close the program.");

	getchar();
	getchar();
}