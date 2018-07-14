#include <iostream>
#include <armadillo>
#include <time.h>
#include <cmath>

int main()
{
	int Nx = 10, Ny = 10, Nz = 10, Nxy, i, j, k;
	double a = 3.2094, c = 5.2108, rc = 6.5, A, B;
	int x0 = 0, y0 = 0, z0 = 0;
	arma::cube xla, yla, zla;
	xla.set_size(Ny, Nx, Nz);
	yla.set_size(Ny, Nx, Nz);
	zla.set_size(Ny, Nx, Nz);
	
	Nxy=Nx*Ny;
	
	struct timeval start, end;
	gettimeofday(&start, NULL);
	
	for (k = 0; k < Nz; k++) {
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Ny; j++)  
			{
				A = 1 - std::pow(-1, j - Ny/2+1);
				B = 1 - std::pow(-1, k - Nz/2+1);
				xla(j, i, k) = (i - Nx/2+1) * a + A * a / 4 + B * a / 4 + x0;
				yla(j, i, k) = (j - Ny/2+1) * a * sqrt(3) / 2 + B * a * sqrt(3) / 12 + y0;
				zla(j,i, k) = (k - Nz/2+1) * c / 2 + z0;
			}
		}
	} //xla, yla, zla
	
////	xla.save("xla.dat", arma::raw_ascii);
////	yla.save("yla.dat", arma::raw_ascii);
////	zla.save("zla.dat", arma::raw_ascii);

	// type 1:
	arma::cube rla;
	rla.set_size(Ny, Nx, Nz);
	rla = arma::sqrt(arma::square(xla, 2) + arma::square(yla, 2) + arma::square(zla, 2));
	arma::uvec ind_nb1 = arma::find ((rla <= rc)&&(rla>1.e-2));
	indz1=ind_nb1/Nxy;
	indx1=(ind_nb1-indz1*Nxy)/Ny;
	indy1=ind_nb1-indz1*Nxy-indx1*Ny;
	indx1 -= Nx/2;
	indy1 -= Ny/2;
	indz1 -= Nz/2;
	
	arma::umat indy_1_1 = indy_1.row(0), indx_1_1 = indx_1.row(0), indz_1_1 = indz_1.row(0); //��һ��Ϊ������
	arma::umat indy1 = indy_1_1 - (Ny/2 -1); //���⣺����֮��С������������˷�Χ Ӧ���Ƕ�����������Ͳ���!!!

//	//rla.print("����rla = \n");
//	//ind_nb1.print("ind_nb1= \n");
//	//indy1.print("indy1:");
//	//std::cout << indy1.size() << std::endl;

//	//type 2:
//	rla = arma::sqrt(pow(xla-a/2, 2) + pow(yla-a*sqrt(3)/2, 2) + pow(zla, 2));
//	arma::uvec ind_nb2 = arma::find(rla <= rc);
//	arma::umat indy_2 = arma::ind2sub(size(yla), ind_nb2), indx_2 = arma::ind2sub(size(xla), ind_nb2), indz_2 = arma::ind2sub(size(zla), ind_nb2);
//	arma::umat indy_2_2 = indy_2.row(0), indx_2_2 = indx_2.row(0), indz_2_2 = indz_2.row(0);

//	//indy_2_2.print("�����: \n");

//	//type 3:
//	rla = arma::sqrt(pow(xla - a / 2, 2) + pow(yla - a * sqrt(3) / 6, 2) + pow(zla - c / 2, 2));
//	arma::uvec ind_nb3 = arma::find(rla <= rc);
//	arma::umat indy_3 = arma::ind2sub(size(yla), ind_nb3), indx_3 = arma::ind2sub(size(xla), ind_nb3), indz_3 = arma::ind2sub(size(zla), ind_nb3);
//	arma::umat indy_3_3 = indy_3.row(0), indx_3_3 = indx_3.row(0), indz_3_3 = indz_3.row(0);
//	//indy_3_3.print("�����: \n");

//	//type 4:
//	rla = arma::sqrt(pow(xla - a, 2) + pow(yla - a * sqrt(3) * 2 / 3, 2) + pow(zla - c / 2, 2));
//	arma::uvec ind_nb4 = arma::find(rla <= rc);
//	arma::umat indy_4 = arma::ind2sub(size(yla), ind_nb4), indx_4 = arma::ind2sub(size(xla), ind_nb4), indz_4 = arma::ind2sub(size(zla), ind_nb4);
//	arma::umat indy_4_4 = indy_4.row(0), indx_4_4 = indx_4.row(0), indz_4_4 = indz_4.row(0);
//	//indy_4_4.print("�����: \n");

	gettimeofday(&end, NULL);
	std::cout<< ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6<<std::endl;

	return 0;
}