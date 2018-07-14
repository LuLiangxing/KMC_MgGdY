#include <iostream>
#include <armadillo>
#include <vector>
#include <fstream>

int main()
{
	int Nx = 10, Ny = 10, Nz = 10;
	double a = 3.2094, c = 5.2108, rc = 6.5;
	int x0 = 0, y0 = 0, z0 = 0;
	arma::cube xla, yla, zla;
	xla.set_size(Ny, Nx, Nz);
	yla.set_size(Ny, Nx, Nz);
	zla.set_size(Ny, Nx, Nz);
	
	for (int k = 0; k < Nz; k++){
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				double A = 1 - pow(-1, i);
				double B = 1 - pow(-1, k);
				xla(i, j, k) = (j - (Ny/2-1)) * a + A * a / 4 + B * a / 4 + x0;
			}
			
		}
	} //xla

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				double B = 1 - pow(-1, k);
				yla(i, j, k) = (i - (Nx/2-1)) * a * sqrt(3) / 2 + B * a * sqrt(3) / 12 + y0;
			}

		}
	} //yla

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				zla(i, j, k) = (k - (Nz/2-1)) * c / 2 + z0;
			}

		}
	} //zla
	
	//xla.print("æÿ’Ûxla = \n");
	//yla.print("æÿ’Ûyla = \n");
	//zla.print("æÿ’Ûzla = \n");

	arma::icube ix, iy, iz;
	ix.set_size(Ny, Nx, Nz);
	iy.set_size(Ny, Nx, Nz);
	iz.set_size(Ny, Nx, Nz);

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				ix(i, j, k) = (j - 4);
			}

		}
	} 

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				iy(i, j, k) = (i - 4) ;
			}

		}
	} 

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				iz(i, j, k) = (k - 4) ;
			}

		}
	}

	//ix.print("ix="), iy.print("iy="), iz.print("iz=");

	//tape 1:
	arma::cube rla;
	rla.set_size(Ny, Nx, Nz);
	rla = arma::sqrt(pow(xla, 2) + pow(yla, 2) + pow(zla, 2));
	arma::uvec ind_nb1 = arma::find (rla <= rc && rla>0);
	arma::mat rnb1 = rla.elem(ind_nb1);
	arma::imat indy1 = iy.elem(ind_nb1), indx1 = ix.elem(ind_nb1), indz1 = iz.elem(ind_nb1);
	//arma::umat indz_1 = ind_nb1 / (Nx*Ny), indx_1 = (ind_nb1 - indz_1 * Nx*Ny) / Ny, indy_1 = ind_nb1 - indz_1 * Nx*Ny - indx_1 * Ny;
	//arma::mat indz_1_1 = indz_1;

	//rla.print("æÿ’Ûrla = \n");
	//ind_nb1.print("ind_nb1= \n");
	//indy1.print("indy1= \n"), indx1.print("indx1= \n"), indz1.print("indz1= \n");
	//rnb1.print("rnb1:");
	//std::cout << indy_1.size() << std::endl;

	//tape 2:
	rla = arma::sqrt(pow(xla-a/2, 2) + pow(yla-a*sqrt(3)/2, 2) + pow(zla, 2));
	arma::uvec ind_nb2 = arma::find(rla <= rc && rla>0);
	arma::mat rnb2 = rla.elem(ind_nb2);
	arma::imat indy2 = iy.elem(ind_nb2)-1, indx2 = ix.elem(ind_nb2), indz2 = iz.elem(ind_nb2);
	//arma::umat indy_2 = arma::ind2sub(size(yla), ind_nb2), indx_2 = arma::ind2sub(size(xla), ind_nb2), indz_2 = arma::ind2sub(size(zla), ind_nb2);
	//arma::umat indy_2_2 = indy_2.row(0), indx_2_2 = indx_2.row(0), indz_2_2 = indz_2.row(0);
	//indy2.print("indy2= \n"), indx2.print("indx2= \n"), indz2.print("indz2= \n");

	//tape 3:
	rla = arma::sqrt(pow(xla - a / 2, 2) + pow(yla - a * sqrt(3) / 6, 2) + pow(zla - c / 2, 2));
	arma::uvec ind_nb3 = arma::find(rla <= rc && rla>0);
	arma::mat rnb3 = rla.elem(ind_nb3);
	arma::imat indy3 = iy.elem(ind_nb3), indx3 = ix.elem(ind_nb3), indz3 = iz.elem(ind_nb3)-1;
	//arma::umat indy_3 = arma::ind2sub(size(yla), ind_nb3), indx_3 = arma::ind2sub(size(xla), ind_nb3), indz_3 = arma::ind2sub(size(zla), ind_nb3);
	//arma::umat indy_3_3 = indy_3.row(0), indx_3_3 = indx_3.row(0), indz_3_3 = indz_3.row(0);
	//indy3.print("indy3= \n"), indx3.print("indx3= \n"), indz3.print("indz3= \n");

	//tape 4:
	rla = arma::sqrt(pow(xla - a, 2) + pow(yla - a * sqrt(3) * 2 / 3, 2) + pow(zla - c / 2, 2));
	arma::uvec ind_nb4 = arma::find(rla <= rc && rla>0);
	arma::mat rnb4 = rla.elem(ind_nb4);
	arma::imat indy4 = iy.elem(ind_nb4)-1, indx4 = ix.elem(ind_nb4), indz4 = iz.elem(ind_nb4)-1;
	//arma::umat indy_4 = arma::ind2sub(size(yla), ind_nb4), indx_4 = arma::ind2sub(size(xla), ind_nb4), indz_4 = arma::ind2sub(size(zla), ind_nb4);
	//arma::umat indy_4_4 = indy_4.row(0), indx_4_4 = indx_4.row(0), indz_4_4 = indz_4.row(0);
	//indy4.print("indy4= \n"), indx4.print("indx4= \n"), indz4.print("indz4= \n");

	ix += (Nx/2 - 1), iy += (Ny/2 - 1), iz += (Nz/2 - 1);
	//ix.print("ix="), iy.print("iy="), iz.print("iz=");

	int numy, numx, numz;
	int A = indx1.n_elem;
	arma::imat nbs_x, nbs_y, nbs_z;
	//arma::uvec nbs;
	arma::umat nbs;
	//int s = 1;
	for (int s = 0; s < Ny*Nx*Nz; s++) 
	{
	    //int i, j, k;
	    //k = s / (Nx * Ny), i = (s - k * (Nx*Ny)) / Ny, j = s - k * (Nx*Ny) - i * Ny;
		numy = iy(s), numx = ix(s), numz = iz(s);

		if (( numy - 2*(numy/2) == 0) && (numz - 2*(numz/2) == 0)) //tape1
		{
			nbs_x = numx + indx1, nbs_y = numy + indy1, nbs_z = numz + indz1;
		}
		else if ((numy - 2*(numy/2) == 1) && (numz - 2*(numz/2) == 0)) //tape2
		{
			nbs_x = numx + indx2, nbs_y = numy + indy2, nbs_z = numz + indz2;
		}
		else if ((numy - 2 * (numy / 2) == 0) && (numz - 2 * (numz / 2) == 1)) //tape3
		{
			nbs_x = numx + indx3, nbs_y = numy + indy3, nbs_z = numz + indz3;
		}
		else if ((numy - 2 * (numy / 2) == 1) && (numz - 2 * (numz / 2) == 1)) //tape4
		{
			nbs_x = numx + indx4, nbs_y = numy + indy4, nbs_z = numz + indz4;
		}
		//nbs_x.print("nbs_x="), nbs_y.print("nbs_y"), nbs_z.print("nbs_z");

		nbs.set_size(Ny*Nx*Nz, A);
		//nbs.set_size(1, A);
		for (int a = 0; a < A; a++) 
		{
			int bx = nbs_x(a), by = nbs_y(a), bz = nbs_z(a);
			if (bx < 0) { nbs_x(a) = (bx+Nx) % Nx ; }
			else { nbs_x(a) = bx % Nx ; }
			if (by < 0) { nbs_y(a) = (by+Ny) % Ny ; }
			else { nbs_y(a) = by % Ny ; }
			if (bz < 0) { nbs_z(a) = (bz+Nz) % Nz ; }
			else { nbs_z(a) = bz % Nz ; }
			
			//int cx = nbs_x(a), cy = nbs_y(a), cz = nbs_z(a);
			//nbs(s,a) = arma::find(ix == cx && iy == cy && iz == cz);
			nbs(s, a) = nbs_z(a) * Nx * Ny + nbs_x(a) * Ny + nbs_y(a);
			
		}
	}

	//nbs_x.print("nbs_x="), nbs_y.print("nbs_y"), nbs_z.print("nbs_z");
	//nbs.print("nbs");
	//nbs.save("nbs.text", arma::arma_ascii);
	//arma::mat test;
	//test.load("test.txt", arma::auto_detect);
	//test.print("test=");
	arma::mat FMg, FMgMg, FMgY, FY, FYY, r, ro, Ro_Mg, Ro_Y;
	FMg.load("FMg.txt", arma::auto_detect), FMgMg.load("FMgMg.txt", arma::auto_detect), FMgY.load("FMgY.txt", arma::auto_detect);
	FY.load("FY.txt", arma::auto_detect), FYY.load("FYY.txt", arma::auto_detect), r.load("r.txt", arma::auto_detect);
	ro.load("ro.txt", arma::auto_detect), Ro_Mg.load("Ro_Mg.txt", arma::auto_detect), Ro_Y.load("Ro_Y.txt", arma::auto_detect);
	//FMg.print("FMg=");Ro_Y.print("Ro_Y=");
	FMgMg = FMgMg / 3.6; FYY = FYY / 3.6; FMgY = FMgY / 3.6;
	//FMgMg.print("FMgMg=");
	arma::mat roMgnb1, roYnb1, FMgMgnb1, FMgYnb1, FYYnb1;
	arma::interp1(r, Ro_Mg, rnb1, roMgnb1);
	arma::interp1(r, Ro_Y, rnb1, roYnb1);
	arma::interp1(r, FMgMg, rnb1, FMgMgnb1);
	arma::interp1(r, FMgY, rnb1, FMgYnb1);
	arma::interp1(r, FYY, rnb1, FYYnb1);
	//roMgnb1.print("roMgnb1="); FYYnb1.print("FYYnb1=");

	arma::mat roMgnb2, roYnb2, FMgMgnb2, FMgYnb2, FYYnb2;
	arma::interp1(r, Ro_Mg, rnb2, roMgnb2);
	arma::interp1(r, Ro_Y, rnb2, roYnb2);
	arma::interp1(r, FMgMg, rnb2, FMgMgnb2);
	arma::interp1(r, FMgY, rnb2, FMgYnb2);
	arma::interp1(r, FYY, rnb2, FYYnb2);
	//roMgnb2.print("roMgnb2="); FYYnb2.print("FYYnb2=");

	arma::mat roMgnb3, roYnb3, FMgMgnb3, FMgYnb3, FYYnb3;
	arma::interp1(r, Ro_Mg, rnb3, roMgnb3);
	arma::interp1(r, Ro_Y, rnb3, roYnb3);
	arma::interp1(r, FMgMg, rnb3, FMgMgnb3);
	arma::interp1(r, FMgY, rnb3, FMgYnb3);
	arma::interp1(r, FYY, rnb3, FYYnb3);
	//roMgnb3.print("roMgnb3="); FYYnb3.print("FYYnb3=");

	arma::mat roMgnb4, roYnb4, FMgMgnb4, FMgYnb4, FYYnb4;
	arma::interp1(r, Ro_Mg, rnb4, roMgnb4);
	arma::interp1(r, Ro_Y, rnb4, roYnb4);
	arma::interp1(r, FMgMg, rnb4, FMgMgnb4);
	arma::interp1(r, FMgY, rnb4, FMgYnb4);
	arma::interp1(r, FYY, rnb4, FYYnb4);
	//roMgnb4.print("roMgnb4="); FYYnb4.print("FYYnb4=");

	arma::mat F, elem;
	elem.ones(Ny*Nx*Nz, 1);
	F.zeros(1, Ny*Nx*Nz);
	elem(444,0) = 2;
	//elem.print("elem=");
	//int s = 3; 
	for (int s = 0; s < Ny*Nx*Nz; s++) 
	{
		numy = iy(s), numx = ix(s), numz = iz(s);
		arma::mat elemMg(1, A), elemY(1, A);
		for (int a = 0; a < A; a++)
		{
			int B = nbs(s, a);
			if (elem(B) == 1)
			{
				elemMg(0, a) = 1;
			}
			else { elemMg(0, a) = 0; }
		}
		for (int a = 0; a < A; a++)
		{
			int B = nbs(s, a);
			if (elem(B) == 2)
			{
				elemY(0, a) = 1;
			}
			else { elemY(0, a) = 0; }
		}
		//elemMg.print("elemMg=");
		//elemY.print("elemY=");
		if ((numy - 2 * (numy / 2) == 0) && (numz - 2 * (numz / 2) == 0))//tape1
		{
			arma::mat rola1 = elemMg * roMgnb1 + elemY * roYnb1;
			arma::mat F1, F2;
			//rola1.print("rola1=");
			if (elem(s) == 1) 
			{
				arma::interp1(ro, FMg, rola1, F1);
				//F1.print("F1=");
				F1 = F1 + (elemMg * FMgMgnb1)*0.5 + elemY * FMgYnb1;
				//F1.print("F1=");
				F(0, s) = F1(0);
				//F.print("F=");
			}
			else if (elem(s) == 2) { arma::interp1(ro, FY, rola1, F2); F2 = F2 + (elemMg * FMgYnb1)*0.5 + elemY * FYYnb1; F(0, s) = F2(0);}
		}
		else if ((numy - 2 * (numy / 2) == 1) && (numz - 2 * (numz / 2) == 0))//tape2
		{
			arma::mat rola2 = elemMg * roMgnb2 + elemY * roYnb2;
			arma::mat F1, F2;
			//rola2.print("rola2=");
			if (elem(s) == 1)
			{
				arma::interp1(ro, FMg, rola2, F1);
				F1 = F1 + (elemMg * FMgMgnb2)*0.5 + elemY * FMgYnb2;
				F(0, s) = F1(0);
			}
			else if (elem(s) == 2) { arma::interp1(ro, FY, rola2, F2); F2 = F2 + (elemMg * FMgYnb2)*0.5 + elemY * FYYnb2; F(0, s) = F2(0); }
		}
		else if ((numy - 2 * (numy / 2) == 0) && (numz - 2 * (numz / 2) == 1))//tape3
		{
			arma::mat rola3 = elemMg * roMgnb3 + elemY * roYnb3;
			arma::mat F1, F2;
			//rola3.print("rola3=");
			if (elem(s) == 1)
			{
				arma::interp1(ro, FMg, rola3, F1);
				F1 = F1 + (elemMg * FMgMgnb3)*0.5 + elemY * FMgYnb3;
				F(0, s) = F1(0);
			}
			else if (elem(s) == 2) { arma::interp1(ro, FY, rola3, F2); F2 = F2 + (elemMg * FMgYnb3)*0.5 + elemY * FYYnb3; F(0, s) = F2(0); }
		}
		else if ((numy - 2 * (numy / 2) == 1) && (numz - 2 * (numz / 2) == 1))//tape4
		{
			arma::mat rola4 = elemMg * roMgnb4 + elemY * roYnb4;
			arma::mat F1, F2;
			//rola4.print("rola4=");
			if (elem(s) == 1)
			{
				arma::interp1(ro, FMg, rola4, F1);
				F1 = F1 + (elemMg * FMgMgnb4)*0.5 + elemY * FMgYnb4;
				F(0, s) = F1(0);
			}
			else if (elem(s) == 2) { arma::interp1(ro, FY, rola4, F2); F2 = F2 + (elemMg * FMgYnb4)*0.5 + elemY * FYYnb4; F(0, s) = F2(0); }
		}
		//F.print("F=");
	}
	//F.print("F=");
	double Ftot=0;
	for (int s = 0; s < Ny*Nx*Nz; s++) 
	{
		Ftot += F(s);
	}
	arma::cout << Ftot << arma::endl;

	return 0;
}