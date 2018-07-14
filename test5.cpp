#include <iostream>
#include <armadillo>
#include <vector>

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
				xla(i, j, k) = (j - 4) * a + A * a / 4 + B * a / 4 + x0;
			}
			
		}
	} //xla

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				double B = 1 - pow(-1, k);
				yla(i, j, k) = (i - 4) * a * sqrt(3) / 2 + B * a * sqrt(3) / 12 + y0;
			}

		}
	} //yla

	for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Nx; j++) {
			for (int i = 0; i < Ny; i++) {
				zla(i, j, k) = (k - 4) * c / 2 + z0;
			}

		}
	} //zla
	
	//xla.print("����xla = \n");
	//yla.print("����yla = \n");
	//zla.print("����zla = \n");

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

	//rla.print("����rla = \n");
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

	ix += 4, iy += 4, iz += 4;
	//ix.print("ix="), iy.print("iy="), iz.print("iz=");

	int numy, numx, numz;
	arma::imat nbs_x, nbs_y, nbs_z;
	//arma::vec nbs;
	arma::umat nbs;//����Ϊimat������
	//int s = 1;
	for (int s = 0; s < Ny*Nx*Nz; s++) 
	{
	    //int i, j, k;
	    //k = s / (Nx * Ny), i = (s - k * (Nx*Ny)) / Ny, j = s - k * (Nx*Ny) - i * Ny;
		numy = iy(s), numx = ix(s), numz = iz(s);//�ӡ�arma::s64��ת������int�������ܶ�ʧ����

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

		//int B = nbs_x.n_elem;//�ӡ�const arma::uword��ת������int�������ܶ�ʧ����
		int A = 56;
		nbs.set_size(Ny*Nx*Nz, A);
		//nbs.set_size(1, A);
		for (int a = 0; a < A; a++) 
		{
			int bx = nbs_x(a), by = nbs_y(a), bz = nbs_z(a);//�ӡ�arma::s64��ת������int�������ܶ�ʧ����
			if (bx < 0) { nbs_x(a) = (bx+Nx) % Nx ; }
			else { nbs_x(a) = bx % Nx ; }
			if (by < 0) { nbs_y(a) = (by+Ny) % Ny ; }
			else { nbs_y(a) = by % Ny ; }
			if (bz < 0) { nbs_z(a) = (bz+Nz) % Nz ; }
			else { nbs_z(a) = bz % Nz ; }
			
			int cx = nbs_x(a), cy = nbs_y(a), cz = nbs_z(a);
			//nbs(s,a) = arma::find(ix == cx && iy == cy && iz == cz);
			nbs(s, a) = nbs_z(a) * Nx * Ny + nbs_x(a) * Ny + nbs_y(a);
			
		}
	}
	//arma::uvec nbs = arma::find(ix == 9 && iy == 9 && iz == 8);
	//nbs_x.print("nbs_x="), nbs_y.print("nbs_y"), nbs_z.print("nbs_z");
//	nbs.print("nbs");
	
//	arma::umat xx(Nx*Ny*Nz,56);
//	for (int i=0;i<Nx*Ny*Nz*56;i++)
//		xx.at(i) = (nbs.at(i));
	
	
//	xx.print("xx");
	xla(nbs(arma::span(0),arma::span(0,15))).print("xla(nbs)");
//	xla(nbs(arma::span(0),arma::span(0,9))).print("xla(nbs(0,0:9))");

//	arma::ivec xnb = (nbs.row(0));
//	arma::urowvec xxnb = arma::urowvec(xnb, 3, false) ;

//	arma::uword __aux_urowvec_1 [] = {2, 1, 3} ;
//	_aux_urowvec_1 = urowvec(__aux_urowvec_1, 3, false) ;
	
	
	
//	nbs.save("nbs.dat", arma::raw_ascii);
	
	return 0;
}