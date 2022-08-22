// headers C
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <string>
#include <unistd.h>
// headers Lorene
#include "et_rot_mag.h"
#include "eos.h"
#include "utilitaires.h"
//#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"
#include "metric.h"	    
#include <tensor.h>
#include <tbl.h>
namespace Lorene{
// Local prototype (for drawings only)
Cmp raccord_c1(const Cmp& uu, int l1) ; 
}

using namespace Lorene ;

// Definition de f_j et de M_j
Cmp f_j(const Cmp& x, double a_j){ 
  Cmp resu(x.get_mp() ) ;
  resu = a_j ;
  return resu ;
}
Cmp M_j(const Cmp& x, double a_j){ 
  return - a_j*x ;
}
double r(Et_magnetisation  star,double theta)
{
  const Map& mp=star.get_mp() ;
  const Mg3d& mg = *(mp.get_mg()) ;
  int type_t = mg.get_type_t() ; 
 #ifndef NDEBUG
  int type_p = mg.get_type_p() ; 
 #endif
  int nt = mg.get_nt(0) ;        
  assert( (type_t == SYM) || (type_t == NONSYM) ) ; 
  assert( (type_p == SYM) || (type_p == NONSYM) ) ; 
  int k = 0 ; 
  int j = (type_t == SYM ? nt-1 : nt / 2); 
  int l = star.l_surf()(k, j) ; 
  double xi = star.xi_surf()(k, j) ; 
  return mp.val_r(l, xi, theta, 0); 
}

double energy(Et_magnetisation &star, int num)
{
  double r=star.ray_pole();
  double wyn=0;
  for (int i = 0; i < num; i++)
  {
    wyn+=star.get_ener()().val_point(r*i/num,0,0)*4*M_PI*r/num*pow(r*i/num,2);
  }
  return wyn*Lorene::Unites::rhonuc_si*pow(Lorene::Unites::c_si,2)*pow(Lorene::Unites::r_unit,3);
}

double energy_domain(Et_magnetisation &star, int num)
{
  Tenseur energy_dens=star.get_a_car()*star.get_ener()*star.get_bbb();
  energy_dens.set_std_base();
  Tbl wyn=energy_dens().integrale_domains();
  return wyn(num)*Lorene::Unites::rhonuc_si*pow(Lorene::Unites::c_si,2)*pow(Lorene::Unites::r_unit,3);
}

double energy(Et_magnetisation &star)
{
  Tenseur energy_dens=star.get_a_car()*star.get_ener()*star.get_bbb();
  energy_dens.set_std_base();
  double *wyn=new double(energy_dens().integrale());
  return (*wyn)*Lorene::Unites::rhonuc_si*pow(Lorene::Unites::c_si,2)*pow(Lorene::Unites::r_unit,3);
}
void wy(ofstream& out,Et_magnetisation  star,int nt)
{
  //double theta[nt];
  out<<star.MagMom()/1.e32<<" ";
  out<<star.Magn()(0)(0,0,0,0)*Lorene::Unites_mag::mag_unit/1.e13<<" ";
  for (int i=0;i<nt;i++)
  {
    out<<r(star,i*M_PI/nt)*10<<" ";
  }
  out<<endl;
}
int core(int i,int num)
{
  ofstream out;
  ofstream fout;
  out.open("temp"+to_string(num)+".txt",ios_base::app);
  fout.open("surf"+to_string(num)+".txt",ios_base::app);
  using namespace Unites_mag ;
  char blabla[120] ;
  int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
  delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
  nt, np ; 
  double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, thres_adapt, aexp_mass, relax, relax_poisson, precis_adapt ;  
  double Q0, a_j0, Q_ini, a_j_ini, magmom_wanted ;
  int mer_mag, mer_change_mag, mer_fix_mag, mag_in_eos, mer_magmom, use_magnetisation ;

  //------------------------------------------------------------------
  //	    Parameters of the computation 
  //------------------------------------------------------------------
  ifstream fich("parrot.d") ;
  fich.getline(blabla, 120) ;
  fich >> relat_i ; fich.getline(blabla, 120) ;
  bool relat = (relat_i == 1) ; 
  fich >> ent_c ; fich.getline(blabla, 120) ;
  ent_c+= -i*0.01;
  fich >> freq_si ; fich.getline(blabla, 120) ;
  fich >> fact_omega ; fich.getline(blabla, 120) ;
  fich >> mbar_wanted ; fich.getline(blabla, 120) ;
  mbar_wanted *= msol ; 
  fich >> magmom_wanted ; fich.getline(blabla, 120) ;
  fich.getline(blabla, 120) ;
  fich >> Q0 ; fich.getline(blabla,120) ;
  fich >> a_j0 ; fich.getline(blabla,120) ;
  fich >> Q_ini; fich.getline(blabla,120) ;
  fich >> a_j_ini ; fich.getline(blabla,120) ;
  fich >> mer_mag ; fich.getline(blabla,120) ;
  fich >> mer_change_mag ; fich.getline(blabla,120) ;
  fich >> mer_fix_mag ; fich.getline(blabla,120);
  fich >> mag_in_eos ; fich.getline(blabla,120);
  fich >> use_magnetisation ; fich.getline(blabla,120);
  fich.getline(blabla,120);
  fich >> mer_max ; fich.getline(blabla, 120) ;
  fich >> precis ; fich.getline(blabla, 120) ;
  fich >> mer_rot ; fich.getline(blabla, 120) ;
  fich >> freq_ini_si ; fich.getline(blabla, 120) ;
  fich >> mer_change_omega ; fich.getline(blabla, 120) ;
  fich >> mer_fix_omega ; fich.getline(blabla, 120) ;
  fich >> delta_mer_kep ; fich.getline(blabla, 120) ;
  fich >> thres_adapt ; fich.getline(blabla, 120) ;
  fich >> mer_mass ; fich.getline(blabla, 120) ;
  fich >> mer_magmom ; fich.getline(blabla, 120) ;
  fich >> aexp_mass ; fich.getline(blabla, 120) ;
  fich >> relax ; fich.getline(blabla, 120) ;
  fich >> mermax_poisson ; fich.getline(blabla, 120) ;
  fich >> relax_poisson ; fich.getline(blabla, 120) ;
  fich >> precis_adapt ; fich.getline(blabla, 120) ;
  fich >> graph ; fich.getline(blabla, 120) ;
  fich.getline(blabla, 120) ;
  fich >> nz ; fich.getline(blabla, 120) ;
  fich >> nzet; fich.getline(blabla, 120) ;
  fich >> nzadapt; fich.getline(blabla, 120) ;
  fich >> nt; fich.getline(blabla, 120) ;
  fich >> np; fich.getline(blabla, 120) ;
  if (i==0)
  {
    ofstream out1;
    out1.open("results.txt");
    out1<<a_j0<<endl;
    out1<<"R_circ"<< setw(15) << "R_eq" << setw(15) << "R_pole"<< setw(15) 
  << "M_g"<< setw(15) << "M_b"<< setw(15) << "MM(1e32)"  << setw(15) << "Cent"<<setw(15)<<"RPolar"<<setw(15)<<"Teq"<<setw(15)<<"grv2_err"<<setw(15)<<"grv3_err"<<setw(15);
  out1<<"number_den"<<setw(15)<<" Ttoal energy "<<setw(15)<<" central pressure "<<setw(15)<<" central energy density "<<setw(15)<<" quadrupole momentum "<<endl;
    out1.close();
  }
  
  int* nr = new int[nz];
  int* nt_tab = new int[nz];
  int* np_tab = new int[nz];
  double* bornes = new double[nz+1];
     
  fich.getline(blabla, 120);
  for (int l=0; l<nz; l++) 
  {
	  fich >> nr[l]; 
	  fich >> bornes[l]; fich.getline(blabla, 120) ;
	  np_tab[l] = np ; 
	  nt_tab[l] = nt ; 
  }
  bornes[nz] = __infinity ;

  Tbl ent_limit(nzet) ;
  ent_limit.set_etat_qcq() ;
  ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
  for (int l=0; l<nzet-1; l++) 
  {
    fich >> ent_limit.set(l) ; fich.getline(blabla, 120) ;
  }


  fich.close();
      
  // Particular cases
  // ----------------

  // Initial frequency = final frequency
  if ( freq_ini_si < 0 ) 
  {
	  freq_ini_si = freq_si ; 
	  mer_change_omega = mer_rot ; 
	  mer_fix_omega = mer_rot + 1 ;  
  }

    
  //-----------------------------------------------------------------------
  //		Equation of state
  //-----------------------------------------------------------------------

  fich.open("par_eos.d") ;
  //unique_ptr<Eos>  peos(Eos::eos_from_file(fich) );
  Eos* peos = Eos::eos_from_file(fich) ;
  Eos& eos = *peos ;

  fich.close() ;

  //-----------------------------------------------------------------------
  //		Construction of the multi-grid and the mapping
  //-----------------------------------------------------------------------

  // Rescale of bornes in the case where there more than 1 domain inside
  //   the star

  for (int l=0; l<nzet-1; l++) 
  {
    bornes[l+1] = bornes[nzet] * sqrt(1 - ent_limit(l) / ent_c) ;
  }

  // Type of r sampling :
  int* type_r = new int[nz];
  type_r[0] = RARE ; 
  for (int l=1; l<nz-1; l++)
  {
	  type_r[l] = FIN ; 
  }
  type_r[nz-1] = UNSURR ; 
    
  // Type of sampling in theta and phi :
  int type_t = SYM ; 
  int type_p = SYM ; 
    
  Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

  Map_et mp(mg, bornes) ;
   
  // Cleaning
  // --------
      
  delete [] nr ; 
  delete [] nt_tab ; 
  delete [] np_tab ; 
  delete [] type_r ; 
  delete [] bornes ; 

      
  cout << endl 
	<< "==========================================================" << endl
	<< "                    Physical parameters                   " << endl
	<< "=========================================================="
	<< endl ; 
  cout << endl ;

  cout << endl << "Equation of state : " 
	<< endl << "=================   " << endl ;
  cout << eos << endl ; 

  cout << "Central enthalpy : " << ent_c << " c^2" << endl ; 
  cout << "Rotation frequency : " << freq_si << " Hz" << endl ; 
  if ( abs(mer_mass) < mer_max ) {
	  cout << "Required Baryon mass [M_sol] : " 
	  << mbar_wanted / msol << endl ; 
  }
  if ( mer_magmom < mer_max ) {
	  cout << "Required magnetic moment [SI] : " 
	  << magmom_wanted << endl ; 
  }
  cout << endl 
	<< "==========================================================" << endl
	<< "               Computational parameters                   " << endl
	<< "=========================================================="
	<< endl << endl ; 

  cout << "Maximum number of steps in the main iteration : " 
	<< mer_max << endl ; 
  cout << "Relaxation factor in the main iteration  : " 
	<< relax << endl ; 
  cout << "Threshold on the enthalpy relative change for ending the computation : " 
	<< precis << endl ; 
  cout << "Maximum number of steps in Map_et::poisson : " 
	<< mermax_poisson << endl ; 
  cout << "Relaxation factor in Map_et::poisson : " 
	<< relax_poisson << endl ; 
  cout << "Step from which the baryon mass is forced to converge : " 
	<< mer_mass << endl ; 
  cout << "Step from which the magnetic moment is forced to converge : " 
	<< mer_magmom << endl ; 
  cout << "Exponent for the increase factor of the central enthalpy/CFA : " 
	<< aexp_mass << endl ; 
  cout << 
  "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping"
  << endl << thres_adapt << endl ; 


  cout << endl << "Multi-grid : " 
	<< endl << "==========" << endl << mg << endl ; 
  cout << "Mapping : " 
	<< endl << "=======" << endl << mp << endl ; 
  

  //-----------------------------------------------------------------------
  //		Construction of the star
  //-----------------------------------------------------------------------
    
  Et_magnetisation star(mp, nzet, relat, eos, (use_magnetisation!=0), 
			  (mag_in_eos!=0) ) ; 

    
  if ( star.is_relativistic() ) {
	  cout << "========================" << endl ;
	  cout << "Relativistic computation" << endl ;
	  cout << "========================" << endl ;
  }
  else {
	  cout << "=====================" << endl ;
	  cout << "Newtonian computation" << endl ;
	  cout << "=====================" << endl ;
  }
  
  //-----------------------------------------------------------------------
  //		Initialization of the enthalpy field
  //-----------------------------------------------------------------------


  const Coord& r = mp.r ;
  double ray0 = mp.val_r(nzet-1, 1., 0., 0.) ;  
  Cmp ent0(mp) ; 
  ent0 = ent_c * ( 1 - r*r / (ray0*ray0) ) ; 
  ent0.annule(nz-1) ; 
  ent0.std_base_scal() ; 
  star.set_enthalpy(ent0) ;  
    
  // Initialization of (n,e,p) from H
  star.equation_of_state() ; 

  // Initialization of (E,S,U,etc...) (quantities relative to the Eulerian obs)
  star.hydro_euler() ; 
  // Initialization of eulerian electromagnetic quantities.
  star.MHD_comput() ;

  /*
  cout << endl << "Initial star : " 
  << endl << "============   " << endl ;

  cout << star << endl ; 
  */

  //-----------------------------------------------------------------------
  //		Computation of the rotating equilibrium
  //-----------------------------------------------------------------------
  double omega = 2 * M_PI * freq_si / f_unit ; 
  double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 

  Itbl icontrol(11) ;
  icontrol.set_etat_qcq() ; 
  icontrol.set(0) = mer_max ; 
  icontrol.set(1) = mer_rot ; 
  icontrol.set(2) = mer_change_omega ; 
  icontrol.set(3) = mer_fix_omega ; 
  icontrol.set(4) = mer_mass ; 
  icontrol.set(5) = mermax_poisson ;  
  icontrol.set(6) = delta_mer_kep ; 
  icontrol.set(7) = mer_mag ;
  icontrol.set(8) = mer_change_mag ;
  icontrol.set(9)= mer_fix_mag ;
  icontrol.set(10) = mer_magmom ;
    
  Tbl control(8) ; 
  control.set_etat_qcq() ; 
  control.set(0) = precis ; 
  control.set(1) = omega_ini ; 
  control.set(2) = relax ; 
  control.set(3) = relax_poisson ; 
  control.set(4) = thres_adapt ; 
  control.set(5) = precis_adapt ; 
  control.set(6) = Q_ini ;
  control.set(7) = a_j_ini ;

  Tbl diff(1) ;     
  star.equilibrium_mag(ent_c, omega, fact_omega, nzadapt, ent_limit, 
	  icontrol, control, mbar_wanted, magmom_wanted,
	  aexp_mass, diff, Q0, a_j0, &f_j, &M_j) ;
  
  cout << endl << "Final star : " 
	<< endl << "==========   " << endl ;

  cout.precision(10) ; 
  cout << star << endl ; 
   
  const Base_vect_spher& bspher = mp.get_bvect_spher() ;
  Sym_tensor gij(mp, COV, bspher) ;
  gij.set_etat_zero() ;
  gij.set(1,1) = star.get_a_car()() ;
  gij.set(2,2) = star.get_a_car()() ;
  gij.set(3,3) = star.get_b_car()() ;
  Metric gam(gij) ;

  Scalar fac = sqrt(star.get_a_car()()) ;
  fac.std_spectral_base() ;
  Vector Bmag(mp, CON, bspher) ;
  Bmag.set(1) = Scalar(star.Magn()(0)) / fac ;
  Bmag.set(1).dec_dzpuis(2) ;
  Bmag.set(2) = Scalar(star.Magn()(1)) / fac ;
  Bmag.set(2).dec_dzpuis(2) ;
  Bmag.set(3) = 0 ;

  Scalar norm_B_Euler 
    = sqrt( Scalar(star.get_a_car()())*( Bmag(1)*Bmag(1) + Bmag(2)*Bmag(2) ) ) ;
  norm_B_Euler.std_spectral_base() ;
  Tbl maxB = 0.5*(max(abs(Bmag(1))) + max(abs(Bmag(2)))) ;
  Scalar divB = Bmag.divergence(gam) ;
  
  
  cout << "div(B) / max(B) in each domain :  " 
	<< max(abs(divB)) / maxB ; 
  cout << "Maximal magnetic field in the Eulerian frame: " 
	<< max(max(norm_B_Euler))*mag_unit << " [T]" << endl ;
  

  // Saveguard of the whole configuration
  // ------------------------------------
  /*
  FILE* fresu = fopen("resu.d", "w") ;
    
  star.get_mp().get_mg()->sauve(fresu) ;	// writing of the grid
  star.get_mp().sauve(fresu) ;                // writing of the mapping
  star.get_eos().sauve(fresu) ;  		// writing of the EOS
  star.sauve(fresu) ;                         // writing of the star
    
  fclose(fresu) ;
    */  
  out <<star.r_circ()*10 <<setw(15) << star.ray_eq ()*10<< setw(15);
  out<< star.ray_pole ()*10 << setw(15) << star.mass_g ()/msol << setw(15);
  out<< star.mass_b ()/msol << setw(15) << star.MagMom()/1.e32 << setw(15) ;
  out<< star.Magn()(0)(0,0,0,0)*mag_unit/1.e13 <<setw(15);
  out<<star.Magn()(0).va.val_point(star.l_surf()(0,0),star.xi_surf()(0,0),0.,0.)*mag_unit/1.e13;
  out<<setw(15)<<star.Magn()(1).va.val_point(star.l_surf()(0,0),star.xi_surf()(0,0),M_PI/2.,0.)*mag_unit/1.e13 <<setw(15)<<star.grv2()<<setw(15)<<star.grv3()<<setw(15) ;
  out<<star.get_nbar()()(0,0,0,0)<<setw(15)<<energy_domain(star,0)+energy_domain(star,1)<<setw(15)<<star.get_press()()(0,0,0,0)<<setw(15)<<star.get_ener()()(0,0,0,0)<<setw(15)<<star.mom_quad()<<endl;    
  ofstream k;
  k.open("rot_par"+to_string(num)+".txt",ios_base::app);
  k<<star.mass_b()/msol<<setw(15)<<star.Magn()(0)(0,0,0,0)*mag_unit/1.e13<<setw(15)<<star.get_omega_c()*f_unit/(2*M_PI)<<setw(15)<<star.angu_mom()<<endl;
  k.close();
  wy(fout,star,nt);
  // Drawings
  // --------
  /*
  if (graph == 1) {

	char title[80] ;

	// Cmp defining the surface of the star (via the enthalpy field)
	Cmp surf = star.get_ent()() ; 
	Cmp surf_ext(mp) ; 
	surf_ext = - 0.2 * surf(0, 0, 0, 0) ; 
	surf_ext.annule(0, star.get_nzet()-1) ; 
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ; 

	int nzdes = star.get_nzet() ; 

	des_coupe_y(star.get_At(), 0., nzdes, "A\\dt\\u Potential", &surf,
	  2.) ; 
  rename("pgplot.png","potential.png");
	des_coupe_y(star.get_Aphi(), 0., nzdes, "Magnetic field", &surf) ; 
  rename("pgplot.png","mag_field.png");
	des_coupe_y(norm_B_Euler, 0., nzdes, "Magnetic field norm", &surf) ; 
  rename("pgplot.png","mag_norm.png");
	des_coupe_y(star.get_ent()(), 0., nzdes, "Enthalpy", &surf) ; 
  rename("pgplot.png","enthalpy.png");

	strcpy(title, "Gravitational potential \\gn") ;   

	des_coupe_y(star.get_logn()(), 0., nzdes, title, &surf) ; 
	rename("pgplot.png","grav.png");
	strcpy(title, "Azimuthal shift N\\u\\gf") ; 
	des_coupe_y(star.get_nphi()(), 0., nzdes, title, &surf) ; 
	rename("pgplot.png","azi_shift.png");
	strcpy(title, "Metric potential \\gz") ; 
	des_coupe_y(star.get_dzeta()(), 0., nzdes, title, &surf) ; 
  rename("pgplot.png","mag_pot.png");
	
	strcpy(title, "Metric potential (NB-1) r sin\\gh") ; 
	des_coupe_y(star.get_tggg()(), 0., nzdes, title, &surf) ; 
	rename("pgplot.png","metric_potential.png");
	strcpy(title, "A\\u2\\dK\\uij\\dK\\dij\\u") ; 
	des_coupe_y(star.get_ak_car()(), 0., nzdes, title, &surf) ; 
  rename("pgplot.png","A.png");

  }*/
    // Cleaning
    // --------
  delete peos ;
  out.close();
  fout.close();
	//exit(EXIT_SUCCESS) ;
	return EXIT_SUCCESS;
}

void worker(int number,int thread_num,int max_it)
{
  sleep(number*1.5);
  for (int i=number;i<max_it;i+=thread_num)
  {
    try
    {
      core(i,number);
    }
    catch(const std::invalid_argument& e)
    {
      std::cerr << e.what() << '\n';
    }
    
  }
}

int main(int argc, char *argv[])
{
  int thread_num=atoi(argv[1]);
  int max_it=atoi(argv[2]);
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  worker(rank,thread_num,max_it);
  MPI_Finalize();
  return 0;
}

