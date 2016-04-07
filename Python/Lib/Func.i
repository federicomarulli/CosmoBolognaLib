%module cblFunc

%include "stl.i"
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include <std_shared_ptr.i>

%{
#include "Func.h"
#include "Data.h"
#include "Data1D.h"
#include "Data2D.h"
#include "Field3D.h"

  class EnumTypes {
  public:
    enum {_Omega_matter_LCDM_, _Omega_matter_, _Omega_baryon_, _Omega_neutrinos_, _massless_neutrinos_, _massive_neutrinos_, _Omega_DE_, _Omega_radiation_, _H0_, _scalar_amp_, _n_spec_, _w0_, _wa_, _fNL_, _sigma8_};
    enum {_1D_, _2D_}; 
    enum {_linear_, _logarithmic_};
    enum {_IdentityPrior_, _GaussianPrior_, _PoissonPrior_, _FunctionPrior_};
    enum {_1D_data_, _2D_data_, _1D_collection_data}; 
    enum {_angular_lin_, _angular_log_, _comoving_lin_, _comoving_log_, _comovingCartesian_linlin_, _comovingCartesian_linlog_, _comovingCartesian_loglin_, _comovingCartesian_loglog_, _comovingPolar_linlin_, _comovingPolar_linlog_, _comovingPolar_loglin_, _comovingPolar_loglog_};
    enum {_1D_monopole_, _1D_projected_, _1D_deprojected_, _1D_multipoles_, _1D_wedges_, _1D_filtered_, _1D_angular_, _2D_Cartesian_, _2D_polar_};
    enum {_Poisson_, _Jackknife_, _Bootstrap_};
    enum {_comoving_theta_, _comoving_side_};
    enum {_angular_connected_, _angular_reduced_, _comoving_connected_, _comoving_reduced_};
    enum {_GenericObject_, _RandomObject_, _Mock_, _Halo_, _Galaxy_, _Cluster_, _Void_}; 
    enum {_X_, _Y_, _Z_, _RA_, _Dec_, _Redshift_, _Dc_, _Weight_, _Mass_, _Richness_, _Magnitude_, _Vx_, Vy, Vz, _Region_, _Radius_, _Generic_};
    enum {_createRandom_box_, _createRandom_shuffle_, _createRandom_cone_, _createRandom_mock_, _createRandom_VIPERS_};
  };
  
  string cosmobl::par::DirCosmo = "~/CosmoBolognaLib/";
  string cosmobl::par::DirLoc = "./";
      
  static const double yotta = 1.e24;
  static const double zetta = 1.e21;
  static const double exa = 1.e18;
  static const double peta = 1.e15;
  static const double tera = 1.e12;
  static const double giga = 1.e9;
  static const double mega = 1.e6;
  static const double kilo = 1.e3;
  static const double ecto = 1.e2;
  static const double deca = 10.;
  static const double deci = 1.e-1;
  static const double centi = 1.e-2;
  static const double milli = 1.e-3;
  static const double micro = 1.e-6;
  static const double nano = 1.e-9;
  static const double pico = 1.e-12;
  static const double femto = 1.e-15;
  static const double atto = 1.e-18;
  static const double zepto = 1.e-21;
  static const double yocto = 1.e-24;
  static const char fINT[] = "%i";
  static const char fLONG[] = "%lli";
  static const char fDP0[] = "%1.0f"; 
  static const char fDP1[] = "%2.1f"; 
  static const char fDP2[] = "%3.2f"; 
  static const char fDP3[] = "%4.3f"; 
  static const char fDP4[] = "%5.4f"; 
  static const char fDP5[] = "%6.5f"; 
  static const char fDP6[] = "%7.6f"; 
  static const char ee3[] = "%4.3e";
  static const double pi = 3.1415926535897932;
  static const double ee = 2.7182818284590452;
  static const double hbar = 1.054571726e-3;
  static const double cc = 299792.458;      
  static const double kB = 1.3806488e23;
  static const double sSB = 5.670373e-8;
  static const double el = 1.602176565e-19;
  static const double alpha = 7.2973525698e-3;
  static const double epsilon0 = 8.854187817e-12;
  static const double mu0 = 12.566370614e-7; 
  static const double NAv = 6.02214129e23; 
  static const double GN = 6.6738480e-11;
  static const double gn = 9.80665;
  static const double lP = 1.616199e-35; 
  static const double MP = 2.17651e-8;
  static const double Msol = 1.98855e30;
  static const double me = 9.10938291e-31; 
  static const double mn = 1.674927351e-27; 
  static const double mp = 1.672621777e-27;
  static const double au = 149597870700; 
  static const double pc = 3.0856775814671916e16;
  static const double TCMB = 2.72548;
  static const double yr = 31557600;      
  static const string col_default = "\033[0m";
  static const string col_red = "\033[0;31m";
  static const string col_green = "\033[0;32m";
  static const string col_blue = "\033[0;34m";

  %}

%include "Func.h"
%include "Data.h"
%include "Data1D.h"
%include "Data2D.h"
%include "Field3D.h"

static const double yotta = 1.e24;
static const double zetta = 1.e21;
static const double exa = 1.e18;
static const double peta = 1.e15;
static const double tera = 1.e12;
static const double giga = 1.e9;
static const double mega = 1.e6;
static const double kilo = 1.e3;
static const double ecto = 1.e2;
static const double deca = 10.;
static const double deci = 1.e-1;
static const double centi = 1.e-2;
static const double milli = 1.e-3;
static const double micro = 1.e-6;
static const double nano = 1.e-9;
static const double pico = 1.e-12;
static const double femto = 1.e-15;
static const double atto = 1.e-18;
static const double zepto = 1.e-21;
static const double yocto = 1.e-24;
static const char fINT[] = "%i";
static const char fLONG[] = "%lli";
static const char fDP0[] = "%1.0f"; 
static const char fDP1[] = "%2.1f"; 
static const char fDP2[] = "%3.2f"; 
static const char fDP3[] = "%4.3f"; 
static const char fDP4[] = "%5.4f"; 
static const char fDP5[] = "%6.5f"; 
static const char fDP6[] = "%7.6f"; 
static const char ee3[] = "%4.3e";
static const double pi = 3.1415926535897932;
static const double ee = 2.7182818284590452;
static const double hbar = 1.054571726e-3;
static const double cc = 299792.458;      
static const double kB = 1.3806488e23;
static const double sSB = 5.670373e-8;
static const double el = 1.602176565e-19;
static const double alpha = 7.2973525698e-3;
static const double epsilon0 = 8.854187817e-12;
static const double mu0 = 12.566370614e-7; 
static const double NAv = 6.02214129e23; 
static const double GN = 6.6738480e-11;
static const double gn = 9.80665;
static const double lP = 1.616199e-35; 
static const double MP = 2.17651e-8;
static const double Msol = 1.98855e30;
static const double me = 9.10938291e-31; 
static const double mn = 1.674927351e-27; 
static const double mp = 1.672621777e-27;
static const double au = 149597870700; 
static const double pc = 3.0856775814671916e16;
static const double TCMB = 2.72548;
static const double yr = 31557600;      
static const string col_default = "\033[0m";
static const string col_red = "\033[0;31m";
static const string col_green = "\033[0;32m";
static const string col_blue = "\033[0;34m";

%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<std::string>;
%template(DoubleVectorVector) std::vector<std::vector<double> >;

class EnumTypes {
 public:
  enum {_Omega_matter_LCDM_, _Omega_matter_, _Omega_baryon_, _Omega_neutrinos_, _massless_neutrinos_, _massive_neutrinos_, _Omega_DE_, _Omega_radiation_, _H0_, _scalar_amp_, _n_spec_, _w0_, _wa_, _fNL_, _sigma8_};
  enum {_1D_, _2D_}; 
  enum {_linear_, _logarithmic_};
  enum {_1D_data_, _2D_data_, _1D_collection_data}; 
  enum {_IdentityPrior_, _GaussianPrior_, _PoissonPrior_, _FunctionPrior_};
  enum {_angular_lin_, _angular_log_, _comoving_lin_, _comoving_log_, _comovingCartesian_linlin_, _comovingCartesian_linlog_, _comovingCartesian_loglin_, _comovingCartesian_loglog_, _comovingPolar_linlin_, _comovingPolar_linlog_, _comovingPolar_loglin_, _comovingPolar_loglog_};
  enum {_1D_monopole_, _1D_projected_, _1D_deprojected_, _1D_multipoles_, _1D_wedges_, _1D_filtered_, _1D_angular_, _2D_Cartesian_, _2D_polar_};
  enum {_Poisson_, _Jackknife_, _Bootstrap_};
  enum {_comoving_theta_, _comoving_side_};
  enum {_angular_connected_, _angular_reduced_, _comoving_connected_, _comoving_reduced_};
  enum {_GenericObject_, _RandomObject_, _Mock_, _Halo_, _Galaxy_, _Cluster_, _Void_}; 
  enum {_X_, _Y_, _Z_, _RA_, _Dec_, _Redshift_, _Dc_, _Weight_, _Mass_, _Richness_, _Magnitude_, _Vx_, Vy, Vz, _Region_, _Radius_, _Generic_};
  enum {_createRandom_box_, _createRandom_shuffle_, _createRandom_cone_, _createRandom_mock_, _createRandom_VIPERS_};
};
  
