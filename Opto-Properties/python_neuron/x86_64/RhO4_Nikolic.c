/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__RhO4_Nikolic
#define _nrn_initial _nrn_initial__RhO4_Nikolic
#define nrn_cur _nrn_cur__RhO4_Nikolic
#define _nrn_current _nrn_current__RhO4_Nikolic
#define nrn_jacob _nrn_jacob__RhO4_Nikolic
#define nrn_state _nrn_state__RhO4_Nikolic
#define _net_receive _net_receive__RhO4_Nikolic 
#define kin kin__RhO4_Nikolic 
#define rates rates__RhO4_Nikolic 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define phiOn _p[0]
#define phiOn_columnindex 0
#define Dt_delay _p[1]
#define Dt_delay_columnindex 1
#define Dt_on _p[2]
#define Dt_on_columnindex 2
#define Dt_off _p[3]
#define Dt_off_columnindex 3
#define nPulses _p[4]
#define nPulses_columnindex 4
#define phi_m _p[5]
#define phi_m_columnindex 5
#define E _p[6]
#define E_columnindex 6
#define gam _p[7]
#define gam_columnindex 7
#define g0 _p[8]
#define g0_columnindex 8
#define v0 _p[9]
#define v0_columnindex 9
#define v1 _p[10]
#define v1_columnindex 10
#define k1 _p[11]
#define k1_columnindex 11
#define k2 _p[12]
#define k2_columnindex 12
#define kf _p[13]
#define kf_columnindex 13
#define kb _p[14]
#define kb_columnindex 14
#define Gf0 _p[15]
#define Gf0_columnindex 15
#define Gb0 _p[16]
#define Gb0_columnindex 16
#define p _p[17]
#define p_columnindex 17
#define q _p[18]
#define q_columnindex 18
#define Gd1 _p[19]
#define Gd1_columnindex 19
#define Gd2 _p[20]
#define Gd2_columnindex 20
#define Gr0 _p[21]
#define Gr0_columnindex 21
#define phi _p[22]
#define phi_columnindex 22
#define i _p[23]
#define i_columnindex 23
#define C1 _p[24]
#define C1_columnindex 24
#define O1 _p[25]
#define O1_columnindex 25
#define O2 _p[26]
#define O2_columnindex 26
#define C2 _p[27]
#define C2_columnindex 27
#define Ga1 _p[28]
#define Ga1_columnindex 28
#define Ga2 _p[29]
#define Ga2_columnindex 29
#define Gf _p[30]
#define Gf_columnindex 30
#define Gb _p[31]
#define Gb_columnindex 31
#define h1 _p[32]
#define h1_columnindex 32
#define h2 _p[33]
#define h2_columnindex 33
#define fphi _p[34]
#define fphi_columnindex 34
#define fv _p[35]
#define fv_columnindex 35
#define tally _p[36]
#define tally_columnindex 36
#define DC1 _p[37]
#define DC1_columnindex 37
#define DO1 _p[38]
#define DO1_columnindex 38
#define DO2 _p[39]
#define DO2_columnindex 39
#define DC2 _p[40]
#define DC2_columnindex 40
#define v _p[41]
#define v_columnindex 41
#define _g _p[42]
#define _g_columnindex 42
#define _tsav _p[43]
#define _tsav_columnindex 43
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_rates(void*);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "rates", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "Dt_off", 0, 1e+09,
 "Dt_on", 0, 1e+09,
 "Dt_delay", 0, 1e+09,
 "nPulses", 0, 1000,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Dt_delay", "ms",
 "Dt_on", "ms",
 "Dt_off", "ms",
 "nPulses", "1",
 "E", "mV",
 "gam", "1",
 "g0", "pS",
 "v0", "mV",
 "v1", "mV",
 "k1", "/ms",
 "k2", "/ms",
 "kf", "/ms",
 "kb", "/ms",
 "Gf0", "/ms",
 "Gb0", "/ms",
 "p", "1",
 "q", "1",
 "Gd1", "/ms",
 "Gd2", "/ms",
 "Gr0", "/ms",
 "i", "nA",
 0,0
};
 static double C20 = 0;
 static double C10 = 0;
 static double O20 = 0;
 static double O10 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"RhO4_Nikolic",
 "phiOn",
 "Dt_delay",
 "Dt_on",
 "Dt_off",
 "nPulses",
 "phi_m",
 "E",
 "gam",
 "g0",
 "v0",
 "v1",
 "k1",
 "k2",
 "kf",
 "kb",
 "Gf0",
 "Gb0",
 "p",
 "q",
 "Gd1",
 "Gd2",
 "Gr0",
 0,
 "phi",
 "i",
 0,
 "C1",
 "O1",
 "O2",
 "C2",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 44, _prop);
 	/*initialize range parameters*/
 	phiOn = 1e+18;
 	Dt_delay = 25;
 	Dt_on = 100;
 	Dt_off = 50;
 	nPulses = 1;
 	phi_m = 1e+16;
 	E = 0;
 	gam = 0.05;
 	g0 = 75000;
 	v0 = 43;
 	v1 = 17;
 	k1 = 0.05;
 	k2 = 0.015;
 	kf = 0.03;
 	kb = 0.0115;
 	Gf0 = 0.01;
 	Gb0 = 0.015;
 	p = 0.7;
 	q = 0.47;
 	Gd1 = 0.11;
 	Gd2 = 0.025;
 	Gr0 = 0.0004;
  }
 	_prop->param = _p;
 	_prop->param_size = 44;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _RhO4_Nikolic_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 3,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 44, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 RhO4_Nikolic /home/gjgpb9/OpsinStudy/Opto-Properties/python_neuron/modfiles/RhO4_Nikolic.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Nikolic 4-state rhodopsin model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 extern double *_nrn_thread_getelm(SparseObj*, int, int);
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[4], _dlist1[4]; static double *_temp1;
 static int kin();
 
static int kin (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<4;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ phi ) ;
   /* ~ C1 <-> O1 ( Ga1 , Gd1 )*/
 f_flux =  Ga1 * C1 ;
 b_flux =  Gd1 * O1 ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  Ga1 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 3 ,1)  -= _term;
 _term =  Gd1 ;
 _MATELM1( 1 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ O1 <-> O2 ( Gf , Gb )*/
 f_flux =  Gf * O1 ;
 b_flux =  Gb * O2 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  Gf ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  Gb ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ O2 <-> C2 ( Gd2 , Ga2 )*/
 f_flux =  Gd2 * O2 ;
 b_flux =  Ga2 * C2 ;
 _RHS1( 2) -= (f_flux - b_flux);
 
 _term =  Gd2 ;
 _MATELM1( 2 ,2)  += _term;
 _term =  Ga2 ;
 _MATELM1( 2 ,0)  -= _term;
 /*REACTION*/
  /* ~ C2 <-> C1 ( Gr0 , 0.0 )*/
 f_flux =  Gr0 * C2 ;
 b_flux =  0.0 * C1 ;
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  Gr0 ;
 _MATELM1( 1 ,0)  -= _term;
 _term =  0.0 ;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
   /* C1 + O1 + O2 + C2 = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= C2 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= O2 ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= O1 ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= C1 ;
 /*CONSERVATION*/
   } return _reset;
 }
 
static int  rates ( _threadargsprotocomma_ double _lphi ) {
   if ( _lphi > 0.0 ) {
     h1 = 1.0 / ( 1.0 + pow ( phi_m , p ) / pow ( _lphi , p ) ) ;
     h2 = 1.0 / ( 1.0 + pow ( phi_m , q ) / pow ( _lphi , q ) ) ;
     }
   else {
     h1 = 0.0 ;
     h2 = 0.0 ;
     }
   Ga1 = k1 * h1 ;
   Ga2 = k2 * h1 ;
   Gf = Gf0 + kf * h2 ;
   Gb = Gb0 + kb * h2 ;
    return 0; }
 
static double _hoc_rates(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 1.0 ) {
     phi = phiOn ;
     net_send ( _tqitem, _args, _pnt, t +  Dt_on , 0.0 ) ;
     }
   else {
     phi = 0.0 ;
     if ( tally > 0.0 ) {
       net_send ( _tqitem, _args, _pnt, t +  Dt_off , 1.0 ) ;
       tally = tally - 1.0 ;
       }
     }
   } }
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<4;_i++) _p[_dlist1[_i]] = 0.0;}
 rates ( _threadargscomma_ phi ) ;
 /* ~ C1 <-> O1 ( Ga1 , Gd1 )*/
 f_flux =  Ga1 * C1 ;
 b_flux =  Gd1 * O1 ;
 DC1 -= (f_flux - b_flux);
 DO1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ O1 <-> O2 ( Gf , Gb )*/
 f_flux =  Gf * O1 ;
 b_flux =  Gb * O2 ;
 DO1 -= (f_flux - b_flux);
 DO2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ O2 <-> C2 ( Gd2 , Ga2 )*/
 f_flux =  Gd2 * O2 ;
 b_flux =  Ga2 * C2 ;
 DO2 -= (f_flux - b_flux);
 DC2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> C1 ( Gr0 , 0.0 )*/
 f_flux =  Gr0 * C2 ;
 b_flux =  0.0 * C1 ;
 DC2 -= (f_flux - b_flux);
 DC1 += (f_flux - b_flux);
 
 /*REACTION*/
   /* C1 + O1 + O2 + C2 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<4;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ phi ) ;
 /* ~ C1 <-> O1 ( Ga1 , Gd1 )*/
 _term =  Ga1 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 3 ,1)  -= _term;
 _term =  Gd1 ;
 _MATELM1( 1 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ O1 <-> O2 ( Gf , Gb )*/
 _term =  Gf ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  Gb ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ O2 <-> C2 ( Gd2 , Ga2 )*/
 _term =  Gd2 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 0 ,2)  -= _term;
 _term =  Ga2 ;
 _MATELM1( 2 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ C2 <-> C1 ( Gr0 , 0.0 )*/
 _term =  Gr0 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 1 ,0)  -= _term;
 _term =  0.0 ;
 _MATELM1( 0 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
   /* C1 + O1 + O2 + C2 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 4, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  C2 = C20;
  C1 = C10;
  O2 = O20;
  O1 = O10;
 {
   C1 = 1.0 ;
   O1 = 0.0 ;
   O2 = 0.0 ;
   C2 = 0.0 ;
   phi = 0.0 ;
   rates ( _threadargscomma_ phi ) ;
   i = 0.0 ;
   tally = nPulses ;
   if ( tally > 0.0 ) {
     net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  Dt_delay , 1.0 ) ;
     tally = tally - 1.0 ;
     }
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   fphi = O1 + gam * O2 ;
   fv = ( 1.0 - exp ( - ( v - E ) / v0 ) ) / ( ( v - E ) / v1 ) ;
   i = g0 * fphi * fv * ( v - E ) * ( 1e-6 ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {  sparse_thread(&_thread[_spth1]._pvoid, 4, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = C2_columnindex;  _dlist1[0] = DC2_columnindex;
 _slist1[1] = C1_columnindex;  _dlist1[1] = DC1_columnindex;
 _slist1[2] = O2_columnindex;  _dlist1[2] = DO2_columnindex;
 _slist1[3] = O1_columnindex;  _dlist1[3] = DO1_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/gjgpb9/OpsinStudy/Opto-Properties/python_neuron/modfiles/RhO4_Nikolic.mod";
static const char* nmodl_file_text = 
  "TITLE Nikolic 4-state rhodopsin model\n"
  "\n"
  "NEURON {\n"
  "    POINT_PROCESS RhO4_Nikolic\n"
  "    NONSPECIFIC_CURRENT i\n"
  "    RANGE i, E, gam, v0, v1, g0 :, fphi, fv\n"
  "    RANGE k1, k2, kf, kb, Gf0, Gb0, Gd1, Gd2, Gr0, p, q\n"
  "    RANGE phiOn, phi, phi_m, Dt_delay, Dt_on, Dt_off, nPulses\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (nA) = (nanoamp)\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)\n"
  "    (pS) = (picosiemens)\n"
  "    :(photons)  =(1)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER { : Initialise parameters to defaults. These may be changed through hoc files\n"
  "\n"
  ": Illumination\n"
  "    phiOn   = 1e18 :ph/s/mm^2 :(mW/mm2)                                      _______\n"
  "    Dt_delay    = 25    (ms)    <0, 1e9>        : delay before ON phase         |  ON   |  OFF\n"
  "    Dt_on     = 100   (ms)    <0, 1e9>        : duration of ON phase  <-Dt_delay->|<-Dt_on->|<-Dt_off->\n"
  "    Dt_off    = 50    (ms)    <0, 1e9>        : duration of OFF phase ________|       |________\n"
  "    nPulses = 1     (1)     <0, 1e3>        : num pulses to deliver         <-- one pulse -->\n"
  "\n"
  ": Illumination constants   \n"
  "    phi_m    = 1e16  : (photons/sec mm2)     : Hill Constant\n"
  ":   lambda  = 470   : (nm)\n"
  "\n"
  ": Conductance\n"
  "    E       = 0         (mV)                : Channel reversal potential\n"
  "    gam     = 0.05      (1)                 : Ratio of open-state conductances\n"
  "    g0      = 75000     (pS)                : gbar : defined in the main prog as HR_expression*Area\n"
  "\n"
  ": Inward rectifier conductance    \n"
  "    v0      = 43        (mV)\n"
  "    v1      = 17       (mV)    :4.1 as original\n"
  "\n"
  ": State transition rate parameters (/ms)\n"
  "    k1      = 0.05      (/ms)\n"
  "    k2      = 0.015     (/ms)\n"
  "    kf      = 0.03      (/ms)\n"
  "    kb      = 0.0115    (/ms)\n"
  "    Gf0     = 0.01      (/ms)\n"
  "    Gb0     = 0.015     (/ms)\n"
  "    p       = 0.7       (1)                 : Hill Coefficient Ga{1,2}\n"
  "    q       = 0.47      (1)                 : Hill Coefficient G{f,b}\n"
  "    Gd1     = 0.11      (/ms)\n"
  "    Gd2     = 0.025     (/ms)\n"
  "    Gr0     = 0.0004    (/ms)\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "    phi     :(photons/s mm2)                : Instantaneous flux\n"
  "\n"
  "    Ga1     (/ms)  : C1 -> O1\n"
  "    Ga2     (/ms)  : C2 -> O2\n"
  "    Gf      (/ms)  : O1 -> O2\n"
  "    Gb      (/ms)  : O2 -> O1\n"
  "    h1      (1)\n"
  "    h2      (1)\n"
  "\n"
  "    fphi    (1)                             : Fractional conductance\n"
  "    fv      (1)                             : Fractional conductance\n"
  "    :g      (1) \n"
  "    v       (mV)\n"
  "    tally                                   : Number of pulses left to deliver\n"
  "    i       (nA)\n"
  "}\n"
  "\n"
  "\n"
  "STATE { C1 O1 O2 C2 }\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE kin METHOD sparse\n"
  "    fphi    = O1+gam*O2                     : Light dependency op=[0:1] \n"
  "    fv      = (1-exp(-(v-E)/v0))/((v-E)/v1) : Voltage dependency (equal 1 at Vm=-70mV)\n"
  "    i       = g0*fphi*fv*(v-E)*(1e-6)       : Photocurrent\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL {   : Initialise variables\n"
  ": State occupancy proportions - starts in C1\n"
  "    C1      = 1     : Closed state 1 (Ground state)\n"
  "    O1      = 0     : Open state 1 (High conductivity)\n"
  "    O2      = 0     : Open state 2 (Low conductivity)\n"
  "    C2      = 0     : Closed state 2\n"
  "\n"
  "    phi     = 0\n"
  "    rates(phi) : necessary?\n"
  "    i       = 0\n"
  "    \n"
  "    tally   = nPulses : Set the tally to the number of pulses required \n"
  "    if (tally > 0) {\n"
  "        net_send(Dt_delay, 1)\n"
  "        tally = tally - 1\n"
  "    }\n"
  "}\n"
  "\n"
  "\n"
  "KINETIC kin {\n"
  "    rates(phi)\n"
  "    ~ C1 <-> O1 (Ga1, Gd1)\n"
  "    ~ O1 <-> O2 (Gf,  Gb)\n"
  "    ~ O2 <-> C2 (Gd2, Ga2)\n"
  "    ~ C2 <-> C1 (Gr0, 0)\n"
  "    :O1=0\n"
  "    :O2=0\n"
  "    CONSERVE C1 + O1 + O2 + C2 = 1\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(phi) { : Define equations for calculating transition rates\n"
  "\n"
  "    :p1 = pow(phi,p)\n"
  "    :p2 = pow(phi,q)\n"
  "    :d1 = pow(phi_m,p)\n"
  "    :d2 = pow(phi_m,q)\n"
  "    \n"
  "    if (phi > 0) {\n"
  "        h1 = 1/(1 + pow(phi_m,p)/pow(phi,p)) : p1/(p1 + d1)\n"
  "        h2 = 1/(1 + pow(phi_m,q)/pow(phi,q)) : p2/(p2 + d2)\n"
  "    } else {\n"
  "        h1 = 0\n"
  "        h2 = 0\n"
  "    }\n"
  "    \n"
  "        :Ga1 = k1*(phi/phi0)\n"
  "    Ga1 = k1 * h1             :Ga1 = k1 * pow(phi,p)/(pow(phi,p) + pow(phi_m,p))\n"
  "        :Ga2 = k2*(phi/phi0)\n"
  "    Ga2 = k2 * h1             :Ga2 = k2 * pow(phi,p)/(pow(phi,p) + pow(phi_m,p))\n"
  "        :Gf = Gf0 + kf*log(1+(phi/phi0))\n"
  "    Gf = Gf0 + kf * h2      :Gf = Gf0 + kf * pow(phi,q)/(pow(phi,q) + pow(phi_m,q))\n"
  "        :Gb = Gb0 + kb*log(1+(phi/phi0))\n"
  "    Gb = Gb0 + kb * h2      :Gb = Gb0 + kb * pow(phi,q)/(pow(phi,q) + pow(phi_m,q))\n"
  "    \n"
  "}\n"
  "\n"
  "\n"
  ": Add functions for calculating transition rates outside of simulations\n"
  ":FUNCTION A1(phi) {\n"
  ":    A1 = a10*(phi/phi0)\n"
  ":}\n"
  ": ...\n"
  "\n"
  "\n"
  "NET_RECEIVE (w) {\n"
  "    if (flag == 1) { : ignore any but self-events with flag == 1. This may not be necessary... see e.g. nrn/src/nrnoc/intfire1.mod\n"
  "        phi = phiOn\n"
  "        net_send(Dt_on, 0) : Schedule the next off phase\n"
  "    } else { : Turn the light off\n"
  "        phi = 0\n"
  "        if (tally > 0) { : Schedule the next on phase\n"
  "            net_send(Dt_off, 1)\n"
  "            tally = tally - 1\n"
  "        }\n"
  "    }\n"
  "}\n"
  ;
#endif
