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
 
#define nrn_init _nrn_init__RhO3
#define _nrn_initial _nrn_initial__RhO3
#define nrn_cur _nrn_cur__RhO3
#define _nrn_current _nrn_current__RhO3
#define nrn_jacob _nrn_jacob__RhO3
#define nrn_state _nrn_state__RhO3
#define _net_receive _net_receive__RhO3 
#define kin kin__RhO3 
#define rates rates__RhO3 
 
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
#define g0 _p[7]
#define g0_columnindex 7
#define v0 _p[8]
#define v0_columnindex 8
#define v1 _p[9]
#define v1_columnindex 9
#define k_a _p[10]
#define k_a_columnindex 10
#define p _p[11]
#define p_columnindex 11
#define k_r _p[12]
#define k_r_columnindex 12
#define q _p[13]
#define q_columnindex 13
#define Gd _p[14]
#define Gd_columnindex 14
#define Gr0 _p[15]
#define Gr0_columnindex 15
#define phi _p[16]
#define phi_columnindex 16
#define i _p[17]
#define i_columnindex 17
#define C _p[18]
#define C_columnindex 18
#define O _p[19]
#define O_columnindex 19
#define D _p[20]
#define D_columnindex 20
#define Ga _p[21]
#define Ga_columnindex 21
#define Gr _p[22]
#define Gr_columnindex 22
#define fphi _p[23]
#define fphi_columnindex 23
#define fv _p[24]
#define fv_columnindex 24
#define tally _p[25]
#define tally_columnindex 25
#define DC _p[26]
#define DC_columnindex 26
#define DO _p[27]
#define DO_columnindex 27
#define DD _p[28]
#define DD_columnindex 28
#define v _p[29]
#define v_columnindex 29
#define _g _p[30]
#define _g_columnindex 30
#define _tsav _p[31]
#define _tsav_columnindex 31
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
 "Gr0", 0, 1e+09,
 "Gd", 0, 1e+09,
 "nPulses", 0, 1000,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Dt_delay", "ms",
 "Dt_on", "ms",
 "Dt_off", "ms",
 "nPulses", "1",
 "E", "mV",
 "g0", "pS",
 "v0", "mV",
 "v1", "mV",
 "k_a", "/ms",
 "p", "1",
 "k_r", "/ms",
 "q", "1",
 "Gd", "/ms",
 "Gr0", "/ms",
 "i", "nA",
 0,0
};
 static double C0 = 0;
 static double D0 = 0;
 static double O0 = 0;
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
"RhO3",
 "phiOn",
 "Dt_delay",
 "Dt_on",
 "Dt_off",
 "nPulses",
 "phi_m",
 "E",
 "g0",
 "v0",
 "v1",
 "k_a",
 "p",
 "k_r",
 "q",
 "Gd",
 "Gr0",
 0,
 "phi",
 "i",
 0,
 "C",
 "O",
 "D",
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
 	_p = nrn_prop_data_alloc(_mechtype, 32, _prop);
 	/*initialize range parameters*/
 	phiOn = 1.0002e+16;
 	Dt_delay = 25;
 	Dt_on = 100;
 	Dt_off = 50;
 	nPulses = 1;
 	phi_m = 7.7e+17;
 	E = 0;
 	g0 = 27700;
 	v0 = 43;
 	v1 = 17.1;
 	k_a = 93.25;
 	p = 1;
 	k_r = 0.01;
 	q = 1;
 	Gd = 0.278;
 	Gr0 = 2e-05;
  }
 	_prop->param = _p;
 	_prop->param_size = 32;
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

 void _RhO3_reg() {
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
  hoc_register_prop_size(_mechtype, 32, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 RhO3 /home/gjgpb9/OpsinStudy/Opto-Properties/python_neuron/modfiles/RhO3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Nikolic 3-state rhodopsin model";

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
 static int _slist1[3], _dlist1[3]; static double *_temp1;
 static int kin();
 
static int kin (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<3;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ phi ) ;
   /* ~ C <-> O ( Ga , 0.0 )*/
 f_flux =  Ga * C ;
 b_flux =  0.0 * O ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  Ga ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 2 ,1)  -= _term;
 _term =  0.0 ;
 _MATELM1( 1 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ O <-> D ( Gd , 0.0 )*/
 f_flux =  Gd * O ;
 b_flux =  0.0 * D ;
 _RHS1( 2) -= (f_flux - b_flux);
 
 _term =  Gd ;
 _MATELM1( 2 ,2)  += _term;
 _term =  0.0 ;
 _MATELM1( 2 ,0)  -= _term;
 /*REACTION*/
  /* ~ D <-> C ( Gr , 0.0 )*/
 f_flux =  Gr * D ;
 b_flux =  0.0 * C ;
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  Gr ;
 _MATELM1( 1 ,0)  -= _term;
 _term =  0.0 ;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
   /* C + O + D = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= D ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= O ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= C ;
 /*CONSERVATION*/
   } return _reset;
 }
 
static int  rates ( _threadargsprotocomma_ double _lphi ) {
   if ( _lphi > 0.0 ) {
     Ga = k_a * 1.0 / ( 1.0 + pow ( phi_m , p ) / pow ( _lphi , p ) ) ;
     Gr = Gr0 + k_r * 1.0 / ( 1.0 + pow ( phi_m , q ) / pow ( _lphi , q ) ) ;
     }
   else {
     Ga = 0.0 ;
     Gr = Gr0 ;
     }
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
 {int _i; for(_i=0;_i<3;_i++) _p[_dlist1[_i]] = 0.0;}
 rates ( _threadargscomma_ phi ) ;
 /* ~ C <-> O ( Ga , 0.0 )*/
 f_flux =  Ga * C ;
 b_flux =  0.0 * O ;
 DC -= (f_flux - b_flux);
 DO += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ O <-> D ( Gd , 0.0 )*/
 f_flux =  Gd * O ;
 b_flux =  0.0 * D ;
 DO -= (f_flux - b_flux);
 DD += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ D <-> C ( Gr , 0.0 )*/
 f_flux =  Gr * D ;
 b_flux =  0.0 * C ;
 DD -= (f_flux - b_flux);
 DC += (f_flux - b_flux);
 
 /*REACTION*/
   /* C + O + D = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<3;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ phi ) ;
 /* ~ C <-> O ( Ga , 0.0 )*/
 _term =  Ga ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 2 ,1)  -= _term;
 _term =  0.0 ;
 _MATELM1( 1 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ O <-> D ( Gd , 0.0 )*/
 _term =  Gd ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 0 ,2)  -= _term;
 _term =  0.0 ;
 _MATELM1( 2 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ D <-> C ( Gr , 0.0 )*/
 _term =  Gr ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 1 ,0)  -= _term;
 _term =  0.0 ;
 _MATELM1( 0 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
   /* C + O + D = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 3;}
 
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
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 3, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
  C = C0;
  D = D0;
  O = O0;
 {
   C = 1.0 ;
   O = 0.0 ;
   D = 0.0 ;
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
   fphi = O ;
   fv = ( 1.0 - exp ( - ( v - E ) / v0 ) ) * v1 / ( v - E ) ;
   i = 0.001 * g0 * fphi * fv * ( v - E ) ;
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
 {  sparse_thread(&_thread[_spth1]._pvoid, 3, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
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
 _slist1[0] = D_columnindex;  _dlist1[0] = DD_columnindex;
 _slist1[1] = C_columnindex;  _dlist1[1] = DC_columnindex;
 _slist1[2] = O_columnindex;  _dlist1[2] = DO_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/gjgpb9/OpsinStudy/Opto-Properties/python_neuron/modfiles/RhO3.mod";
static const char* nmodl_file_text = 
  "TITLE Nikolic 3-state rhodopsin model\n"
  "\n"
  "NEURON  {\n"
  "    POINT_PROCESS RhO3\n"
  "    NONSPECIFIC_CURRENT i\n"
  "    RANGE i, E, v0, v1, g0 :, fphi, fv\n"
  "    RANGE k_a, k_r, Gd, Gr0, p, q\n"
  "    RANGE phiOn, phi, phi_m, Dt_delay, Dt_on, Dt_off, nPulses\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (nA) = (nanoamp)\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)\n"
  "    (pS) = (picosiemens)\n"
  "    :(photons) = (1)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER { : Initialise parameters to defaults. These may be changed through hoc files\n"
  "\n"
  ": Illumination\n"
  "    phiOn   = 1.0002e16      :(photons/sec mm2)  : Flux [Irradiance :(mW/mm^2)]    : 1e18 as original  _______\n"
  "    Dt_delay    = 25        (ms)    <0, 1e9>    : delay before ON phase             |  ON   |  OFF\n"
  "    Dt_on     = 100       (ms)    <0, 1e9>    : duration of each ON phase <-Dt_delay->|<-Dt_on->|<-Dt_off->\n"
  "    Dt_off    = 50        (ms)    <0, 1e9>    : duration of each OFF phase________|       |________\n"
  "    nPulses = 1         (1)     <0, 1e3>    : num pulses to deliver             <-- one pulse -->\n"
  "\n"
  ": Illumination constants   \n"
  ":   lambda  = 470       :(nm)\n"
  "    phi_m   = 7.7e17      :(photons/sec mm2)  : Hill Constant  1e16 as original\n"
  "  \n"
  ": Conductance    \n"
  "    E       = 0         (mV)                : Channel reversal potential\n"
  "    g0      = 27700      (pS)                : defined in the main prog as HR_expression*Area\n"
  "\n"
  ": Inward rectifier conductance              fv = (1-exp(-(v-E)/v0))*v1/(v-E) := 1 at Vm=-70mV\n"
  "    v0      = 43        (mV)\n"
  "    v1      = 17.1      (mV)  : 4.1 is original one\n"
  "    \n"
  ": State transition rate parameters (/ms)\n"
  "    k_a     = 93.25      (/ms)               : Quantum efficiency * number of photons absorbed by a RhO molecule  per unit time 0.5*1.2e-14 /1.1   0.28 as original\n"
  "    p       = 1       (1)                 : Hill Coefficient   :0.4 as original\n"
  "    k_r     = 0.01      (/ms)    :0.28 as original\n"
  "    q       = 1       (1)                 : Hill Coefficient   :0.4 as original\n"
  "    Gd      = 0.278    (/ms) <0, 1e9>   : 0.0909 original : @ 1mW mm^-2\n"
  "    Gr0     = 0.00002    (/ms) <0, 1e9>  :0.0002 original   : tau_r,dark = 5-10s p405 Nikolic et al. 2009\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "    phi             : \"flux\" should now be \"intensity\". Er is normalised to be dimensionless\n"
  "    Ga      (/ms)\n"
  "    Gr      (/ms)\n"
  "\n"
  "    fphi    (1)\n"
  "    fv      (1)\n"
  "    :g      (pS)\n"
  "    v       (mV) \n"
  "    i       (nA)\n"
  "    tally           : how many more pulses to deliver\n"
  "}\n"
  "\n"
  "\n"
  "STATE { C O D }\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE kin METHOD sparse\n"
  "    fphi    = O                             : Light dependency on open state O=[0:1] \n"
  "    fv      = (1-exp(-(v-E)/v0))*v1/(v-E)   : Voltage dependency (equal 1 at Vm=-70mV)\n"
  "    i       = 0.001*g0*fphi*fv*(v-E):*(1e-6)       : Photocurrent\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL {   : Initialise variables to fully dark-adapted state\n"
  ": State occupancy proportions - starts in C\n"
  "    C       = 1: 1\n"
  "    O       = 0: 0\n"
  "    D       = 0: 0\n"
  "    \n"
  "    phi     = 0\n"
  "    rates(phi) : necessary?\n"
  "    i       = 0\n"
  "    \n"
  "    tally   = nPulses\n"
  "    if (tally > 0) {\n"
  "        net_send(Dt_delay, 1)   : Schedule an on event at t+Dt_delay\n"
  "        tally = tally - 1\n"
  "    }\n"
  "}\n"
  "\n"
  ":DERIVATIVE deriv {\n"
  ":   rates(flux)\n"
  ":   C' = (Gr * D) - (Ga * C)\n"
  ":   O' = (Ga * C) - (Gd * D)\n"
  ":   D' = (Gd * O) - (Gr * D)\n"
  ":}\n"
  "KINETIC kin {\n"
  "    rates(phi)\n"
  "    ~ C <-> O (Ga, 0)\n"
  "    ~ O <-> D (Gd, 0)\n"
  "    ~ D <-> C (Gr, 0)\n"
  "    CONSERVE C + O + D = 1\n"
  "}\n"
  "PROCEDURE rates(phi) {          : Define equations for calculating transition rates\n"
  "    if (phi>0) {\n"
  "        :Ga = phi * k           : k = quantum efficiency\n"
  "        Ga = k_a * 1/(1+pow(phi_m,p)/pow(phi,p))    : pow(phi,p)/(pow(phi,p) + pow(phi_m,p))\n"
  "        Gr = Gr0 + k_r * 1/(1+pow(phi_m,q)/pow(phi,q))    : pow(phi,q)/(pow(phi,q) + pow(phi_m,q))\n"
  "    } else {\n"
  "        Ga = 0\n"
  "        Gr = Gr0\n"
  "    }\n"
  "}  \n"
  ": NET_RECEIVE (next (ms), phiNext) {\n"
  ":   if (flag == 1) {            : Switch light on\n"
  ":        phi = phiNext\n"
  ":        net_send(next, 0)       : Schedule an off event at t+Dt_on\n"
  ":    } else {                    : Switch light off\n"
  ":        phi = 0\n"
  ":        if (tally > 0) {        : More pulses to deliver\n"
  ":            net_send(next, 1)   : Schedule an on event at t+Dt_off\n"
  ":            tally = tally - 1\n"
  ":        }\n"
  ":    }\n"
  ":}\n"
  "NET_RECEIVE (w) {\n"
  "    if (flag == 1) {            : Switch light on\n"
  "        phi = phiOn\n"
  "        net_send(Dt_on, 0)        : Schedule an off event at t+Dt_on\n"
  "    } else {                    : Switch light off\n"
  "        phi = 0\n"
  "        if (tally > 0) {        : More pulses to deliver\n"
  "            net_send(Dt_off, 1)   : Schedule an on event at t+Dt_off\n"
  "            tally = tally - 1\n"
  "        }\n"
  "    }\n"
  "}\n"
  ;
#endif
