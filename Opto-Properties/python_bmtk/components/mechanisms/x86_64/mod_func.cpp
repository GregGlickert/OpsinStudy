#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _AMPA_NMDA_STP_reg(void);
extern void _cadad_reg(void);
extern void _cal2_reg(void);
extern void _can_mig_reg(void);
extern void _ChR2H134R_william_event_reg(void);
extern void _ChR2H134R_william_reg(void);
extern void _ChR2H134R_william_original_reg(void);
extern void _GABA_AB_STP_reg(void);
extern void _GABA_A_STP_reg(void);
extern void _gap_reg(void);
extern void _Gfluct_reg(void);
extern void _h_kole_reg(void);
extern void _imCA3_reg(void);
extern void _kap_BS_reg(void);
extern void _kBK_reg(void);
extern void _kdmc_BS_reg(void);
extern void _kdr_BS_reg(void);
extern void _kdrCA3_reg(void);
extern void _kdrinter_reg(void);
extern void _leak_reg(void);
extern void _nainter_reg(void);
extern void _napCA3_reg(void);
extern void _natCA3_reg(void);
extern void _nax_BS_reg(void);
extern void _RhO3_reg(void);
extern void _RhO4_Nikolic_reg(void);
extern void _vecevent_coreneuron_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"modfiles//AMPA_NMDA_STP.mod\"");
    fprintf(stderr, " \"modfiles//cadad.mod\"");
    fprintf(stderr, " \"modfiles//cal2.mod\"");
    fprintf(stderr, " \"modfiles//can_mig.mod\"");
    fprintf(stderr, " \"modfiles//ChR2H134R_william_event.mod\"");
    fprintf(stderr, " \"modfiles//ChR2H134R_william.mod\"");
    fprintf(stderr, " \"modfiles//ChR2H134R_william_original.mod\"");
    fprintf(stderr, " \"modfiles//GABA_AB_STP.mod\"");
    fprintf(stderr, " \"modfiles//GABA_A_STP.mod\"");
    fprintf(stderr, " \"modfiles//gap.mod\"");
    fprintf(stderr, " \"modfiles//Gfluct.mod\"");
    fprintf(stderr, " \"modfiles//h_kole.mod\"");
    fprintf(stderr, " \"modfiles//imCA3.mod\"");
    fprintf(stderr, " \"modfiles//kap_BS.mod\"");
    fprintf(stderr, " \"modfiles//kBK.mod\"");
    fprintf(stderr, " \"modfiles//kdmc_BS.mod\"");
    fprintf(stderr, " \"modfiles//kdr_BS.mod\"");
    fprintf(stderr, " \"modfiles//kdrCA3.mod\"");
    fprintf(stderr, " \"modfiles//kdrinter.mod\"");
    fprintf(stderr, " \"modfiles//leak.mod\"");
    fprintf(stderr, " \"modfiles//nainter.mod\"");
    fprintf(stderr, " \"modfiles//napCA3.mod\"");
    fprintf(stderr, " \"modfiles//natCA3.mod\"");
    fprintf(stderr, " \"modfiles//nax_BS.mod\"");
    fprintf(stderr, " \"modfiles//RhO3.mod\"");
    fprintf(stderr, " \"modfiles//RhO4_Nikolic.mod\"");
    fprintf(stderr, " \"modfiles//vecevent_coreneuron.mod\"");
    fprintf(stderr, "\n");
  }
  _AMPA_NMDA_STP_reg();
  _cadad_reg();
  _cal2_reg();
  _can_mig_reg();
  _ChR2H134R_william_event_reg();
  _ChR2H134R_william_reg();
  _ChR2H134R_william_original_reg();
  _GABA_AB_STP_reg();
  _GABA_A_STP_reg();
  _gap_reg();
  _Gfluct_reg();
  _h_kole_reg();
  _imCA3_reg();
  _kap_BS_reg();
  _kBK_reg();
  _kdmc_BS_reg();
  _kdr_BS_reg();
  _kdrCA3_reg();
  _kdrinter_reg();
  _leak_reg();
  _nainter_reg();
  _napCA3_reg();
  _natCA3_reg();
  _nax_BS_reg();
  _RhO3_reg();
  _RhO4_Nikolic_reg();
  _vecevent_coreneuron_reg();
}

#if defined(__cplusplus)
}
#endif
