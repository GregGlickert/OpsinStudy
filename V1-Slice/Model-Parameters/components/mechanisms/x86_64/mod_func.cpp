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
extern void _exp2syn_stp_reg(void);
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
extern void _vecevent_coreneuron_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"modfiles_core//AMPA_NMDA_STP.mod\"");
    fprintf(stderr, " \"modfiles_core//cadad.mod\"");
    fprintf(stderr, " \"modfiles_core//cal2.mod\"");
    fprintf(stderr, " \"modfiles_core//can_mig.mod\"");
    fprintf(stderr, " \"modfiles_core//ChR2H134R_william_event.mod\"");
    fprintf(stderr, " \"modfiles_core//exp2syn_stp.mod\"");
    fprintf(stderr, " \"modfiles_core//GABA_A_STP.mod\"");
    fprintf(stderr, " \"modfiles_core//gap.mod\"");
    fprintf(stderr, " \"modfiles_core//Gfluct.mod\"");
    fprintf(stderr, " \"modfiles_core//h_kole.mod\"");
    fprintf(stderr, " \"modfiles_core//imCA3.mod\"");
    fprintf(stderr, " \"modfiles_core//kap_BS.mod\"");
    fprintf(stderr, " \"modfiles_core//kBK.mod\"");
    fprintf(stderr, " \"modfiles_core//kdmc_BS.mod\"");
    fprintf(stderr, " \"modfiles_core//kdr_BS.mod\"");
    fprintf(stderr, " \"modfiles_core//kdrCA3.mod\"");
    fprintf(stderr, " \"modfiles_core//kdrinter.mod\"");
    fprintf(stderr, " \"modfiles_core//leak.mod\"");
    fprintf(stderr, " \"modfiles_core//nainter.mod\"");
    fprintf(stderr, " \"modfiles_core//napCA3.mod\"");
    fprintf(stderr, " \"modfiles_core//natCA3.mod\"");
    fprintf(stderr, " \"modfiles_core//nax_BS.mod\"");
    fprintf(stderr, " \"modfiles_core//vecevent_coreneuron.mod\"");
    fprintf(stderr, "\n");
  }
  _AMPA_NMDA_STP_reg();
  _cadad_reg();
  _cal2_reg();
  _can_mig_reg();
  _ChR2H134R_william_event_reg();
  _exp2syn_stp_reg();
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
  _vecevent_coreneuron_reg();
}

#if defined(__cplusplus)
}
#endif
