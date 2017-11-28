#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _CaP_reg(void);
extern void _CaP_nmodl_reg(void);
extern void _CaT_reg(void);
extern void _CurrentClampExt_reg(void);
extern void _K2_reg(void);
extern void _K2_nmodl_reg(void);
extern void _KA_reg(void);
extern void _KA_nmodl_reg(void);
extern void _KC_reg(void);
extern void _KC_nmodl_reg(void);
extern void _KMnew2_reg(void);
extern void _KMnew2_nmodl_reg(void);
extern void _Kdr_reg(void);
extern void _Kdr_nmodl_reg(void);
extern void _Kh1_reg(void);
extern void _Kh1_nmodl_reg(void);
extern void _Kh2_reg(void);
extern void _Kh2_nmodl_reg(void);
extern void _LeakConductance_reg(void);
extern void _NaF_reg(void);
extern void _NaF_nmodl_reg(void);
extern void _NaP_reg(void);
extern void _NaP_nmodl_reg(void);
extern void _cad_reg(void);
extern void _cad_nmodl_altered_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," CaP.mod");
    fprintf(stderr," CaP_nmodl.mod");
    fprintf(stderr," CaT.mod");
    fprintf(stderr," CurrentClampExt.mod");
    fprintf(stderr," K2.mod");
    fprintf(stderr," K2_nmodl.mod");
    fprintf(stderr," KA.mod");
    fprintf(stderr," KA_nmodl.mod");
    fprintf(stderr," KC.mod");
    fprintf(stderr," KC_nmodl.mod");
    fprintf(stderr," KMnew2.mod");
    fprintf(stderr," KMnew2_nmodl.mod");
    fprintf(stderr," Kdr.mod");
    fprintf(stderr," Kdr_nmodl.mod");
    fprintf(stderr," Kh1.mod");
    fprintf(stderr," Kh1_nmodl.mod");
    fprintf(stderr," Kh2.mod");
    fprintf(stderr," Kh2_nmodl.mod");
    fprintf(stderr," LeakConductance.mod");
    fprintf(stderr," NaF.mod");
    fprintf(stderr," NaF_nmodl.mod");
    fprintf(stderr," NaP.mod");
    fprintf(stderr," NaP_nmodl.mod");
    fprintf(stderr," cad.mod");
    fprintf(stderr," cad_nmodl_altered.mod");
    fprintf(stderr, "\n");
  }
  _CaP_reg();
  _CaP_nmodl_reg();
  _CaT_reg();
  _CurrentClampExt_reg();
  _K2_reg();
  _K2_nmodl_reg();
  _KA_reg();
  _KA_nmodl_reg();
  _KC_reg();
  _KC_nmodl_reg();
  _KMnew2_reg();
  _KMnew2_nmodl_reg();
  _Kdr_reg();
  _Kdr_nmodl_reg();
  _Kh1_reg();
  _Kh1_nmodl_reg();
  _Kh2_reg();
  _Kh2_nmodl_reg();
  _LeakConductance_reg();
  _NaF_reg();
  _NaF_nmodl_reg();
  _NaP_reg();
  _NaP_nmodl_reg();
  _cad_reg();
  _cad_nmodl_altered_reg();
}
