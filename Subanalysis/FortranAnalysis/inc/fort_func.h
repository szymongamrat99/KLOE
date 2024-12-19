// Author: Szymon Gamrat
// Date of last update: 26.05.2024

#ifndef FORT_FUNC_H
#define FORT_FUNC_H

#ifdef __cplusplus
extern "C"
{
#endif

  void clearstruct_(); // Setting all elements of INTERF to 0
  void GetKslEvent_(int *ntmc, int *mother, int *vtxmc, int *pidmc, float *xvmc, float *yvmc,
                    float *zvmc, float *pxmc, float *pymc, float *pzmc, int *nvtxmc, float *ipmc, float *KchMC, float *KneMC, float *DtMC, float *DlMC, bool *truth, bool *truthreg, bool *truthsemi, bool *truththree, bool *truthomega, bool *truthelse,
                    bool *truthdouble);

  void GetTdiff_(float *Kch, float *Kne, float *PhiVtx, float *PhiP, float *Broots, float *DL,
                 float *DT);

  void find_kchrec_(int *findKS, int *findKL, int *last_vtx, int *findClose, float *Bx, float *By,
                    float *Bz, int *qualv, int *nv, int *ntv, int *IV, float *CurV, float *PhiV, float *CotV, float *xv, float *yv, float *zv, int *vtaken, float *KchRec, float *trk1, float *trk2, float *cosTrk);

  void cor_ip_boost_(float *Kchrec, float *Bpx, float *Bpy, float *Bpz, float *Bx, float *By, float *Bz, float *Broots, float *trk1, float *trk2, float *KchBoost, float *ip_closest, float *ip_plane, float *chdist, float *Qmiss);

  int cos_pipi_kaonCM_ (float *Bpx, float *Bpy, float *Bpz, float *Broots, float *ip_closest, float *KchLAB, float *Ppich1, float *Ppich2, float *costrk);

  void find_neuvtx_(float *Bpx, float *Bpy, float *Bpz, float *Broots, float *KchBoost, int *ncl, float *enecl, int *ncll, float *xcl, float *ycl, float *zcl, float *ip, float *tcl, float *cldist, float *KneRecLor, float *trc, int *nclwrong, int *ncllwrong, float *KneRec, float *minv4gam, float *pi0, int *g4taken, float *trcv, float *PgamRecTaken, int *gpairtaken, float *test, float *g4vtxerr);

  void test_array_(double *array, float *numbers, double *array1);
  void summary_();

#ifdef __cplusplus
}
#endif

#endif // !FORT_FUNC_H