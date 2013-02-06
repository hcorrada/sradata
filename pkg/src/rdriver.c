#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>
#include <string.h>

#include "sra.h"

#define SRATABLE(argname) (((sra_tbl_holder *) R_ExternalPtrAddr(argname))->tbl)
#define ARRAYINDEX(i, j, k, ni, nj) (k*ni*nj + j*ni + i)

typedef struct sra_tbl_holder {
  const SRAMgr *sra;
  const SRATable *tbl;
} sra_tbl_holder;

void sra_finalizer(SEXP sra_tbl) {
  sra_tbl_holder *h = (sra_tbl_holder *) R_ExternalPtrAddr(sra_tbl);

  if (! h) return;
  closeTable(h->sra, h->tbl);
  Free(h);
}

SEXP sra_open(SEXP filePath) {
  int rc;
  sra_tbl_holder *h = (sra_tbl_holder *) Calloc(1, sra_tbl_holder);
  const char *table_to_read;

  table_to_read = CHAR(STRING_ELT(filePath, 0));
  rc = openTable(table_to_read, &(h->sra), &(h->tbl));
  if (rc != 0) {
    error("Error opening SRA table %s\n", table_to_read);
  }
  
  SEXP e_ptr = R_MakeExternalPtr(h, install("sra_handle"), R_NilValue);
  R_RegisterCFinalizerEx(e_ptr, sra_finalizer, TRUE);

  return e_ptr;
}

SEXP sra_get_info(SEXP sraTable) {
  SEXP res;
  SEXP maxSpotId;
  SEXP nLanes, lanes, nTiles;
  SEXP starts, rlens;

  const SRATable *tbl;
  rc_t rc;
  
  int MAX_NLANES = 8;
  int MAX_NTILES = 999;
  
  tbl = SRATABLE(sraTable);

  PROTECT(maxSpotId = allocVector(INTSXP, 1));
  PROTECT(nLanes = allocVector(INTSXP, 1));
  PROTECT(lanes = allocVector(INTSXP, MAX_NLANES));
  PROTECT(nTiles = allocVector(INTSXP, MAX_NLANES));
  PROTECT(starts = allocVector(INTSXP, MAX_NLANES * MAX_NTILES));
  PROTECT(rlens = allocVector(INTSXP, MAX_NLANES));

  rc = getInfo2(tbl, INTEGER(maxSpotId),
		INTEGER(nLanes),
		INTEGER(lanes),
		INTEGER(nTiles),
		INTEGER(starts),
		INTEGER(rlens));
  if (rc != 0) {
    error("Error getting table info\n");
  }

  PROTECT(res = allocVector(VECSXP, 6));
  SET_VECTOR_ELT(res, 0, maxSpotId);
  SET_VECTOR_ELT(res, 1, nLanes);
  SET_VECTOR_ELT(res, 2, lanes);
  SET_VECTOR_ELT(res, 3, nTiles);
  SET_VECTOR_ELT(res, 4, starts);
  SET_VECTOR_ELT(res, 5, rlens);
  UNPROTECT(7);
  return res;
}

SEXP sra_getBases(SEXP sraTable, SEXP startSpotId, SEXP endSpotId, SEXP rlen) {
  SEXP reads;
  spotid_t startId, endId, id;
  size_t nreads, ncycles;
  char *tmpchar;
  int i;

  const SRATable *tbl;
  SRAColumn const *read;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;

  tbl = SRATABLE(sraTable);
  startId = *INTEGER(startSpotId);
  endId = *INTEGER(endSpotId);
  nreads = endId - startId + 1;
  ncycles = *INTEGER(rlen);
  
  tmpchar = (char *) R_alloc(ncycles+1, sizeof(char));
  PROTECT(reads = allocVector(STRSXP, nreads));

  rc = SRATableOpenColumnRead(tbl, &read, "READ", insdc_fasta_t);
  if (rc != 0) {
    error("Error opening column READ (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
  }

  for (i=0, id=startId; id<=endId; id++, i++) {
    rc = SRAColumnRead(read, id, &col_data, &off, &sz);
    sz2 = sz/8;
    if (rc != 0 || off!=0 || (sz2 != ncycles)) {
      printf("Error reading spot %d in column READ (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }
    
    memcpy((void *) tmpchar, col_data, sz2);
    tmpchar[sz2] = '\0';
    SET_STRING_ELT(reads, i, mkChar(tmpchar));
  }

  SRAColumnRelease(read);
  UNPROTECT(1);
  return reads;
}

SEXP sra_getIds(SEXP sraTable, SEXP startSpotId, SEXP endSpotId) {
  SEXP names;
  spotid_t startId, endId, id;
  size_t nreads;
  char *tmpchar;
  int i;
  int MAX_NAMELEN=512;

  const SRATable *tbl;
  SRAColumn const *name;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;

  tbl = SRATABLE(sraTable);
  startId = *INTEGER(startSpotId);
  endId = *INTEGER(endSpotId);
  nreads = endId - startId + 1;
  
  tmpchar = (char *) R_alloc(MAX_NAMELEN, sizeof(char));
  PROTECT(names = allocVector(STRSXP, nreads));

  rc = SRATableOpenColumnRead(tbl, &name, "NAME", vdb_ascii_t);
  if (rc != 0) {
    error("Error opening column NAME (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
  }

  for (i=0, id=startId; id<=endId; id++, i++) {
    rc = SRAColumnRead(name, id, &col_data, &off, &sz);
    sz2 = sz/8;
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column NAME (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }
    
    memcpy((void *) tmpchar, col_data, sz2);
    tmpchar[sz2] = '\0';
    SET_STRING_ELT(names, i, mkChar(tmpchar));
  }

  SRAColumnRelease(name);
  UNPROTECT(1);
  return names;
}

SEXP sra_getQuals(SEXP sraTable, SEXP startSpotId, SEXP endSpotId, SEXP rlen) {
  SEXP quals;
  spotid_t startId, endId, id;
  size_t nreads, ncycles;
  uint8_t *buffer;
  char *tmpchar;

  int i, j;

  const SRATable *tbl;
  SRAColumn const *qual;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;

  tbl = SRATABLE(sraTable);
  startId = *INTEGER(startSpotId);
  endId = *INTEGER(endSpotId);
  nreads = endId - startId + 1;
  ncycles = *INTEGER(rlen);

  buffer = (uint8_t *) R_alloc(ncycles+1, sizeof(uint8_t));
  tmpchar = (char *) R_alloc(ncycles+1, sizeof(char));
  PROTECT(quals = allocVector(STRSXP, nreads));

  rc = SRATableOpenColumnRead(tbl, &qual, "QUALITY", ncbi_qual1_t);
  if (rc != 0) {
    error("Error opening column QUALITY (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
  }

  for (i=0, id=startId; id<=endId; id++, i++) {
    rc = SRAColumnRead(qual, id, &col_data, &off, &sz);
    sz2 = sz/8;

    if (rc != 0 || off!=0 || (sz2 != ncycles)) {
      printf("Error reading spot %d in column QUALITY (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }
    
    memcpy((void *) buffer, col_data, sz2);
    for (j=0; j<ncycles; j++) {
      tmpchar[j] = (char) (buffer[j] + 33);
    }

    tmpchar[sz2] = '\0';
    SET_STRING_ELT(quals, i, mkChar(tmpchar));
  }

  SRAColumnRelease(qual);
  UNPROTECT(1);
  return quals;
}

SEXP sra_getReadInfo(SEXP sraTable, SEXP startSpotId, SEXP endSpotId) {
  SEXP res;
  spotid_t startId, endId, id;
  size_t nreads;
  int i;

  const SRATable *tbl;
  SRAColumn const *name;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;

  SEXP pLane, pTile, pX, pY;
  int *lane, *tile, *x, *y;

  tbl = SRATABLE(sraTable);
  startId = *INTEGER(startSpotId);
  endId = *INTEGER(endSpotId);
  nreads = endId - startId + 1;
  
  PROTECT(pLane = allocVector(INTSXP, nreads));
  PROTECT(pTile = allocVector(INTSXP, nreads));
  PROTECT(pX = allocVector(INTSXP, nreads));
  PROTECT(pY = allocVector(INTSXP, nreads));

  lane = INTEGER(pLane);
  tile = INTEGER(pTile);
  x = INTEGER(pX);
  y = INTEGER(pY);

  rc = SRATableOpenColumnRead(tbl, &name, "NAME", vdb_ascii_t);
  if (rc != 0) {
    error("Error opening column NAME (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
  }

  for (i=0, id=startId; id<=endId; id++, i++) {
    rc = SRAColumnRead(name, id, &col_data, &off, &sz);
    sz2 = sz/8;
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column NAME (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }
    
    const SRATableData *tdata = SRATableGetTableData(tbl);
    lane[i] = tdata->coord.lane;
    tile[i] = tdata->coord.tile;
    x[i] = tdata->coord.x;
    y[i] = tdata->coord.y;
  }

  SRAColumnRelease(name);

  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, pLane);
  SET_VECTOR_ELT(res, 1, pTile);
  SET_VECTOR_ELT(res, 2, pX);
  SET_VECTOR_ELT(res, 3, pY);
  
  UNPROTECT(5);
  return res;
}

SEXP sra_getIntensity(SEXP sraTable, SEXP startSpotId, SEXP endSpotId, SEXP rlen) {
  SEXP pInts;
  double *ints;

  spotid_t startId, endId, id;
  size_t nreads, ncycles;
  float *buffer;

  const SRATable *tbl;
  SRAColumn const *intensity;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;
  int j;

  tbl = SRATABLE(sraTable);
  startId = *INTEGER(startSpotId);
  endId = *INTEGER(endSpotId);
  nreads = endId - startId + 1;
  ncycles = *INTEGER(rlen);

  buffer = (float *) R_alloc(4 * ncycles, sizeof(float));
  PROTECT(pInts = allocVector(REALSXP, nreads * 4 * ncycles));
  ints = REAL(pInts);

  rc = SRATableOpenColumnRead(tbl, &intensity, "INTENSITY", ncbi_fsamp4_t);
  if (rc != 0) {
    error("Error opening column INTENSITY (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
  }

  for (id=startId; id<=endId; id++, ints += 4*ncycles) {
    rc = SRAColumnRead(intensity, id, &col_data, &off, &sz);
    sz2 = sz/8;

    if (rc != 0 || off!=0 || (sz2 != 4 * ncycles * sizeof(float))) {
      printf("Error reading spot %d in column INTENSITY (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }
    
    memcpy((void *) buffer, col_data, sz2);
    for (j=0; j<4*ncycles; j++) {
      ints[j] = (double) buffer[j];
    }
  }

  SRAColumnRelease(intensity);
  UNPROTECT(1);
  return pInts;
}

