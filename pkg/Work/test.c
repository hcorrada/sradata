#include <sra.h>
#include <stdio.h>
#include <stdlib.h>

int main (int *argc, char **argv) {
  char *tab;
  char *accession;

  SRAMgr const *sra;
  SRATable const *tbl;
  spotid_t nspots;
  rc_t rc;

  tab = argv[1];
  accession = argv[2];

  rc = openTable(tab, &sra, &tbl, &nspots);
  if (rc != 0) {
    return rc;
  }
  printf("There are %d spots\n", nspots);

  //  rc = getReads(tbl, accession, 10);
  rc = getBases(tbl, 1, 10);
  if (rc != 0) {
    closeTable(sra, tbl);
    return 1;
  }

  printf("\n");
  rc = getQuals(tbl, 1, 10);
  if (rc != 0) {
    closeTable(sra, tbl);
    return 1;
  }

  getInfo(tbl, 1, 10);
  getRanges(tbl);
  closeTable(sra, tbl);
  return 0;
}
