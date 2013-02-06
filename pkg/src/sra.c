#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sra.h"


int openTable (const char *table_to_read, const SRAMgr **sra, const SRATable **tbl) {
  rc_t rc;

  rc = SRAMgrMakeRead(sra);
  if (rc != 0) {
    printf("failed initizalizing sra mgr (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  rc = SRAMgrOpenTableRead(*sra, tbl, table_to_read);
  if (rc != 0) {
    printf("failed opening table (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    SRAMgrRelease(*sra);
    return SRAErrMake(rc);
  }

  return SRAErrMake(rc);
}

int closeTable(const SRAMgr *sra, const SRATable *tbl) {
  SRATableRelease(tbl);
  SRAMgrRelease(sra);
  return 0;
}

int getBases(const SRATable *tbl, const spotid_t minSpotId, const spotid_t maxSpotId) {
  spotid_t id;
  SRAColumn const *read;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;
  char buffer[512 * 1024];

  rc = SRATableOpenColumnRead(tbl, &read, "READ", insdc_fasta_t);
  if (rc != 0) {
    printf("Error opening column READ (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  for (id=minSpotId; id<maxSpotId; id++) {
    rc = SRAColumnRead(read, id, &col_data, &off, &sz);
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column READ (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }
    
    sz2 = sz/8;
    memcpy((void *) buffer, col_data, sz2);
    buffer[sz2] = '\0';
    printf("%s\n", buffer);
  }
  SRAColumnRelease(read);
}

int getQuals(const SRATable *tbl, const spotid_t minSpotId, const spotid_t maxSpotId) {
  spotid_t id;
  SRAColumn const *quals;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;
  uint8_t buffer[512 * 1024];
  int i;

  rc = SRATableOpenColumnRead(tbl, &quals, "QUALITY", ncbi_qual1_t);
  if (rc != 0) {
    printf("Error opening column QUALITY (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  for (id=minSpotId; id<maxSpotId; id++) {
    rc = SRAColumnRead(quals, id, &col_data, &off, &sz);
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column QUALITY (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }

    if (sz != 0) {
      sz2 = sz/8;
      memcpy((void *) buffer, col_data, sz2);
      for (i=0; i<sz2; i++) {
	printf("%hu ", buffer[i]);
      }
      printf("\n");
    } else {
      printf("Empty quality string for spot %d\n", id);
    }
  }
  SRAColumnRelease(quals);
}

int getInfo(const SRATable *tbl, const spotid_t minSpotId, const spotid_t maxSpotId) {
  spotid_t id;
  SRAColumn const *name;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  size_t sz2;
  char buffer[512 * 1024];
  int i;

  rc = SRATableOpenColumnRead(tbl, &name, "NAME", vdb_ascii_t);
  if (rc != 0) {
    printf("Error opening column NAME (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  for (id=minSpotId; id<maxSpotId; id++) {
    rc = SRAColumnRead(name, id, &col_data, &off, &sz);
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column NAME (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }

    if (sz != 0) {
      sz2 = sz/8;
      memcpy((void *) buffer, col_data, sz2);
      printf("%s\n", buffer);
    } else {
      printf("Empty quality string for spot %d\n", id);
    }

    const SRATableData *tdata = SRATableGetTableData(tbl);
    printf("%d\t%d\t%d\t%d\n", tdata->coord.lane, tdata->coord.tile, tdata->coord.x, tdata->coord.y);
  }
  SRAColumnRelease(name);
}

int getInfo2(const SRATable *tbl, spotid_t *maxSpotId, int *nLanes, int *lanes, int *nTiles, int *starts, int *rlens) {
  rc_t rc;
  spotid_t id;

  SRAColumn const *name;
  SRAColumn const *read;

  const void *col_data;
  bitsz_t off, sz;
  int lane, tile;
  int curNTiles;
  int k;

  rc = SRATableMaxSpotId(tbl, maxSpotId);
  if (rc != 0) {
    printf("Failed getting max spot id (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  rc = SRATableOpenColumnRead(tbl, &name, "NAME", vdb_ascii_t);
  if (rc != 0) {
    printf("Error opening column NAME (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  rc = SRATableOpenColumnRead(tbl, &read, "READ", insdc_fasta_t);
  if (rc != 0) {
    printf("Error opening column READ (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  *nLanes = 0;
  lane = -1; tile = -1;
  k=0;

  for (id=1; id<(*maxSpotId); id++) {
    rc = SRAColumnRead(name, id, &col_data, &off, &sz);
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column NAME (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }

    const SRATableData *tdata = SRATableGetTableData(tbl);

    // check if new lane
    if (tdata->coord.lane != lane) {
      // save the number of tiles in lane just finished, if there is one
      if (lane != -1) {
	nTiles[*nLanes-1] = curNTiles;
      }

      // bump the number of lanes
      (*nLanes)++;
      k++;

      // save the new lane index
      lanes[(*nLanes)-1] = tdata->coord.lane;
      starts[k-1] = id;

      // get the read length
      rc = SRAColumnRead(read, id, &col_data, &off, &sz);
      if (rc != 0 || off!=0) {
	printf("Error reading spot %d in column READ (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
	sz = 0;
	continue;
      }
      rlens[(*nLanes)-1] = sz/8;

      // initialize the tile counter
      curNTiles = 1;

      // keep track of what you've seen
      lane = tdata->coord.lane;
      tile = tdata->coord.tile;
    }

    if (tdata->coord.tile != tile) {
      curNTiles++;
      tile = tdata->coord.tile;

      k++;
      starts[k-1] = id;
    }
  }
  nTiles[*nLanes-1] = curNTiles;

  SRAColumnRelease(name);
  return SRAErrMake(rc);
}

int getRanges(const SRATable *tbl) {
  spotid_t maxSpotId, id;

  SRAColumn const *name;
  rc_t rc;
  const void *col_data;
  bitsz_t off, sz;
  int lane, tile;

  rc = SRATableMaxSpotId(tbl, &maxSpotId);
  if (rc != 0) {
    printf("Failed getting max spot id (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  rc = SRATableOpenColumnRead(tbl, &name, "NAME", vdb_ascii_t);
  if (rc != 0) {
    printf("Error opening column NAME (%s)\n", SRAErrToEnglish(SRAErrMake(rc), NULL));
    return SRAErrMake(rc);
  }

  lane = -1; tile=-1;
  for (id=1; id<maxSpotId; id++) {
    rc = SRAColumnRead(name, id, &col_data, &off, &sz);
    if (rc != 0 || off!=0) {
      printf("Error reading spot %d in column NAME (%s)\n", id, SRAErrToEnglish(SRAErrMake(rc), NULL));
      continue;
    }

    const SRATableData *tdata = SRATableGetTableData(tbl);
    if (tdata->coord.lane != lane || tdata->coord.tile != tile) {
      if (id > 1) {
	printf("%d\n", id-1);
      }

      lane=tdata->coord.lane;
      tile=tdata->coord.tile;
      printf("%d %d: %d-", lane, tile, id);
    }
  }
  printf("%d\n", id);
  SRAColumnRelease(name);
}
