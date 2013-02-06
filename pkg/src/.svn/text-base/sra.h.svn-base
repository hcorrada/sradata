#ifndef _SRA_H
#define _SRA_H

#include <klib/defs.h>
#include <vdb/types.h>
#include <klib/rc.h>
#include <sra/sradb.h>
#include <sra/sradb-priv.h>
#include <sra/illumina.h>

int openTable (const char *table_to_read, const SRAMgr **sra, const SRATable **tbl);
int closeTable (const SRAMgr *sra, const SRATable *tbl);
int getInfo2 (const SRATable *tbl, spotid_t *maxSpotId, int *nLanes,
	      int *lanes, int *nTiles, int *starts, int *rlens);

#endif
