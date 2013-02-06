/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#include <stdio.h>
#include <stdlib.h> 

/* be sure to add <path-to-sdk>/itf to your include path */
#include <sra/sradb.h>

int main ( int argc, char *argv [] ) {
    rc_t            rc;
    SRAMgr const    *sra;
    SRATable const  *tbl;
    spotid_t        max;
    const char      *table_to_read = "myTable";

    SRAColumn const *read;
    SRAColumn const *qual;
    const void      *col_data;
    bitsz_t         off, sz;
    spotid_t        id;
    

    if (argc == 2)
        table_to_read = argv[1];
    
    printf("initializing manager\n");

    rc = SRAMgrMakeRead(&sra);

    if (rc != 0) {
      /*              fprintf(stderr,"failed initializing sra mgr (%s)\n",SRAErrToEnglish(rc,NULL));*/
      fprintf(stderr, "failed initializing sra mgr\n");
       /*
       * return SRAErrMake(rc);
       */
      return rc;
    }

    printf("opening table\n");

    rc = SRAMgrOpenTableRead (sra, &tbl, table_to_read);

    if (rc != 0) {
      /*        fprintf(stderr,"failed opening table (%s)\n",SRAErrToEnglish(SRAErrMake(rc),NULL));*/
      fprintf(stderr, "failed opening table\n");
      SRAMgrRelease(sra);
      
        return rc;
    }

    printf("getting max spot id\n");

    rc = SRATableMaxSpotId(tbl, &max);

    if (rc != 0) {
        fprintf(stderr,"failed getting max spot id\n");
        SRATableRelease(tbl);
        SRAMgrRelease(sra);

	/*        return SRAErrMake(rc);*/
	return rc;
    }
    
    printf("opening READ column\n");

    rc = SRATableOpenColumnRead(tbl, &read, "READ", NULL /*"INDSC:fasta" */);
    
    if (rc != 0) {
      /*        fprintf(stderr,"failed opening READ column (%s)",SRAErrToEnglish(SRAErrMake(rc),NULL));*/
      fprintf(stderr, "failed opening READ column\n");
        SRATableRelease(tbl);
        SRAMgrRelease(sra);

        return rc;
    }

    printf("opening QUALITY column\n"); 

    rc = SRATableOpenColumnRead(tbl, &qual, "QUALITY", NULL /* "NCBI:qual1" */);

    if (rc != 0) {
      /*        fprintf(stderr,"failed opening QUALITY column (%s)",SRAErrToEnglish(SRAErrMake(rc),NULL));*/
      fprintf(stderr, "failed opening QUALITY column\n");
        SRAColumnRelease(read);
        SRATableRelease(tbl);
        SRAMgrRelease(sra);

	/*        return SRAErrMake(rc);*/
	return rc;
    }

    for (id = 1; id <= max; ++id) {

        printf("reading READ column\n");
        rc = SRAColumnRead(read, id, &col_data, &off, &sz);

        if (rc) {
	  /*            fprintf(stderr,"failed reading READ column (%s)",SRAErrToEnglish(SRAErrMake(rc),NULL));*/
	  fprintf(stderr, "failed reading READ column\n");
            continue;
        }

        printf("off: %lu, sz: %lu, bases: %lu\n", off, sz, sz / 8 );
        

        printf("reading QUALITY column\n");

        rc = SRAColumnRead(qual, id, &col_data, &off, &sz);

        if (rc) {
	  /*            fprintf(stderr,"failed reading QUALITY column (%s)",SRAErrToEnglish(SRAErrMake(rc),NULL));*/
	  fprintf(stderr, "failed reading QUALITY column\n");
            continue;
        }

        printf("off: %lu, sz: %lu, bases: %lu\n", off, sz, sz / 8);

    } /* end for-loop on spot-id's */


    /*
     * NOTE - It is VERY important to release the toolkit objects once their use is no 
     * longer required
     */
    SRAColumnRelease(qual);
    SRAColumnRelease(read);
    SRATableRelease(tbl);
    SRAMgrRelease(sra);

    /*    return SRAErrMake(rc);*/
    return rc;
}
