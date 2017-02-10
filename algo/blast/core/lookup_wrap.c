/* $Id: lookup_wrap.c,v 1.21 2006/04/27 19:33:39 madden Exp $
 * ===========================================================================
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
 * Author: Ilya Dondoshansky
 *
 */

/** @file lookup_wrap.c
 * Wrapper for different flavors of lookup tables allowing a uniform interface in the code.
 * The wrapper (LookupTableWrap) contains an unsigned byte specifying the type of lookup 
 * table as well as a void pointer pointing to the actual lookup table.  Examples of different 
 * types of lookup tables are those for protein queries, the "standard" nucleotide one, the 
 * megablast lookup table, etc.
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: lookup_wrap.c,v 1.21 2006/04/27 19:33:39 madden Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <algo/blast/core/lookup_wrap.h>
#include <algo/blast/core/lookup_util.h>
#include <algo/blast/core/blast_lookup.h>
#include <algo/blast/core/mb_lookup.h>
#include <algo/blast/core/phi_lookup.h>
#include <algo/blast/core/blast_rps.h>

Int2 LookupTableWrapInit(BLAST_SequenceBlk* query, 
        const LookupTableOptions* lookup_options,	
        BlastSeqLoc* lookup_segments, BlastScoreBlk* sbp, 
        LookupTableWrap** lookup_wrap_ptr, const BlastRPSInfo *rps_info)
{
   Int2 status = 0;
   Int4 num_table_entries;
   LookupTableWrap* lookup_wrap;
   const Int4 kNucEntriesCutoff = 8500;  /* probably machine dependent */

   /* Construct the lookup table. */
   *lookup_wrap_ptr = lookup_wrap = 
      (LookupTableWrap*) calloc(1, sizeof(LookupTableWrap));
   lookup_wrap->lut_type = lookup_options->lut_type;

   switch ( lookup_options->lut_type ) {
   case AA_LOOKUP_TABLE:
       {
       Int4** matrix = NULL;
       Boolean has_pssm = FALSE;
       if (sbp->psi_matrix && sbp->psi_matrix->pssm) {
           matrix = sbp->psi_matrix->pssm->data;
           has_pssm = TRUE;
       } else {
           matrix = sbp->matrix->data;
       }
       BlastAaLookupNew(lookup_options, (BlastLookupTable* *)
                        &lookup_wrap->lut);
       ((BlastLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
       BlastAaLookupIndexQuery( (BlastLookupTable*) lookup_wrap->lut, matrix, 
                                 query, lookup_segments);
       _BlastAaLookupFinalize((BlastLookupTable*) lookup_wrap->lut);
       }
      break;
   case NA_LOOKUP_TABLE:
   case MB_LOOKUP_TABLE:
      /* choose either a standard or a megablast lookup table,
         depending on the query size and word size. Megablast
         tables are used for large queries and word sizes, and
         standard tables are used otherwise. For word size 11 
         the standard lookup table is especially efficient, 
         so the cutoff is larger in that case. Discontiguous
         megablast must always use a megablast table */

      num_table_entries = EstimateNumTableEntries(lookup_segments);

      if (lookup_options->mb_template_length > 0 ||
          (lookup_options->word_size == 11 && 
           num_table_entries > 2*kNucEntriesCutoff) ||
          (lookup_options->word_size >= 10 && 
           lookup_options->word_size != 11 &&
           num_table_entries > kNucEntriesCutoff) ) {

         lookup_wrap->lut_type = MB_LOOKUP_TABLE;
         MB_LookupTableNew(query, lookup_segments, 
                           (BlastMBLookupTable* *) &(lookup_wrap->lut), 
                           lookup_options, num_table_entries);
      }
      else {
         lookup_wrap->lut_type = NA_LOOKUP_TABLE;
         LookupTableNew(lookup_options, 
                        (BlastLookupTable* *) &(lookup_wrap->lut), 
                        num_table_entries, FALSE);
	    
         BlastNaLookupIndexQuery((BlastLookupTable*) lookup_wrap->lut, query,
                                 lookup_segments);
         _BlastAaLookupFinalize((BlastLookupTable*) lookup_wrap->lut);
      }
      break;
   case PHI_AA_LOOKUP: case PHI_NA_LOOKUP:
       {
           Blast_Message* error_msg = NULL;
           const Boolean kIsDna = (lookup_options->lut_type == PHI_NA_LOOKUP);
           status = SPHIPatternSearchBlkNew(lookup_options->phi_pattern, kIsDna, sbp,
                             (SPHIPatternSearchBlk* *) &(lookup_wrap->lut),
                             &error_msg);
           /** @todo FIXME: this error message must be passed further up!!! */
           if (error_msg) {
               Blast_MessagePost(error_msg);
               Blast_MessageFree(error_msg);
           }
           break;
       }
   case RPS_LOOKUP_TABLE:
      status = RPSLookupTableNew(rps_info, (BlastRPSLookupTable* *)(&lookup_wrap->lut));
      break;
      
   default:
      {
         /* FIXME - emit error condition here */
      }
   } /* end switch */

   return status;
}

LookupTableWrap* LookupTableWrapFree(LookupTableWrap* lookup)
{
   if (!lookup)
       return NULL;

   if (lookup->lut_type == MB_LOOKUP_TABLE) {
      lookup->lut = (void*) 
         MBLookupTableDestruct((BlastMBLookupTable*)lookup->lut);
   } else if (lookup->lut_type == PHI_AA_LOOKUP || 
              lookup->lut_type == PHI_NA_LOOKUP) {
       lookup->lut = (void*)
           SPHIPatternSearchBlkFree((SPHIPatternSearchBlk*)lookup->lut);
   } else if (lookup->lut_type == RPS_LOOKUP_TABLE) {
      lookup->lut = (void*) 
         RPSLookupTableDestruct((BlastRPSLookupTable*)lookup->lut);
   } else {
      lookup->lut = (void*) 
         LookupTableDestruct((BlastLookupTable*)lookup->lut);
   }
   sfree(lookup);
   return NULL;
}

Int4 GetOffsetArraySize(LookupTableWrap* lookup)
{
   Int4 offset_array_size;

   switch (lookup->lut_type) {
   case MB_LOOKUP_TABLE:
      offset_array_size = OFFSET_ARRAY_SIZE + 
         ((BlastMBLookupTable*)lookup->lut)->longest_chain;
      break;
   case AA_LOOKUP_TABLE: case NA_LOOKUP_TABLE:
      offset_array_size = OFFSET_ARRAY_SIZE + 
         ((BlastLookupTable*)lookup->lut)->longest_chain;
      break;
   default:
      offset_array_size = OFFSET_ARRAY_SIZE;
      break;
   }
   return offset_array_size;
}
