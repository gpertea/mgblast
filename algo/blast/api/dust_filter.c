static char const rcsid[] = "$Id: dust_filter.c,v 1.6 2006/01/06 15:04:38 madden Exp $";

/*
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
 * File Name:  $RCSfile: dust_filter.c,v $
 *
 * Author: Tom Madden
 *
 */

/** @file dust_filter.c
 * Dust filtering API for the new BLAST code
 */

/* Prototypes of functions defined below */
#include <algo/blast/api/dust_filter.h>
#include <algo/blast/api/blast_api.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/api/seqsrc_readdb.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_inline.h>
#include <algo/blast/core/blast_dust.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

static Int2
s_GetFilteringLocations(BLAST_SequenceBlk* query_blk, BlastQueryInfo* query_info, SeqLoc* query_seqloc, const SBlastFilterOptions* filter_options, BlastMaskLoc** filter_maskloc)
{
    Int4 context = 0; /* loop variable. */
    const Boolean kIsNucl = TRUE;
    Boolean no_forward_strand = (query_info->first_context > 0);  /* filtering needed on reverse strand. */
    SeqLoc* slp_var = query_seqloc;

    ASSERT(query_info && query_blk && filter_maskloc && query_seqloc);

    *filter_maskloc = BlastMaskLocNew(query_info->last_context+1);

    for (context = query_info->first_context;
         context <= query_info->last_context && slp_var; ++context) {
      
        Boolean reverse = BlastIsReverseStrand(kIsNucl, context);
        Int4 query_length = query_info->contexts[context].query_length;
        

        /* For each query, check if forward strand is present */
        if (query_length <= 0)
        {
            if (kIsNucl && (context & 1) == 0)  /* Needed only for blastn, or does this not apply FIXME */
               no_forward_strand = TRUE;  /* No plus strand, we cannot simply infer locations by going from plus to minus */
            continue;
        }
        else if (!reverse)  /* This is a plus strand, safe to set no_forward_strand to FALSE as clearly there is one. */
               no_forward_strand = FALSE;

        if (!reverse || no_forward_strand)
        {
            BlastSeqLoc *filter_slp = NULL;   /* Used to hold combined SeqLoc's */
            Int4 filter_index = context;
            Int4 context_offset = query_info->contexts[context].query_offset;
            Uint1* buffer = &query_blk->sequence[context_offset];
            SDustOptions* dust_options = filter_options->dustOptions;

            if (BlastIsReverseStrand(kIsNucl, context) == TRUE)
            {  /* Reverse this as it's on minus strand. */
                  BlastSeqLoc* filter_slp_rev = NULL;
                  SeqBufferDust(buffer, query_length, 0, dust_options->level, 
                       dust_options->window, dust_options->linker, &filter_slp);

                  /* Reverse this relative to the part of the query being searched, leave it up to 
                    BlastMaskLocToSeqLoc to put it into the context of the entire query sequence (and
                    not just that part being searched). */
                  filter_slp_rev = BlastSeqLocReverse(filter_slp, query_length);
                  filter_slp = BlastSeqLocFree(filter_slp);
                  filter_slp = filter_slp_rev;
            }
            else
            {
                   SeqBufferDust(buffer, query_length, 0, dust_options->level, 
                       dust_options->window, dust_options->linker, &filter_slp);
            }

             (*filter_maskloc)->seqloc_array[filter_index] = filter_slp;
        }

        if (slp_var->choice == SEQLOC_WHOLE) 
        {
            if (BlastIsReverseStrand(kIsNucl, context) == TRUE)
                slp_var = slp_var->next;
        }
        else
        {
            slp_var = slp_var->next;
        }
    }

    return 0;
}

Int2
Blast_FindDustSeqLoc(SeqLoc* query_seqloc,
                             const SBlastOptions* options,
                             SeqLoc* *mask_loc)
{
    Int2 status = 0;
    BlastMaskLoc* filter_loc = NULL;
    BLAST_SequenceBlk* query_blk = NULL;
    BlastQueryInfo* query_info = NULL;
    QuerySetUpOptions* qsup_options = options->query_options;
    const EBlastProgramType kProgram = eBlastTypeBlastn;

    *mask_loc = NULL;

    /* If dust filtering not requested, return success. */
    if (qsup_options->filtering_options == NULL || qsup_options->filtering_options->dustOptions == NULL)
        return 0;

    BLAST_SetUpQuery(kProgram, query_seqloc, qsup_options, NULL, &query_info, &query_blk);

    status = s_GetFilteringLocations(query_blk, query_info, query_seqloc, qsup_options->filtering_options, &filter_loc);

    query_info = BlastQueryInfoFree(query_info);
    query_blk = BlastSequenceBlkFree(query_blk);

    if (filter_loc)
    {
            *mask_loc = BlastMaskLocToSeqLoc(kProgram, filter_loc, query_seqloc);
            filter_loc = BlastMaskLocFree(filter_loc);
    }
    
    return status;
}

/* @} */

