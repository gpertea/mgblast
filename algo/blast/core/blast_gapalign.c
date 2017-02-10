/* $Id: blast_gapalign.c,v 1.170 2006/04/12 20:28:44 camacho Exp $
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
 *
 */

/** @file blast_gapalign.c
 * Functions to perform gapped alignment
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_gapalign.c,v 1.170 2006/04/12 20:28:44 camacho Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_gapalign.h>
#include <algo/blast/core/blast_util.h> /* for NCBI2NA_UNPACK_BASE macros */
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/greedy_align.h>
#include "blast_gapalign_priv.h"
#include "blast_hits_priv.h"
#include "blast_itree.h"

static Int2 s_BlastDynProgNtGappedAlignment(BLAST_SequenceBlk* query_blk, 
   BLAST_SequenceBlk* subject_blk, BlastGapAlignStruct* gap_align, 
   const BlastScoringParameters* score_params, BlastInitHSP* init_hsp);
static Int4 s_BlastAlignPackedNucl(Uint1* B, Uint1* A, Int4 N, Int4 M, 
   Int4* pej, Int4* pei, BlastGapAlignStruct* gap_align,
   const BlastScoringParameters* score_params, Boolean reverse_sequence);

static Int2 s_BlastProtGappedAlignment(EBlastProgramType program, 
   BLAST_SequenceBlk* query_in, BLAST_SequenceBlk* subject_in,
   BlastGapAlignStruct* gap_align,
   const BlastScoringParameters* score_params, BlastInitHSP* init_hsp);

/** Lower bound for scores. Divide by two to prevent underflows. */
#define MININT INT4_MIN/2

/** Minimal size of a chunk for state array allocation. */
#define	CHUNKSIZE	2097152

/** Retrieve the state structure corresponding to a given length
 * @param head Pointer to the first element of the state structures 
 *        array [in]
 * @param length The length for which the state structure has to be 
 *        found [in]
 * @return The found or created state structure
 */
static GapStateArrayStruct*
s_GapGetState(GapStateArrayStruct** head, Int4 length)

{
   GapStateArrayStruct*	retval,* var,* last;
   Int4	chunksize = MAX(CHUNKSIZE, length + length/3);

   length += length/3;	/* Add on about 30% so the end will get reused. */
   retval = NULL;
   if (*head == NULL) {
      retval = (GapStateArrayStruct*) 
         malloc(sizeof(GapStateArrayStruct));
      retval->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
      retval->length = chunksize;
      retval->used = 0;
      retval->next = NULL;
      *head = retval;
   } else {
      var = *head;
      last = *head;
      while (var) {
         if (length < (var->length - var->used)) {
            retval = var;
            break;
         } else if (var->used == 0) { 
            /* If it's empty and too small, replace. */
            sfree(var->state_array);
            var->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
            var->length = chunksize;
            retval = var;
            break;
         }
         last = var;
         var = var->next;
      }
      
      if (var == NULL)
      {
         retval = (GapStateArrayStruct*) malloc(sizeof(GapStateArrayStruct));
         retval->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
         retval->length = chunksize;
         retval->used = 0;
         retval->next = NULL;
         last->next = retval;
      }
   }

#ifdef ERR_POST_EX_DEFINED
   if (retval->state_array == NULL)
      ErrPostEx(SEV_ERROR, 0, 0, "state array is NULL");
#endif
		
   return retval;

}

/** Remove a state from a state structure */
static Boolean
s_GapPurgeState(GapStateArrayStruct* state_struct)
{
   while (state_struct)
   {
      /*
	memset(state_struct->state_array, 0, state_struct->used);
      */
      state_struct->used = 0;
      state_struct = state_struct->next;
   }
   
   return TRUE;
}

/** Deallocate the memory for greedy gapped alignment */
static SGreedyAlignMem* 
s_BlastGreedyAlignsFree(SGreedyAlignMem* gamp)
{
   if (gamp->last_seq2_off) {
      sfree(gamp->last_seq2_off[0]);
      sfree(gamp->last_seq2_off);
   } else {
      if (gamp->last_seq2_off_affine) {
         sfree(gamp->last_seq2_off_affine[0]);
         sfree(gamp->last_seq2_off_affine);
      }
      sfree(gamp->diag_bounds);
   }
   sfree(gamp->max_score);
   if (gamp->space)
      MBSpaceFree(gamp->space);
   sfree(gamp);
   return NULL;
}

/** Allocate memory for the greedy gapped alignment algorithm
 * @param score_params Parameters related to scoring [in]
 * @param ext_params Options and parameters related to the extension [in]
 * @param max_dbseq_length The length of the longest sequence in the 
 *        database [in]
 * @return The allocated SGreedyAlignMem structure
 */
static SGreedyAlignMem* 
s_BlastGreedyAlignMemAlloc(const BlastScoringParameters* score_params,
		       const BlastExtensionParameters* ext_params,
		       Int4 max_dbseq_length)
{
   SGreedyAlignMem* gamp;
   Int4 max_d, max_d_1, Xdrop, d_diff, max_cost, gd, i;
   Int4 reward, penalty, gap_open, gap_extend;
   Int4 Mis_cost, GE_cost;
   Boolean do_traceback;
   
   if (score_params == NULL || ext_params == NULL) 
      return NULL;
   
   do_traceback = 
      (ext_params->options->ePrelimGapExt != eGreedyExt);

   if (score_params->reward % 2 == 1) {
      reward = 2*score_params->reward;
      penalty = -2*score_params->penalty;
      Xdrop = 2*MAX(ext_params->gap_x_dropoff,
                    ext_params->gap_x_dropoff_final);
      gap_open = 2*score_params->gap_open;
      gap_extend = 2*score_params->gap_extend;
   } else {
      reward = score_params->reward;
      penalty = -score_params->penalty;
      Xdrop = MAX(ext_params->gap_x_dropoff,
                  ext_params->gap_x_dropoff_final);
      gap_open = score_params->gap_open;
      gap_extend = score_params->gap_extend;
   }

   if (gap_open == 0 && gap_extend == 0)
      gap_extend = reward / 2 + penalty;

   max_d = MIN(GREEDY_MAX_COST,
               max_dbseq_length / GREEDY_MAX_COST_FRACTION + 1);

   gamp = (SGreedyAlignMem*) calloc(1, sizeof(SGreedyAlignMem));

   if (score_params->gap_open==0 && score_params->gap_extend==0) {
      d_diff = (Xdrop+reward/2)/(penalty+reward)+1;
   
      gamp->last_seq2_off = (Int4**) malloc((max_d + 2) * sizeof(Int4*));
      if (gamp->last_seq2_off == NULL) {
         sfree(gamp);
         return NULL;
      }
      gamp->last_seq2_off[0] = 
         (Int4*) malloc((max_d + max_d + 6) * sizeof(Int4) * 2);
      if (gamp->last_seq2_off[0] == NULL) {
#ifdef ERR_POST_EX_DEFINED
	 ErrPostEx(SEV_WARNING, 0, 0, 
              "Failed to allocate %ld bytes for greedy alignment", 
              (max_d + max_d + 6) * sizeof(Int4) * 2);
#endif
         s_BlastGreedyAlignsFree(gamp);
         return NULL;
      }

      gamp->last_seq2_off[1] = gamp->last_seq2_off[0] + max_d + max_d + 6;
      gamp->last_seq2_off_affine = NULL;
      gamp->diag_bounds = NULL;
   } else {
      gamp->last_seq2_off = NULL;
      Mis_cost = reward + penalty;
      GE_cost = gap_extend+reward/2;
      max_d_1 = max_d;
      max_d *= GE_cost;
      max_cost = MAX(Mis_cost, gap_open+GE_cost);
      gd = BLAST_Gdb3(&Mis_cost, &gap_open, &GE_cost);
      d_diff = (Xdrop+reward/2)/gd+1;
      gamp->diag_bounds = (Int4*) calloc(2*(max_d+1+max_cost), sizeof(Int4));
      gamp->last_seq2_off_affine = (SGreedyOffset**) 
	 malloc((MAX(max_d, max_cost) + 2) * sizeof(SGreedyOffset*));
      if (!gamp->diag_bounds || !gamp->last_seq2_off_affine) {
         s_BlastGreedyAlignsFree(gamp);
         return NULL;
      }
      gamp->last_seq2_off_affine[0] = (SGreedyOffset*)
	 calloc((2*max_d_1 + 6) , sizeof(SGreedyOffset) * (max_cost+1));
      for (i = 1; i <= max_cost; i++)
	 gamp->last_seq2_off_affine[i] = 
	    gamp->last_seq2_off_affine[i-1] + 2*max_d_1 + 6;
      if (!gamp->last_seq2_off_affine || !gamp->last_seq2_off_affine[0]) {
         s_BlastGreedyAlignsFree(gamp);
         return NULL;
      }
   }
   gamp->max_score = (Int4*) malloc(sizeof(Int4) * (max_d + 1 + d_diff));

   if (do_traceback)
      gamp->space = MBSpaceNew(0);
   if (!gamp->max_score || (do_traceback && !gamp->space))
      /* Failure in one of the memory allocations */
      s_BlastGreedyAlignsFree(gamp);

   return gamp;
}

/* Documented in blast_gapalign.h */
BlastGapAlignStruct* 
BLAST_GapAlignStructFree(BlastGapAlignStruct* gap_align)
{
   if (!gap_align)
      return NULL;

   GapEditScriptDelete(gap_align->edit_script);
   GapPrelimEditBlockFree(gap_align->fwd_prelim_tback);
   GapPrelimEditBlockFree(gap_align->rev_prelim_tback);
   if (gap_align->greedy_align_mem)
      s_BlastGreedyAlignsFree(gap_align->greedy_align_mem);
   GapStateFree(gap_align->state_struct);
   sfree(gap_align->dp_mem);

   sfree(gap_align);
   return NULL;
}

/* Documented in blast_gapalign.h */
Int2
BLAST_GapAlignStructNew(const BlastScoringParameters* score_params, 
   const BlastExtensionParameters* ext_params, 
   Uint4 max_subject_length,
   BlastScoreBlk* sbp, BlastGapAlignStruct** gap_align_ptr)
{
   Int2 status = 0;
   BlastGapAlignStruct* gap_align;

   if (!gap_align_ptr || !sbp || !score_params || !ext_params)
      return -1;

   gap_align = (BlastGapAlignStruct*) calloc(1, sizeof(BlastGapAlignStruct));

   *gap_align_ptr = gap_align;

   gap_align->sbp = sbp;

   gap_align->gap_x_dropoff = ext_params->gap_x_dropoff;

   if (ext_params->options->ePrelimGapExt == eDynProgExt) {
      /* allocate structures for ordinary dynamic programming */
      gap_align->dp_mem_alloc = 1000;
      gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
      if (!gap_align->dp_mem)
         gap_align = BLAST_GapAlignStructFree(gap_align);
   }
   else {
      /* allocate structures for greedy dynamic programming */
      max_subject_length = MIN(max_subject_length, MAX_DBSEQ_LEN);
      gap_align->greedy_align_mem = 
         s_BlastGreedyAlignMemAlloc(score_params, ext_params, 
                                    max_subject_length);
      if (!gap_align->greedy_align_mem)
         gap_align = BLAST_GapAlignStructFree(gap_align);
   }

   if (!gap_align)
      return -1;

   gap_align->positionBased = (sbp->psi_matrix != NULL);

   gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
   gap_align->rev_prelim_tback = GapPrelimEditBlockNew();

   return status;
}

/** Values for the editing script operations in traceback */
enum {
    SCRIPT_SUB           = eGapAlignSub,     /**< Substitution */
    SCRIPT_GAP_IN_A      = eGapAlignDel,     /**< Deletion */
    SCRIPT_GAP_IN_B      = eGapAlignIns,     /**< Insertion */
    SCRIPT_OP_MASK       = 0x07, /**< Mask for edit script operations */

    SCRIPT_EXTEND_GAP_A  = 0x10, /**< continue a gap in A */
    SCRIPT_EXTEND_GAP_B  = 0x40, /**< continue a gap in B */
};

/** Low level function to perform dynamic programming gapped extension 
 * with traceback.
 * @param A The query sequence [in]
 * @param B The subject sequence [in]
 * @param M Maximal extension length in query [in]
 * @param N Maximal extension length in subject [in]
 * @param a_offset Resulting starting offset in query [out]
 * @param b_offset Resulting starting offset in subject [out]
 * @param edit_block Structure to hold traceback generated [out]
 * @param gap_align Structure holding various information and allocated 
 *        memory for the gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param query_offset The starting offset in query [in]
 * @param reversed Has the sequence been reversed? Used for psi-blast [in]
 * @param reverse_sequence Do reverse the sequence [in]
 * @return The best alignment score found.
*/
Int4
ALIGN_EX(Uint1* A, Uint1* B, Int4 M, Int4 N, Int4* a_offset, 
	Int4* b_offset, GapPrelimEditBlock *edit_block, 
        BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params, Int4 query_offset, 
        Boolean reversed, Boolean reverse_sequence)
	
{ 
    /* See Blast_SemiGappedAlign for more general comments on 
       what this code is doing; comments in this function
       only apply to the traceback computations */

    Int4 i; 
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    Uint1* b_ptr;
  
    BlastGapDP* score_array;

    Int4 gap_open;
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;
    Int4 best_score;
  
    Int4** matrix = NULL;
    Int4** pssm = NULL;
    Int4* matrix_row = NULL;
  
    Int4 score;
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;
  
    GapStateArrayStruct* state_struct;
    Uint1* edit_script_row;
    Uint1** edit_script;
    Int4 *edit_start_offset;
    Int4 edit_script_num_rows;
    Int4 orig_b_index;
    Uint1 script, next_script, script_row, script_col;
    Int4 num_extra_cells;

    matrix = gap_align->sbp->matrix->data;
    if (gap_align->positionBased) {
        pssm = gap_align->sbp->psi_matrix->pssm->data;
    }

    *a_offset = 0;
    *b_offset = 0;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    x_dropoff = gap_align->gap_x_dropoff;
  
    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;
  
    if(N <= 0 || M <= 0) 
        return 0;
  
    /* Initialize traceback information. edit_script[] is
       a 2-D array which is filled in row by row as the 
       dynamic programming takes place */

    s_GapPurgeState(gap_align->state_struct);

    /* Conceptually, traceback requires a byte array of size
       MxN. To save on memory, each row of the array is allocated
       separately. edit_script[i] points to the storage reserved
       for row i, and edit_start_offset[i] gives the offset into
       the B sequence corresponding to element 0 of edit_script[i].
       
       Also make the number of edit script rows grow dynamically */

    edit_script_num_rows = 100;
    edit_script = (Uint1**) malloc(sizeof(Uint1*) * edit_script_num_rows);
    edit_start_offset = (Int4*) malloc(sizeof(Int4) * edit_script_num_rows);

    /* allocate storage for the first row of the traceback
       array. Because row elements correspond to gaps in A,
       the alignment can only go x_dropoff/gap_extend positions
       at most before failing the X dropoff criterion */

    if (gap_extend > 0)
        num_extra_cells = x_dropoff / gap_extend + 3;
    else
        num_extra_cells = N + 3;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                                      2 * gap_align->dp_mem_alloc);
        sfree(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                                  sizeof(BlastGapDP));
    }

    state_struct = s_GapGetState(&gap_align->state_struct, num_extra_cells);

    edit_script[0] = state_struct->state_array;
    edit_start_offset[0] = 0;
    edit_script_row = state_struct->state_array;

    score = -gap_open_extend;
    score_array = gap_align->dp_mem;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
  
    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff) 
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score -= gap_extend;
        edit_script_row[i] = SCRIPT_GAP_IN_A;
    }
    state_struct->used = i + 1;
  
    b_size = i;
    best_score = 0;
    first_b_index = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;
  
    for (a_index = 1; a_index <= M; a_index++) {

        /* Set up the next row of the edit script; this involves
           allocating memory for the row, then pointing to it.
           It is not necessary to allocate space for offsets less
           than first_b_index (since the inner loop will never 
           look at them); 
           
           It is unknown at this point how far to the right the 
           current traceback row will extend; all that's known for
           sure is that the previous row fails the X-dropoff test
           after b_size cells, and that the current row can go at
           most num_extra_cells beyond that before failing the test */

        if (gap_extend > 0)
            state_struct = s_GapGetState(&gap_align->state_struct, 
                           b_size - first_b_index + num_extra_cells);
        else
            state_struct = s_GapGetState(&gap_align->state_struct, 
                                        N + 3 - first_b_index);

        if (a_index == edit_script_num_rows) {
            edit_script_num_rows = edit_script_num_rows * 2;
            edit_script = (Uint1 **)realloc(edit_script, 
                                            edit_script_num_rows *
                                            sizeof(Uint1 *));
            edit_start_offset = (Int4 *)realloc(edit_start_offset, 
                                                edit_script_num_rows *
                                                sizeof(Int4));
        }

        edit_script[a_index] = state_struct->state_array + 
                                state_struct->used + 1;
        edit_start_offset[a_index] = first_b_index;

        /* the inner loop assumes the current traceback
           row begins at offset zero of B */

        edit_script_row = edit_script[a_index] - first_b_index;
        orig_b_index = first_b_index;

        if (!(gap_align->positionBased)) {
            if(reverse_sequence)
                matrix_row = matrix[ A[ M - a_index ] ];
            else
                matrix_row = matrix[ A[ a_index ] ];
        }
        else {
            if(reversed || reverse_sequence)
                matrix_row = pssm[M - a_index];
            else
                matrix_row = pssm[a_index + query_offset];
        }

        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;

        for (b_index = first_b_index; b_index < b_size; b_index++) {

            b_ptr += b_increment;
            score_gap_col = score_array[b_index].best_gap;
            next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

            /* script, script_row and script_col contain the
               actions specified by the dynamic programming.
               when the inner loop has finished, 'script' con-
               tains all of the actions to perform, and is
               written to edit_script[a_index][b_index]. Otherwise,
               this inner loop is exactly the same as the one
               in Blast_SemiGappedAlign() */

            script = SCRIPT_SUB;
            script_col = SCRIPT_EXTEND_GAP_B;
            script_row = SCRIPT_EXTEND_GAP_A;

            if (score < score_gap_col) {
                script = SCRIPT_GAP_IN_B;
                score = score_gap_col;
            }
            if (score < score_gap_row) {
                script = SCRIPT_GAP_IN_A;
                score = score_gap_row;
            }

            if (best_score - score > x_dropoff) {

                if (first_b_index == b_index)
                    first_b_index++;
                else
                    score_array[b_index].best = MININT;
            }
            else {
                last_b_index = b_index;
                if (score > best_score) {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                if (score_gap_col < (score - gap_open_extend)) {
                    score_array[b_index].best_gap = score - gap_open_extend;
                }
                else {
                    score_array[b_index].best_gap = score_gap_col;
                    script += script_col;
                }

                if (score_gap_row < (score - gap_open_extend)) 
                    score_gap_row = score - gap_open_extend;
                else
                    script += script_row;

                score_array[b_index].best = score;
            }

            score = next_score;
            edit_script_row[b_index] = script;
        }
  
        if (first_b_index == b_size)
            break;

        if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                                          2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                                               gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }


        if (last_b_index < b_size - 1) {
            b_size = last_b_index + 1;
        }
        else {
            while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {

                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                edit_script_row[b_size] = SCRIPT_GAP_IN_A;
                b_size++;
            }
        }

        /* update the memory allocator to reflect the exact number
           of traceback cells this row needed */

        state_struct->used += MAX(b_index, b_size) - orig_b_index + 1;

        if (b_size <= N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }
    
    /* Pick the optimal path through the now complete
       edit_script[][]. This is equivalent to flattening 
       the 2-D array into a 1-D list of actions. */

    a_index = *a_offset;
    b_index = *b_offset;
    script = SCRIPT_SUB;

    while (a_index > 0 || b_index > 0) {
        /* Retrieve the next action to perform. Rows of
           the traceback array do not necessarily start
           at offset zero of B, so a correction is needed
           to point to the correct position */

        next_script = 
            edit_script[a_index][b_index - edit_start_offset[a_index]];

        switch(script) {
        case SCRIPT_GAP_IN_A:
            script = next_script & SCRIPT_OP_MASK;
            if (next_script & SCRIPT_EXTEND_GAP_A)
                script = SCRIPT_GAP_IN_A;
            break;

        case SCRIPT_GAP_IN_B:
            script = next_script & SCRIPT_OP_MASK;
            if (next_script & SCRIPT_EXTEND_GAP_B)
                script = SCRIPT_GAP_IN_B;
            break;

        default:
            script = next_script & SCRIPT_OP_MASK;
            break;
        }

        if (script == SCRIPT_GAP_IN_A) {
            b_index--;
        }
        else if (script == SCRIPT_GAP_IN_B) {
            a_index--;
        }
        else {
            a_index--;
            b_index--;
        }
        GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
    }

    sfree(edit_start_offset);
    sfree(edit_script);
    return best_score;
}

/** Low level function to perform gapped extension in one direction with 
 * or without traceback.
 * @param A The query sequence [in]
 * @param B The subject sequence [in]
 * @param M Maximal extension length in query [in]
 * @param N Maximal extension length in subject [in]
 * @param a_offset Resulting starting offset in query [out]
 * @param b_offset Resulting starting offset in subject [out]
 * @param score_only Only find the score, without saving traceback [in]
 * @param edit_block Structure to hold generated traceback [out]
 * @param gap_align Structure holding various information and allocated 
 *        memory for the gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param query_offset The starting offset in query [in]
 * @param reversed Has the sequence been reversed? Used for psi-blast [in]
 * @param reverse_sequence Do reverse the sequence [in]
 * @return The best alignment score found.
 */
Int4 
Blast_SemiGappedAlign(Uint1* A, Uint1* B, Int4 M, Int4 N,
   Int4* a_offset, Int4* b_offset, Boolean score_only, 
   GapPrelimEditBlock *edit_block, BlastGapAlignStruct* gap_align, 
   const BlastScoringParameters* score_params, 
   Int4 query_offset, Boolean reversed, Boolean reverse_sequence)
{
    Int4 i;                     /* sequence pointers and indices */
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    Uint1* b_ptr;
  
    BlastGapDP* score_array;

    Int4 gap_open;              /* alignment penalty variables */
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;
  
    Int4** matrix = NULL;       /* pointers to the score matrix */
    Int4** pssm = NULL;
    Int4* matrix_row = NULL;
  
    Int4 score;                 /* score tracking variables */
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;
    Int4 best_score;
    Int4 num_extra_cells;
  
    if (!score_only) {
        return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align, 
                      score_params, query_offset, reversed, reverse_sequence);
    }
    
    /* do initialization and sanity-checking */

    matrix = gap_align->sbp->matrix->data;
    if (gap_align->positionBased) {
        pssm = gap_align->sbp->psi_matrix->pssm->data;
    }
    *a_offset = 0;
    *b_offset = 0;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    x_dropoff = gap_align->gap_x_dropoff;
  
    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;
  
    if(N <= 0 || M <= 0) 
        return 0;
  
    /* Allocate and fill in the auxiliary bookeeping structures.
       Since A and B could be very large, maintain a window
       of auxiliary structures only large enough to contain to current
       set of DP computations. The initial window size is determined
       by the number of cells needed to fail the x-dropoff test */

    if (gap_extend > 0)
        num_extra_cells = x_dropoff / gap_extend + 3;
    else
        num_extra_cells = N + 3;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                                      2 * gap_align->dp_mem_alloc);
        sfree(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                                  sizeof(BlastGapDP));
    }

    score_array = gap_align->dp_mem;
    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
  
    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff) 
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score -= gap_extend;
    }
  
    /* The inner loop below examines letters of B from 
       index 'first_b_index' to 'b_size' */

    b_size = i;
    best_score = 0;
    first_b_index = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;
  
    for (a_index = 1; a_index <= M; a_index++) {
        /* pick out the row of the score matrix 
           appropriate for A[a_index] */

        if (!(gap_align->positionBased)) {
            if(reverse_sequence)
                matrix_row = matrix[ A[ M - a_index ] ];
            else
                matrix_row = matrix[ A[ a_index ] ];
        }
        else {
            if(reversed || reverse_sequence)
                matrix_row = pssm[M - a_index];
            else 
                matrix_row = pssm[a_index + query_offset];
        }

        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        /* initialize running-score variables */
        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;

        for (b_index = first_b_index; b_index < b_size; b_index++) {

            b_ptr += b_increment;
            score_gap_col = score_array[b_index].best_gap;
            next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

            if (score < score_gap_col)
                score = score_gap_col;

            if (score < score_gap_row)
                score = score_gap_row;

            if (best_score - score > x_dropoff) {

                /* the current best score failed the X-dropoff
                   criterion. Note that this does not stop the
                   inner loop, only forces future iterations to
                   skip this column of B. 

                   Also, if the very first letter of B that was
                   tested failed the X dropoff criterion, make
                   sure future inner loops start one letter to 
                   the right */

                if (b_index == first_b_index)
                    first_b_index++;
                else
                    score_array[b_index].best = MININT;
            }
            else {
                last_b_index = b_index;
                if (score > best_score) {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }

                /* If starting a gap at this position will improve
                   the best row, or column, score, update them to 
                   reflect that. */

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                score_array[b_index].best_gap = MAX(score - gap_open_extend,
                                                    score_gap_col);
                score_gap_row = MAX(score - gap_open_extend, score_gap_row);
                score_array[b_index].best = score;
            }

            score = next_score;
        }

        /* Finish aligning if the best scores for all positions
           of B will fail the X-dropoff test, i.e. the inner loop 
           bounds have converged to each other */

        if (first_b_index == b_size)
            break;

        /* enlarge the window for score data if necessary */

        if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                                          2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                                               gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }

        if (last_b_index < b_size - 1) {
            /* This row failed the X-dropoff test earlier than
               the last row did; just shorten the loop bounds
               before doing the next row */

            b_size = last_b_index + 1;
        }
        else {
            /* The inner loop finished without failing the X-dropoff
               test; initialize extra bookkeeping structures until
               the X dropoff test fails or we run out of letters in B. 
               The next inner loop will have larger bounds */

            while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size++;
            }
        }

        if (b_size <= N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }
    
    return best_score;
}

/** Editing script operations for out-of-frame traceback */
enum {
    SCRIPT_AHEAD_ONE_FRAME       = 1, /**< Shift 1 frame in sequence A
                                         (gap 2 nucleotides) */ 
    SCRIPT_AHEAD_TWO_FRAMES      = 2, /**< Shift 2 frames in sequence A
                                         (gap 1 nucleotide) */
    SCRIPT_NEXT_IN_FRAME         = 3, /**< Shift to next base (substitution) */
    SCRIPT_NEXT_PLUS_ONE_FRAME   = 4, /**< Shift to next base plus 1 frame
                                         (gap 1 nucleotide in sequence B) */
    SCRIPT_NEXT_PLUS_TWO_FRAMES  = 5, /**< Shift to next base plus 2 frames
                                         (gap 2 nucleotides in sequence B) */
    SCRIPT_OOF_OPEN_GAP          = 0x08 /**< Opening a gap */
};

/** Low level function to perform gapped extension with out-of-frame
 * gapping with traceback.  
 * @param A The query sequence [in]
 * @param B The subject sequence [in]
 * @param M Maximal extension length in query [in]
 * @param N Maximal extension length in subject [in]
 * @param a_offset Resulting starting offset in query [out]
 * @param b_offset Resulting starting offset in subject [out]
 * @param edit_block Structure to hold generated traceback [out]
 * @param gap_align Structure holding various information and allocated 
 *        memory for the gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param query_offset The starting offset in query [in]
 * @param reversed Has the sequence been reversed? Used for psi-blast [in]
 * @return The best alignment score found.
 */
static Int4 
s_OutOfFrameAlignWithTraceback(Uint1* A, Uint1* B, Int4 M, Int4 N,
   Int4* a_offset, Int4* b_offset, GapPrelimEditBlock *edit_block,
   BlastGapAlignStruct* gap_align, const BlastScoringParameters* score_params,
   Int4 query_offset, Boolean reversed)
{
    Int4 i, increment;          /* sequence pointers and indices */
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index;
  
    BlastGapDP* score_array;

    Int4 gap_open;              /* alignment penalty variables */
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 shift_penalty;
    Int4 x_dropoff;
  
    Int4** matrix = NULL;       /* pointers to the score matrix */
    Int4** pssm = NULL;
    Int4* matrix_row = NULL;
  
    Int4 score;                 /* score tracking variables */
    Int4 score_row1; 
    Int4 score_row2; 
    Int4 score_row3;
    Int4 score_gap_col; 
    Int4 score_col1; 
    Int4 score_col2; 
    Int4 score_col3;
    Int4 score_other_frame1; 
    Int4 score_other_frame2; 
    Int4 best_score;
  
    GapStateArrayStruct* state_struct;
    Uint1** edit_script;
    Uint1* edit_script_row;
    Int4 *edit_start_offset;
    Int4 edit_script_num_rows;
    Int4 orig_b_index;
    Int1 script, next_script;
    Int4 num_extra_cells;

    /* do initialization and sanity-checking */

    matrix = gap_align->sbp->matrix->data;
    if (gap_align->positionBased) {
        pssm = gap_align->sbp->psi_matrix->pssm->data;
    }
    *a_offset = 0;
    *b_offset = -2;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    shift_penalty = score_params->shift_pen;
    x_dropoff = gap_align->gap_x_dropoff;
  
    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;
  
    if(N <= 0 || M <= 0) 
        return 0;
  
    /* Initialize traceback information. edit_script[] is
       a 2-D array which is filled in row by row as the 
       dynamic programming takes place */

    s_GapPurgeState(gap_align->state_struct);

    /* Conceptually, traceback requires a byte array of size
       MxN. To save on memory, each row of the array is allocated
       separately. edit_script[i] points to the storage reserved
       for row i, and edit_start_offset[i] gives the offset into
       the B sequence corresponding to element 0 of edit_script[i] 
       
       Also make the number of edit script rows grow dynamically */

    edit_script_num_rows = 100;
    edit_script = (Uint1**) malloc(sizeof(Uint1*) * edit_script_num_rows);
    edit_start_offset = (Int4*) malloc(sizeof(Int4) * edit_script_num_rows);

    /* allocate storage for the first row of the traceback
       array. Because row elements correspond to gaps in A,
       the alignment can only go at most x_dropoff/gap_extend 
       positions, in all three frames, before failing the 
       X dropoff criterion */

    if (gap_extend > 0)
        num_extra_cells = CODON_LENGTH * (x_dropoff / gap_extend + 5);
    else
        num_extra_cells = N + 5;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                                      2 * gap_align->dp_mem_alloc);
        sfree(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                                  sizeof(BlastGapDP));
    }

    state_struct = s_GapGetState(&gap_align->state_struct, num_extra_cells);

    edit_script[0] = state_struct->state_array;
    edit_start_offset[0] = 0;
    edit_script_row = state_struct->state_array;

    score_array = gap_align->dp_mem;
    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
  
    for (i = 3; i <= N + 2; i += 3) {
        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        edit_script_row[i] = SCRIPT_GAP_IN_B;

        score_array[i-1].best = MININT;
        score_array[i-1].best_gap = MININT;
        edit_script_row[i-1] = SCRIPT_GAP_IN_B;

        score_array[i-2].best = MININT;
        score_array[i-2].best_gap = MININT;
        edit_script_row[i-2] = SCRIPT_GAP_IN_B;

        if (score < -x_dropoff) 
            break;
        score -= gap_extend;
    }
  
    /* The inner loop below examines letters of B from 
       index 'first_b_index' to 'b_size' */

    score_array[i].best = MININT;
    score_array[i].best_gap = MININT;
    b_size = i - 2;
    state_struct->used = b_size + 1;

    best_score = 0;
    first_b_index = 0;
    if (reversed) {
        increment = -1;
    }
    else {
        /* Allow for a backwards frame shift */
        B -= 2;
        increment = 1;
    }
  
    for (a_index = 1; a_index <= M; a_index++) {

        /* Set up the next row of the edit script; this involves
           allocating memory for the row, then pointing to it.
           It is not necessary to allocate space for offsets less
           than first_b_index (since the inner loop will never 
           look at them). 

           It is unknown at this point how far to the right the 
           current traceback row will extend; all that's known for
           sure is that the previous row fails the X-dropoff test
           after b_size cells, and that the current row can go at
           most num_extra_cells beyond that before failing the test */

        if (gap_extend > 0)
            state_struct = s_GapGetState(&gap_align->state_struct, 
                           b_size - first_b_index + num_extra_cells);
        else
            state_struct = s_GapGetState(&gap_align->state_struct, 
                                        N + 5 - first_b_index);

        if (a_index == edit_script_num_rows) {
            edit_script_num_rows = edit_script_num_rows * 2;
            edit_script = (Uint1 **)realloc(edit_script, 
                                            edit_script_num_rows *
                                            sizeof(Uint1 *));
            edit_start_offset = (Int4 *)realloc(edit_start_offset, 
                                                edit_script_num_rows *
                                                sizeof(Int4));
        }

        edit_script[a_index] = state_struct->state_array + 
                                state_struct->used + 1;
        edit_start_offset[a_index] = first_b_index;

        /* the inner loop assumes the current traceback
           row begins at offset zero of B */

        edit_script_row = edit_script[a_index] - first_b_index;
        orig_b_index = first_b_index;

        if (!(gap_align->positionBased)) {
            matrix_row = matrix[ A[ a_index * increment ] ];
        }
        else {
            if(reversed)
                matrix_row = pssm[M - a_index];
            else 
                matrix_row = pssm[a_index + query_offset];
        }

        score = MININT;
        score_row1 = MININT; 
        score_row2 = MININT; 
        score_row3 = MININT;
        score_gap_col = MININT; 
        score_col1 = MININT; 
        score_col2 = MININT; 
        score_col3 = MININT;
        score_other_frame1 = MININT; 
        score_other_frame2 = MININT; 
        last_b_index = first_b_index;
        b_index = first_b_index;

        /* The inner loop is identical to that of OOF_SEMI_G_ALIGN,
           except that traceback operations are sprinkled throughout. */

        while (b_index < b_size) {

            /* FRAME 0 */

            score = MAX(score_other_frame1, score_other_frame2) - shift_penalty;
            score = MAX(score, score_col1);
            if (score == score_col1) {
                script = SCRIPT_NEXT_IN_FRAME;
            }
            else if (score + shift_penalty == score_other_frame1) {
                if (score_other_frame1 == score_col2)
                    script = SCRIPT_AHEAD_TWO_FRAMES;
                else
                    script = SCRIPT_NEXT_PLUS_TWO_FRAMES;
            }
            else {
                if (score_other_frame2 == score_col3)
                    script = SCRIPT_AHEAD_ONE_FRAME;
                else
                    script = SCRIPT_NEXT_PLUS_ONE_FRAME;
            }
            score += matrix_row[ B[ b_index * increment ] ];

            score_other_frame1 = MAX(score_col1, score_array[b_index].best);
            score_col1 = score_array[b_index].best;
            score_gap_col = score_array[b_index].best_gap;

            if (score < MAX(score_gap_col, score_row1)) {
                if (score_gap_col > score_row1) {
                    score = score_gap_col;
                    script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_A;
                }
                else {
                    score = score_row1;
                    script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;
                }

                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    score_array[b_index].best_gap = score_gap_col - gap_extend;
                    score_row1 -= gap_extend;
                }
            }
            else {
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    if (score > best_score) {
                        best_score = score;
                        *a_offset = a_index;
                        *b_offset = b_index;
                    }

                    score -= gap_open_extend;
                    score_row1 -= gap_extend;
                    if (score > score_row1)
                        score_row1 = score;
                    else
                        script |= SCRIPT_EXTEND_GAP_B;

                    score_gap_col -= gap_extend;
                    if (score < score_gap_col) {
                        score_array[b_index].best_gap = score_gap_col;
                        script |= SCRIPT_EXTEND_GAP_A;
                    }
                    else {
                        score_array[b_index].best_gap = score;
                    }
                }
            }

            edit_script_row[b_index++] = script;
            if (b_index >= b_size) {
                score = score_row1;
                score_row1 = score_row2;
                score_row2 = score_row3;
                score_row3 = score;
                break;
            }

            /* FRAME 1 */

            score = MAX(score_other_frame1, score_other_frame2) - shift_penalty;
            score = MAX(score, score_col2);
            if (score == score_col2) {
                script = SCRIPT_NEXT_IN_FRAME;
            }
            else if (score + shift_penalty == score_other_frame1) {
                if (score_other_frame1 == score_col1)
                    script = SCRIPT_AHEAD_ONE_FRAME;
                else
                    script = SCRIPT_NEXT_PLUS_ONE_FRAME;
            }
            else {
                if (score_other_frame2 == score_col3)
                    script = SCRIPT_AHEAD_TWO_FRAMES;
                else
                    script = SCRIPT_NEXT_PLUS_TWO_FRAMES;
            }
            score += matrix_row[ B[ b_index * increment ] ];
            score_other_frame2 = MAX(score_col2, score_array[b_index].best);
            score_col2 = score_array[b_index].best;
            score_gap_col = score_array[b_index].best_gap;

            if (score < MAX(score_gap_col, score_row2)) {
                score = MAX(score_gap_col, score_row2);
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    if (score == score_gap_col)
                        script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_A;
                    else 
                        script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;

                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    score_array[b_index].best_gap = score_gap_col - gap_extend;
                    score_row2 -= gap_extend;
                }
            }
            else {
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    if (score > best_score) {
                        best_score = score;
                        *a_offset = a_index;
                        *b_offset = b_index;
                    }
                    score -= gap_open_extend;
                    score_row2 -= gap_extend;
                    if (score > score_row2)
                        score_row2 = score;
                    else
                        script |= SCRIPT_EXTEND_GAP_B;

                    score_gap_col -= gap_extend;
                    if (score < score_gap_col) {
                        score_array[b_index].best_gap = score_gap_col;
                        script |= SCRIPT_EXTEND_GAP_A;
                    }
                    else {
                        score_array[b_index].best_gap = score;
                    }
                }
            }

            edit_script_row[b_index++] = script;
            if (b_index >= b_size) {
                score = score_row2;
                score_row2 = score_row1;
                score_row1 = score_row3;
                score_row3 = score;
                break;
            }

            /* FRAME 2 */

            score = MAX(score_other_frame1, score_other_frame2) - shift_penalty;
            score = MAX(score, score_col3);
            if (score == score_col3) {
                script = SCRIPT_NEXT_IN_FRAME;
            }
            else if (score + shift_penalty == score_other_frame1) {
                if (score_other_frame1 == score_col1)
                    script = SCRIPT_AHEAD_TWO_FRAMES;
                else
                    script = SCRIPT_NEXT_PLUS_TWO_FRAMES;
            }
            else {
                if (score_other_frame2 == score_col2)
                    script = SCRIPT_AHEAD_ONE_FRAME;
                else
                    script = SCRIPT_NEXT_PLUS_ONE_FRAME;
            }
            score += matrix_row[ B[ b_index * increment ] ];
            score_other_frame1 = score_other_frame2;
            score_other_frame2 = MAX(score_col3, score_array[b_index].best);
            score_col3 = score_array[b_index].best;
            score_gap_col = score_array[b_index].best_gap;

            if (score < MAX(score_gap_col, score_row3)) {
                score = MAX(score_gap_col, score_row3);
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    if (score == score_gap_col)
                        script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_A;
                    else 
                        script = SCRIPT_OOF_OPEN_GAP | SCRIPT_GAP_IN_B;

                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    score_array[b_index].best_gap = score_gap_col - gap_extend;
                    score_row3 -= gap_extend;
                }
            }
            else {
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    if (score > best_score) {
                        best_score = score;
                        *a_offset = a_index;
                        *b_offset = b_index;
                    }
                    score -= gap_open_extend;
                    score_row3 -= gap_extend;
                    if (score > score_row3)
                        score_row3 = score;
                    else
                        script |= SCRIPT_EXTEND_GAP_B;

                    score_gap_col -= gap_extend;
                    if (score < score_gap_col) {
                        score_array[b_index].best_gap = score_gap_col;
                        script |= SCRIPT_EXTEND_GAP_A;
                    }
                    else {
                        score_array[b_index].best_gap = score;
                    }
                }
            }
            edit_script_row[b_index++] = script;
        }
  
        /* Finish aligning if the best scores for all positions
           of B will fail the X-dropoff test, i.e. the inner loop 
           bounds have converged to each other */

        if (first_b_index == b_size)
            break;

        /* Enlarge the window for score data if necessary */

        if (last_b_index + num_extra_cells + 5 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                                          2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                                               gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }

        if (last_b_index < b_size - 1) {
            /* This row failed the X-dropoff test earlier than
               the last row did; just shorten the loop bounds
               before doing the next row */

            b_size = last_b_index + 1;
        }
        else {
            /* The inner loop finished without failing the X-dropoff
               test; initialize extra bookkeeping structures until
               the X dropoff test fails or we run out of letters in B. 
               The next inner loop will have larger bounds. 
             
               Keep initializing extra structures until the X dropoff
               test fails in all frames for this row */

            score = MAX(score_row1, score_row2);
            score = MAX(score, score_row3);
            while (score >= (best_score - x_dropoff) && b_size < N + 1) {

                score_array[b_size].best = score_row1;
                score_array[b_size].best_gap = score_row1 - gap_open_extend;
                score_row1 -= gap_extend;
                edit_script_row[b_size] = SCRIPT_OOF_OPEN_GAP | 
                                          SCRIPT_GAP_IN_B;

                score_array[b_size+1].best = score_row2;
                score_array[b_size+1].best_gap = score_row2 - gap_open_extend;
                score_row2 -= gap_extend;
                edit_script_row[b_size+1] = SCRIPT_OOF_OPEN_GAP | 
                                            SCRIPT_GAP_IN_B;

                score_array[b_size+2].best = score_row3;
                score_array[b_size+2].best_gap = score_row3 - gap_open_extend;
                score_row3 -= gap_extend;
                edit_script_row[b_size+2] = SCRIPT_OOF_OPEN_GAP | 
                                            SCRIPT_GAP_IN_B;
                b_size += 3;
                score -= gap_extend;
            }
        }

        /* update the memory allocator to reflect the exact number
           of traceback cells this row needed */

        b_size = MIN(b_size, N + 1);
        state_struct->used += MAX(b_index, b_size) - orig_b_index + 1;

        /* chop off the best score in each frame */

        last_b_index = MIN(b_size + 4, N + 3);
        while (b_size < last_b_index) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }

    a_index = *a_offset;
    b_index = *b_offset;
    script = SCRIPT_AHEAD_ONE_FRAME;

    while (a_index > 0 || b_index > 0) {
        /* Retrieve the next action to perform. Rows of
           the traceback array do not necessarily start
           at offset zero of B, so a correction is needed
           to point to the correct position */

        next_script = 
            edit_script[a_index][b_index - edit_start_offset[a_index]];

        switch (script) {
        case SCRIPT_GAP_IN_A:
            script = next_script & SCRIPT_OP_MASK;
            if (next_script & (SCRIPT_OOF_OPEN_GAP | SCRIPT_EXTEND_GAP_A))
                script = SCRIPT_GAP_IN_A;
            break;
        case SCRIPT_GAP_IN_B:
            script = next_script & SCRIPT_OP_MASK;
            if (next_script & (SCRIPT_OOF_OPEN_GAP | SCRIPT_EXTEND_GAP_B))
                script = SCRIPT_GAP_IN_B;
            break;
        default:
            script = next_script & SCRIPT_OP_MASK;
            break;
        }

        if (script == SCRIPT_GAP_IN_B) {
            b_index -= 3;
        }
        else {
            b_index -= script;
            a_index--;
        }

        GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
    }

    sfree(edit_start_offset);
    sfree(edit_script);

    if (!reversed)
        *b_offset -= 2;
    return best_score;
}

/** Low level function to perform gapped extension with out-of-frame
 * gapping with or without traceback.
 * @param A The query sequence [in]
 * @param B The subject sequence [in]
 * @param M Maximal extension length in query [in]
 * @param N Maximal extension length in subject [in]
 * @param a_offset Resulting starting offset in query [out]
 * @param b_offset Resulting starting offset in subject [out]
 * @param score_only Only find the score, without saving traceback [in]
 * @param edit_block Structure to hold generated traceback [out]
 * @param gap_align Structure holding various information and allocated 
 *        memory for the gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param query_offset The starting offset in query [in]
 * @param reversed Has the sequence been reversed? Used for psi-blast [in]
 * @return The best alignment score found.
 */
static Int4 
s_OutOfFrameGappedAlign(Uint1* A, Uint1* B, Int4 M, Int4 N,
   Int4* a_offset, Int4* b_offset, Boolean score_only,
   GapPrelimEditBlock *edit_block, BlastGapAlignStruct* gap_align, 
   const BlastScoringParameters* score_params,
   Int4 query_offset, Boolean reversed)
{
    Int4 i, increment;          /* sequence pointers and indices */
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index;
  
    Int4 gap_open;              /* alignment penalty variables */
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 shift_penalty;
    Int4 x_dropoff;
  
    BlastGapDP* score_array;
    Int4 num_extra_cells;

    Int4** matrix = NULL;       /* pointers to the score matrix */
    Int4** pssm = NULL;
    Int4* matrix_row = NULL;
  
    Int4 score;                 /* score tracking variables */
    Int4 score_row1; 
    Int4 score_row2; 
    Int4 score_row3;
    Int4 score_gap_col; 
    Int4 score_col1; 
    Int4 score_col2; 
    Int4 score_col3;
    Int4 score_other_frame1; 
    Int4 score_other_frame2; 
    Int4 best_score;
  
    if (!score_only) {
        return s_OutOfFrameAlignWithTraceback(A, B, M, N, a_offset, b_offset, 
                                              edit_block, gap_align, 
                                              score_params, query_offset, 
                                              reversed);
    }
    
    /* do initialization and sanity-checking */

    matrix = gap_align->sbp->matrix->data;
    if (gap_align->positionBased) {
        pssm = gap_align->sbp->psi_matrix->pssm->data;
    }
    *a_offset = 0;
    *b_offset = -2;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    shift_penalty = score_params->shift_pen;
    x_dropoff = gap_align->gap_x_dropoff;
  
    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;
  
    if(N <= 0 || M <= 0) 
        return 0;
  
    /* Allocate and fill in the auxiliary bookeeping structures.
       Since A and B could be very large, maintain a window
       of auxiliary structures only large enough to contain to current
       set of DP computations. The initial window size is determined
       by the number of cells needed to fail the x-dropoff test */

    if (gap_extend > 0)
        num_extra_cells = CODON_LENGTH * (x_dropoff / gap_extend + 5);
    else
        num_extra_cells = N + 5;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                                      2 * gap_align->dp_mem_alloc);
        sfree(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                                  sizeof(BlastGapDP));
    }

    score_array = gap_align->dp_mem;
    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
  
    for (i = 3; i <= N + 2; i += 3) {
        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score_array[i-1].best = MININT;
        score_array[i-1].best_gap = MININT;
        score_array[i-2].best = MININT;
        score_array[i-2].best_gap = MININT;
        if (score < -x_dropoff) 
            break;
        score -= gap_extend;
    }
  
    /* The inner loop below examines letters of B from 
       index 'first_b_index' to 'b_size' */

    b_size = i - 2;
    score_array[i].best = MININT;
    score_array[i].best_gap = MININT;

    best_score = 0;
    first_b_index = 0;
    if (reversed) {
        increment = -1;
    }
    else {
        /* Allow for a backwards frame shift */
        B -= 2;
        increment = 1;
    }
  
    for (a_index = 1; a_index <= M; a_index++) {

        /* pick out the row of the score matrix 
           appropriate for A[a_index] */

        if (!(gap_align->positionBased)) {
            matrix_row = matrix[ A[ a_index * increment ] ];
        }
        else {
            if(reversed)
                matrix_row = pssm[M - a_index];
            else 
                matrix_row = pssm[a_index + query_offset];
        }

        /* initialize running-score variables */
        score = MININT;
        score_row1 = MININT; 
        score_row2 = MININT; 
        score_row3 = MININT;
        score_gap_col = MININT; 
        score_col1 = MININT; 
        score_col2 = MININT; 
        score_col3 = MININT;
        score_other_frame1 = MININT; 
        score_other_frame2 = MININT; 
        last_b_index = first_b_index;
        b_index = first_b_index;

        while (b_index < b_size) {

            /* FRAME 0 */

            /* Pick the best score among all frames */
            score = MAX(score_other_frame1, score_other_frame2) - shift_penalty;
            score = MAX(score, score_col1) + 
                                matrix_row[ B[ b_index * increment ] ];
            score_other_frame1 = MAX(score_col1, score_array[b_index].best);
            score_col1 = score_array[b_index].best;
            score_gap_col = score_array[b_index].best_gap;

            /* Use the row and column scores if they improve
               the score overall */

            if (score < MAX(score_gap_col, score_row1)) {
                score = MAX(score_gap_col, score_row1);
                if (best_score - score > x_dropoff) {

                   /* the current best score failed the X-dropoff
                      criterion. Note that this does not stop the
                      inner loop, only forces future iterations to
                      skip this column of B. 
   
                      Also, if the very first letter of B that was
                      tested failed the X dropoff criterion, make
                      sure future inner loops start one letter to 
                      the right */

                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    /* update the row and column running scores */
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    score_array[b_index].best_gap = score_gap_col - gap_extend;
                    score_row1 -= gap_extend;
                }
            }
            else {
                if (best_score - score > x_dropoff) {

                   /* the current best score failed the X-dropoff
                      criterion. */

                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    /* The current best score exceeds the
                       row and column scores, and thus may
                       improve on the current optimal score */

                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    if (score > best_score) {
                        best_score = score;
                        *a_offset = a_index;
                        *b_offset = b_index;
                    }

                    /* compute the best scores that include gaps
                       or gap extensions */

                    score -= gap_open_extend;
                    score_row1 -= gap_extend;
                    score_row1 = MAX(score, score_row1);
                    score_array[b_index].best_gap = MAX(score, 
                                                  score_gap_col - gap_extend);
                }
            }

            /* If this was the last letter of B checked, rotate
               the row scores so that code beyond the inner loop
               works correctly */

            if (++b_index >= b_size) {
                score = score_row1;
                score_row1 = score_row2;
                score_row2 = score_row3;
                score_row3 = score;
                break;
            }

            /* FRAME 1 */

            /* This code, and that for Frame 2, are essentially the
               same as the preceeding code. The only real difference
               is the updating of the other_frame best scores */

            score = MAX(score_other_frame1, score_other_frame2) - shift_penalty;
            score = MAX(score, score_col2) + 
                                matrix_row[ B[ b_index * increment ] ];
            score_other_frame2 = MAX(score_col2, score_array[b_index].best);
            score_col2 = score_array[b_index].best;
            score_gap_col = score_array[b_index].best_gap;

            if (score < MAX(score_gap_col, score_row2)) {
                score = MAX(score_gap_col, score_row2);
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    score_array[b_index].best_gap = score_gap_col - gap_extend;
                    score_row2 -= gap_extend;
                }
            }
            else {
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    if (score > best_score) {
                        best_score = score;
                        *a_offset = a_index;
                        *b_offset = b_index;
                    }
                    score -= gap_open_extend;
                    score_row2 -= gap_extend;
                    score_row2 = MAX(score, score_row2);
                    score_array[b_index].best_gap = MAX(score, 
                                                  score_gap_col - gap_extend);
                }
            }

            if (++b_index >= b_size) {
                score = score_row2;
                score_row2 = score_row1;
                score_row1 = score_row3;
                score_row3 = score;
                break;
            }

            /* FRAME 2 */

            score = MAX(score_other_frame1, score_other_frame2) - shift_penalty;
            score = MAX(score, score_col3) + 
                                matrix_row[ B[ b_index * increment ] ];
            score_other_frame1 = score_other_frame2;
            score_other_frame2 = MAX(score_col3, score_array[b_index].best);
            score_col3 = score_array[b_index].best;
            score_gap_col = score_array[b_index].best_gap;

            if (score < MAX(score_gap_col, score_row3)) {
                score = MAX(score_gap_col, score_row3);
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    score_array[b_index].best_gap = score_gap_col - gap_extend;
                    score_row3 -= gap_extend;
                }
            }
            else {
                if (best_score - score > x_dropoff) {
                    if (first_b_index == b_index) 
                        first_b_index = b_index + 1;
                    else
                        score_array[b_index].best = MININT;
                }
                else {
                    last_b_index = b_index;
                    score_array[b_index].best = score;
                    if (score > best_score) {
                        best_score = score;
                        *a_offset = a_index;
                        *b_offset = b_index;
                    }
                    score -= gap_open_extend;
                    score_row3 -= gap_extend;
                    score_row3 = MAX(score, score_row3);
                    score_array[b_index].best_gap = MAX(score, 
                                                  score_gap_col - gap_extend);
                }
            }
            b_index++;
        }
  
        /* Finish aligning if the best scores for all positions
           of B will fail the X-dropoff test, i.e. the inner loop 
           bounds have converged to each other */

        if (first_b_index == b_size)
            break;

        /* Enlarge the window for score data, if necessary */

        if (b_size + num_extra_cells + 5 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(b_size + num_extra_cells + 100,
                                          2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                                               gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }

        if (last_b_index < b_size - 1) {
            /* This row failed the X-dropoff test earlier than
               the last row did; just shorten the loop bounds
               before doing the next row */

            b_size = last_b_index + 1;
        }
        else {
            /* The inner loop finished without failing the X-dropoff
               test; initialize extra bookkeeping structures until
               the X dropoff test fails or we run out of letters in B. 
               The next inner loop will have larger bounds. 
             
               Keep initializing extra structures until the X dropoff
               test fails in all frames for this row */

            score = MAX(score_row1, score_row2);
            score = MAX(score, score_row3);
            while (score >= (best_score - x_dropoff) && b_size < N + 1) {

                score_array[b_size].best = score_row1;
                score_array[b_size].best_gap = score_row1 - gap_open_extend;
                score_row1 -= gap_extend;

                score_array[b_size+1].best = score_row2;
                score_array[b_size+1].best_gap = score_row2 - gap_open_extend;
                score_row2 -= gap_extend;

                score_array[b_size+2].best = score_row3;
                score_array[b_size+2].best_gap = score_row3 - gap_open_extend;
                score_row3 -= gap_extend;

                b_size += 3;
                score -= gap_extend;
            }
        }

        /* chop off the best score in each frame */
        b_size = MIN(b_size, N + 1);
        last_b_index = MIN(b_size + 4, N + 3);
        while (b_size < last_b_index) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }

    if (!reversed) {
        /* The sequence was shifted, so length should be adjusted as well */
        *b_offset -= 2;
    }
    return best_score;
}

/** Find the HSP offsets relative to the individual query sequence instead of
 * the concatenated sequence.
 * @param query Query sequence block [in]
 * @param query_info Query information structure, including context offsets 
 *                   array [in]
 * @param init_hsp Initial HSP [in] [out]
 * @param query_out Optional: query sequence block with modified sequence 
 *                  pointer and sequence length [out]
 * @param init_hsp_out Optional: Pointer to initial HSP structure to hold
 *                     adjusted offsets; the input init_hsp will be modified if
 *                     this parameter is NULL [out]
 * @param context_out Which context this HSP belongs to? [out]
 */
static void 
s_GetRelativeCoordinates(const BLAST_SequenceBlk* query, 
                       const BlastQueryInfo* query_info, 
                       BlastInitHSP* init_hsp, BLAST_SequenceBlk* query_out,
                       BlastInitHSP* init_hsp_out, Int4* context_out)
{
   Int4 context = 0;
   Int4 query_start = 0, query_length = 0;

   context = BSearchContextInfo(init_hsp->offsets.qs_offsets.q_off, query_info);

   if (query && query->oof_sequence) {
       /* Out-of-frame blastx case: all frames of the same parity are mixed
          together in a special sequence. */
       Int4 query_end_pt = 0;
       Int4 mixed_frame_context = context - context % CODON_LENGTH;
       
       query_start = query_info->contexts[mixed_frame_context].query_offset;
       query_end_pt =
           query_info->contexts[mixed_frame_context+CODON_LENGTH-1].query_offset +
           query_info->contexts[mixed_frame_context+CODON_LENGTH-1].query_length;
       query_length =
           query_end_pt - query_start;
   } else {
       query_start  = query_info->contexts[context].query_offset;
       query_length = query_info->contexts[context].query_length;
   }
   
   if (query && query_out) {
      if (query->oof_sequence) {
         query_out->sequence = NULL;
         query_out->oof_sequence = query->oof_sequence + query_start;
      } else {
         query_out->sequence = query->sequence + query_start;
         query_out->oof_sequence = NULL;
      }
      query_out->length = query_length;
   }

   if (init_hsp_out) {
      init_hsp_out->offsets.qs_offsets.q_off = 
          init_hsp->offsets.qs_offsets.q_off - query_start;
      init_hsp_out->offsets.qs_offsets.s_off = 
          init_hsp->offsets.qs_offsets.s_off;
      if (init_hsp->ungapped_data) {
         if (!init_hsp_out->ungapped_data) {
            init_hsp_out->ungapped_data = (BlastUngappedData*) 
               BlastMemDup(init_hsp->ungapped_data, sizeof(BlastUngappedData));
         } else {
            memcpy(init_hsp_out->ungapped_data, init_hsp->ungapped_data,
                   sizeof(BlastUngappedData));
         }
         init_hsp_out->ungapped_data->q_start = 
            init_hsp->ungapped_data->q_start - query_start;
      }
   } else {
      init_hsp->offsets.qs_offsets.q_off -= query_start;
      if (init_hsp->ungapped_data)
         init_hsp->ungapped_data->q_start -= query_start;
   }

   *context_out = context;
}


/** Convert the initial list of traceback actions from a non-OOF
 *  gapped alignment into a blast edit script. Note that this routine
 *  assumes the input edit blocks have not been reversed or rearranged
 *  by calling code
 *  @param rev_prelim_tback Traceback from extension to the left [in]
 *  @param fwd_prelim_tback Traceback from extension to the right [in]
 *  @return Pointer to the resulting edit script, or NULL if there
 *          are no traceback actions specified
 */
GapEditScript*
Blast_PrelimEditBlockToGapEditScript (GapPrelimEditBlock* rev_prelim_tback,
                                      GapPrelimEditBlock* fwd_prelim_tback)
{
   Boolean merge_ops = FALSE;
   GapEditScript* esp;
   GapPrelimEditScript *op;
   Int4 i;
   Int4 index=0;
   Int4 size = 0;

   if (rev_prelim_tback == NULL || fwd_prelim_tback == NULL)
      return NULL;

   /* The fwd_prelim_tback script will get reversed here as the traceback started from the highest scoring point
     and worked backwards. The rev_prelim_tback script does NOT get reversed.  Since it was reversed when the 
     traceback was produced it's already "forward" */

   if (fwd_prelim_tback->num_ops > 0 && rev_prelim_tback->num_ops > 0 &&
       fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].op_type == 
         rev_prelim_tback->edit_ops[(rev_prelim_tback->num_ops)-1].op_type)
     merge_ops = TRUE;
 
   size = fwd_prelim_tback->num_ops+rev_prelim_tback->num_ops;
   if (merge_ops)
     size--;

   esp = GapEditScriptNew(size);

   index = 0;
   for (i=0; i < rev_prelim_tback->num_ops; i++) {
      op = rev_prelim_tback->edit_ops + i;
      esp->op_type[index] = op->op_type;
      esp->num[index] = op->num;
      index++;
   }

   if (fwd_prelim_tback->num_ops == 0)
      return esp;

   if (merge_ops)
       esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;

   /* If we merge, then we skip the first one. */
   if (merge_ops)
      i = fwd_prelim_tback->num_ops - 2;
   else
      i = fwd_prelim_tback->num_ops - 1;

   for (; i >= 0; i--) {
      op = fwd_prelim_tback->edit_ops + i;
      esp->op_type[index] = op->op_type;
      esp->num[index] = op->num;
      index++;
   }

   return esp;
}

/** Fills the BlastGapAlignStruct structure with the results of a gapped 
 * extension.
 * @param gap_align the initialized gapped alignment structure [in] [out]
 * @param q_start The starting offset in query [in]
 * @param s_start The starting offset in subject [in]
 * @param q_end The ending offset in query [in]
 * @param s_end The ending offset in subject [in]
 * @param score The alignment score [in]
 * @param esp The edit script containing the traceback information [in]
 */
static Int2
s_BlastGapAlignStructFill(BlastGapAlignStruct* gap_align, Int4 q_start, 
   Int4 s_start, Int4 q_end, Int4 s_end, Int4 score, GapEditScript* esp)
{
   gap_align->query_start = q_start;
   gap_align->query_stop = q_end;
   gap_align->subject_start = s_start;
   gap_align->subject_stop = s_end;
   gap_align->score = score;

   if (esp)
      gap_align->edit_script = esp;

   return 0;
}

Int2 
BLAST_GreedyGappedAlignment(Uint1* query, Uint1* subject, 
   Int4 query_length, Int4 subject_length, BlastGapAlignStruct* gap_align,
   const BlastScoringParameters* score_params, 
   Int4 q_off, Int4 s_off, Boolean compressed_subject, Boolean do_traceback)
{
   Uint1* q;
   Uint1* s;
   Int4 score;
   Int4 X;
   Int4 q_avail, s_avail;
   Int4 q_ext_l, q_ext_r, s_ext_l, s_ext_r;
   GapPrelimEditBlock *fwd_prelim_tback = NULL;
   GapPrelimEditBlock *rev_prelim_tback = NULL;
   Uint1 rem;
   GapEditScript* esp = NULL;
   
   q_avail = query_length - q_off;
   s_avail = subject_length - s_off;
   
   q = query + q_off;
   if (!compressed_subject) {
      s = subject + s_off;
      rem = 4; /* Special value to indicate that sequence is already
                  uncompressed */
   } else {
      s = subject + s_off/4;
      rem = s_off % 4;
   }

   X = gap_align->gap_x_dropoff;
   
   if (do_traceback) {
      fwd_prelim_tback = gap_align->fwd_prelim_tback;
      rev_prelim_tback = gap_align->rev_prelim_tback;
      GapPrelimEditBlockReset(fwd_prelim_tback);
      GapPrelimEditBlockReset(rev_prelim_tback);
   }
   
   /* extend to the right */
   score = BLAST_AffineGreedyAlign(q, q_avail, s, s_avail, FALSE, X,
              score_params->reward, -score_params->penalty, 
              score_params->gap_open, score_params->gap_extend,
              &q_ext_r, &s_ext_r, gap_align->greedy_align_mem, 
              fwd_prelim_tback, rem);

   if (compressed_subject)
      rem = 0;

   /* extend to the left */
   score += BLAST_AffineGreedyAlign(query, q_off, 
               subject, s_off, TRUE, X, 
               score_params->reward, -score_params->penalty, 
               score_params->gap_open, score_params->gap_extend, 
               &q_ext_l, &s_ext_l, gap_align->greedy_align_mem, 
               rev_prelim_tback, rem);

   /* In basic case the greedy algorithm returns number of 
      differences, hence we need to convert it to score */
   if (score_params->gap_open==0 && score_params->gap_extend==0) {
      score = 
         (q_ext_r + s_ext_r + q_ext_l + s_ext_l)*score_params->reward/2 - 
         score*(score_params->reward - score_params->penalty);
   } else if (score_params->reward % 2 == 1) {
      score /= 2;
   }

   if (do_traceback) {
      esp = Blast_PrelimEditBlockToGapEditScript(rev_prelim_tback, 
                                             fwd_prelim_tback);
   }
   
   s_BlastGapAlignStructFill(gap_align, q_off-q_ext_l, s_off-s_ext_l, 
                             q_off+q_ext_r, s_off+s_ext_r, score, esp);
   return 0;
}

/** Performs dynamic programming style gapped extension, given an initial 
 * HSP (offset pair), two sequence blocks and scoring and extension options.
 * Note: traceback is not done in this function.
 * @param query_blk The query sequence [in]
 * @param subject_blk The subject sequence [in]
 * @param gap_align The auxiliary structure for gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param init_hsp The initial HSP that needs to be extended [in]
*/
static Int2 
s_BlastDynProgNtGappedAlignment(BLAST_SequenceBlk* query_blk, 
   BLAST_SequenceBlk* subject_blk, BlastGapAlignStruct* gap_align, 
   const BlastScoringParameters* score_params, BlastInitHSP* init_hsp)
{
   Boolean found_start, found_end;
   Int4 q_length=0, s_length=0, score_right, score_left, 
      private_q_start, private_s_start;
   Uint1 offset_adjustment;
   Uint1* query,* subject;
   
   found_start = FALSE;
   found_end = FALSE;

   query = query_blk->sequence;
   subject = subject_blk->sequence;
   score_left = 0;
   /* If subject offset is not at the start of a full byte, 
      s_BlastAlignPackedNucl won't work, so shift the alignment start
      to the left */
   offset_adjustment = 
       COMPRESSION_RATIO - 
       (init_hsp->offsets.qs_offsets.s_off % COMPRESSION_RATIO);
   q_length = init_hsp->offsets.qs_offsets.q_off + offset_adjustment;
   s_length = init_hsp->offsets.qs_offsets.s_off + offset_adjustment;
   if (q_length != 0 && s_length != 0) {
      found_start = TRUE;
      score_left = s_BlastAlignPackedNucl(query, subject, q_length, s_length, 
                      &private_q_start, &private_s_start, gap_align, 
                      score_params, TRUE);
      if (score_left < 0) 
         return -1;
      gap_align->query_start = q_length - private_q_start;
      gap_align->subject_start = s_length - private_s_start;
   }

   score_right = 0;
   if (q_length < query_blk->length && 
       s_length < subject_blk->length)
   {
      found_end = TRUE;
      score_right = s_BlastAlignPackedNucl(query+q_length-1, 
         subject+(s_length+3)/COMPRESSION_RATIO - 1, 
         query_blk->length-q_length, 
         subject_blk->length-s_length, &(gap_align->query_stop),
         &(gap_align->subject_stop), gap_align, score_params, FALSE);
      if (score_right < 0) 
         return -1;
      gap_align->query_stop += q_length;
      gap_align->subject_stop += s_length;
   }

   if (found_start == FALSE) {
      /* Start never found */
      gap_align->query_start = init_hsp->offsets.qs_offsets.q_off;
      gap_align->subject_start = init_hsp->offsets.qs_offsets.s_off;
   }

   if (found_end == FALSE) {
      gap_align->query_stop = init_hsp->offsets.qs_offsets.q_off - 1;
      gap_align->subject_stop = init_hsp->offsets.qs_offsets.s_off - 1;
   }

   gap_align->score = score_right+score_left;

   return 0;
}

/** Aligns two nucleotide sequences, one (A) should be packed in the
 * same way as the BLAST databases, the other (B) should contain one
 * basepair/byte. Traceback is not done in this function.
 * @param B The query sequence [in]
 * @param A The subject sequence [in]
 * @param N Maximal extension length in query [in]
 * @param M Maximal extension length in subject [in]
 * @param b_offset Resulting starting offset in query [out]
 * @param a_offset Resulting starting offset in subject [out]
 * @param gap_align The auxiliary structure for gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param reverse_sequence Reverse the sequence.
 * @return The best alignment score found.
*/
static Int4 
s_BlastAlignPackedNucl(Uint1* B, Uint1* A, Int4 N, Int4 M, 
	Int4* b_offset, Int4* a_offset, 
        BlastGapAlignStruct* gap_align,
        const BlastScoringParameters* score_params, 
        Boolean reverse_sequence)
{ 
    Int4 i;                     /* sequence pointers and indices */
    Int4 a_index, a_base_pair;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    Uint1* b_ptr;
  
    BlastGapDP* score_array;
    Int4 num_extra_cells;

    Int4 gap_open;              /* alignment penalty variables */
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;
  
    Int4* *matrix;              /* pointers to the score matrix */
    Int4* matrix_row;
  
    Int4 score;                 /* score tracking variables */
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;
    Int4 best_score;
  
    /* do initialization and sanity-checking */

    matrix = gap_align->sbp->matrix->data;
    *a_offset = 0;
    *b_offset = 0;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    x_dropoff = gap_align->gap_x_dropoff;
  
    /* The computations below assume that alignment will
       never be declined */

    ASSERT(score_params->decline_align >= INT2_MAX);

    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;
  
    if(N <= 0 || M <= 0) 
        return 0;
  
    /* Allocate and fill in the auxiliary bookeeping structures.
       Since A and B could be very large, maintain a window
       of auxiliary structures only large enough to contain the current
       set of DP computations. The initial window size is determined
       by the number of cells needed to fail the x-dropoff test */

    if (gap_extend > 0)
        num_extra_cells = x_dropoff / gap_extend + 3;
    else
        num_extra_cells = N + 3;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                                      2 * gap_align->dp_mem_alloc);
        sfree(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                                  sizeof(BlastGapDP));
    }

    score_array = gap_align->dp_mem;
    score = -gap_open_extend;
    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
  
    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff) 
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend; 
        score -= gap_extend;
    }
  
    /* The inner loop below examines letters of B from 
       index 'first_b_index' to 'b_size' */

    b_size = i;
    best_score = 0;
    first_b_index = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;
  
    for (a_index = 1; a_index <= M; a_index++) {

        /* pick out the row of the score matrix 
           appropriate for A[a_index] */

        if(reverse_sequence) {
            a_base_pair = NCBI2NA_UNPACK_BASE(A[(M-a_index)/4], 
                                               ((a_index-1)%4));
            matrix_row = matrix[a_base_pair];
        } 
        else {
            a_base_pair = NCBI2NA_UNPACK_BASE(A[1+((a_index-1)/4)], 
                                               (3-((a_index-1)%4)));
            matrix_row = matrix[a_base_pair];
        }

        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        /* initialize running-score variables */
        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;

        for (b_index = first_b_index; b_index < b_size; b_index++) {

            b_ptr += b_increment;
            score_gap_col = score_array[b_index].best_gap;
            next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
            
            if (score < score_gap_col)
                score = score_gap_col;

            if (score < score_gap_row)
                score = score_gap_row;

            if (best_score - score > x_dropoff) {

                /* the current best score failed the X-dropoff
                   criterion. Note that this does not stop the
                   inner loop, only forces future iterations to
                   skip this column of B. 

                   Also, if the very first letter of B that was
                   tested failed the X dropoff criterion, make
                   sure future inner loops start one letter to 
                   the right */

                if (b_index == first_b_index)
                    first_b_index++;
                else
                    score_array[b_index].best = MININT;
            }
            else {
                last_b_index = b_index;
                if (score > best_score) {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }

                /* If starting a gap at this position will improve
                   the best row or column score, update them to 
                   reflect that. */

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                score_array[b_index].best_gap = MAX(score - gap_open_extend,
                                                    score_gap_col);
                score_gap_row = MAX(score - gap_open_extend, score_gap_row);

                score_array[b_index].best = score;
            }

            score = next_score;
        }
  
        /* Finish aligning if the best scores for all positions
           of B will fail the X-dropoff test, i.e. the inner loop 
           bounds have converged to each other */

        if (first_b_index == b_size)
            break;

        if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                                          2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                                               gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }

        if (last_b_index < b_size - 1) {
            /* This row failed the X-dropoff test earlier than
               the last row did; just shorten the loop bounds
               before doing the next row */

            b_size = last_b_index + 1;
        }
        else {
            /* The inner loop finished without failing the X-dropoff
               test; initialize extra bookkeeping structures until
               the X dropoff test fails or we run out of letters in B. 
               The next inner loop will have larger bounds */

            while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size++;
            }
        }

        if (b_size <= N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }
    
    return best_score;
}

Int4 
BlastGetStartForGappedAlignment (Uint1* query, Uint1* subject,
   const BlastScoreBlk* sbp, Uint4 q_start, Uint4 q_length, 
   Uint4 s_start, Uint4 s_length)
{
    Int4 index1, max_offset, score, max_score, hsp_end;
    Uint1* query_var,* subject_var;
    Boolean positionBased = (sbp->psi_matrix != NULL);
    
    if (q_length <= HSP_MAX_WINDOW) {
        max_offset = q_start + q_length/2;
        return max_offset;
    }

    hsp_end = q_start + HSP_MAX_WINDOW;
    query_var = query + q_start;
    subject_var = subject + s_start;
    score=0;
    for (index1=q_start; index1<hsp_end; index1++) {
        if (!(positionBased))
            score += sbp->matrix->data[*query_var][*subject_var];
        else
            score += sbp->psi_matrix->pssm->data[index1][*subject_var];
        query_var++; subject_var++;
    }
    max_score = score;
    max_offset = hsp_end - 1;
    hsp_end = q_start + MIN(q_length, s_length);
    for (index1=q_start + HSP_MAX_WINDOW; index1<hsp_end; index1++) {
        if (!(positionBased)) {
            score -= sbp->matrix->data[*(query_var-HSP_MAX_WINDOW)][*(subject_var-HSP_MAX_WINDOW)];
            score += sbp->matrix->data[*query_var][*subject_var];
        } else {
            score -= sbp->psi_matrix->pssm->data[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
            score += sbp->psi_matrix->pssm->data[index1][*subject_var];
        }
        if (score > max_score) {
            max_score = score;
            max_offset = index1;
        }
        query_var++; subject_var++;
    }
    if (max_score > 0)
       max_offset -= HSP_MAX_WINDOW/2;
    else 
       max_offset = q_start;

    return max_offset;
}

Int2 BLAST_GetGappedScore (EBlastProgramType program_number, 
        BLAST_SequenceBlk* query, BlastQueryInfo* query_info, 
        BLAST_SequenceBlk* subject, 
        BlastGapAlignStruct* gap_align,
        const BlastScoringParameters* score_params,
        const BlastExtensionParameters* ext_params,
        const BlastHitSavingParameters* hit_params,
        BlastInitHitList* init_hitlist,
        BlastHSPList** hsp_list_ptr, BlastGappedStats* gapped_stats)

{
   Int4 index;
   BlastInitHSP* init_hsp = NULL;
   BlastInitHSP* init_hsp_array;
   Int4 q_start, s_start, q_end, s_end;
   Boolean is_prot;
   Boolean is_greedy;
   Int4 max_offset;
   Int2 status = 0;
   BlastHSPList* hsp_list = NULL;
   const BlastHitSavingOptions* hit_options = hit_params->options;
   BLAST_SequenceBlk query_tmp;
   Int4 context;
   BlastIntervalTree *tree;
   Int4 score;
   Int4 **rpsblast_pssms = NULL;   /* Pointer to concatenated PSSMs in
                                       RPS-BLAST database */
   const int kHspNumMax = BlastHspNumMax(TRUE, hit_options);

   if (!query || !subject || !gap_align || !score_params || !ext_params ||
       !hit_params || !init_hitlist || !hsp_list_ptr)
      return 1;

   if (init_hitlist->total == 0)
      return 0;

   is_prot = (program_number != eBlastTypeBlastn &&
              program_number != eBlastTypePhiBlastn);
   is_greedy = (ext_params->options->ePrelimGapExt != eDynProgExt);

   if (program_number == eBlastTypeRpsTblastn || 
       program_number == eBlastTypeRpsBlast) {
       rpsblast_pssms = gap_align->sbp->psi_matrix->pssm->data;
   }

   ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist));

   if (*hsp_list_ptr == NULL)
      *hsp_list_ptr = hsp_list = Blast_HSPListNew(kHspNumMax);
   else 
      hsp_list = *hsp_list_ptr;

   init_hsp_array = init_hitlist->init_hsp_array;

   /* Initialize the interval tree with the maximum possible
      query and subject offsets. For query sequences this is always
      query->length, and for subject sequences it is subject->length
      except for out-of-frame tblastn. In that case, subject->length
      is the length of one strand of the subject sequence, but the
      subject offsets of each ungapped alignment are wrt all six subject
      frames concatenated together. Finally, add 1 to the maximum, 
      to account for HSPs that include the end of a sequence */

   if (program_number == eBlastTypeTblastn &&
       score_params->options->is_ooframe) {
      tree = Blast_IntervalTreeInit(0, query->length+1,
                                    0, 2*(subject->length + CODON_LENGTH)+1);
   }
   else {
      tree = Blast_IntervalTreeInit(0, query->length+1,
                                    0, subject->length+1);
   }

   for (index=0; index<init_hitlist->total; index++)
   {
      BlastHSP tmp_hsp;
      init_hsp = &init_hsp_array[index];

      /* Now adjust the initial HSP's coordinates. */
      s_GetRelativeCoordinates(query, query_info, init_hsp, &query_tmp, 
                             NULL, &context);

      if (rpsblast_pssms)
         gap_align->sbp->psi_matrix->pssm->data = rpsblast_pssms + 
             query_info->contexts[context].query_offset;

      if (!init_hsp->ungapped_data) {
         q_start = q_end = init_hsp->offsets.qs_offsets.q_off;
         s_start = s_end = init_hsp->offsets.qs_offsets.s_off;
         score = INT4_MIN;
      } else {
         q_start = init_hsp->ungapped_data->q_start;
         q_end = q_start + init_hsp->ungapped_data->length;
         s_start = init_hsp->ungapped_data->s_start;
         s_end = s_start + init_hsp->ungapped_data->length;
         score = init_hsp->ungapped_data->score;
      }

      tmp_hsp.score = score;
      tmp_hsp.context = context;
      tmp_hsp.query.offset = q_start;
      tmp_hsp.query.end = q_end;
      tmp_hsp.query.frame = query_info->contexts[context].frame;
      tmp_hsp.subject.offset = s_start;
      tmp_hsp.subject.end = s_end;
      tmp_hsp.subject.frame = subject->frame;

      if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
                                    hit_options->min_diag_separation))
      {
         BlastHSP* new_hsp;

         if (gapped_stats) {
            ++gapped_stats->extensions;
         }
 
         if(is_prot && !score_params->options->is_ooframe) {
            max_offset = 
               BlastGetStartForGappedAlignment(query_tmp.sequence, 
                  subject->sequence, gap_align->sbp,
                  init_hsp->ungapped_data->q_start,
                  init_hsp->ungapped_data->length,
                  init_hsp->ungapped_data->s_start,
                  init_hsp->ungapped_data->length);
            init_hsp->offsets.qs_offsets.s_off += 
                max_offset - init_hsp->offsets.qs_offsets.q_off;
            init_hsp->offsets.qs_offsets.q_off = max_offset;
         }

         if (is_prot) {
            status =  s_BlastProtGappedAlignment(program_number, &query_tmp, 
                         subject, gap_align, score_params, init_hsp);
         } else if (is_greedy) {
            status = BLAST_GreedyGappedAlignment(
                         query_tmp.sequence, subject->sequence, 
                         query_tmp.length, subject->length, 
                         gap_align, score_params, 
                         init_hsp->offsets.qs_offsets.q_off, 
                         init_hsp->offsets.qs_offsets.s_off, 
                         (Boolean) TRUE, 
                         (Boolean) (ext_params->options->ePrelimGapExt == 
                                    eGreedyWithTracebackExt));
         } else {
            /*  Start the gapped alignment on the fourth character of the
             *  eight character word that seeded the alignment; the start
             *  is included in the leftward extension. */
            init_hsp->offsets.qs_offsets.s_off += 3;
            init_hsp->offsets.qs_offsets.q_off += 3;
            status = s_BlastDynProgNtGappedAlignment(&query_tmp, subject, 
                         gap_align, score_params, init_hsp);
         }

         if (status) {
             if (rpsblast_pssms) {
                 gap_align->sbp->psi_matrix->pssm->data = rpsblast_pssms;
             }
            return status;
         }
        
         if (gap_align->score >= hit_params->cutoff_score) {
             Int2 query_frame = 0;
             /* For mixed-frame search, the query frame is determined 
                from the offset, not only from context. */
             if (score_params->options->is_ooframe && 
                 program_number == eBlastTypeBlastx) {
                 query_frame = gap_align->query_start % CODON_LENGTH + 1;
                 if ((context % NUM_FRAMES) >= CODON_LENGTH)
                     query_frame = -query_frame;
             } else {
                 query_frame = query_info->contexts[context].frame;
             }

             Blast_HSPInit(gap_align->query_start, 
                           gap_align->query_stop, gap_align->subject_start, 
                           gap_align->subject_stop, 
                           init_hsp->offsets.qs_offsets.q_off, 
                           init_hsp->offsets.qs_offsets.s_off, context, 
                           query_frame, subject->frame, gap_align->score, 
                           &(gap_align->edit_script), &new_hsp);
             Blast_HSPListSaveHSP(hsp_list, new_hsp);
             BlastIntervalTreeAddHSP(new_hsp, tree, query_info, 
                                     eQueryAndSubject);
         }
         else {
            /* a greedy alignment may have traceback associated with it;
               free that traceback if the alignment will not be used */
            if (is_greedy) {
               gap_align->edit_script = GapEditScriptDelete(
                                             gap_align->edit_script);
            }
         }
      }
   }   

   tree = Blast_IntervalTreeFree(tree);
   if (rpsblast_pssms) {
       gap_align->sbp->psi_matrix->pssm->data = rpsblast_pssms;
   }

   *hsp_list_ptr = hsp_list;
   return status;
}

/** Out-of-frame gapped alignment wrapper function.
 * @param query Query sequence [in]
 * @param subject Subject sequence [in]
 * @param q_off Offset in query [in]
 * @param s_off Offset in subject [in]
 * @param private_q_start Extent of alignment in query [out]
 * @param private_s_start Extent of alignment in subject [out]
 * @param score_only Return score only, without traceback [in]
 * @param edit_block Structure to hold generated traceback [out]
 * @param gap_align Gapped alignment information and preallocated 
 *                  memory [in] [out]
 * @param score_params Scoring parameters [in]
 * @param psi_offset Starting position in PSI-BLAST matrix [in]
 * @param reversed Direction of the extension [in]
 * @param switch_seq Sequences need to be switched for blastx, 
 *                   but not for tblastn [in]
 */
static Int4 
s_OutOfFrameSemiGappedAlignWrap(Uint1* query, Uint1* subject, Int4 q_off, 
   Int4 s_off, Int4* private_q_start, Int4* private_s_start, 
   Boolean score_only, GapPrelimEditBlock *edit_block, 
   BlastGapAlignStruct* gap_align, 
   const BlastScoringParameters* score_params, Int4 psi_offset, 
   Boolean reversed, Boolean switch_seq)
{
   if (switch_seq) {
      return s_OutOfFrameGappedAlign(subject, query, s_off, q_off,
                private_s_start, private_q_start, score_only, edit_block,
                gap_align, score_params, psi_offset, reversed);
   } else {
      return s_OutOfFrameGappedAlign(query, subject, q_off, s_off,
                private_q_start, private_s_start, score_only, edit_block,
                gap_align, score_params, psi_offset, reversed);
   }
}

/** Maximal subject length after which the offsets are adjusted to a 
 * subsequence.
 */
#define MAX_SUBJECT_OFFSET 90000
/** Approximate upper bound on a number of gaps in an HSP, needed to determine
 * the length of the subject subsequence to be retrieved for alignment with
 * traceback. 
 */
#define MAX_TOTAL_GAPS 3000

void 
AdjustSubjectRange(Int4* subject_offset_ptr, Int4* subject_length_ptr, 
		   Int4 query_offset, Int4 query_length, Int4* start_shift)
{
   Int4 s_offset;
   Int4 subject_length = *subject_length_ptr;
   Int4 max_extension_left, max_extension_right;
   
   /* If subject sequence is not too long, leave everything as is */
   if (subject_length < MAX_SUBJECT_OFFSET) {
      *start_shift = 0;
      return;
   }

   s_offset = *subject_offset_ptr;
   /* Maximal extension length is the remaining length in the query, plus 
      an estimate of a maximal total number of gaps. */
   max_extension_left = query_offset + MAX_TOTAL_GAPS;
   max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;

   if (s_offset <= max_extension_left) {
      *start_shift = 0;
   } else {
      *start_shift = s_offset - max_extension_left;
      *subject_offset_ptr = max_extension_left;
   }

   *subject_length_ptr = 
      MIN(subject_length, s_offset + max_extension_right) - *start_shift;
}

/** Performs gapped extension for protein sequences, given two
 * sequence blocks, scoring and extension options, and an initial HSP 
 * with information from the previously performed ungapped extension
 * @param program BLAST program [in]
 * @param query_blk The query sequence block [in]
 * @param subject_blk The subject sequence block [in]
 * @param gap_align The auxiliary structure for gapped alignment [in]
 * @param score_params Parameters related to scoring [in]
 * @param init_hsp The initial HSP information [in]
 */
static Int2 
s_BlastProtGappedAlignment(EBlastProgramType program, 
   BLAST_SequenceBlk* query_blk, BLAST_SequenceBlk* subject_blk, 
   BlastGapAlignStruct* gap_align,
   const BlastScoringParameters* score_params, BlastInitHSP* init_hsp)
{
   Boolean found_start, found_end;
   Int4 q_length=0, s_length=0, score_right, score_left;
   Int4 private_q_start, private_s_start;
   Uint1* query=NULL,* subject=NULL;
   Boolean switch_seq = FALSE;
   Int4 query_length = query_blk->length;
   Int4 subject_length = subject_blk->length;
   Int4 subject_shift = 0;
   BlastScoringOptions *score_options = score_params->options;
    
   if (gap_align == NULL)
      return -1;
   
   if (score_options->is_ooframe) {
      q_length = init_hsp->offsets.qs_offsets.q_off;
      /* For negative subject frames, make subject offset relative to the part
         of the mixed-frame sequence corresponding to the reverse strand. */
      if (program == eBlastTypeTblastn && subject_blk->frame < 0)
          init_hsp->offsets.qs_offsets.s_off -= subject_length + 1;
      
      s_length = init_hsp->offsets.qs_offsets.s_off;

      if (program == eBlastTypeBlastx) {
         subject = subject_blk->sequence + s_length;
         query = query_blk->oof_sequence + CODON_LENGTH + q_length;
         query_length -= CODON_LENGTH - 1;
         switch_seq = TRUE;
      } else if (program == eBlastTypeTblastn) {
         subject = subject_blk->oof_sequence + CODON_LENGTH + s_length;
         query = query_blk->sequence + q_length;
         subject_length -= CODON_LENGTH - 1;
      }
   } else {
      q_length = init_hsp->offsets.qs_offsets.q_off + 1;
      s_length = init_hsp->offsets.qs_offsets.s_off + 1;
      query = query_blk->sequence;
      subject = subject_blk->sequence;
   }

   AdjustSubjectRange(&s_length, &subject_length, q_length, query_length, 
                      &subject_shift);

   found_start = FALSE;
   found_end = FALSE;
    
   /* Looking for "left" score */
   score_left = 0;
   if (q_length != 0 && s_length != 0) {
      found_start = TRUE;
      if(score_options->is_ooframe) {
         score_left = 
             s_OutOfFrameSemiGappedAlignWrap(query, subject, q_length, s_length,
                 &private_q_start, &private_s_start, TRUE, NULL, 
                 gap_align, score_params, q_length, TRUE, switch_seq);
      } else {
         score_left = 
             Blast_SemiGappedAlign(query, subject+subject_shift, q_length, s_length,
                                   &private_q_start, &private_s_start, TRUE, 
                                   NULL, gap_align, score_params, 
                                   init_hsp->offsets.qs_offsets.q_off, 
                                   FALSE, TRUE);
      }
        
      gap_align->query_start = q_length - private_q_start;
      gap_align->subject_start = s_length - private_s_start + subject_shift;
      
   }

   score_right = 0;
   if (q_length < query_length && s_length < subject_length) {
      found_end = TRUE;
      if(score_options->is_ooframe) {
         score_right = 
             s_OutOfFrameSemiGappedAlignWrap(query-1, subject-1, 
                 query_length-q_length+1, subject_length-s_length+1,
                 &(gap_align->query_stop), &(gap_align->subject_stop), 
                 TRUE, NULL, gap_align, score_params, q_length, FALSE, 
                 switch_seq);
         gap_align->query_stop += q_length;
         gap_align->subject_stop += s_length + subject_shift;
      } else {
         score_right = 
             Blast_SemiGappedAlign(query+init_hsp->offsets.qs_offsets.q_off, 
                                   subject+init_hsp->offsets.qs_offsets.s_off, 
                                   query_length-q_length, 
                                   subject_length-s_length, 
                                   &(gap_align->query_stop), 
                                   &(gap_align->subject_stop), 
                                   TRUE, NULL, gap_align, score_params, 
                                   init_hsp->offsets.qs_offsets.q_off, FALSE, 
                                   FALSE);
         /* Make end offsets point to the byte after the end of the 
            alignment */
         gap_align->query_stop += init_hsp->offsets.qs_offsets.q_off + 1;
         gap_align->subject_stop += init_hsp->offsets.qs_offsets.s_off + 1;
      }
   }
   
   if (found_start == FALSE) {  /* impossible for in-frame */
      gap_align->query_start = init_hsp->offsets.qs_offsets.q_off;
      gap_align->subject_start = init_hsp->offsets.qs_offsets.s_off;
   }
   if (found_end == FALSE) {    /* impossible for out-of-frame */
      gap_align->query_stop = q_length;
      gap_align->subject_stop = s_length;
   }
   
   gap_align->score = score_right+score_left;

   return 0;
}

/** Converts OOF traceback from a gapped alignment to a GapEditScript.
 * This function is for out-of-frame gapping.
 * @param rev_prelim_tback Traceback operations for left extension [in]
 * @param fwd_prelim_tback Traceback operations for right extension [in]
 * @param nucl_align_length Length of alignment on the translated
 *                          nucleotide sequence [in]
 * @param edit_script_ptr The resulting edit block [out]
 * @return Always zero
 */
static Int2
s_BlastOOFTracebackToGapEditScript(GapPrelimEditBlock *rev_prelim_tback,
                                   GapPrelimEditBlock *fwd_prelim_tback,
                                   Int4 nucl_align_length,
                                   GapEditScript** edit_script_ptr)
{
    GapEditScript* e_script;
    EGapAlignOpType last_op;
    Int4 last_num;
    GapPrelimEditBlock *tmp_prelim_tback;
    Int4 i, num_nuc;
    int extra_needed=0;
    
    /* prepend a substitution, since the input sequences were
       shifted prior to the OOF alignment */

    tmp_prelim_tback = GapPrelimEditBlockNew();
    last_op = eGapAlignSub;
    last_num = 1;

    /* Propagate the extra substitution through the edit script. */

    for (i = 0; i < rev_prelim_tback->num_ops; i++) {

        EGapAlignOpType next_op = rev_prelim_tback->edit_ops[i].op_type;
        Int4 next_num = rev_prelim_tback->edit_ops[i].num;

        if (next_op == last_op) {

            /* merge consecutive operations that are identical */

            last_num += next_num;
        }
        else if (next_op == eGapAlignIns || next_op == eGapAlignDel) {

            /* Propagating the initial extra substitution requires
               that edit operations which are not in-frame gaps must
               'flow through' an in-frame gap; the last non-gap
               operation before the gap will appear after the gap. 
               This could mean either one or two edit blocks
               getting created */

            if (last_num > 1) {
                GapPrelimEditBlockAdd(tmp_prelim_tback, last_op, last_num - 1);
            }
            GapPrelimEditBlockAdd(tmp_prelim_tback, next_op, next_num);
            last_num = 1;          /* last_op unchanged; one 'orphan' now */
        }
        else {

            /* write out the last edit operation and
               move forward in the script */

            GapPrelimEditBlockAdd(tmp_prelim_tback, last_op, last_num);
            last_op = next_op;
            last_num = next_num;
        }
    }

    /* Handle the crossing over from the left extension to
       the right extension. Begin by writing out the final 
       group of ops from the left extension, except for the 
       one op that must be merged */

    last_num--;
    if (last_num > 0) {
        GapPrelimEditBlockAdd(tmp_prelim_tback, last_op, last_num);
    }

    if (last_op != eGapAlignSub) {

        /* the last traceback letter from the left extension must
           be folded into the right extension. */

        for (i = fwd_prelim_tback->num_ops - 1; i >= 0; i--) {

            GapPrelimEditScript *new_script = fwd_prelim_tback->edit_ops + i;

            if (new_script->op_type == eGapAlignIns ||
                new_script->op_type == eGapAlignDel) {

                /* do not merge with in-frame gaps; just write out
                   the gaps and 'flow through' the last op */

                GapPrelimEditBlockAdd(tmp_prelim_tback, new_script->op_type,
                                      new_script->num);
            }
            else {

                /* merge with the first letter of new_script */

                last_op += new_script->op_type - eGapAlignSub;
                GapPrelimEditBlockAdd(tmp_prelim_tback, last_op, 1);

                /* borrow one letter from new_script to complete the
                   merge, and skip new_script if it only had one letter */
                new_script->num--;
                if (new_script->num == 0)
                    i--;
                break;
            }
        }
        fwd_prelim_tback->num_ops = i + 1;
    }

    /* form the final edit block */
    e_script = Blast_PrelimEditBlockToGapEditScript(tmp_prelim_tback, fwd_prelim_tback);
    GapPrelimEditBlockFree(tmp_prelim_tback);

    /* postprocess the edit script */
    num_nuc = 0;
    for (i=0; i<e_script->size; i++)
    {
        int total_actions=0;
        /* Count the number of nucleotides in the next 
           traceback operation and delete any traceback operations
           that would make the alignment too long. This check is
           not needed for in-frame alignment because the
           traceback is never changed after it is computed */

        last_op = e_script->op_type[i];

        if (last_op == eGapAlignIns)
            last_op = eGapAlignSub;

        total_actions = last_op * e_script->num[i];

        if (num_nuc + total_actions >= nucl_align_length) {
            e_script->num[i] = (nucl_align_length - num_nuc + 
                             last_op - 1) / last_op;
            break; /* We delete the rest of the script. */
        }
        else {
            num_nuc += total_actions;;
        }
    }
    e_script->size = MIN(i+1, e_script->size);  /* If we broke out early then we truncate the edit script. */

    extra_needed = 0;
    for (i=0; i<e_script->size; i++)
    {
        if (e_script->op_type[i] % 3 != 0 && e_script->num[i] > 1) {
           extra_needed += e_script->num[i] - 1;
        }
    }

    if (extra_needed)
    {
        GapEditScript* new_esp = GapEditScriptNew(extra_needed+e_script->size);
        int new_esp_i=0;
        for (i=0; i<e_script->size; i++)
        {
           /* Only one frame shift operation at a time is
           allowed in a single edit script element. */
           new_esp->num[new_esp_i] = e_script->num[i];
           new_esp->op_type[new_esp_i] = e_script->op_type[i];
           new_esp_i++;
           last_op = e_script->op_type[i];
           if (last_op % 3 != 0 && e_script->num[i] > 1) {
               Int4 num_ops = e_script->num[i];
               int esp_index=0;
               new_esp->num[new_esp_i-1] = 1;
               for (esp_index = 1; esp_index < num_ops; esp_index++) {
                   new_esp->num[new_esp_i] = 1;
                   new_esp->op_type[new_esp_i] = last_op;
                   new_esp_i++;
               }
           }
       }
       e_script = GapEditScriptDelete(e_script);
       e_script = new_esp;
    }
    *edit_script_ptr = e_script;

    /* finally, add one to the size of any block of substitutions
       that follows a frame shift op */
    last_op = e_script->op_type[0];
    for (i=1; i<e_script->size; i++)
    {
        if (e_script->op_type[i] == eGapAlignSub && (last_op % 3) != 0)
            (e_script->num[i])++;

        last_op = e_script->op_type[i];
    }

    return 0;
}

Int2 BLAST_GappedAlignmentWithTraceback(EBlastProgramType program, Uint1* query,
        Uint1* subject, BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length)
{
    Boolean found_start, found_end;
    Int4 score_right, score_left, private_q_length, private_s_length;
    Int4 q_length, s_length;
    Boolean is_ooframe = score_params->options->is_ooframe;
    Int2 status = 0;
    Boolean switch_seq = FALSE;
    GapPrelimEditBlock *fwd_prelim_tback;
    GapPrelimEditBlock *rev_prelim_tback;
    
    if (gap_align == NULL)
        return -1;
    
    fwd_prelim_tback = gap_align->fwd_prelim_tback;
    rev_prelim_tback = gap_align->rev_prelim_tback;
    GapPrelimEditBlockReset(fwd_prelim_tback);
    GapPrelimEditBlockReset(rev_prelim_tback);

    found_start = FALSE;
    found_end = FALSE;
    
    q_length = query_length;
    s_length = subject_length;
    if (is_ooframe) {
       /* The mixed frame sequence is shifted to the 3rd position, so its 
          maximal available length for extension is less by 2 than its 
          total length. */
       switch_seq = (program == eBlastTypeBlastx);
       if (switch_seq) {
          q_length -= CODON_LENGTH - 1;
       } else {
          s_length -= CODON_LENGTH - 1;
       }
    }

    score_left = 0;
    found_start = TRUE;
        
    if(is_ooframe) {
       /* NB: Left extension does not include starting point corresponding
          to the offset pair; the right extension does. */
       score_left =
          s_OutOfFrameSemiGappedAlignWrap(query+q_start, subject+s_start, 
             q_start, s_start, &private_q_length, &private_s_length, 
             FALSE, rev_prelim_tback, gap_align, score_params, q_start, 
             TRUE, switch_seq);
       gap_align->query_start = q_start - private_q_length;
       gap_align->subject_start = s_start - private_s_length;
    } else {        
       /* NB: The left extension includes the starting point 
          [q_start,s_start]; the right extension does not. */
       score_left = 
          Blast_SemiGappedAlign(query, subject, q_start+1, s_start+1,
             &private_q_length, &private_s_length, FALSE, rev_prelim_tback,
             gap_align, score_params, q_start, FALSE, TRUE);
       gap_align->query_start = q_start - private_q_length + 1;
       gap_align->subject_start = s_start - private_s_length + 1;
    }
    
    score_right = 0;

    if ((q_start < q_length) && (s_start < s_length)) {
       found_end = TRUE;
       if(is_ooframe) {
          score_right = 
             s_OutOfFrameSemiGappedAlignWrap(query+q_start-1, subject+s_start-1,
                q_length-q_start, s_length-s_start, 
                &private_q_length, &private_s_length, FALSE, fwd_prelim_tback,
                gap_align, score_params, q_start, FALSE, switch_seq);
        } else {
            score_right = 
               Blast_SemiGappedAlign(query+q_start, subject+s_start, 
                  q_length-q_start-1, s_length-s_start-1, &private_q_length, 
                  &private_s_length, FALSE, fwd_prelim_tback, gap_align, 
                  score_params, q_start, FALSE, FALSE);
        }

        gap_align->query_stop = q_start + private_q_length + 1;
        gap_align->subject_stop = s_start + private_s_length + 1;
    }
    
    if (found_start == FALSE) {	/* Start never found */
        gap_align->query_start = q_start;
        gap_align->subject_start = s_start;
    }
    
    if (found_end == FALSE) {
        gap_align->query_stop = q_start - 1;
        gap_align->subject_stop = s_start - 1;
    }

    if(is_ooframe) {
        Int4 nucl_align_length;
        if (program == eBlastTypeBlastx) {
            nucl_align_length = gap_align->query_stop - 
                                gap_align->query_start;
        }
        else {
            nucl_align_length = gap_align->subject_stop - 
                                gap_align->subject_start;
        }
        status = s_BlastOOFTracebackToGapEditScript(rev_prelim_tback, 
                       fwd_prelim_tback, nucl_align_length, 
                       &gap_align->edit_script);
    } else {
        gap_align->edit_script = 
            Blast_PrelimEditBlockToGapEditScript(rev_prelim_tback,
                                             fwd_prelim_tback);
    }

    gap_align->score = score_right + score_left;
    return status;
}

Int2 BLAST_GetUngappedHSPList(BlastInitHitList* init_hitlist,
                              BlastQueryInfo* query_info, 
                              BLAST_SequenceBlk* subject, 
                              const BlastHitSavingOptions* hit_options, 
                              BlastHSPList** hsp_list_ptr)
{
   BlastHSPList* hsp_list = NULL;
   Int4 index;
   BlastInitHSP* init_hsp;
   Int4 context;
   const int kHspNumMax = BlastHspNumMax(FALSE, hit_options);

   /* The BlastHSPList structure can be allocated and passed from outside */
   if (*hsp_list_ptr != NULL)
      hsp_list = *hsp_list_ptr;

   if (!init_hitlist) {
      if (!hsp_list)
         *hsp_list_ptr = NULL;
      else
         hsp_list->hspcnt = 0;
      return 0;
   }

   for (index = 0; index < init_hitlist->total; ++index) {
      BlastHSP* new_hsp;
      BlastUngappedData* ungapped_data=NULL;
      init_hsp = &init_hitlist->init_hsp_array[index];
      if (!init_hsp->ungapped_data) 
         continue;

      /* Adjust the initial HSP's coordinates in case of concatenated 
         multiple queries/strands/frames */
      s_GetRelativeCoordinates(NULL, query_info, init_hsp, NULL, 
                             NULL, &context);
      if (!hsp_list) {
         hsp_list = Blast_HSPListNew(kHspNumMax);
         *hsp_list_ptr = hsp_list;
      }
      ungapped_data = init_hsp->ungapped_data;
      Blast_HSPInit(ungapped_data->q_start, 
                    ungapped_data->length+ungapped_data->q_start,
                    ungapped_data->s_start, 
                    ungapped_data->length+ungapped_data->s_start,
                    init_hsp->offsets.qs_offsets.q_off, 
                    init_hsp->offsets.qs_offsets.s_off, 
                    context, query_info->contexts[context].frame, 
                    subject->frame, ungapped_data->score, NULL, &new_hsp);
      Blast_HSPListSaveHSP(hsp_list, new_hsp);
   }

   /* Sort the HSP array by score */
   Blast_HSPListSortByScore(hsp_list);

   return 0;
}

