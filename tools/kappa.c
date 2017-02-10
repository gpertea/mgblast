static char const rcsid[] = "$Id: kappa.c,v 6.83 2006/01/30 15:53:09 madden Exp $";

/* $Id: kappa.c,v 6.83 2006/01/30 15:53:09 madden Exp $ 
*   ==========================================================================
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
* ===========================================================================*/

/*****************************************************************************

File name: kappa.c

Authors: Alejandro Schaffer, Mike Gertz

Contents: Utilities for doing Smith-Waterman alignments and adjusting
    the scoring system for each match in blastpgp

Functions in this file call functions in the composition adjustment
module to perform adjustment and recalculation of alignments.  This
file is intended to be (thin) glue code that:

- adapts C toolkit data structures by the composition adjustment module; 
- calls functions in the composition adjustment module;
- converts the output of the composition adjustment module into the
  format needed by the C toolkit.

This file is very similar to the blast_kappa.c file, which is glue
code the performs the same function for the C++ toolkit.  The sole
external symbol from this file should be RedoAlignmentCore. Static
functions in this file have a 'Kappa_' prefix.  Please adhere to this
convention to avoid a name clash with functions in blast_kappa.c (the
name clash can matter in debuggers and search engines.)

 $Revision: 6.83 $

 Revision 6.75  2005/11/07 15:28:56  coulouri
 From Mike Gertz:

 Record the composition mode used for each aligment; formerly the
 composition mode was recorded for each subject sequence.  For blastp,
 the composition is only done once per subject sequence so the former
 behavior was correct, but it was incorrect for tblastn.

 Fixed a bug when the query frame is checked in an intermediate
 containment test.  For tblastn, this effectively disabled the test.
 This bug can effect tblastn alignments when SmithWaterman is not used
 in the rare occasion that omitting the intermediate containment test
 causes a different alignment to be found.

 Added a sort key to location_compare windows to sort first by query
 index.  For tblastn alignments when query concatination is used, not
 having this sort key may on rare occassion cause a different alignment
 to be reported.

 Cast a few Nlm_FloatHi values to Int4 to suppress compiler warnings.

 Revision 6.74  2005/09/29 17:40:18  coulouri
 from mike gertz:
     - Many changes to support query concatenation; including changing
       the API of RedoAlignementCore.
     - In RedoAlignmentCore, only call BlastLinkHsps if
       search->pbp->do_sum_stats is true.
     - In RedoAlignmentCore, in the non-Smith-Waterman case, compute a
       separate composition for each HSP, even if multiple HSPs are in
       the same window.

 Revision 6.73  2005/09/13 16:23:01  coulouri
 correct include paths for windows

 Revision 6.72  2005/09/13 14:55:06  coulouri
 From Mike Gertz: use libblastcompadj

 Revision 6.71  2005/09/08 19:48:36  coulouri
 From Mike Gertz:
 Changed the behavior of a SWheap so that space for only a small number
 (100) of hits is allocated initially, but the SWheap will grow to hold
 as many hits as are needed.

 Revision 6.70  2005/09/08 14:01:49  coulouri
 From Mike Gertz:
   - Fixed a bug in which the extent of the subject match was being
     calculated incorrectly setting forbidden ranges for a
     Smith-Waterman alignment.

 Revision 6.69  2005/08/31 20:33:25  coulouri
     New code to use the value of options->kappa_expect_value as the cutoff
     e-value within RedoAlignmentCore. The code saves and restores the
     value of search->pbp->cutoff_e.

 Revision 6.68  2005/08/30 18:20:46  coulouri
 From Mike Gertz:
     Changes to routines used in mode > 1 of compositional adjustment of
     scoring matrices:
       - Removed REscaleInitialScores as it is now redundant; rescaling
         is now done as part of the call s_ScatterScores.
       - Refactored Kappa_AdjustBXZ into the routines
         Kappa_SetNonstandardAaScores and Kappa_AverageScores.
       - Kappa_SetNonstandardAaScores sets scores for U, * and -,
         unlike its predecessor Kappa_AdjustBXZ.
       - Fixed a bug in the formula used to compute the score of (B,B)
         matches.
       - In s_RoundScoreMatrix, check whether the unrounded value is <
         BLAST_SCORE_MIN and if so set the rounded value to
         BLAST_SCORE_MIN; do not read the incoming value of the "matrix"
         parameter.

     Changes to mode 1 (traditional) composition-based statistics:
       - scaleMatrix renamed s_ScaleMatrx with slightly different
         parameter list.
       - s_ScaleMatrix does not scale substitution scores for the
         nonstandard X, U, - and * characters (but does scale B and Z).
         The previous version scaled scores for X and *.

     Made minor changes to Kappa_CompositionBasedStats and
     Kappa_AdjustSearch to call the new routines described above.

 Revision 6.67  2005/08/05 12:04:53  coulouri
 From Mike Gertz:
 - Move setting gap_align->translate2 to Kappa_RecordInitialSearch;
   for tblastn and some option settings it would not otherwise get set.
 - Remove a now redundant setting of gap_align->translate2 in the
   HitToGapAlign routine.
 - Fixes to comments.

 Revision 6.66  2005/08/02 14:40:29  coulouri
 From Mike Gertz:
 - Fixes to comments
 - Added enum constant eGapChar; renamed eStarChar to eStopChar.
 - Change the integer type of some variables to suppress warnings about
   assigning a wider type to a narrower type.
 - Made the routines BLbasicSmithWatermanScoreOnly and BLSmithWatermanFindStart
   static.
 - Renamed s_ScatterFreqRatios -> s_ScatterScores; renamed parameters.
 - In NewAlignmentUsingXdrop, use the translate2 field from gap_align
   to set the same field in the edit block.
 - Changed the Kappa_WindowInfo datatype to hold a list of HSPs in the
   window; added logic in several places, notably WindowsFromHSPs, to
   generate and maintain these lists.
 - Refactored Kappa_AdjustSearch.  Use a more sophisticated rule,
   implemented in the new Kappa_GetSubjectComposition routine, to
   determine the subject sequence composition for tblastn.
 - Removed unused parameters from several routines.
 - In RedoAlignmentCore, delete NRrecord to avoid a memory leak.

 Revision 6.65  2005/07/28 14:34:06  coulouri
 From Mike Gertz:
 - Made minor fixes to whitespace and formatting.
 - Fixed the comment to SWheap to accurately reflect recent changes in
   the way SWheapRecordCompare is used.
 - Changed some #define'd constants to enums or "static const"
   variables.  Renamed these constants appropriately.
 - Made the array alphaConvert and the routine scalePosMatrix static.
 - Made sure that gap_align->translate2 is set for tblastn.
 - Fixed a bug Kappa_SequenceGetTranslatedWindow: the wrong formula had
   been used to compute num_nucleotides.

 Revision 6.64  2005/07/26 13:07:20  coulouri
 - Changed #include "NRdefs.h" to #include <NRdefs.h> to be consistent
   with toolbox conventions.

 - Removed the unused kScoreMatrixRange constant

 - Extended comments for some routines.  Fixed spelling and other
   errors in comments.

 - Fixed white space and formatting in several locations.

 - Renamed functions:
   o permuteLetterProbs -> s_GatherLetterProbs
   o permuteFreqRatios  -> s_ScatterFreqRatios
   o adjustBXZ -> Kappa_AdjustBXZ
   o roundScoreMatrix -> s_RoundScoreMatrix

 - Eliminated the use of a temporary variable in REscaleInitialScores.

 - In sScatterFreqRatios, replaced the NRrecord parameter by a pointer
   to the frequency ratio matrix contained in the record.

 Revision 6.63  2005/07/14 20:19:45  coulouri
    - In Kappa_SequenceGetWindow, change all selenocysteine (U)
      residues in the subject sequence to X.
    - In scaleMatrix, do not compute the log of frequency ratios that are
      zero, set the matrix element to BLAST_SCORE_MIN instead.
    - Removed the startMatrix parameter to scaleMatrix, as it is no
      longer used.

 Revision 6.61  2005/06/08 19:31:31  papadopo
 From Michael Gertz:
 1. The use of the score for comparing collections of alignments
    was removed in March 2004; revert to the previous rule
 2. Added a new field "bestScore" to the SWheapRecord structure
 3. Various additional changes to support the use of score as a
    key in a SWheapRecord
 4. The comparison function is now used to determine whether new HSPs
    may be added to the heap.  Previously, the evalue only was used
 5. Removed a complex test that sometimes terminated computation of
    Smith-Waterman alignments for a given subject sequence. It is more
    consistent with all other modes of operation to omit this test
 6. Removed the "capacity" parameter to SWheapInitialize and used the
    heapThreshold parameter to set the capacity for the heap

 Revision 6.60  2005/05/18 21:27:33  papadopo
 make fillResidueProbability unconditional in Kappa_AdjustSearch

 Revision 6.59  2005/05/16 18:10:13  kans
 include Mode_condition.h to get prototype for chooseMode

 Revision 6.58  2005/05/16 17:46:48  papadopo
 From Mike Gertz/Alejandro Schaffer: Added support for compositional
 adjustment of score matrices, both conditionally and unconditionally.
 The adjustParameters argument to RedoAlignmentCore is no longer Boolean
 because there are now 4 allowed options:
       0 no adjustment
       1 composition-based statistics
       2 conditional compositional adjustment
       3 unconditional compositional adjustment

 Revision 6.57  2005/04/13 17:18:02  coulouri
 changes to WindowsFromHSPs for tblastn composition-based statistics

 Revision 6.56  2005/01/24 15:52:54  camacho
 doxygen fixes from Mike Gertz

 Revision 6.55  2005/01/20 16:28:59  camacho
 doxygen fixes

 Revision 6.54  2005/01/18 14:54:13  camacho
 Change in tie-breakers for score comparison, suggestion by Mike Gertz

 Revision 6.53  2004/11/24 15:42:33  coulouri
 do not dereference null pointer

 Revision 6.52  2004/11/23 21:29:02  camacho
 Changes from Mike Gertz:
 - If no alignments are found, I no longer create an empty HSP list and
   process it; I just move on to the next match.
 - For tblastn, culling by containment now occurs before link_hsps as
   it should. This can not have any effect on the current output, as
   tblastn + RedoAlignmentCore is not enabled.
 - I changed some Nlm_Mallocs, etc. to MemNew.
 - Other cosmetic changes, doxygen friendly comments

 Revision 6.51  2004/11/01 20:43:48  camacho
 Added error handling for missing PSSM

 Revision 6.50  2004/10/27 21:00:02  camacho
 Change the order of elements in Score* returned by GetScoreSetFromBlastHsp to
 be consistent with score ordering in other contexts of BLAST.

 Revision 6.49  2004/09/09 13:47:47  papadopo
 from Michael Gertz: remove unnecessary check for the presence of Selenocysteine in translated sequences

 Revision 6.48  2004/08/23 19:37:38  papadopo
 From Michael Gertz: fix memory leak in getStartFreqRatios

 Revision 6.47  2004/08/23 17:12:21  papadopo
 From Michael Gertz: Backported some changes to
 getStartFreqRatios and computeScaledStandardMatrix
 from algo/blast/core/blast_kappa.c code. Also a few
 straightforward changes to getStartFreqRatios,
 because the code has diverged slightly

 Revision 6.46  2004/07/28 17:19:54  kans
 HeapSort callbacks are int LIBCALLBACK for compilation on the PC

 Revision 6.45  2004/07/27 19:59:00  papadopo
 Changes by Michael Gertz to allow RedoAlignmentCore to be used
 for tblastn searches. The current version of the code uses two heuristics:
 1) a heuristic that skips doing re-alignment for an HSP if the old
    alignment of the HSP is contained in a higher-scoring HSP that
    has already been re-aligned.
 2) a heuristic that attempts to ensure that all alignments in the list
    of redone alignments has distinct ends.
 The original code also used these heuristics, but the code that
 implements the heuristics has been rewritten. As a result,
 RedoAlignmentCore will occasionally report better-scoring alignments.
 This is still a heuristic; which HSPs are reported still depends on
 the order in which they are re-aligned.

 Revision 6.44  2004/06/23 14:53:58  camacho
 Use renamed FreqRatios structure from posit.h

 Revision 6.43  2004/06/22 14:16:56  camacho
 Use SFreqRatios structure from algo/blast/core to obtain underlying matrices' frequency ratios.

 Revision 6.42  2004/06/18 15:50:25  madden
 For secondary alignments with Smith-Waterman on, use the E-value from the X-drop alignment computed by ALIGN that has more flexibility in which positions are allowed in the dynamic program.

 Add a sort of alignments for the same query-subject pair because the use of X-drop alignments occasionally reorders such alignments.

 Changes from Mike Gertz, submitted by Alejandro Schaffer.

 Revision 6.41  2004/06/14 21:11:05  papadopo
 From Michael Gertz:
 - Added several casts where casts occur in blast_kappa.c.  These casts
      should have no real effect; the log of blast_kappa.c indicates that
      they suppress compiler warnings.
 - Changed the type of one variable that holds a score from
      Nlm_FloatHi to Int4.
 - moved the definition Kappa_ForbiddenRanges and relevant
      routines earlier in the file.
 - fixed some comments.
 - made a few (~5) changes in whitespace.

 Revision 6.40  2004/06/03 16:10:50  dondosha
 Fix in Kappa_SearchParametersNew: allocate correct number of rows for matrices

 Revision 6.39  2004/03/31 18:12:13  papadopo
 Mike Gertz' refactoring of RedoAlignmentCore

 Revision 6.38  2004/02/06 19:25:16  camacho
 Alejandro Schaffer's corrections to RedoAlignmentCore pointed out by
 Mike Gertz:
 1. Corrected the rule for assigning newSW->isfirstAlignment
 2. Changed upper bound on assignment to forbiddenRanges
 3. Assigned a value to newScore earlier
 4. Eliminated use of skipThis
 5. Restored value of search->gapAlign at the end

 Revision 6.37  2004/01/27 20:31:52  madden
 remove extra setting of kbp_gap

 Revision 6.36  2004/01/06 17:48:44  dondosha
 Do not free Karlin block in RedoAlignmentCore, because its pointer is passed outside

 Revision 6.35  2003/12/01 19:15:27  madden
 Add one byte to filteredMatchingSequence to prevent ABR/ABW

 Revision 6.34  2003/10/22 20:37:19  madden
 Set kbp to rescaled values, use upper-case for SCALING_FACTOR define

 Revision 6.33  2003/10/02 19:59:34  kans
 BlastGetGapAlgnTbck needed FALSE instead of NULL in two parameters - Mac compiler complaint

 Revision 6.32  2003/10/02 19:31:24  madden
 In RedoAlignmentCore, call procedure BlastGetGapAlgnTbck instead of ALIGN
 to redo alignments; this allows the endpoints of the alignment to change.
 Because  BlastGetGapAlgnTbck returns a list of SeqAlign's while ALIGN
 passes back a single-alignment, the post-processing of the redone
 alignments is changed including the addition of the procedure
 concatenateListOfSeqaligns.

 Revision 6.30  2003/05/30 17:25:36  coulouri
 add rcsid

 Revision 6.29  2003/05/13 16:02:53  coulouri
 make ErrPostEx(SEV_FATAL, ...) exit with nonzero status

 Revision 6.28  2002/12/19 14:40:35  kans
 changed C++-style comment to C-style

 Revision 6.27  2002/12/10 22:58:42  bealer
 Keep mappings to sequences from readdb so that "num_ident" code does not
 segfault with multiple databases.

 Revision 6.26  2002/11/06 20:31:14  dondosha
 Added recalculation of the number of identities when computing seqalign from SWResults

 Revision 6.25  2002/10/16 13:33:58  madden
 Corrected initialization of initialUngappedLambda in RedoAlignmentCore

 Revision 6.24  2002/09/03 13:55:14  kans
 changed NULL to 0 for Mac compiler

 Revision 6.23  2002/08/29 15:47:49  camacho
 Changed RedoAlignmentCore to work without readdb facility

 Revision 6.22  2002/08/20 15:43:08  camacho
 Fixed memory leak in getStartFreqRatios

 Revision 6.21  2002/05/23 20:14:21  madden
 Correction for last checkin to cover SmithWaterman type alignment

 Revision 6.20  2001/12/28 18:02:33  dondosha
 Keep score and scoreThisAlign for each local alignment, so as to allow tie-breaking by score

 Revision 6.19  2001/07/26 12:52:25  madden
 Fix memory leaks

 Revision 6.18  2001/07/09 15:12:47  shavirin
 Functions BLbasicSmithWatermanScoreOnly() and BLSmithWatermanFindStart()
 used to calculate Smith-waterman alignments on low level become external.

 Revision 6.17  2001/05/25 19:40:46  vakatov
 Nested comment typo fixed

 Revision 6.16  2001/04/13 20:47:36  madden
 Eliminated use of PRO_K_MULTIPLIER in adjusting E-values Added allocateStartFreqs and freeStartFreqs and getStartFreqRatios to enable the score matrix scaling to work entirely with frequency ratios and round to integers only at the very end of the scaling calculation.

 Revision 6.15  2001/03/20 15:07:34  madden
 Fix from AS for (near) exact matches

 Revision 6.14  2001/01/04 13:48:44  madden
 Correction for 3-parameter gapping

 Revision 6.13  2000/12/05 19:31:33  madden
 Avoid duplicate insertion into heap when ((doThis == FALSE) && (currentState = SWPurging))

 Revision 6.12  2000/10/16 21:08:05  madden
 segResult takes BioseqPtr as argument, produced from readdb_get_bioseq

 Revision 6.11  2000/10/13 19:32:58  madden
 Fix for bug that caused crash

 Revision 6.10  2000/10/10 21:46:03  shavirin
 Added support for BLOSUM50, BLOSUM90, PAM250 with -t T

 Revision 6.9  2000/10/10 19:45:51  shavirin
 Checked for NULL hsp_array in the function RedoAlignmentCore().

 Revision 6.8  2000/08/18 21:28:24  madden
 support for BLOSUM62_20A and BLOSUM62_20B, prevent LambdaRatio from getting above 1.0

 Revision 6.7  2000/08/09 21:10:15  shavirin
 Added paramter discontinuous to the function newConvertSWalignsToSeqAligns()

 Revision 6.6  2000/08/08 21:45:04  shavirin
 Removed initialization of decline_aline parameter to INT2_MAX.

 Revision 6.5  2000/08/03 22:20:10  shavirin
 Added creation of the default posFreqs array if it is empty in
 RedoAlignmentCore(). Fixed some memory leaks.

 Revision 6.4  2000/08/02 18:00:34  shavirin
 Fixed memory leak in RedoAlignmentCore. Fixed a problem specific
 to BLOSUM45 and BLOSUM80 in the procedure scalePosMatrix().

 Revision 6.3  2000/07/26 19:34:13  kans
 removed unix-only headers

 Revision 6.2  2000/07/26 17:06:31  lewisg
 add LIBCALLs

 Revision 6.1  2000/07/25 17:40:05  shavirin
 Initial revision.


******************************************************************************/

#include<assert.h>
#include<string.h>

#include <ncbi.h>
#include <blast.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blastpri.h>
#include <txalign.h>
#include <simutil.h>
#include <posit.h>
#include <gapxdrop.h>
#include <fcntl.h>
#include <profiles.h>

#include <algo/blast/composition_adjustment/nlm_linear_algebra.h>
#include <algo/blast/composition_adjustment/matrix_frequency_data.h>
#include <algo/blast/composition_adjustment/composition_adjustment.h>
#include <algo/blast/composition_adjustment/compo_heap.h>
#include <algo/blast/composition_adjustment/smith_waterman.h>
#include <algo/blast/composition_adjustment/redo_alignment.h>


/**
 * Create a score set from the data in an HSP.
 *
 * @param hsp               the HSP of interest [in]
 * @param logK              Karlin-Altschul statistical parameter [in]
 * @param lambda            Karlin-Altschul statistical parameter [in]
 * @param scoreDivisor      the value by which reported scores are to be
 *                          divided [in]
 */
static ScorePtr
Kappa_GetScoreSetFromBlastHsp(
  const BLAST_HSPPtr hsp,
  Nlm_FloatHi lambda,
  Nlm_FloatHi logK,
  Nlm_FloatHi scoreDivisor)
{
  ScorePtr      score_set = NULL;       /* the new score set */
  Int4          score;          /* the score, scaled using scoreDivisor */
  Nlm_FloatHi   bit_score;      /* the integer-valued score, in bits */
  Nlm_FloatHi   evalue;         /* the e-value, with numbers too close to zero
                                   set to zero */

  score = Nlm_Nint(((Nlm_FloatHi) hsp->score) / scoreDivisor);
  if (score > 0) {
    MakeBlastScore(&score_set, "score", 0.0, score);
  }

  evalue = (hsp->evalue >= 1.0e-180) ? hsp->evalue : 0;
  if (hsp->num > 1) {
    MakeBlastScore(&score_set, "sum_n", 0.0, hsp->num);
    MakeBlastScore(&score_set, "sum_e", evalue, 0);
    if( hsp->ordering_method == BLAST_SMALL_GAPS) {
      MakeBlastScore(&score_set, "small_gap", 0.0, 1);
    }
  } else {
    MakeBlastScore(&score_set, "e_value", evalue, 0);
  }

  /* Compute the bit score using the newly computed scaled score. */
  bit_score = (score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;
  if (bit_score >= 0.) {
    MakeBlastScore(&score_set, "bit_score", bit_score, 0);
  }

  if (hsp->ordering_method > 3) {
    /* In new tblastn this means splice junction was found */
    MakeBlastScore(&score_set, "splice_junction", 0.0, 1);
  }

  if (hsp->num_ident > 0) {
    MakeBlastScore(&score_set, "num_ident", 0.0, hsp->num_ident);
  }

  MakeBlastScore(&score_set, "comp_adjustment_method",0.0,
                 hsp->comp_adjustment_method);
  return score_set;
}


/**
 * Converts a collection of alignments represented by a BLAST_HitList
 * into a collection of instances of SeqAlign. Returns the new
 * collection.
 *
 * @param hitlist           the hitlist to convert to SeqAligns [in]
 * @param subject_id        the identifier of the subject sequence [in]
 * @param query_id          the identifier of the query sequence [in]
 * @param queryOrigin       origin of this specific query in the
 *                          concatenated query [in]
 * @param queryLength       length of this specific query [in]
 * @param lambda            Karlin-Altschul statistical parameter [in]
 * @param logK              Karlin-Altschul statistical parameter [in]
 * @param scoreDivisor      value by which reported scores are to be
 *                          divided [in]
 * @param bestEvalue)       the best (smallest) e-value of all the alignments
 *                          in the hitlist [out]
 */
static SeqAlignPtr
Kappa_SeqAlignsFromHitlist(
  BLAST_HitListPtr hitlist,
  SeqIdPtr subject_id,
  SeqIdPtr query_id,
  Int4 queryOrigin,
  Int4 queryLength,
  Nlm_FloatHi lambda,
  Nlm_FloatHi logK,
  Nlm_FloatHi scoreDivisor)
{
  SeqAlignPtr aligns = NULL;  /* list of SeqAligns to be returned */
  Int4        hsp_index;

  for( hsp_index = hitlist->hspcnt - 1; hsp_index >= 0; hsp_index-- ) {
    /* iterate in reverse order over all HSPs in the hitlist */
    BLAST_HSPPtr hsp;           /* HSP for this iteration */
    SeqAlignPtr seqAlign;       /* the new SeqAlign */

    hsp      = hitlist->hsp_array[hsp_index];

    /* Shift the edit block to query local coordinates */
    hsp->gap_info->start1           -= queryOrigin;
    hsp->gap_info->length1          -= queryOrigin;
    hsp->gap_info->original_length1  = queryLength;

    seqAlign = GapXEditBlockToSeqAlign(hsp->gap_info, subject_id, query_id);
    seqAlign->score =
        Kappa_GetScoreSetFromBlastHsp(hsp, lambda, logK, scoreDivisor);

    seqAlign->next  = aligns;
    aligns          = seqAlign;
  }
  return aligns;
}


/**
 * Converts a list of objects of type BlastCompo_Alignment to an
 * new object of type BLAST_HitList and returns the result. Conversion
 * in this direction is lossless.  The list passed to this routine is
 * freed to ensure that there is no aliasing of fields between the
 * list of BlastCompo_Alignments and the new hitlist.
 *
 * @param search            general search information
 * @param alignments        a list of distinct alignments
 */
static BLAST_HitListPtr
Kappa_SortedHitlistFromAligns(
  BlastSearchBlkPtr search,
  BlastCompo_Alignment ** alignments)
{
  /* The context of the query is always zero in current
   * implementations of RedoAlignmentCore. */
  const Int2 context = 0;

  BLAST_HitListPtr hitlist;             /* the new hitlist */
  BlastCompo_Alignment * align;      /* represents the current
                                         * alignment in loops */
  if(search->current_hitlist != NULL) {
    search->current_hitlist->hspcnt_max = search->current_hitlist->hspcnt;
    BlastHitListPurge(search->current_hitlist);
  } else {
    search->current_hitlist =  BlastHitListNew(search);
  }
  hitlist = search->current_hitlist;
  hitlist->do_not_reallocate = FALSE;

  for(align = *alignments; align != NULL; align = align->next) {
    BLAST_HSPPtr hsp;           /* the new HSP for this alignment */
    /* queryExtent and matchExtent represent the extent of the
       alignment in the query and subject sequences respectively. */
    Int4 queryExtent = align->queryEnd - align->queryStart;
    Int4 matchExtent = align->matchEnd - align->matchStart;

    BlastSaveCurrentHsp(search, align->score, align->queryStart,
                        align->matchStart, matchExtent, context);

    hsp = hitlist->hsp_array[hitlist->hspcnt - 1];

    hsp->query.length   = queryExtent;
    hsp->query.end      = align->queryEnd;

    hsp->subject.length = matchExtent;
    hsp->subject.end    = align->matchEnd;

    hsp->subject.frame  = (Int2) align->frame;
    hsp->evalue         = 0.0;  /* E-values are computed after the full
                                   hitlist has been created. */
    hsp->gap_info       = (GapXEditBlockPtr) align->context;
    switch (align->matrix_adjust_rule) {
    case eDontAdjustMatrix:
        hsp->comp_adjustment_method = eNoCompositionBasedStats;
      break;
    case eCompoScaleOldMatrix:
        hsp->comp_adjustment_method = eCompositionBasedStats;
      break;
    default:
        hsp->comp_adjustment_method = eCompositionMatrixAdjust;
      break;
    }
    /* Break the aliasing between align->context and hsp->gap_info. */
    align->context = NULL;
  }
  BlastCompo_AlignmentsFree(alignments, NULL);
  if (hitlist->hspcnt > 0) {
    HeapSort(hitlist->hsp_array, hitlist->hspcnt, sizeof(BLAST_HSPPtr),
             score_compare_hsps);
  }
  return hitlist;
}


/**
 * Calculate expect values for a list of alignments and remove those 
 * that don't meet an evalue cutoff from the list.
 * 
 * @param *pbestScore      the largest score observed in the list
 * @param *pbestEvalue     the best (smallest) evalue observed in the
 *                         list
 * @param full_subject_length   the untranslated length of the subject 
 *                              sequence
 * @param search           contains search parameters used to compute 
 *                         evalues and the evalue cutoff
 * @param do_link_hsps     use hsp linking when computing evalues
 */
static void
Kappa_HitlistEvaluateAndPurge(int * pbestScore, double *pbestEvalue,
                             int full_subject_length,
                             BlastSearchBlkPtr search,
                             int do_link_hsps)
{
  BLAST_HitListPtr hitlist; /* a hitlist containing the newly-computed
                                   * alignments */
  Nlm_FloatHi bestEvalue;   /* best evalue among alignments in the
                                     hitlist */
  Int4 bestScore;           /* best score among alignments in the
                               hitlist */
  Int4 hsp_index;           /* index of the current HSP */
  
  hitlist = search->current_hitlist;
  search->subject->length  = full_subject_length;
  if (do_link_hsps) {
    BlastLinkHsps(search);
  } else {
    BlastGetNonSumStatsEvalue(search);
  }
  BlastReapHitlistByEvalue(search);
  /* Find the evalue of the best alignment in the list -- the list
   * is typically sorted by score and not by evalue, so a search is
   * necessary. */
  bestEvalue = DBL_MAX;
  bestScore  = 0;
  for (hsp_index = 0;  hsp_index < hitlist->hspcnt;  hsp_index++) {
    if( hitlist->hsp_array[hsp_index]->evalue < bestEvalue ) {
      bestEvalue = hitlist->hsp_array[hsp_index]->evalue;
    }
    if( hitlist->hsp_array[hsp_index]->score > bestScore ) {
      bestScore = hitlist->hsp_array[hsp_index]->score;
    }
  }
  *pbestScore  = bestScore;
  *pbestEvalue = bestEvalue;
}


/**
 * A callback routine: compute lambda for the given score frequencies.
 * (@sa calc_lambda_type).
 */
static double
Kappa_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
{
  int i, n;
  double avg;
  BLAST_ScoreFreq freq;

  n = max_score - min_score + 1;
  avg = 0.0;
  for (i = 0;  i < n;  i++) {
    avg += (min_score + i) * probs[i];
  }        
  freq.obs_min = min_score;
  freq.obs_max = max_score;
  freq.sprob = &probs[-min_score];
  freq.score_avg = avg;
  
  return impalaKarlinLambdaNR(&freq, lambda0);
}


/**
 * Get a matrix of the frequency ratios that underlie the score
 * matrix being used on this pass. The returned matrix is
 * position-specific, so if we are in the first pass, use query to
 * convert the 20x20 standard matrix into a position-specific
 * variant. matrixName is the name of the underlying 20x20 score
 * matrix used. numPositions is the length of the query;
 * startNumerator is the matrix of frequency ratios as stored in
 * posit.h. It needs to be divided by the frequency of the second
 * character to get the intended ratio
 *
 * @param returnRatios      an allocated matrix to hold the frequency
 *                          ratios [out]
 * @param search            general search information [in]
 * @param query             the query sequence [in]
 * @param matrixName        name of the underlying matrix [in]
 * @param startNumerator    matrix of frequency ratios [in]
 * @param numPositions      length of the query [in]
 * @param positionSpecific  is this a position-specific search? [in]
 */
static void
Kappa_GetStartFreqRatios(Nlm_FloatHi ** returnRatios,
                   BlastSearchBlkPtr search,
                   Uint1Ptr query,
                   const char *matrixName,
                   Nlm_FloatHi **startNumerator,
                   Int4 numPositions,
                   Boolean positionSpecific)
{
   Int4 i,j;
   FreqRatios * stdFreqRatios = NULL;
   /* a small cutoff used to determine whether it is necessary
    * to reverse the multiplication done in posit.c */
   const double posEpsilon = 0.0001;

   stdFreqRatios = PSIMatrixFrequencyRatiosNew(matrixName);
   if (positionSpecific) {
    for(i = 0; i < numPositions; i++) {
       for(j = 0; j < PROTEIN_ALPHABET; j++) {
         returnRatios[i][j] = stdFreqRatios->data[query[i]][j];
       }
     }
   } else {
     for (i = 0; i < PROTEIN_ALPHABET; i++) {
       for (j = 0; j < PROTEIN_ALPHABET; j++) {
         returnRatios[i][j] = stdFreqRatios->data[i][j];
       }
     }
   }
   stdFreqRatios = PSIMatrixFrequencyRatiosFree(stdFreqRatios);

   if (positionSpecific) {
     Nlm_FloatHi *standardProb; /*probabilities of each letter*/
     BLAST_ResFreqPtr stdrfp; /* gets standard frequencies in prob field */

     stdrfp = BlastResFreqNew(search->sbp);
     BlastResFreqStdComp(search->sbp,stdrfp);
     standardProb = &stdrfp->prob[0];

     /*reverse multiplication done in posit.c*/
     for(i = 0; i < numPositions; i++) {
       for(j = 0; j < PROTEIN_ALPHABET; j++) {
         if ((standardProb[query[i]] > posEpsilon) &&
             (standardProb[j] > posEpsilon) &&
             (j != eStopChar) && (j != Xchar) &&
             (startNumerator[i][j] > posEpsilon))
         {
           returnRatios[i][j] = startNumerator[i][j]/standardProb[j];
         }
       }
     }
     stdrfp = BlastResFreqDestruct(stdrfp);
   }
}


/** SCALING_FACTOR is a multiplicative factor used to get more bits of
 * precision in the integer matrix scores. It cannot be arbitrarily
 * large because we do not want total alignment scores to exceed
 * -(INT4_MIN) */
#define SCALING_FACTOR 32


/**
 * produce a scaled-up version of the position-specific matrix
 * starting from posFreqs
 *
 * @param fillPosMatrix     is the matrix to be filled
 * @param nonposMatrix      is the underlying position-independent matrix,
 *                          used to fill positions where frequencies are
 *                          irrelevant
 * @param matrixName        name of the postion-independent matrix
 * @param posFreq           frequecy ratios for the position-dependent
 *                          matrix
 * @param query             query sequence data
 * @param queryLength       length of the query
 * @param sbp               stores various parameters of the search
 */
static void
Kappa_ScalePosMatrix(BLAST_Score **fillPosMatrix, BLAST_Score **nonposMatrix,
                     const Char *matrixName, Nlm_FloatHi **posFreqs,
                     Uint1 *query, Int4 queryLength, BLAST_ScoreBlkPtr sbp,
                     Nlm_FloatHi localScalingFactor)
{

     posSearchItems *posSearch; /*used to pass data into scaling routines*/
     compactSearchItems *compactSearch; /*used to pass data into scaling routines*/
     Int4 i; /*loop index*/
     BLAST_ResFreqPtr stdrfp; /* gets standard frequencies in prob field */
     Int4 a; /*index over characters*/


     posSearch = (posSearchItems *) MemNew (1 * sizeof(posSearchItems));
     compactSearch = (compactSearchItems *) MemNew (1 * sizeof(compactSearchItems));
     posSearch->posMatrix = (BLAST_Score **) MemNew((queryLength + 1) * sizeof(BLAST_Score *));
     posSearch->posPrivateMatrix = fillPosMatrix;
     posSearch->posFreqs = posFreqs;
     posSearch->stdFreqRatios = PSIMatrixFrequencyRatiosNew(matrixName);
     for(i = 0; i <= queryLength; i++) 
       posSearch->posMatrix[i] = (BLAST_Score *) MemNew(PROTEIN_ALPHABET * sizeof(BLAST_Score));

     compactSearch->query = (Uint1Ptr) query;
     compactSearch->qlength = queryLength;
     compactSearch->alphabetSize = PROTEIN_ALPHABET;
     compactSearch->gapped_calculation = TRUE;
     compactSearch->matrix = nonposMatrix;
     compactSearch->lambda =  sbp->kbp_gap_std[0]->Lambda;
     compactSearch->kbp_std = sbp->kbp_std;
     compactSearch->kbp_psi = sbp->kbp_psi;
     compactSearch->kbp_gap_psi = sbp->kbp_gap_psi;
     compactSearch->kbp_gap_std = sbp->kbp_gap_std;
     compactSearch->lambda_ideal = sbp->kbp_ideal->Lambda;
     compactSearch->K_ideal = sbp->kbp_ideal->K;

     stdrfp = BlastResFreqNew(sbp);
     BlastResFreqStdComp(sbp,stdrfp); 
     compactSearch->standardProb = MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));
     for(a = 0; a < compactSearch->alphabetSize; a++)
       compactSearch->standardProb[a] = stdrfp->prob[a];
     stdrfp = BlastResFreqDestruct(stdrfp);

     posFreqsToMatrix(posSearch,compactSearch);
     impalaScaling(posSearch, compactSearch, localScalingFactor, FALSE);

     for(i = 0; i <= queryLength; i++)
       MemFree(posSearch->posMatrix[i]);

     MemFree(compactSearch->standardProb);
     MemFree(posSearch->posMatrix);
     PSIMatrixFrequencyRatiosFree(posSearch->stdFreqRatios);
     MemFree(posSearch);
     MemFree(compactSearch);
}


/** Convert an array of HSPs into a list of BlastCompo_Alignment
 * objects.
 *
 * @param search              search parameters
 * @param hsp_array           existing array of HSPs
 * @param hspcnt              length of hsp_array
 * @param localScalingFactor  factor by which scores are scaled in the
 *                            new list of alignments
 * @return                    the new list of alignments.
 */
static BlastCompo_Alignment * 
Kappa_ResultHspToDistinctAlign(BlastSearchBlkPtr search,
                               BLASTResultHsp hsp_array[], Int4 hspcnt,
                               double localScalingFactor)
{
  BlastCompo_Alignment *aligns;  /* list of alignments to be returned */ 
  BlastCompo_Alignment **tail;   /* next location to store a pointer to
                                       a new alignment */
  int hsp_index;                    /* loop index */

  aligns = NULL;
  tail = &aligns;
  for (hsp_index = 0;  hsp_index < hspcnt;  hsp_index++) {
    BlastCompo_Alignment *new_align;  /* a new alignment */
    int queryIndex, queryEnd, matchEnd;
    BLASTResultHsp * hsp = &hsp_array[hsp_index];
    queryEnd = hsp->query_offset + hsp->query_length;
    matchEnd = hsp->subject_offset + hsp->subject_length;
    if(search->mult_queries != NULL) {
      queryIndex = GetQueryNum(search->mult_queries,
                               hsp->query_offset, queryEnd - 1, 0);
    } else {
      queryIndex = 0;
    }
    new_align =
      BlastCompo_AlignmentNew((int) (hsp->score * localScalingFactor),
                                 eDontAdjustMatrix,
                                 hsp->query_offset, queryEnd,
                                 queryIndex, hsp->subject_offset,
                                 matchEnd, hsp->subject_frame, hsp);
    if (new_align == NULL) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0, 
                "Failed to allocate a new alignment");
    }
    *tail = new_align;
    tail = &new_align->next;
  }
  return aligns;
}


/**
 * Redo a S-W alignment using an x-drop alignment.  The result will
 * usually be the same as the S-W alignment. The call to ALIGN
 * attempts to force the endpoints of the alignment to match the
 * optimal endpoints determined by the Smith-Waterman algorithm.
 * ALIGN is used, so that if the data structures for storing BLAST
 * alignments are changed, the code will not break
 *
 * @param query         the query data
 * @param queryStart    start of the alignment in the query sequence
 * @param queryEnd      end of the alignment in the query sequence,
 *                      as computed by the Smith-Waterman algorithm
 * @param subject       the subject (database) sequence
 * @param matchStart    start of the alignment in the subject sequence
 * @param matchEnd      end of the alignment in the query sequence,
 *                      as computed by the Smith-Waterman algorithm
 * @param frame         translation frame of the database sequence (0
 *                      if not translated)
 * @param gap_align     parameters for a gapped alignment
 * @param score         score computed by the Smith-Waterman algorithm
 */
static void
Kappa_SWFindFinalEndsUsingXdrop(
  BlastCompo_SequenceData * query,
  Int4 queryStart,
  Int4 queryEnd,
  BlastCompo_SequenceData * subject,
  Int4 matchStart,
  Int4 matchEnd,
  Int4 frame,
  GapAlignBlkPtr gap_align,
  Int4 score)
{
  Int4 XdropAlignScore;         /* alignment score obtained using X-dropoff
                                 * method rather than Smith-Waterman */
  Int4 doublingCount = 0;       /* number of times X-dropoff had to be
                                 * doubled */
  Int4 *alignScript, *dummy;    /* the alignment script that will be
                                   generated below by the ALIGN
                                   routine. */
  GapXEditBlockPtr editBlock;   /* traceback info for this alignment */
  /* Extent of the alignment as computed by an x-drop alignment
   * (usually the same as (queryEnd - queryStart) and (matchEnd -
   * matchStart)) */
  Int4 queryExtent, matchExtent;

  gap_align->query_start = queryStart;
  gap_align->subject_start = matchStart;
  do {
    alignScript =
      (Int4 *) MemNew((subject->length + query->length + 3) * sizeof(Int4));

    XdropAlignScore =
      ALIGN(&query->data[queryStart - 1], &subject->data[matchStart - 1],
            queryEnd - queryStart + 1, matchEnd - matchStart + 1,
            alignScript, &queryExtent, &matchExtent, &dummy,
            gap_align, queryStart - 1, FALSE);
    gap_align->score = XdropAlignScore;
    gap_align->query_stop = gap_align->query_start + queryExtent - 1;
    gap_align->subject_stop = gap_align->subject_start + matchExtent - 1;

    gap_align->x_parameter *= 2;
    doublingCount++;
    if((XdropAlignScore < score) && (doublingCount < 3)) {
      MemFree(alignScript);
    }
  } while((XdropAlignScore < score) && (doublingCount < 3));
  editBlock =
    TracebackToGapXEditBlock(NULL, NULL, queryExtent, matchExtent,
                             alignScript, queryStart, matchStart);
  alignScript = MemFree(alignScript);
  editBlock->length1          = query->length;
  editBlock->length2          = subject->length;
  editBlock->discontinuous    = gap_align->discontinuous;
  editBlock->translate2       = gap_align->translate2;
  editBlock->frame2           = (Int2) frame;

  gap_align->edit_block = editBlock;
}


typedef struct Kappa_SequenceLocalData {
  Uint1         prog_number;    /**< identifies the type of blast search being
                                     performed. The type of search determines
                                     how sequence data should be obtained. */
  CharPtr       genetic_code;   /**< genetic code for translated searches */
  BioseqPtr     bsp_db;         /**< An object that represents this sequence
                                   in low level routines. */
  Boolean       bioseq_needs_unlock;    /**< true if the bsp_db must be
                                           disposed of by a call to
                                           BioseqUnlock, rather than
                                           BioseqFree */
  ReadDBFILEPtr rdfp;           /**< A pointer to a database from which
                                     sequences may be obtained */
} Kappa_SequenceLocalData;


/**
 * Initialize a new matching sequence, obtaining information about the
 * sequence from the search.
 *
 * @param self              object to be initialized
 * @param search            search information
 * @param subject_index     index of the matching sequence in the database
 */
static void
Kappa_MatchingSequenceInitialize(
  BlastCompo_MatchingSequence * self,
  BlastSearchBlkPtr search,
  Int4 subject_index)
{
  Kappa_SequenceLocalData * local_data;
  
  local_data = self->local_data = MemNew(sizeof(Kappa_SequenceLocalData));
  self->index = subject_index;
  
  local_data->prog_number  = search->prog_number;
  local_data->rdfp         = search->rdfp;
  if(local_data->prog_number ==  blast_type_tblastn) {
    local_data ->genetic_code = search->db_genetic_code;
  } else {
    /* the sequence will not be translated */
    local_data ->genetic_code = NULL;
  }
  if(local_data->rdfp) {
    local_data->rdfp->parameters |= READDB_KEEP_HDR_AND_SEQ;
    local_data->bsp_db = readdb_get_bioseq(local_data->rdfp, self->index );
    local_data->bioseq_needs_unlock = FALSE;
  } else {
    local_data->bsp_db = BioseqLockById(search->subject_info->sip);
    local_data->bioseq_needs_unlock = TRUE;
  }
  self->length = local_data->bsp_db->length;
}


/** Release the resources associated with a matching sequence. */
static void
Kappa_MatchingSequenceRelease(BlastCompo_MatchingSequence * self)
{
  Kappa_SequenceLocalData * local_data = self->local_data;
  if(local_data->bsp_db) {
    if(local_data->bioseq_needs_unlock) {
      BioseqUnlock(local_data->bsp_db);
      local_data->bsp_db = NULL;
    } else {
      local_data->bsp_db = BioseqFree(local_data->bsp_db);
    }
  } 
  MemFree(local_data);
  self->local_data = NULL;
}


/** Character to use for masked residues */
#define BLASTP_MASK_RESIDUE eXchar
/** Default instructions and mask residue for SEG filtering */
#define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"


/**
 * Obtain a string of translated data
 *
 * @param self          the sequence from which to obtain the data [in]
 * @param range         the range, in amino acid coordinates, of data
 *                      to get; includes the translation frame [in]
 * @param seqData       the resulting data [out]
 */
static void
Kappa_SequenceGetTranslatedRange(const BlastCompo_MatchingSequence * self,
                                 const BlastCompo_SequenceRange * range,
                                 BlastCompo_SequenceData * seqData )
{
  Int4 i;
  Int4 nucleotide_start;        /* position of the first nucleotide to be
                                 * translated */
  Int4 num_nucleotides;         /* number of nucleotides to translate */
  Int2 translation_frame = (Int2) range->context;
  Kappa_SequenceLocalData * local_data = self->local_data;

  nucleotide_start   = 3 * range->begin;
  num_nucleotides    =
    3 * (range->end - range->begin) + ABS(translation_frame) - 1;
  { /* scope of nucleotide_data */
    Uint1Ptr nucleotide_data = MemNew((num_nucleotides + 1) * sizeof(Uint1));
    { /* scope of spp */
      SeqPortPtr spp;       /* a SeqPort used to read the sequence data */
      Uint1 strand;         /* a flag indicating which strand to read */
      Uint1 nucleotide;     /* an individual nucleotide */

      if (translation_frame >= 0) {
          strand = Seq_strand_plus;
      } else {
          strand = Seq_strand_minus;
      }
      spp    = SeqPortNew(local_data->bsp_db, FIRST_RESIDUE, LAST_RESIDUE,
                          strand, Seq_code_ncbi4na);
      SeqPortSeek(spp, nucleotide_start, SEEK_SET);

      for(i = 0;
          i < num_nucleotides &&
            (nucleotide = SeqPortGetResidue(spp)) != SEQPORT_EOF;
          i++ ) {
        /* for all nucleotides in the translation range */
        nucleotide_data[i] = nucleotide;
      }

      spp = SeqPortFree(spp);
    }  /* end scope of spp */
    seqData->buffer =
      GetTranslation(nucleotide_data, num_nucleotides, translation_frame,
                     &seqData->length, local_data->genetic_code);
    seqData->data = seqData->buffer + 1; /* This is a protein sequence,
                                            so the first byte is nil. */
    MemFree(nucleotide_data);
  } /* end scope of nucleotide_data */

#ifndef KAPPA_NO_SEG_SEQUENCE_TBLASTN
  { /* scope of variables used for SEG filtering */
    BioseqPtr bsp_temp; /* a Bioseq for the translated sequence */
    ObjectIdPtr oip;    /* a unique ObjectId for the translated sequence */
    SeqLocPtr seg_slp;  /* a SeqLoc for SEG filtering */
    bsp_temp     = BlastMakeTempProteinBioseq(seqData->data, seqData->length,
                                              Seq_code_ncbistdaa);
    bsp_temp->id = SeqIdSetFree(bsp_temp->id);

    oip = (ObjectIdPtr) UniqueLocalId();
    ValNodeAddPointer(&(bsp_temp->id), SEQID_LOCAL, oip);
    SeqMgrAddToBioseqIndex(bsp_temp);

    seg_slp = BlastBioseqFilter(bsp_temp, BLASTP_MASK_INSTRUCTIONS);
    if (seg_slp) {
      HackSeqLocId(seg_slp, local_data->bsp_db->id);
      BlastMaskTheResidues(seqData->data, seqData->length,
                           BLASTP_MASK_RESIDUE, seg_slp, FALSE, 0);
      seg_slp = SeqLocSetFree(seg_slp);
    }

    SeqMgrDeleteFromBioseqIndex(bsp_temp);
    bsp_temp->id = SeqIdSetFree(bsp_temp->id);
    bsp_temp     = BioseqFree(bsp_temp);
  } /* end scope of variables used for SEG filtering */
#endif
}


/**
 * Obtain the sequence data that lies within the given range.
 *
 * @param self          sequence information [in]
 * @param range         the range, in amino acid coordinates, of data
 *                      to get [in]
 * @param seqData       the sequence data obtained [out]
 * 
 * @returns   always 0 (posts a fatal error on failure rather than
 *            returning an error code.)
 */
static int
Kappa_SequenceGetRange(
  const BlastCompo_MatchingSequence * self,
  const BlastCompo_SequenceRange * range,
  BlastCompo_SequenceData * seqData )
{
  Kappa_SequenceLocalData * local_data = self->local_data;
  if(local_data->prog_number ==  blast_type_tblastn) {
    /* The sequence must be translated. */
    Kappa_SequenceGetTranslatedRange(self, range, seqData);
  } else {
    /* The sequence does not need to be translated. */
    /* Obtain the entire sequence (necessary for SEG filtering.) */
    if(local_data->rdfp != NULL) {
      Uint1Ptr origData;        /* data obtained from readdb_get_sequence;
                                 * this data cannot be modified, so we copy
                                 * it. */
      Int4       idx;
      seqData->length    = readdb_get_sequence(local_data->rdfp, self->index,
                                               (Uint1Ptr PNTR) & origData );
      seqData->buffer    = MemNew((seqData->length + 2) * sizeof(Uint1));
      seqData->buffer[0] = '\0';
      seqData->data      = &seqData->buffer[1];
      for( idx = 0; idx < seqData->length; idx++ ) {
        /* Copy the sequence data, replacing occurrences of amino acid
         * number 24 (Selenocysteine) with number 21 (Undetermined or
         * atypical). */
        if(origData[idx] != eSelenocysteine) {
            seqData->data[idx] = origData[idx];
        } else {
          seqData->data[idx] = eXchar;
          fprintf(stderr, "Selenocysteine (U) at position %ld"
                  " replaced by X\n",
                  (long) idx + 1);
        }
      }
      seqData->data[seqData->length] = '\0';
    } else { /* self->rdfp is NULL */
      SeqPortPtr spp = NULL;      /* a SeqPort used to read the
                                     sequence data */
      Uint1      residue;         /* an individual residue */
      Int4       idx;

      seqData->length    = local_data->bsp_db->length;
      seqData->buffer    = MemNew((seqData->length + 2) * sizeof(Uint1));
      seqData->buffer[0] = '\0';
      seqData->data      = seqData->buffer + 1;

      spp =
        SeqPortNew(local_data->bsp_db, FIRST_RESIDUE, LAST_RESIDUE,
                   Seq_strand_unknown, Seq_code_ncbistdaa);

      idx = 0;
      while((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF) {
        if(IS_residue(residue)) {
          /* Replace occurrences of amino acid number 24
             (Selenocysteine) with number 21 (Undetermined or
             atypical). */
          if(residue == eSelenocysteine) {
            residue = eXchar;
            fprintf(stderr, "Selenocysteine (U) at position %ld"
                    " replaced by X\n",
                    (long) idx + 1);
          }
          seqData->data[idx++] = residue;
        }
      }
      seqData->data[idx] = 0;    /* terminate the string */
      spp = SeqPortFree(spp);
    } /* end else self->rdfp is NULL */

#ifndef KAPPA_BLASTP_NO_SEG_SEQUENCE
    {
      SeqLocPtr seg_slp;  /*pointer to structure for SEG filtering*/

      seg_slp =
        BlastBioseqFilter(local_data->bsp_db, BLASTP_MASK_INSTRUCTIONS);
      if (seg_slp) {
        BlastMaskTheResidues(seqData->data, seqData->length,
                             BLASTP_MASK_RESIDUE, seg_slp, FALSE, 0);
        seg_slp = SeqLocSetFree(seg_slp);
      }
    }
#endif
    /* Fit the data to the range. */
    seqData ->data    = &seqData->data[range->begin - 1];
    *seqData->data++  = '\0';
    seqData ->length  = range->end - range->begin;
  } /* end else the sequence does not need to be translated */
  return 0;
}


/**
 * Converts an HSP, obtained from a blast search, to a GapAlignBlk that
 * is in a state in which it has all information necessary to redo the
 * computation of a traceback.
 *
 * @param gap_align     the GapAlignBlk to be modified  [out]
 * @param search        general information about the search [in]
 * @param hsp           the HSP to be converted [in]
 * @param queryOrigin   origin of the query participating in this
 *                      alignment within the concatenated query [in]
 * @param subject_range,       the subject_range used to compute the traceback,
 *                             in amino acid coordinates [in]
 * @param query,        the query data [in]
 * @param subject       the subject data [in]
 */
static void
Kappa_HitToGapAlign(
  GapAlignBlkPtr gap_align,
  BlastSearchBlkPtr search,
  BLAST_HSPPtr       hsp,
  Int4 queryOrigin,
  BlastCompo_SequenceRange * subject_range,
  BlastCompo_SequenceData * query,
  BlastCompo_SequenceData * subject)
{
  gap_align->query           = query->data;
  gap_align->query_length    = query->length;
  gap_align->query_frame     = 0;

  gap_align->subject         = subject->data;
  gap_align->subject_length  = subject->length;
  gap_align->subject_frame   = hsp->subject.frame;

  /* Shift the hsp to use query/subject local coordinates. */
  hsp->query.offset         -= queryOrigin;
  hsp->query.gapped_start   -= queryOrigin;
  hsp->query.end            -= queryOrigin;
  hsp->subject.offset       -= subject_range->begin;
  hsp->subject.gapped_start -= subject_range->begin;
  hsp->subject.end          -= subject_range->begin;
  if(CheckStartForGappedAlignment(search, hsp, gap_align->query,
                                   gap_align->subject, search->sbp->matrix)) {
    /* We may use the starting point supplied by the HSP. */
    gap_align->q_start = hsp->query.gapped_start;
    gap_align->s_start = hsp->subject.gapped_start;
  } else {
    /* We must recompute the start for the gapped alignment, as the
       one in the HSP was unacceptable.*/
    gap_align->q_start =
      GetStartForGappedAlignment(search, hsp, gap_align->query,
                                 gap_align->subject, search->sbp->matrix);
    gap_align->s_start =
      (hsp->subject.offset - hsp->query.offset) + gap_align->q_start;
  }
}


/**
 * Reads a GapAlignBlk that has been used to compute a traceback, and
 * return a BlastCompo_Alignment representing the alignment.
 *
 * @param gap_align         the GapAlignBlk
 * @param matrix_adjust_rule    the rule used to compute the scoring matrix
 * @param queryRange        range of the query data in the
 *                          concatenated query
 * @param ccat_query_length length of the concatenated query
 * @param subjectRange      range of the subject data in the full
 *                          database sequence, expressed in amino acid
 *                          coordinates
 * @param subjectLength     original. untranslated length of the
 *                          subject sequence
 */
static BlastCompo_Alignment *
Kappa_NewAlignFromGapAlign(
  GapAlignBlkPtr gap_align,
  EMatrixAdjustRule matrix_adjust_rule,
  BlastCompo_SequenceRange * query_range,
  Int4 ccat_query_length,
  BlastCompo_SequenceRange * subject_range,
  Int4 subjectLength)
{
  int query_index, translation_frame;
  BlastCompo_Alignment * obj; /* the new alignment */

  /* The gap_align is in coordinates that are relative to the
     query/subject subject_range.  Shift to global coordinates. */
  int queryStart = gap_align->query_start   + query_range->begin;
  int queryEnd   = gap_align->query_stop    + query_range->begin   + 1;
  int matchStart = gap_align->subject_start + subject_range->begin;
  int matchEnd   = gap_align->subject_stop  + subject_range->begin + 1;

  if(gap_align->edit_block != NULL) {
    gap_align->edit_block->start1           += query_range->begin;
    gap_align->edit_block->length1          += query_range->begin;
    gap_align->edit_block->start2           += subject_range->begin;
    gap_align->edit_block->length2          += subject_range->begin;
    gap_align->edit_block->original_length1  = ccat_query_length;
    gap_align->edit_block->original_length2  = subjectLength;
  }
  query_index = query_range->context;
  translation_frame = subject_range->context;
  obj = BlastCompo_AlignmentNew(gap_align->score, matrix_adjust_rule,
                                   queryStart, queryEnd, query_index,
                                   matchStart, matchEnd,
                                   translation_frame,
                                   gap_align->edit_block);
  if (obj == NULL) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0,
                "Failed to allocate a new alignment");
  }
  gap_align->edit_block = NULL; /* set to NULL to avoid aliasing */
  
  return obj;
}


/** A callback used when performing SmithWaterman alignments:
 * Calculate the traceback for one alignment by performing an x-drop
 * alignment in the forward direction, possibly increasing the x-drop
 * parameter until the desired score is attained.  
 *
 * The start, end and score of the alignment should be obtained
 * using the Smith-Waterman algorithm before this routine is called.
 *
 * @param *palign          the new alignment
 * @param *pqueryEnd       on entry, the end of the alignment in the
 *                         query, as computed by the Smith-Waterman
 *                         algorithm.  On exit, the end as computed by
 *                         the x-drop algorithm
 * @param *pmatchEnd       like as *pqueryEnd, but for the subject
 *                         sequence
 * @param queryStart       the starting point in the query
 * @param matchStart       the starting point in the subject
 * @param score            the score of the alignment, as computed by
 *                         the Smith-Waterman algorithm
 * @param query            query sequence data
 * @param query_range      range of this query in the concatenated
 *                         query
 * @param ccat_query_length   total length of the concatenated query
 * @param subject          subject sequence data
 * @param subject_range    range of subject_data in the translated
 *                         query, in amino acid coordinates
 * @param full_subject_length   length of the full subject sequence
 * @param gapping_params        parameters used to compute gapped
 *                              alignments
 * @param matrix_adjust_rule    the rule used to compute the scoring matrix
 * @returns 0   (posts a fatal error if it fails)
 * @sa new_xdrop_align_type
 */
static int
Kappa_NewAlignmentUsingXdrop(BlastCompo_Alignment ** pnewAlign,
                             Int4 * pqueryEnd, Int4 *pmatchEnd,
                             Int4 queryStart, Int4 matchStart, Int4 score,
                             BlastCompo_SequenceData * query,
                             BlastCompo_SequenceRange * query_range,
                             Int4 queryLength,
                             BlastCompo_SequenceData * subject,
                             BlastCompo_SequenceRange * subject_range,
                             Int4 subjectLength,
                             BlastCompo_GappingParams * gapping_params,
                             EMatrixAdjustRule matrix_adjust_rule)
{ 
  /* General search parameters */
  BlastSearchBlkPtr search = gapping_params->context;
  /* Specific parameeters used for gapped alignments */
  GapAlignBlkPtr gap_align = search->gap_align;
  /* the translation frame */
  int frame = subject_range->context;
  
  gap_align->x_parameter = gapping_params->x_dropoff;

  Kappa_SWFindFinalEndsUsingXdrop(query,   queryStart, *pqueryEnd,
                                  subject, matchStart, *pmatchEnd,
                                  frame, gap_align, score);
  *pqueryEnd = gap_align->query_stop + 1;
  *pmatchEnd = gap_align->subject_stop + 1;
  *pnewAlign = Kappa_NewAlignFromGapAlign(gap_align, matrix_adjust_rule,
                                          query_range, queryLength,
                                          subject_range, subjectLength);
  gap_align->x_parameter = gapping_params->x_dropoff;

  return 0;
}


/**
 * A callback: calculate the traceback for one alignment by
 * performing an x-drop alignment in both directions
 *
 * @param in_align         the existing alignment, without traceback
 * @param matrix_adjust_rule    the rule used to compute the scoring matrix
 * @param query_data       query sequence data
 * @param query_range      range of this query in the concatenated
 *                         query
 * @param ccat_query_length   total length of the concatenated query
 * @param subject_data     subject sequence data
 * @param subject_range    range of subject_data in the translated
 *                         query, in amino acid coordinates
 * @param full_subject_length   length of the full subject sequence
 * @param gapping_params        parameters used to compute gapped
 *                              alignments
 * @sa redo_one_alignment_type
 */
static BlastCompo_Alignment *
Kappa_RedoOneAlignment(BlastCompo_Alignment * in_align,
                       EMatrixAdjustRule matrix_adjust_rule,
                       BlastCompo_SequenceData * query_data,
                       BlastCompo_SequenceRange * query_range,
                       int ccat_query_length,
                       BlastCompo_SequenceData * subject_data,
                       BlastCompo_SequenceRange * subject_range,
                       int full_subject_length,
                       BlastCompo_GappingParams * gapping_params)
{
  /* New alignment to be returned */
  BlastCompo_Alignment * newAlign;
  /* The HSP to be redone */
  BLASTResultHspPtr hsp = in_align->context;
  /* A copy of hsp as a BLAST_HSP; needed to fix a type mismatch */
  BLAST_HSPPtr temp_hsp = MemNew(sizeof(BLAST_HSP));
  /* General search information */
  BlastSearchBlkPtr search = gapping_params->context;
  
  CopyResultHspToHSP(hsp, temp_hsp);
  Kappa_HitToGapAlign(search->gap_align, search, temp_hsp, query_range->begin, 
                      subject_range, query_data, subject_data);
  search->gap_align->x_parameter = gapping_params->x_dropoff;

  PerformGappedAlignmentWithTraceback(search->gap_align);

  newAlign = 
      Kappa_NewAlignFromGapAlign(search->gap_align, matrix_adjust_rule,
                                 query_range, ccat_query_length,
                                 subject_range, full_subject_length);
  MemFree(temp_hsp);
  
  return newAlign;
}


/**
 * A Kappa_SearchParameters represents the data needed by
 * RedoAlignmentCore to adjust the parameters of a search, including
 * the original value of these parameters
 */
typedef struct Kappa_SearchParameters {
  Int4          gap_open;        /**< a penalty for the existence of a gap */
  Int4          gapExtend;      /**< a penalty for each residue (or
                                      nucleotide) in the gap */
  Int4          gapDecline;     /**< a penalty for declining to align a pair
                                     of residues */
  int           gap_x_dropoff_final;
  BLAST_Score **origMatrix;     /**< The original matrix values */
  /** a matrix.  Each row represents a query and each column
       represents the probabilties of a particular residue */
  GapAlignBlkPtr orig_gap_align;
  Nlm_FloatHi original_expect_value;    /**< expect value on entry */
  BLAST_KarlinBlkPtr kbp_gap_orig;  /**< copy of the original gapped
                                      Karlin-Altschul block corresponding to
                                      the first context */
  BLAST_KarlinBlkPtr * orig_kbp_gap_array;  /**< pointer to the array of gapped
                                              Karlin-Altschul block for all
                                              contexts; needed to restore
                                              the search to its original 
                                              configuration.  */
  Boolean    use_mq_flag;     /**< the initial value of the
                                   search->mult_queries->use_mq flag;
                                   or false if search->mult_queries ==
                                   NULL */
} Kappa_SearchParameters;


/**
 * Release the data associated with a Kappa_SearchParameters and
 * delete the object
 * @param searchParams the object to be deleted [in][out]
 */
static void
Kappa_SearchParametersFree(Kappa_SearchParameters ** searchParams)
{
  /* for convenience, remove one level of indirection from searchParams */
  Kappa_SearchParameters *sp = *searchParams;

  if(sp->kbp_gap_orig) BlastKarlinBlkDestruct(sp->kbp_gap_orig);

  Nlm_Int4MatrixFree(&sp->origMatrix);

  Nlm_Free(*searchParams);
  *searchParams = NULL;
}


/**
 * Create a new instance of Kappa_SearchParameters
 *
 * @param rows              number of rows in the scoring matrix
 * @param compo_adjust_mode  if >0, use composition-based statistics
 * @param numQueries        the number of queries in the concatenated
 *                          query
 * @param positionBased     if true, the search is position-based
 */
static Kappa_SearchParameters *
Kappa_SearchParametersNew(
  Int4 rows,
  ECompoAdjustModes compo_adjust_mode,
  Boolean positionBased)
{
  Kappa_SearchParameters *sp;   /* the new object */
  sp = MemNew(sizeof(Kappa_SearchParameters));

  sp->orig_kbp_gap_array = NULL;
  sp->kbp_gap_orig       = NULL;
  sp->origMatrix         = NULL;
  
  sp->kbp_gap_orig = BlastKarlinBlkCreate();
  if (compo_adjust_mode != eNoCompositionBasedStats) {
    if (positionBased) {
      sp->origMatrix = Nlm_Int4MatrixNew(rows, PROTEIN_ALPHABET);
    } else {
      sp->origMatrix = Nlm_Int4MatrixNew(PROTEIN_ALPHABET, PROTEIN_ALPHABET);  
    }
    if (sp->origMatrix == NULL) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0, 
                "Failed to allocate a scoring matrix");
    }
  }
  /* end if(compo_adjust_mode) */

  return sp;
}


/**
 * Record the initial value of the search parameters that are to be
 * adjusted.
 *
 * @param searchParams       holds the recorded values [out]
 * @param search             the search parameters [in]
 * @param compo_adjust_mode   mode of composition adjustment [in]
 */
static void
Kappa_RecordInitialSearch(Kappa_SearchParameters * searchParams,
                          BlastSearchBlkPtr search,
                          ECompoAdjustModes compo_adjust_mode)
{
  BLAST_KarlinBlkPtr kbp;     /* statistical parameters used to evaluate a
                               * query-subject pair */
  BLAST_Score **matrix;       /* matrix used to score a local
                                 query-subject alignment */

  searchParams->gap_x_dropoff_final   = search->pbp->gap_x_dropoff_final;
  searchParams->original_expect_value = search->pbp->cutoff_e;

  if (search->mult_queries != NULL) {
    searchParams->use_mq_flag = search->mult_queries->use_mq;
  } else {
    searchParams->use_mq_flag = 0;
  }
  searchParams->gap_open   = search->pbp->gap_open;
  searchParams->gapExtend  = search->pbp->gap_extend;
  searchParams->gapDecline = search->pbp->decline_align;

  searchParams->orig_kbp_gap_array   = search->sbp->kbp_gap;
  if(search->positionBased) {
    kbp    = search->sbp->kbp_gap_psi[0];
    matrix = search->sbp->posMatrix;
  } else {
    kbp    = search->sbp->kbp_gap_std[0];
    matrix = search->sbp->matrix;
  }
  searchParams->kbp_gap_orig->Lambda = kbp->Lambda;
  searchParams->kbp_gap_orig->K      = kbp->K;
  searchParams->kbp_gap_orig->logK   = kbp->logK;
  searchParams->kbp_gap_orig->H      = kbp->H;

  searchParams->orig_gap_align = search->gap_align;
  search->gap_align = NULL; /* break aliasing */
  
  if(compo_adjust_mode != eNoCompositionBasedStats) {
    Int4 i, j;                  /* iteration indices */
    int rows;
    if (search->positionBased) {
      rows = search->context[0].query->length;
    } else {
      rows = PROTEIN_ALPHABET;
    }
    for(i = 0; i < rows; i++) {
      for(j = 0; j < PROTEIN_ALPHABET; j++) {
        searchParams->origMatrix[i][j] = matrix[i][j];
      }
    }
  }
}


/**
 * Rescale the search parameters in the search object to obtain more
 * precision.
 *
 * @param search               the search block to be rescaled
 * @param localScalingFactor   the scale factor
 * @param options              certain command-line options to BLAST
 */
static void
Kappa_RescaleSearch(BlastSearchBlkPtr search,
                double localScalingFactor,
                BLAST_OptionsBlkPtr options)
{
  BLAST_KarlinBlkPtr kbp;     /* the statistical parameters used to
                               * evaluate alignments of a
                               * query-subject pair */
  if(search->positionBased) {
    kbp = search->sbp->kbp_gap_psi[0];
  } else {
    kbp = search->sbp->kbp_gap_std[0];
  }
  kbp->Lambda /= localScalingFactor;
  kbp->logK = log(kbp->K);
  
  if (search->mult_queries != NULL) {
    search->mult_queries->use_mq = 0;
  }
  search->pbp->cutoff_e = options->kappa_expect_value;
  search->pbp->gap_open =
    Nlm_Nint(search->pbp->gap_open * localScalingFactor);
  search->pbp->gap_extend =
    Nlm_Nint(search->pbp->gap_extend * localScalingFactor);
  if(search->pbp->decline_align < INT2_MAX) {
    search->pbp->decline_align =
      Nlm_Nint(search->pbp->decline_align * localScalingFactor);
  }
  search->gap_align                = GapAlignBlkNew(1, 1);
  search->gap_align->matrix        = search->sbp->matrix;
  search->gap_align->positionBased = search->positionBased;
  if(search->positionBased) {
    search->gap_align->posMatrix   = search->sbp->posMatrix;
  }
  search->gap_align->gap_open      = search->pbp->gap_open;
  search->gap_align->gap_extend    = search->pbp->gap_extend;
  search->gap_align->decline_align = search->pbp->decline_align;
  search->gap_align->translate2 = 
      (Boolean) (search->prog_number == blast_type_tblastn);

  search->gap_align->x_parameter   =
    (Int4) (options->gap_x_dropoff_final * NCBIMATH_LN2 / kbp->Lambda);
}


/**
 * Restore the parameters that were adjusted to their original values
 * @param search            the search to be restored [out]
 * @param matrix            the scoring matrix to be restored [out]
 * @param searchParams      a record of the original values [in]
 * @param SmithWaterman     if true, we have performed a Smith-Waterman
 *                          alignment with these search parameters [in]
 */
static void
Kappa_RestoreSearch(
  BlastSearchBlkPtr search,
  BLAST_Score ** matrix,
  Kappa_SearchParameters * searchParams,
  ECompoAdjustModes compo_adjust_mode)
{
  BLAST_KarlinBlkPtr kbp;     /* statistical parameters used to
                                 evaluate the significance of
                                 alignment of a query-subject
                                   pair */
  Int4 i, j;                  /* iteration indices */
  search->pbp->gap_x_dropoff_final = searchParams->gap_x_dropoff_final;
  search->pbp->cutoff_e      = searchParams->original_expect_value;
  search->pbp->gap_open      = searchParams->gap_open;
  search->pbp->gap_extend    = searchParams->gapExtend;
  search->pbp->decline_align = searchParams->gapDecline;
  if (search->mult_queries != NULL) {
      search->mult_queries->use_mq = searchParams->use_mq_flag;
  }
  GapAlignBlkDelete(search->gap_align);
  search->gap_align    = searchParams->orig_gap_align;
  search->sbp->kbp_gap = searchParams->orig_kbp_gap_array;
  
  if(search->positionBased) {
    kbp = search->sbp->kbp_gap_psi[0];
  } else {
    kbp = search->sbp->kbp_gap_std[0];
  }
  kbp->Lambda = searchParams->kbp_gap_orig->Lambda;
  kbp->K      = searchParams->kbp_gap_orig->K;
  kbp->logK   = searchParams->kbp_gap_orig->logK;
  kbp->H      = searchParams->kbp_gap_orig->H;
  
  if(compo_adjust_mode != eNoCompositionBasedStats) {
    int rows;    /* the number of rows in matrix */
    if (search->positionBased) {
      rows = search->context[0].query->length;
    } else {
      rows = PROTEIN_ALPHABET;
    }
    for(i = 0; i < rows; i++) {
      for(j = 0; j < PROTEIN_ALPHABET; j++) {
        matrix[i][j] = searchParams->origMatrix[i][j];
      }
    }
  }
}


/** Initialize an object of type Blast_MatrixInfo.
 * @param self                 object to initialize
 * @param search               search data
 * @param localScalingFactor   amount by which this search is scaled
 * @param matrixName           name of the matrix */
static void
Kappa_MatrixInfoInit(Blast_MatrixInfo * self,
                 double localScalingFactor,
                 BlastSearchBlkPtr search,
                 const char * matrixName)
{
  Uint1Ptr query;             /* the query sequence */
  Int4 queryLength;           /* the length of the query sequence */
  Nlm_FloatHi initialUngappedLambda;

  query       = search->context[0].query->sequence;
  queryLength = search->context[0].query->length;

  if(self->positionBased) {
    if(search->sbp->posFreqs == NULL) {
      search->sbp->posFreqs =
        allocatePosFreqs(queryLength, PROTEIN_ALPHABET);
    }
    Kappa_GetStartFreqRatios(self->startFreqRatios, search, query, matrixName,
                             search->sbp->posFreqs, queryLength, TRUE);
    Kappa_ScalePosMatrix(self->startMatrix, search->sbp->matrix, matrixName,
                         search->sbp->posFreqs, query, queryLength,
                         search->sbp, localScalingFactor);
    initialUngappedLambda = search->sbp->kbp_psi[0]->Lambda;
  } else {
    Kappa_GetStartFreqRatios(self->startFreqRatios, search, query,
                             matrixName, NULL, PROTEIN_ALPHABET,
                             FALSE);
    initialUngappedLambda = search->sbp->kbp_ideal->Lambda;
  }
  self->ungappedLambda = initialUngappedLambda / localScalingFactor;
  if(!search->positionBased) {
    FreqRatios * freqRatios;  /* frequency ratios for the matrix */

    freqRatios = PSIMatrixFrequencyRatiosNew(matrixName);
    if (freqRatios == NULL) {
      ErrPostEx(SEV_FATAL, 1, 0, "blastpgp: Cannot adjust parameters "
                "for matrix %s\n", matrixName);
    }
    Blast_Int4MatrixFromFreq(self->startMatrix, PROTEIN_ALPHABET,
                             freqRatios->data, self->ungappedLambda);
    freqRatios = PSIMatrixFrequencyRatiosFree(freqRatios);
  }
  self->matrixName = strdup(matrixName);
}


/**
 * Record information about all queries in the concatenated query.
 *
 * @param *pquery       a new list of BlastCompo_QueryInfo objects
 * @param *pnumQueries  the lenght of *pquery
 * @param *maxLength    the length of the longest query
 * @param search        a search block, which is used to obtain the 
 *                      query information */
static void
Kappa_GetQueryInfo(BlastCompo_QueryInfo **pquery, int * pnumQueries,
                   int *maxLength, BlastSearchBlkPtr search)
{
  int query_index;
  int numQueries = search->mult_queries ? search->mult_queries->NumQueries : 1;
  BlastCompo_QueryInfo * query =
      Nlm_Calloc(numQueries, sizeof(BlastCompo_QueryInfo));

  *pnumQueries = numQueries;
  *maxLength = 0;
  *pquery = query;
  if (search->mult_queries) {
    for (query_index = 0;  query_index < numQueries;  query_index++) {
      query[query_index].eff_search_space =
        search->mult_queries->SearchSpEff[query_index];
    }
  } else {
    query[0].eff_search_space = search->searchsp_eff;
  }
  if (search->mult_queries != NULL) {
    /* Query concatenation is in use; find each individual query
       within the concatenated query */
    Uint1 * ccat_query = /* the concatenated query */
      search->context[0].query->sequence;
    int length;

    for (query_index = 0;  query_index < numQueries;  query_index++) {
      query[query_index].origin =
        search->mult_queries->QueryStarts[query_index];
      query[query_index].seq.data = &ccat_query[query[query_index].origin];
      length = search->mult_queries->QueryEnds[query_index] -
          query[query_index].origin + 1;
      query[query_index].seq.length = length;
      if (length > *maxLength) {
          *maxLength = length;
      }
    }
  } else {
    /* Query concatenation is not in use; just use the entire query */
    query[0].seq.data   = search->context[0].query->sequence;
    query[0].seq.length = search->context[0].query->length;
    query[0].origin = 0;
    *maxLength = query[0].seq.length;
  }
  for (query_index = 0;  query_index < numQueries;  query_index++) {
    Blast_ReadAaComposition(&query[query_index].composition,
                            query[query_index].seq.data,
                            query[query_index].seq.length);
  }
}


/**
 * Initialize an object that contains the parameters needed to perform
 * a gapped alignment.
 *
 * @param gapping_params    the object to be initialized
 * @param search            holds search parameters
 * @param options           holds BLAST command-line options
 * @param Lambda            the statistical parameter Lambda;
 *                          gives the scale of the scoring system.
 */
static void
Kappa_GappingParamsInit(BlastCompo_GappingParams * gapping_params,
                        BlastSearchBlkPtr search)
{  
  gapping_params->gap_open = search->gap_align->gap_open;
  gapping_params->gap_extend = search->gap_align->gap_extend;
  gapping_params->decline_align = search->gap_align->decline_align;
  gapping_params->x_dropoff = search->gap_align->x_parameter;
  gapping_params->context = search;
}


/** 
 * A callback: Free GapXEditBlock that is pointed to by a (void *).
 */
static void Kappa_FreeEditBlock(void * traceback)
{
    if (traceback != NULL)
        GapXEditBlockDelete((GapXEditBlockPtr)traceback);
}


/**
 * Callbacks used by Blast_RedoOneMatch and
 * Blast_RedoOneMatchSmithWaterman */
static const Blast_RedoAlignCallbacks
redo_align_callbacks = {
    Kappa_CalcLambda, Kappa_SequenceGetRange, Kappa_RedoOneAlignment,
    Kappa_NewAlignmentUsingXdrop, Kappa_FreeEditBlock
};


/** Get a set of alignment parameters, for use by Blast_RedoOneMatch or
 * Blast_RedoOneMatchSmithWaterman.
 * @param search               search parameters
 * @param options              BLAST command line options
 * @param compo_adjust_mode     composition adjustment mode
 * @param localScalingFactor   factor by which scores are scaled */
static Blast_RedoAlignParams *
Kappa_GetAlignParams(BlastSearchBlkPtr search,
                     BLAST_OptionsBlkPtr options,
                     ECompoAdjustModes compo_adjust_mode,
                     double localScalingFactor)
{
  int ccat_query_length;
  int rows;
  int cutoff_s;
  double cutoff_e;
  BlastCompo_GappingParams * gapping_params = NULL;
  Blast_MatrixInfo * scaledMatrixInfo;
  Blast_RedoAlignParams * params;
  int subject_is_translated = search->prog_number == blast_type_tblastn;
  int do_link_hsps = search->pbp->do_sum_stats;

  if(do_link_hsps) {
    cutoff_s = (int) (search->pbp->cutoff_s1 * localScalingFactor);
  } else {
    /* There is no cutoff score; we consider e-values instead */
    cutoff_s = 0;
  }
  cutoff_e = search->pbp->cutoff_e;
  ccat_query_length = search->context[0].query->length;
  rows = search->positionBased ? ccat_query_length : PROTEIN_ALPHABET;
  if (compo_adjust_mode == eNoCompositionBasedStats) {
    /* We cannot compute a value for scaledMatrixInfo if
     * compo_adjust_mode is off.  There are side effects to the
     * computation (specifically a call to updateLambdaK that alters
     * the gapped K and ungapped lambda for position-based searches)
     * that are undesirable and deeply embedded in the code. */
    scaledMatrixInfo = NULL;
  } else {
    scaledMatrixInfo = Blast_MatrixInfoNew(rows, search->positionBased);
    if (scaledMatrixInfo == NULL) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
    }
    Kappa_MatrixInfoInit(scaledMatrixInfo, localScalingFactor, search,
                         options->matrix);
  }
  gapping_params = malloc(sizeof(BlastCompo_GappingParams));
  Kappa_GappingParamsInit(gapping_params, search);
  params =
    Blast_RedoAlignParamsNew(&scaledMatrixInfo, &gapping_params,
                             compo_adjust_mode, search->positionBased,
                             subject_is_translated, ccat_query_length,
                             cutoff_s, cutoff_e, do_link_hsps,
                             &redo_align_callbacks);
  if (params == NULL) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
  }
  return params;
}


/**
 * Convert a Blast_CompoHeap to a flat list of SeqAligns. Note that there may
 * be more than one alignment per element in the heap.  The new list
 * preserves the order of the SeqAligns associated with each
 * HeapRecord.
 *
 * @param self           a Blast_CompoHeap
 */
static SeqAlignPtr
Kappa_CompoHeapToFlatList(BlastCompo_Heap * self)
{
    SeqAlignPtr list = NULL;
    SeqAlignPtr result;

    while (NULL != (result = BlastCompo_HeapPop(self))) {
        SeqAlignPtr oldList = list;
        list = result;
        while (NULL != result->next) {
            result = result->next;
        }
        result->next = oldList;
    }
    return list;
}


/**
 *  Top level routine to recompute alignments for each
 *  match found by the gapped BLAST algorithm
 *
 *  @param search           is the structure with all the information about
 *                          the search
 *  @param options          is used to pass certain command line options
 *                          taken in by BLAST
 *  @param hitlist_count    is the number of old matches
 *  @param compo_adjust_mode determines whether we are to adjust the
 *                          Karlin-Altschul parameters and score matrix
 *  @param SmithWaterman    determines whether the new local alignments
 *                          should be computed by the optimal Smith-Waterman
 *                          algorithm; SmithWaterman false means that
 *                          alignments will be recomputed by the current
 *                          X-drop algorithm as implemented in the procedure
 *                          ALIGN.
 *  @return                 a array of lists of SeqAlign; each element
 *                          in the array is a list of SeqAligns for
 *                          one query in the concatenated query.
 *  It is assumed that at least one of compo_adjust_mode and
 *  SmithWaterman is >0 or true when this procedure is called A linked list
 *  of alignments is returned; the alignments are sorted according to
 *  the lowest E-value of the best alignment for each matching
 *  sequence; alignments for the same matching sequence are in the
 *  list consecutively regardless of the E-value of the secondary
 *  alignments. Ties in sorted order are much rarer than for the
 *  standard BLAST method, but are broken deterministically based on
 *  the index of the matching sequences in the database.
 */
SeqAlignPtr *
RedoAlignmentCore(BlastSearchBlkPtr search,
                  BLAST_OptionsBlkPtr options,
                  Int4 hitlist_count,
                  Int4 adjustParameters,
                  Boolean SmithWaterman)
{
  int status = 0;               /* error status returned by routines */
  Int4 match_index;             /* index over matches */
  SeqAlignPtr * results = NULL; /* an array of lists of SeqAligns to
                                   return */
  Nlm_FloatHi localScalingFactor;       /* the factor by which to
                                         * scale the scoring system in
                                         * order to obtain greater
                                         * precision */
  BLAST_Score      **matrix;    /* score matrix */
  Kappa_SearchParameters *searchParams; /* the values of the search
                                         * parameters that will be
                                         * recorded, altered in the
                                         * search structure in this
                                         * routine, and then restored
                                         * before the routine
                                         * exits. */
  Blast_ForbiddenRanges  forbidden;     /* forbidden ranges for each
                                         * database position (used in
                                         * Smith-Waterman alignments)
                                         */
  BlastCompo_Heap * significantMatches; /* a collection of alignments for each
                                * query sequence with sequences from
                                * the database */
  Blast_CompositionWorkspace
      *NRrecord = NULL;        /* stores all fields needed for
                                * computing a compositionally adjusted
                                * score matrix using Newton's method */
  Int4 query_index;            /* loop index */
  Int4 numQueries;             /* number of queries in the
                                  concatenated query */
  Int4 ccat_query_length;      /* length of the concatenated query, or
                                  of the sole query if query
                                  concatenation is not in use */
  Int4 maxQueryLength;         /* the greatest length among all queries */
  BlastCompo_Alignment * incoming_aligns;  /* existing algnments for a match */
  BlastCompo_QueryInfo * query_info = NULL;
  Blast_RedoAlignParams * redo_align_params;
  double Lambda, logK;
  /* the composition adjustment mode, obtained from adjustParameters
     (we don't make adjustParameters itself an enum so that only
     kappa.c depends on the header that defines ECompoAdjustModes) */
  ECompoAdjustModes compo_adjust_mode;  

  /**** Validate parameters *************/
  if (0 > adjustParameters || eNumCompoAdjustModes <= adjustParameters) {
      /* Unknown composition adjustment mode */
      return NULL;
  } else {
      compo_adjust_mode = (ECompoAdjustModes) adjustParameters;
  }
  if(0 == strcmp(options->matrix, "BLOSUM62_20") &&
     eNoCompositionBasedStats == compo_adjust_mode) {
    /* BLOSUM62_20 only makes sense if compo_adjust_mode is on */
    return NULL;
  }
  if (search->positionBased && search->sbp->posMatrix == NULL) {
    Char* msg = "Cannot perform position-specific search without a PSSM";
    BlastConstructErrorMessage("RedoAlignmentCore", msg, 3,
                               &(search->error_return));
    return NULL;
  }
  if (search->positionBased && compo_adjust_mode > 1) {
    /* Only mode 1 (or 0) works with position-based searches, silently 
     * change the value. */
    compo_adjust_mode = eCompositionBasedStats;
  }
  if (compo_adjust_mode > 1 &&
      !Blast_FrequencyDataIsAvailable(options->matrix)) {
    Char* msg = "Unsupported matrix for compostion-based"
      " matrix adjustment";
    BlastConstructErrorMessage("RedoAlignmentCore", msg, 3,
                               &(search->error_return));
    return NULL;
  }
  /**** End validate parameters *************/

  ccat_query_length = search->context[0].query->length;
  searchParams =
      Kappa_SearchParametersNew(ccat_query_length, compo_adjust_mode, 
                                search->positionBased);
  Kappa_RecordInitialSearch(searchParams, search, compo_adjust_mode);
  if (compo_adjust_mode == eNoCompositionBasedStats) {
      localScalingFactor = 1.0;
  } else {
      localScalingFactor = SCALING_FACTOR;
  }
  if((0 == strcmp(options->matrix, "BLOSUM62_20"))) {
    localScalingFactor /= 10;
  }
  Kappa_RescaleSearch(search, localScalingFactor, options);
  if(search->positionBased) {
    matrix = search->sbp->posMatrix;
  } else {
    matrix = search->sbp->matrix;
  }
  redo_align_params = Kappa_GetAlignParams(search, options, compo_adjust_mode,
                                           localScalingFactor);
  {
      /* Get appropriate values for Lambda and logK */
      BLAST_KarlinBlkPtr kbp;
      if(search->positionBased) {
          kbp    = search->sbp->kbp_gap_psi[0];
      } else {
          kbp    = search->sbp->kbp_gap_std[0];
      }
      Lambda = kbp->Lambda;  
      logK   = kbp->logK;
  }
  Kappa_GetQueryInfo(&query_info, &numQueries, &maxQueryLength, search);
  if(SmithWaterman) {
      status = Blast_ForbiddenRangesInitialize(&forbidden, maxQueryLength);
      if (status != 0) 
          ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
  }
  significantMatches = Nlm_Calloc(numQueries, sizeof(BlastCompo_Heap));
  for (query_index = 0;  query_index < numQueries;  query_index++) {
    status =
      BlastCompo_HeapInitialize(&significantMatches[query_index],
                                options->hitlist_size, options->ethresh);
    if (status != 0) {
        ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
    }
  }
  if (compo_adjust_mode != eNoCompositionBasedStats) {
    NRrecord = Blast_CompositionWorkspaceNew();
    if (NRrecord == NULL) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
    }
    /* We already checke that the frequency data is available for
     * options->matrix, so this call can't fail */
    (void) Blast_CompositionWorkspaceInit(NRrecord, options->matrix);
  }
  for(match_index = 0; match_index < hitlist_count; match_index++) {
    /* for all matching sequences */
    BlastCompo_MatchingSequence matchingSeq; /* the data for a matching
                                         * database sequence */
    BLASTResultHitlistPtr thisMatch;    /* alignment data for the
                                         * current query-subject
                                         * match */
    BlastCompo_Alignment ** alignments; /* array of lists of
                                            * alignments for each
                                            * query to this subject */
    alignments = Nlm_Calloc(numQueries, sizeof(BlastCompo_Alignment *));

    thisMatch = search->result_struct->results[match_index];
    if(thisMatch->hsp_array == NULL) {
      continue;
    }
    if (BlastCompo_EarlyTermination(thisMatch->best_evalue,
                                    significantMatches, numQueries)) {
        break;
    }
    /* Get the sequence for this match */
    Kappa_MatchingSequenceInitialize(&matchingSeq, search,
                                     thisMatch->subject_id);
    incoming_aligns =
      Kappa_ResultHspToDistinctAlign(search, thisMatch->hsp_array,
                                     thisMatch->hspcnt, localScalingFactor);
    if (SmithWaterman) {
      status =
        Blast_RedoOneMatchSmithWaterman(alignments, redo_align_params,
                                        incoming_aligns, thisMatch->hspcnt,
                                        Lambda, logK,
                                        &matchingSeq, query_info, numQueries,
                                        matrix, NRrecord, &forbidden,
                                        significantMatches);
    } else {
      status =
        Blast_RedoOneMatch(alignments, redo_align_params, incoming_aligns,
                           thisMatch->hspcnt, Lambda, &matchingSeq,
                           ccat_query_length, query_info, numQueries,
                           matrix, NRrecord);
    }
    if (status != 0) {
      ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
    }
    for (query_index = 0;  query_index < numQueries;  query_index++) {
      /* Loop over queries */
      if( alignments[query_index] != NULL) { /* alignments were found */
        Nlm_FloatHi bestEvalue;   /* best evalue among alignments in the
                                     hitlist */
        Int4 bestScore;           /* best score among alignments in the
                                     hitlist */        
        BLAST_HitListPtr hitlist; /* a hitlist containing the newly-computed
                                   * alignments */
        Kappa_SortedHitlistFromAligns(search, &alignments[query_index]);
        hitlist = search->current_hitlist;
        if (hitlist->hspcnt > 1 &&
            !(SmithWaterman && search->prog_number == blast_type_blastp)) {
          /* Eliminate HSPs that are contained in a higher-scoring
           * HSP.  With blastp and SmithWaterman alignments, the
           * forbidden ranges rule does not allow one alignment to be
           * contained in another. */
          BLASTCheckHSPInclusion(hitlist->hsp_array, hitlist->hspcnt,
                                 FALSE);
          hitlist->hspcnt =
            HspArrayPurge(hitlist->hsp_array, hitlist->hspcnt, FALSE);
        }
        Kappa_HitlistEvaluateAndPurge(&bestScore, &bestEvalue,
                                      matchingSeq.length,
                                      search, redo_align_params->do_link_hsps);
        if (bestEvalue <= search->pbp->cutoff_e &&
            BlastCompo_HeapWouldInsert(&significantMatches[query_index],
                                       bestEvalue, bestScore,
                                       thisMatch->subject_id)) {
          /* If the best alignment is significant, then create and save
           * a list of SeqAligns. */
          /* SeqAligns for this query-subject pair */
          SeqAlignPtr aligns;
          void * discardedAligns;     
          SeqIdPtr query_id;
          Kappa_SequenceLocalData * local_data = matchingSeq.local_data;

          if (search->mult_queries) {
            query_id = search->mult_queries->FakeBsps[query_index]->id;
          } else {
            query_id = search->query_id;
          }
          aligns =
            Kappa_SeqAlignsFromHitlist(hitlist, local_data->bsp_db->id,
                                       query_id,
                                       query_info[query_index].origin,
                                       query_info[query_index].seq.length,
                                       Lambda, logK, localScalingFactor);
          status =
            BlastCompo_HeapInsert(&significantMatches[query_index],
                                  aligns, bestEvalue, bestScore,
                                  thisMatch->subject_id,
                                  &discardedAligns);
          if (status != 0) {
              ErrPostEx(SEV_FATAL, E_NoMemory, 0, "Failed to allocate memory");
          }
          if (discardedAligns != NULL) {
            SeqAlignSetFree(discardedAligns);
          }
        } /* end if the best alignment is significant */
      } /* end if any alignments were found */
    } /* end loop over queries */
    Kappa_MatchingSequenceRelease(&matchingSeq);
    MemFree(alignments);
    BlastCompo_AlignmentsFree(&incoming_aligns, NULL);
  }
  /* end for all matching sequences */
  results = Nlm_Calloc(numQueries, sizeof(SeqAlignPtr));
  for (query_index = 0;  query_index < numQueries;  query_index++) {
    results[query_index] =
      Kappa_CompoHeapToFlatList(&significantMatches[query_index]);
  }
  /* Clean up */
  free(query_info);
  Blast_RedoAlignParamsFree(&redo_align_params);
  if (search->current_hitlist) {
    search->current_hitlist->hspcnt_max = search->current_hitlist->hspcnt;
    search->current_hitlist = BlastHitListDestruct(search->current_hitlist);
  }
  for (query_index = 0;  query_index < numQueries;  query_index++) {
    BlastCompo_HeapRelease(&significantMatches[query_index]);
  }
  MemFree(significantMatches); significantMatches = NULL;
  if(SmithWaterman) {
      Blast_ForbiddenRangesRelease(&forbidden);
  }
  Kappa_RestoreSearch(search, matrix, searchParams, compo_adjust_mode);
  Kappa_SearchParametersFree(&searchParams);
  if (NULL != NRrecord) {
    Blast_CompositionWorkspaceFree(&NRrecord);
  }
  return (results);
}
