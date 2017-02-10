/* $Id: blast_options_api.h,v 1.11 2006/04/26 12:45:28 madden Exp $
***************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
***************************************************************************
* Author: Ilya Dondoshansky
*
*/

/** @file blast_options_api.h
 * Wrapper API around the core options structures used in BLAST.
 */

#ifndef _BLAST_OPTIONS_API_H_
#define _BLAST_OPTIONS_API_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <algo/blast/core/blast_options.h>
#include <algo/blast/api/blast_returns.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

/** Wrapper structure, containing all BLAST options. */
typedef struct SBlastOptions {
    EBlastProgramType program; /**< Type of BLAST program */
    QuerySetUpOptions* query_options; /**< Options for query sequence set up. */
    LookupTableOptions* lookup_options; /**< Lookup table options. */
    BlastInitialWordOptions* word_options; /**< Options for choosing and
                                              extending initial seeds (words) */
    BlastScoringOptions* score_options; /**< Scoring options */
    BlastExtensionOptions* ext_options; /**< Gapped extension options */
    BlastHitSavingOptions* hit_options; /**< Hit saving options */
    BlastEffectiveLengthsOptions* eff_len_options; /**< Effective lengths 
                                                      options */
    PSIBlastOptions* psi_options; /**< Options specific to PSI BLAST */
    BlastDatabaseOptions* db_options; /**< Structure containing database genetic
                                         code */
    int num_cpus; /**< Number of CPUs to use for preliminary stage of the 
                     search. */
    Boolean believe_query; /**< if TRUE then we are using user Query ID. */
} SBlastOptions;

/** Allocates all core options structures and initializes them with default 
 * values.
 * @param program BLAST program name. [in]
 * @param options Resulting structure [out]
 * @param extra_returns Extra returns structure, to put error message in [out]
 */
Int2 SBlastOptionsNew(const char* program, SBlastOptions** options,
                      Blast_SummaryReturn* extra_returns);

/** Frees the SBlastOptions structure and all core options.
 * @param options Structure to free.
 */
SBlastOptions* SBlastOptionsFree(SBlastOptions* options);

/** Sets the expect value threshold option.
 * @param options Options wrapper structure. [in] [out]
 * @param evalue Value to set. [in]
 */
Int2 SBlastOptionsSetEvalue(SBlastOptions* options, double evalue);

/** Sets the word size option.
 * @param options Options wrapper structure. [in] [out]
 * @param word_size Value to set. [in]
 */
Int2 SBlastOptionsSetWordSize(SBlastOptions* options, Int4 word_size);

/** Sets parameters for discontiguous Mega BLAST.
 * @param options Options wrapper structure. [in] [out]
 * @param template_length Discontiguous word template length value to set. [in]
 * @param template_type Discontiguous word template type value to set. 
 *                      (0 or 1) [in]
 */
Int2 SBlastOptionsSetDiscMbParams(SBlastOptions* options, Int4 template_length,
                                 Int4 template_type);

/** Reset matrix name and gap costs to new values.
 * 
 * @param options Options structure to update. [in] [out]
 * @param matrix_name New matrix name [in]
 * @param gap_open New gap existence cost.  If zero default for matrix is used. [in]
 * @param gap_extend New gap extension cost.  If zero default for matrix is used. [in]
 */
Int2 SBlastOptionsSetMatrixAndGapCosts(SBlastOptions* options, 
                                       const char* matrix_name, 
                                       Int4 gap_open, 
                                       Int4 gap_extend);


/** Reset rewared, penalty and gap costs to new values.
 *  Will suggest and use conservative values if gap_open and gap_extend are zero
 *  and suggest is TRUE.
 * 
 * @param options Options structure to update. [in] [out]
 * @param reward match score [in]
 * @param penalty mismatch score [in]
 * @param gap_open New gap existence cost.  If zero default for reward/penalty is used. [in]
 * @param gap_extend New gap extension cost.  If zero default for reward/penalty is used. [in]
 * @param greedy TRUE specifies that greedy algorithms will be used. [in]
 */
Int2 SBlastOptionsSetRewardPenaltyAndGapCosts(SBlastOptions* options, 
                                       Int4 reward, Int4 penalty,
                                       Int4 gap_open, Int4 gap_extend,
                                       Boolean greedy);

/** Set threshold value.
 * @param options options Options structure to update. [in] [out]
 * @param threshold New value to set, if zero default value for matrix
 *    will be used. [in]
 * @return zero unless error (e.g., threshold is < zero)
 */
Int2 SBlastOptionsSetThreshold(SBlastOptions* options, 
                               Int4 threshold);

/** Set window size for two hit extension.
 * @param options options Options structure to update. [in] [out]
 * @param window_size New value to set, if zero default value for matrix
 *    will be used. [in]
 * @return zero unless error (e.g., window_size is < zero)
 */
Int2 SBlastOptionsSetWindowSize(SBlastOptions* options, 
                               Int4 window_size);

/** Reset database (subject) genetic code option to a new value. 
 * @param options Options structure to update. [in] [out]
 * @param gc New genetic code value. [in]
 * @return Status: non-zero if genetic code string cannotbe found for gc.
 */
Int2 SBlastOptionsSetDbGeneticCode(SBlastOptions* options, Int4 gc);

/** Reset the filter string option.
 * @param options Options structure to update [in] [out]
 * @param str New filter string. [in]
 * @return Status.
 */
Int2 SBlastOptionsSetFilterString(SBlastOptions* options, const char* str);

/** Returns the mask-at-hash option value.
 * @param options The options structure [in]
 * @return Boolean value of the masking at hash option.
 */
Boolean SBlastOptionsGetMaskAtHash(const SBlastOptions* options);


/** sets believe_query flag on SBlastOptions.
 * @param options Object to be modified [in]
 * @param believe_query specifies that query ID was parsed [in]
 * @return zero on success.
 */
Int2 SBlastOptionsSetBelieveQuery(SBlastOptions* options, Boolean believe_query);

/** Gets believe_query field of SBlastOptions.
 * @param options Object to be queried for information [in]
 * @return TRUE if query ID was parsed, otherwise FALSE.
 */
Boolean SBlastOptionsGetBelieveQuery(const SBlastOptions* options);

/* @} */

#ifdef __cplusplus
}
#endif

#endif  /* !_BLAST_OPTIONS_API_H_ */
