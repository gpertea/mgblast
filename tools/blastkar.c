static char const rcsid[] = "$Id: blastkar.c,v 6.113 2006/04/07 13:46:24 madden Exp $";

/* ===========================================================================
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

File name: blastkar.c

Author: Tom Madden

Contents: Functions to calculate BLAST probabilities etc.

Detailed Contents:

	- allocate and deallocate structures used by BLAST to calculate
	probabilities etc.

	- calculate residue frequencies for query and "average" database.

	- read in matrix.

        - calculate sum-p from a collection of HSP's, for both the case
	of a "small" gap and a "large" gap, when give a total score and the
	number of HSP's.

	- calculate expect values for p-values.

	- calculate pseuod-scores from p-values.

****************************************************************************** 
 * $Revision: 6.113 $
 * $Log: blastkar.c,v $
 * Revision 6.113  2006/04/07 13:46:24  madden
 * Improved the comment for NlmKarlinLambdaNR.  Reformatted the
 * function prototype to fit in 80 characters. (from Mike Gertz).
 *
 * Revision 6.112  2005/10/31 14:16:10  madden
 * Add support for reward/penalty of 1/-5, 3/-4, and 3/-2
 *
 * Revision 6.111  2005/10/14 16:32:34  madden
 * Add preliminary support for vecscreen parameters
 *
 * Revision 6.110  2005/10/12 19:21:03  madden
 * Fix bug in s_GetNuclValuesArray
 *
 * Revision 6.109  2005/10/06 12:51:39  madden
 * Add code to produce correct statistics for gapped blastn, changes include:
 *
 * 1.) series of array_of_8 structs with precalculated data
 * 2.) function s_GetNuclValuesArray
 * 3.) function BlastKarlinGetNuclAlphaBeta
 * 4.) function BlastKarlinBlkNuclGappedCalc
 *
 * Revision 6.108  2005/08/09 14:14:46  dondosha
 * From A. Shaffer: added comments to clarify usage of BlastKarlinEtoP and BlastKarlinPtoE
 *
 * Revision 6.107  2005/07/28 14:57:09  coulouri
 * remove dead code
 *
 * Revision 6.106  2005/07/27 17:48:57  coulouri
 * remove hardcoded paths
 *
 * Revision 6.105  2005/04/27 17:20:45  papadopo
 * copy X scores to U scores when building score matrix
 *
 * Revision 6.104  2005/04/20 19:02:15  lavr
 * +<assert.h>
 *
 * Revision 6.103  2005/01/31 17:02:03  dondosha
 * Fixed bug in BlastKarlinLHtoK for blastn with penalty == -reward
 *
 * Revision 6.102  2004/09/28 16:04:19  papadopo
 * From Michael Gertz:
 * 1. Pass the effective database size into BlastSmallGapSumE,
 *     BlastLargeGapSumE and BlastUnevenGapSumE.  The routines use this
 *     value in a simplified formula to compute the e-value of singleton sets.
 * 2. Caused all routines for calculating the significance of multiple
 *     distinct alignments (BlastSmallGapSumE, BlastLargeGapSumE and
 *     BlastUnevenGapSumE) to use
 *
 *        sum_{i in linked_set} (\lambda_i s_i - \ln K_i)
 *
 *     as the weighted sum score, where (\lambda_i, K_i) are taken from
 *     the appropriate query context.
 *
 * Revision 6.101  2004/06/07 20:03:23  coulouri
 * use floating point constants for comparisons with floating point variables
 *
 * Revision 6.100  2004/04/28 14:36:00  madden
 * Changes from Mike Gertz:
 * - I created the new routine BlastGapDecayDivisor that computes a
 *   divisor used to weight the evalue of a collection of distinct
 *   alignments.
 * - I removed  BlastGapDecay and BlastGapDecayInverse which had become
 *   redundant.
 * - I modified the BlastCutoffs routine so that it uses the value
 *   returned by BlastGapDecayDivisor to weight evalues.
 * - I modified BlastSmallGapSumE, BlastLargeGapSumE and
 *   BlastUnevenGapSumE no longer refer to the gap_prob parameter.
 *   Replaced the gap_decay_rate parameter of each of these routines with
 *   a weight_divisor parameter.  Added documentation.
 *
 * Revision 6.99  2004/04/23 13:49:43  madden
 * Cleaned up ifndef in BlastKarlinLHtoK
 *
 * Revision 6.98  2004/04/23 13:19:53  madden
 * Rewrote BlastKarlinLHtoK to do the following and more:
 * 1. fix a bug whereby the wrong formula was used when high score == 1
 *    and low score == -1;
 * 2. fix a methodological error of truncating the first sum
 *    and trying to make it converge quickly by adding terms
 *    of a geometric progression, even though the geometric progression
 *    estimate is not correct in all cases;
 *    the old adjustment code is left in for historical purposes but
 *    #ifdef'd out
 * 3. Eliminate the Boolean bi_modal_score variable.  The old test that
 *    set the value of bi_modal_score would frequently fail to choose the
 *    correct value due to rounding error.
 * 4. changed numerous local variable names to make them more meaningful;
 * 5. added substantial comments to explain what the procedure
 *    is doing and what each variable represents
 *
 * Revision 6.97  2004/04/01 13:43:08  lavr
 * Spell "occurred", "occurrence", and "occurring"
 *
 * Revision 6.96  2004/03/31 17:58:51  papadopo
 * Mike Gertz' changes for length adjustment calculations
 *
 * Revision 6.95  2003/12/12 16:00:34  madden
 * Add gap_decay_rate to BlastCutoffs, remove BlastCutoffs_simple, protection against overflow, removal of defunct _real variables (all from Mike Gertz)
 *
 * Revision 6.94  2003/11/30 03:36:38  camacho
 * Fix compilation error
 *
 * Revision 6.93  2003/11/28 22:39:40  camacho
 * +static keyword to BlastKarlinLtoH
 *
 * Revision 6.92  2003/11/28 15:16:38  camacho
 * Combine newkar.c's contents with blastkar.c
 *
 * Revision 6.91  2003/11/26 19:08:10  madden
 * simplified BlastKarlinLtoH and BlastKarlinLHtoK and provided better protection against overflow, new function NlmKarlinLambdaNR (all from Mike Gertz)
 *
 *
 * Revision 6.90  2003/06/30 20:01:32  dondosha
 * Correction in logic of finding matrices by BLASTMAT environment variable
 *
 * Revision 6.89  2003/05/30 17:25:36  coulouri
 * add rcsid
 *
 * Revision 6.88  2003/05/06 18:54:02  dondosha
 * Set all gap/residue scores in blastn matrix to INT4_MIN/2
 *
 * Revision 6.87  2003/03/07 22:33:25  bealer
 * - Fix UMR when alphabet_type is not set.
 *
 * Revision 6.86  2003/02/27 19:07:56  madden
 * Add functions PrintMatrixMessage and PrintAllowedValuesMessage
 *
 * Revision 6.85  2003/02/26 18:23:49  madden
 * Add functions BlastKarlinkGapBlkFill and BlastKarlinReportAllowedValues, call from BlastKarlinBlkGappedCalcEx
 *
 * Revision 6.84  2002/10/24 22:52:14  dondosha
 * When checking config file for matrices path, allow aa or nt subdirectories too
 *
 * Revision 6.83  2002/07/22 20:10:12  dondosha
 * Correction: previous change did not work for proteins
 *
 * Revision 6.82  2002/07/19 18:34:58  dondosha
 * Ignore bits higher than 4 when computing frequencies - needed for megablast
 *
 * Revision 6.81  2002/05/17 20:30:37  madden
 * Add comments on adding new matrix values
 *
 * Revision 6.80  2002/04/09 18:44:19  madden
 * Do not return if status of BlastScoreBlkMatCreate is non-zero
 *
 * Revision 6.79  2002/02/26 22:25:21  dondosha
 * Return error as soon as it is found that matrix name is not supported
 *
 * Revision 6.78  2001/12/13 14:30:49  madden
 * Add BLASTKAR_SMALL_FLOAT to prevent floating point exception for very small floats
 *
 * Revision 6.77  2001/09/05 20:32:21  dondosha
 * Fixed uninitialized variable bug
 *
 * Revision 6.76  2001/08/23 21:19:05  dondosha
 * Improvements for lambda and K computation when all scores are multiples of a common factor
 *
 * Revision 6.75  2001/02/20 18:31:28  egorov
 * Added protection agains freeing zero pointer
 *
 * Revision 6.74  2001/01/29 16:11:34  madden
 * Added BLOSUM80 values for 25,2
 *
 * Revision 6.73  2000/12/28 16:23:24  madden
 * Function getAlphaBeta from AS
 *
 * Revision 6.72  2000/12/26 17:46:20  madden
 * Add function BlastKarlinGetMatrixValuesEx2 to return alpha and beta
 *
 * Revision 6.71  2000/11/27 16:04:34  dondosha
 * Check if original_matrix not NULL before destructing it
 *
 * Revision 6.70  2000/11/24 22:07:32  shavirin
 * Added new function BlastResFreqFree().
 *
 * Revision 6.69  2000/11/24 21:44:21  shavirin
 * Fixed memory leak in the function BLAST_MatrixDestruct() in case of
 * PSI Blast.
 *
 * Revision 6.68  2000/11/22 15:32:39  madden
 * Remove unneeded line
 *
 * Revision 6.67  2000/11/21 19:06:03  madden
 * Set best value for blosum62_20
 *
 * Revision 6.66  2000/11/21 18:45:56  madden
 * New parameter values from S. Altschul for matrices
 *
 * Revision 6.65  2000/11/03 17:15:13  madden
 * Add 13,2 for blosum80
 *
 * Revision 6.64  2000/10/25 16:39:46  madden
 * Add protection in BlastMatrixToTxMatrix for NULL matrix
 *
 * Revision 6.63  2000/10/23 15:52:18  dondosha
 * Bug in previous revision: INT4_MIN is too small to be a matrix value
 *
 * Revision 6.62  2000/10/20 21:51:09  dondosha
 * Set matrix[15][] values to INT4_MIN for blastn to prevent crossing boundary between strands
 *
 * Revision 6.61  2000/09/27 21:28:40  dondosha
 * Copy square matrix in BLAST_MatrixFill to the original_matrix member of BLAST_Matrix
 *
 * Revision 6.60  2000/09/27 19:29:04  dondosha
 * Check boolean parameter in BLAST_MatrixFill to use position based matrix
 *
 * Revision 6.59  2000/09/15 18:48:16  madden
 * Added new blosum62_20 values
 *
 * Revision 6.58  2000/09/13 16:16:52  dondosha
 * Added LIBCALL to TxMatrix functions declarations
 *
 * Revision 6.57  2000/09/13 15:50:28  dondosha
 * Use SeqMapTable functions for conversion of search matrix to txalign-style matrix
 *
 * Revision 6.56  2000/09/12 16:03:50  dondosha
 * Added functions to create and destroy matrix used in txalign
 *
 * Revision 6.55  2000/09/11 20:49:32  madden
 * New values for BLOSUM62_20 from Stephen Altschul
 *
 * Revision 6.54  2000/08/31 15:45:07  dondosha
 * Added function BlastUnevenGapSumE for sum evalue computation with different gap size on two sequences
 *
 * Revision 6.53  2000/08/28 13:38:31  shavirin
 * Fixed bug in the function BLAST_MatrixFill() detected by Alejandro.
 *
 * Revision 6.52  2000/08/23 18:50:01  madden
 * Changes to support decline-to-align checking
 *
 * Revision 6.51  2000/08/22 12:58:55  madden
 * Add new BLOSUM62_20 values
 *
 * Revision 6.50  2000/08/04 15:49:25  sicotte
 * added BlastScoreBlkMatCreateEx(reward,penalty) and BlastKarlinGetDefaultMatrixValues as external functions
 *
 * Revision 6.49  2000/08/04 15:13:29  dondosha
 * Initialize posFreqs to NULL in BLAST_MatrixFill
 *
 * Revision 6.48  2000/08/04 14:43:58  shavirin
 * Changed values for data corresponding to BLOSUM62_20A and
 * BLOSUM62_20B matrixes.
 *
 * Revision 6.47  2000/08/03 22:23:13  shavirin
 * Added initialization of the posFreqs array in the function BLAST_MatrixFill()
 *
 * Revision 6.46  2000/07/24 18:46:37  shavirin
 * Added stat parameters for the matrixes BLOSUM62_20a and BLOSUM62_20b
 *
 * Revision 6.45  2000/07/24 17:32:42  shavirin
 * Fixed values for size of matrix in the function BlastScoreBlkMaxScoreSet()
 *
 * Revision 6.44  2000/07/20 20:46:30  madden
 * New values for blosum62_20
 *
 * Revision 6.43  2000/06/12 21:37:33  shavirin
 * Adjusted blosum62_values{} and blosum80_values{}.
 *
 * Revision 6.42  2000/06/02 16:20:54  shavirin
 * Fixed minor bug in function BLAST_MatrixFill
 *
 * Revision 6.41  2000/05/26 17:29:54  shavirin
 * Added array of pos frequencies into BLAST_Matrix structure and it's
 * handling.
 *
 * Revision 6.40  2000/04/17 20:41:37  madden
 * Added BLAST_MatrixFetch
 *
 * Revision 6.39  2000/03/27 20:28:19  shavirin
 * Added mulithred-safe protection to the function BlastScoreBlkMatRead()
 *
 * Revision 6.38  2000/02/24 16:39:03  shavirin
 * Added check for existence of matrix file in BlastScoreBlkMatFill().
 *
 * Revision 6.37  2000/01/11 21:23:13  shavirin
 * Increased size of index from Int2 to Int4 in BlastPSIMaxScoreGet() function
 *
 * Revision 6.36  1999/12/22 21:06:34  shavirin
 * Added new function BlastPSIMaxScoreGet().
 *
 * Revision 6.35  1999/12/16 19:17:22  egorov
 * Code cleanup
 *
 * Revision 6.34  1999/11/29 14:17:46  egorov
 * Use LnFactorial instead of log(Factorial)
 *
 * Revision 6.33  1999/11/23 21:38:40  egorov
 * Nlm_Factorial(num) causes overflow if num>170;  return DBL_MAX if num>170.
 * It does not influence program's logic because this large number of HSPs is bogus anyway.
 *
 * Revision 6.32  1999/08/20 14:42:17  madden
 * Changed Robinson frequencies per Stephen A. request
 *
 * Revision 6.31  1999/08/05 13:13:34  madden
 * Improved help messages for permitted matrices and gap values
 *
 * Revision 6.30  1999/07/30 13:25:51  shavirin
 * Fixed bug in the function BLAST_MatrixFill which created wrong matrix size.
 *
 * Revision 6.29  1998/12/31 18:17:04  madden
 * Added strand option
 *
 * Revision 6.28  1998/09/28 12:28:50  madden
 * Protection for a DEC Alpha problem with HUGE_VAL
 *
 * Revision 6.27  1998/09/24 15:26:36  egorov
 * Fix lint complaints
 *
 * Revision 6.26  1998/08/11 12:57:23  madden
 * Changed importance of BlastKarlinBlkGappedCalc error message
 *
 * Revision 6.25  1998/07/17 15:39:57  madden
 * Changes for Effective search space.
 *
 * Revision 6.24  1998/06/12 15:52:50  madden
 * Fixed warnings
 *
 * Revision 6.23  1998/06/02 21:21:20  madden
 * Changes for DNA matrices
 *
 * Revision 6.22  1998/05/03 17:56:45  madden
 * Fixed memory leaks
 *
 * Revision 6.21  1998/04/24 21:51:40  madden
 * Check return value on BlastKarlinBlkCalc
 *
 * Revision 6.20  1998/04/24 19:27:54  madden
 * Added BlastKarlinBlkStandardCalcEx for ideal KarlinBlk
 *
 * Revision 6.19  1998/04/10 20:57:55  madden
 * Added data for BLOSUM80 and BLOSUM45
 *
 * Revision 6.18  1998/04/10 15:05:48  madden
 * Added pref_flags return by value to BlastKarlinGetMatrixValues
 *
 * Revision 6.17  1998/04/02 21:12:29  madden
 * Added infinite values to the arrays returned by BlastKarlinGetMatrixValues
 *
 * Revision 6.16  1998/03/16 17:41:56  madden
 * Fixed leaks
 *
 * Revision 6.15  1998/03/09 17:15:01  madden
 * Added BlastKarlinGetMatrixValues
 *
 * Revision 6.14  1998/02/26 22:34:32  madden
 * Changes for 16 bit windows
 *
 * Revision 6.13  1998/02/06 18:28:05  madden
 * Added functions to produce pseudo-scores from p and e values
 *
 * Revision 6.12  1997/12/15 21:55:44  madden
 * Added new parameters for BLOSUM62_20
 *
 * Revision 6.11  1997/12/11 22:21:02  madden
 * Removed unused variables
 *
 * Revision 6.10  1997/12/10 22:41:58  madden
 * proper casting done
 *
 * Revision 6.9  1997/12/01 22:07:55  madden
 * Added missing values, fixed UMR
 *
 * Revision 6.8  1997/11/14 21:29:58  madden
 * Stephens new parameters
 *
 * Revision 6.7  1997/11/14 17:14:52  madden
 * Added Karlin parameter to matrix structure
 *
 * Revision 6.6  1997/11/07 22:27:29  madden
 * Fix in BLAST_MatrixFill
 *
 * Revision 6.5  1997/11/07 00:48:07  madden
 * Added defintitions and functions for BLAST_Matrix
 *
 * Revision 6.4  1997/10/30 15:41:01  madden
 * Casts and fixes for DEC alpha
 *
 * Revision 6.3  1997/10/22 21:46:59  madden
 * Added BLOSUM90, changed to better BLOSUM62 values
 *
 * Revision 6.2  1997/10/09 19:47:35  madden
 * changes for 1/20th bit BLOSUM62 matrices
 *
 * Revision 6.1  1997/10/06 17:59:17  madden
 * Added checks to BLAST_ScoreBlkDestruct
 *
 * Revision 6.0  1997/08/25 18:52:37  madden
 * Revision changed to 6.0
 *
 * Revision 1.39  1997/07/18 20:56:13  madden
 * Added more BLOSUM50 values
 *
 * Revision 1.38  1997/07/14 20:11:17  madden
 * Removed unused variables
 *
 * Revision 1.37  1997/07/14 15:30:59  madden
 * Changed call to BlastKarlinBlkGappedCalc
 *
 * Revision 1.36  1997/06/23 20:48:37  madden
 * Added support for BLOSUM50 and PAM250 matrices
 *
 * Revision 1.35  1997/06/23 17:53:07  madden
 * BlastKarlinBlkGappedCalc checks for valid values if KarlinBlkPtr NULL
 *
 * Revision 1.34  1997/05/01 15:53:13  madden
 * Addition of extra KarlinBlk's for psi-blast
 *
 * Revision 1.33  1997/04/22  16:36:49  madden
 * Added Robinson amino acid frequencies.
 *
 * Revision 1.32  1997/04/03  19:48:13  madden
 * Changes to use effective database length instead of the length of each
 * sequence in statistical calculations.
 *
 * Revision 1.31  1997/03/07  22:24:43  madden
 * Closed File* for matrix after use.
 *
 * Revision 1.30  1997/02/20  21:33:51  madden
 * Added FindPath call to blastkar.c to look for matrix.
 *
 * Revision 1.29  1997/02/06  23:15:43  madden
 * Added additional gap extension and opening penalties.
 *
 * Revision 1.28  1997/02/05  19:54:59  madden
 * *** empty log message ***
 *
 * Revision 1.27  1997/01/16  22:58:08  madden
 * Secured BlastScoreFreqDestruct against NULL pointers.
 *
 * Revision 1.26  1997/01/13  17:15:29  madden
 * Save matrix name for blastn type matrices.
 *
 * Revision 1.25  1996/12/11  17:59:42  madden
 * Fixed purify leaks.
 *
 * Revision 1.24  1996/12/10  17:30:59  madden
 * Changed statistics for gapped blastp
 *
 * Revision 1.23  1996/12/03  19:13:47  madden
 * Made BlastRes functions non-static.
 *
 * Revision 1.22  1996/11/27  16:39:11  madden
 * Added NULLB to end of array, purify nit.
 *
 * Revision 1.21  1996/11/26  20:38:01  madden
 * Removed unused variable.
 *
 * Revision 1.20  1996/11/26  19:55:25  madden
 * Added check for matrices in standard places.
 *
 * Revision 1.19  1996/11/18  19:32:09  madden
 * Fixed uninitialized variable found by CodeWarrior.
 *
 * Revision 1.18  1996/11/18  14:49:26  madden
 * Rewrote BlastScoreSetAmbigRes to take multiple ambig. chars.
 *
 * Revision 1.17  1996/11/08  21:45:03  madden
 * Fix for blastn matrices.
 *
 * Revision 1.16  1996/11/07  22:31:15  madden
 * Corrected Nucl. matrix for ambiguity characters.
 *
 * Revision 1.15  1996/10/04  20:12:26  madden
 * Fixed memory leaks found by purify.
 *
 * Revision 1.14  1996/10/03  20:49:29  madden
 * Added function BlastKarlinBlkStandardCalc to calculate standard
 * Karlin parameters for blastx and tblastx.
 *
 * Revision 1.13  1996/10/03  18:04:48  madden
 * Each context (or frame) now uses "proper" Karlin parameters.
 *
 * Revision 1.12  1996/10/03  13:07:55  madden
 * Added return statement.
 *
 * Revision 1.11  1996/10/02  20:00:38  madden
 * Calc. different Karlin parameters for each frame.
 *
 * Revision 1.10  1996/10/01  21:24:02  madden
 * Corrected statistics for partial blastn matching.
 *
 * Revision 1.9  1996/09/30  21:56:12  madden
 * Replaced query alphabet of ncbi2na with blastna alphabet.
 *
 * Revision 1.8  1996/09/11  20:08:55  madden
 * *** empty log message ***
 *
 * Revision 1.6  1996/09/11  19:14:09  madden
 * Changes to blastn matrix.
 *
 * Revision 1.5  1996/09/10  19:40:35  madden
 * Added function BlastScoreBlkMatCreate for blastn matrices.
 *
 * Revision 1.4  1996/08/23  19:07:01  madden
 * Adjust changes for NT warning to give correct results.
 *
 * Revision 1.3  1996/08/23  15:32:30  shavirin
 * Fixed a lot of NT compiler warnings about type mismatch
 *
 * Revision 1.2  1996/08/22  20:38:11  madden
 * Changed all doubles and floats to Nlm_FloatHi.
 *
 * Revision 1.1  1996/08/05  19:47:42  madden
 * Initial revision
 *
 * Revision 1.26  1996/06/21  15:23:54  madden
 * Corelibed BlastSumP.
 *
 * Revision 1.25  1996/06/21  15:14:55  madden
 * Made some functions static, added LIBCALL to others.
 *
 * Revision 1.24  1996/06/20  16:52:23  madden
 * Changed "pow" to "Nlm_Powi".
 *
 * Revision 1.23  1996/06/11  15:42:46  madden
 * Replaced strncasempc with toolbox function.
 *
 * Revision 1.22  1996/05/29  12:44:04  madden
 * Added function BlastRepresentativeResidues.
 *
 * Revision 1.21  1996/05/22  20:38:05  madden
 * *** empty log message ***
 *
 * Revision 1.20  1996/05/22  20:21:23  madden
 * removed unused variables.
 *
 * Revision 1.19  1996/05/20  21:18:51  madden
 * Added BLASTMatrixStructure.
 *
 * Revision 1.18  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.17  1996/05/14  21:30:37  madden
 * *** empty log message ***
 *
 * Revision 1.16  1996/04/16  15:33:41  madden
 * Changes for backward-compatability with old blast.
 *
 * Revision 1.15  1996/03/25  16:34:19  madden
 * Changes to mimic old statistics.
 *
 * Revision 1.14  1996/03/01  19:40:37  madden
 * Added factorial correction to SmallGapSum.
 *
 * Revision 1.13  1996/02/28  21:38:08  madden
 * *** empty log message ***
 *
 * Revision 1.12  1996/02/15  23:20:01  madden
 * Added query length to a number of calls.
 *
 * Revision 1.12  1996/02/15  23:20:01  madden
 * Added query length to a number of calls.
 *
 * Revision 1.11  1996/02/09  13:51:08  madden
 * Change to prevent log sign error.
 *
 * Revision 1.10  1996/01/22  22:33:06  madden
 * Changed effective length, fixed (improper) calculation of sum_e.
 *
 * Revision 1.9  1996/01/17  23:18:41  madden
 * Added BlastScoreSetAmbigRes.
 * ci -l blastkar.h
 *
 * Revision 1.8  1996/01/17  17:01:19  madden
 * Fixed BlastSmallGapSumE  and BlastLargeGapSumE.
 *
 * Revision 1.7  1996/01/17  13:46:36  madden
 * Added function BlastKarlinPtoE.
 *
 * Revision 1.6  1996/01/06  18:58:07  madden
 * Added BlastSumP.
 *
 * Revision 1.5  1995/12/28  21:26:16  madden
 * *** empty log message ***
 *
 * Revision 1.4  1995/12/26  23:04:46  madden
 * removed GetBLAST0Matrix, Added "cutoff" functions.
 *
 * Revision 1.3  1995/12/26  20:27:47  madden
 * Replaced BLAST_ScoreMat wiht BLAST_ScorePtr PNTR.
 *
 * Revision 1.2  1995/12/26  14:25:56  madden
 * *** empty log message ***
 *
 * Revision 1.1  1995/12/21  23:07:35  madden
 * Initial revision
 *
 * */
#include <ncbi.h>
#include <assert.h>
#include <ncbimath.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <blastkar.h>
#include <blastpri.h>

static Int2 BlastScoreFreqCalc PROTO((BLAST_ScoreBlkPtr sbp, BLAST_ScoreFreqPtr sfp, BLAST_ResFreqPtr rfp1, BLAST_ResFreqPtr rfp2));

/* OSF1 apparently doesn't like this. */
#if defined(HUGE_VAL) && !defined(OS_UNIX_OSF1)
#define BLASTKAR_HUGE_VAL HUGE_VAL
#else
#define BLASTKAR_HUGE_VAL 1.e30
#endif


/* Allocates and Deallocates the two-dimensional matrix. */
static BLASTMatrixStructurePtr BlastMatrixAllocate PROTO((Int2 alphabet_size));
static BLASTMatrixStructurePtr BlastMatrixDestruct PROTO((BLASTMatrixStructurePtr matrix_struct));

/* performs sump calculation, used by BlastSumPStd */
static Nlm_FloatHi BlastSumPCalc PROTO((int r, Nlm_FloatHi s));

#define COMMENT_CHR	'#'
#define TOKSTR	" \t\n\r"
#define BLAST_MAX_ALPHABET 40 /* ncbistdaa is only 26, this should be enough */

/*
	How many of the first bases are not ambiguous
	(it's four, of course).
*/
#define NUMBER_NON_AMBIG_BP 4
/*
	translates between ncbi4na and blastna.  the first four elements
	of this array match ncbi2na.
*/

Uint1 ncbi4na_to_blastna[] = {
15,/* Gap, 0 */
0, /* A,   1 */
1, /* C,   2 */
6, /* M,   3 */
2, /* G,   4 */
4, /* R,   5 */
9, /* S,   6 */
13, /* V,   7 */
3, /* T,   8 */
8, /* W,   9 */
 5, /* Y,  10 */
 12, /* H,  11 */
 7, /* K,  12 */
 11, /* D,  13 */
 10, /* B,  14 */
 14  /* N,  15 */
};

Uint1 blastna_to_ncbi4na[] = {	 1, /* A, 0 */
				 2, /* C, 1 */
				 4, /* G, 2 */
				 8, /* T, 3 */
				 5, /* R, 4 */
				10, /* Y, 5 */
				 3, /* M, 6 */
				12, /* K, 7 */
				 9, /* W, 8 */
				 6, /* S, 9 */
				14, /* B, 10 */
				13, /* D, 11 */
				11, /* H, 12 */
				 7, /* V, 13 */
				15, /* N, 14 */
				 0, /* Gap, 15 */
			};

/* Used in BlastKarlinBlkGappedCalc */
typedef FloatHi array_of_8[8];

/* Used to temporarily store matrix values for retrieval. */

typedef struct _matrix_info {
	CharPtr		name;			/* name of matrix (e.g., BLOSUM90). */
	array_of_8 	*values;		/* The values (below). */
	Int4		*prefs;			/* Preferences for display. */
	Int4		max_number_values;	/* number of values (e.g., BLOSUM90_VALUES_MAX). */
} MatrixInfo, PNTR MatrixInfoPtr;


/**************************************************************************************

How the statistical parameters for the matrices are stored:
-----------------------------------------------------------
They parameters are stored in a two-dimensional array FloatHi (i.e., 
doubles, which has as it's first dimensions the number of different 
gap existence and extension combinations and as it's second dimension 8.
The eight different columns specify:

1.) gap existence penalty (INT2_MAX denotes infinite).
2.) gap extension penalty (INT2_MAX denotes infinite).
3.) decline to align penalty (INT2_MAX denotes infinite).
4.) Lambda
5.) K
6.) H
7.) alpha
8.) beta

(Items 4-8 are explained in:
Altschul SF, Bundschuh R, Olsen R, Hwa T.
The estimation of statistical parameters for local alignment score distributions.
Nucleic Acids Res. 2001 Jan 15;29(2):351-61.).

Take BLOSUM45 (below) as an example.  Currently (5/17/02) there are
14 different allowed combinations (specified by "#define BLOSUM45_VALUES_MAX 14").
The first row in the array "blosum45_values" has INT2_MAX (i.e., 32767) for gap 
existence, extension, and decline-to-align penalties.  For all practical purposes
this value is large enough to be infinite, so the alignments will be ungapped.
BLAST may also use this value (INT2_MAX) as a signal to skip gapping, so a
different value should not be used if the intent is to have gapless extensions.
The next row is for the gap existence penalty 13 and the extension penalty 3.
The decline-to-align penalty is only supported in a few cases, so it is normally
set to INT2_MAX.


How to add a new matrix to blastkar.c:
--------------------------------------
To add a new matrix to blastkar.c it is necessary to complete 
four steps.  As an example consider adding the matrix
called TESTMATRIX

1.) add a define specifying how many different existence and extensions
penalties are allowed, so it would be necessary to add the line:

#define TESTMATRIX_VALUES_MAX 14

if 14 values were to be allowed.

2.) add a two-dimensional array to contain the statistical parameters:

static Nlm_FloatHi testmatrix_values[TESTMATRIX_VALUES_MAX][8] ={ ...

3.) add a "prefs" array that should hint about the "optimal" 
gap existence and extension penalties:

static Int4 testmatrix_prefs[TESTMATRIX_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
...
};

4.) Go to the function BlastLoadMatrixValues (in this file) and
add two lines before the return at the end of the function: 

        matrix_info = MatrixInfoNew("TESTMATRIX", testmatrix_values, testmatrix_prefs, TESTMATRIX_VALUES_MAX);
        ValNodeAddPointer(&retval, 0, matrix_info);



***************************************************************************************/


	
	

#define BLOSUM45_VALUES_MAX 14
static Nlm_FloatHi  blosum45_values[BLOSUM45_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7},
    {13, 3, (Nlm_FloatHi) INT2_MAX, 0.207, 0.049, 0.14, 1.5, -22},
    {12, 3, (Nlm_FloatHi) INT2_MAX, 0.199, 0.039, 0.11, 1.8, -34},
    {11, 3, (Nlm_FloatHi) INT2_MAX, 0.190, 0.031, 0.095, 2.0, -38},
    {10, 3, (Nlm_FloatHi) INT2_MAX, 0.179, 0.023, 0.075, 2.4, -51},
    {16, 2, (Nlm_FloatHi) INT2_MAX, 0.210, 0.051, 0.14, 1.5, -24},
    {15, 2, (Nlm_FloatHi) INT2_MAX, 0.203, 0.041, 0.12, 1.7, -31},
    {14, 2, (Nlm_FloatHi) INT2_MAX, 0.195, 0.032, 0.10, 1.9, -36},
    {13, 2, (Nlm_FloatHi) INT2_MAX, 0.185, 0.024, 0.084, 2.2, -45},
    {12, 2, (Nlm_FloatHi) INT2_MAX, 0.171, 0.016, 0.061, 2.8, -65},
    {19, 1, (Nlm_FloatHi) INT2_MAX, 0.205, 0.040, 0.11, 1.9, -43},
    {18, 1, (Nlm_FloatHi) INT2_MAX, 0.198, 0.032, 0.10, 2.0, -43},
    {17, 1, (Nlm_FloatHi) INT2_MAX, 0.189, 0.024, 0.079, 2.4, -57},
    {16, 1, (Nlm_FloatHi) INT2_MAX, 0.176, 0.016, 0.063, 2.8, -67},
};

static Int4 blosum45_prefs[BLOSUM45_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_BEST,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};


#define BLOSUM50_VALUES_MAX 16
static Nlm_FloatHi  blosum50_values[BLOSUM50_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.2318, 0.112, 0.3362, 0.6895, -4.0},
    {13, 3, (Nlm_FloatHi) INT2_MAX, 0.212, 0.063, 0.19, 1.1, -16},
    {12, 3, (Nlm_FloatHi) INT2_MAX, 0.206, 0.055, 0.17, 1.2, -18},
    {11, 3, (Nlm_FloatHi) INT2_MAX, 0.197, 0.042, 0.14, 1.4, -25},
    {10, 3, (Nlm_FloatHi) INT2_MAX, 0.186, 0.031, 0.11, 1.7, -34},
    {9, 3, (Nlm_FloatHi) INT2_MAX, 0.172, 0.022, 0.082, 2.1, -48},
    {16, 2, (Nlm_FloatHi) INT2_MAX, 0.215, 0.066, 0.20, 1.05, -15},
    {15, 2, (Nlm_FloatHi) INT2_MAX, 0.210, 0.058, 0.17, 1.2, -20},
    {14, 2, (Nlm_FloatHi) INT2_MAX, 0.202, 0.045, 0.14, 1.4, -27},
    {13, 2, (Nlm_FloatHi) INT2_MAX, 0.193, 0.035, 0.12, 1.6, -32},
    {12, 2, (Nlm_FloatHi) INT2_MAX, 0.181, 0.025, 0.095, 1.9, -41},
    {19, 1, (Nlm_FloatHi) INT2_MAX, 0.212, 0.057, 0.18, 1.2, -21},
    {18, 1, (Nlm_FloatHi) INT2_MAX, 0.207, 0.050, 0.15, 1.4, -28},
    {17, 1, (Nlm_FloatHi) INT2_MAX, 0.198, 0.037, 0.12, 1.6, -33},
    {16, 1, (Nlm_FloatHi) INT2_MAX, 0.186, 0.025, 0.10, 1.9, -42},
    {15, 1, (Nlm_FloatHi) INT2_MAX, 0.171, 0.015, 0.063, 2.7, -76},
};

static Int4 blosum50_prefs[BLOSUM50_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_BEST,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};

#define BLOSUM62_VALUES_MAX 12
static Nlm_FloatHi  blosum62_values[BLOSUM62_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2},
    {11, 2, (Nlm_FloatHi) INT2_MAX, 0.297, 0.082, 0.27, 1.1, -10},
    {10, 2, (Nlm_FloatHi) INT2_MAX, 0.291, 0.075, 0.23, 1.3, -15},
    {9, 2, (Nlm_FloatHi) INT2_MAX, 0.279, 0.058, 0.19, 1.5, -19},
    {8, 2, (Nlm_FloatHi) INT2_MAX, 0.264, 0.045, 0.15, 1.8, -26},
    {7, 2, (Nlm_FloatHi) INT2_MAX, 0.239, 0.027, 0.10, 2.5, -46},
    {6, 2, (Nlm_FloatHi) INT2_MAX, 0.201, 0.012, 0.061, 3.3, -58},
    {13, 1, (Nlm_FloatHi) INT2_MAX, 0.292, 0.071, 0.23, 1.2, -11},
    {12, 1, (Nlm_FloatHi) INT2_MAX, 0.283, 0.059, 0.19, 1.5, -19},
    {11, 1, (Nlm_FloatHi) INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30},
    {10, 1, (Nlm_FloatHi) INT2_MAX, 0.243, 0.024, 0.10, 2.5, -44},
    {9, 1, (Nlm_FloatHi) INT2_MAX, 0.206, 0.010, 0.052, 4.0, -87},
};

static Int4 blosum62_prefs[BLOSUM62_VALUES_MAX] = {
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_BEST,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
};


#define BLOSUM80_VALUES_MAX 10
static Nlm_FloatHi  blosum80_values[BLOSUM80_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.3430, 0.177, 0.6568, 0.5222, -1.6},
    {25, 2, (Nlm_FloatHi) INT2_MAX, 0.342, 0.17, 0.66, 0.52, -1.6},
    {13, 2, (Nlm_FloatHi) INT2_MAX, 0.336, 0.15, 0.57, 0.59, -3},
    {9, 2, (Nlm_FloatHi) INT2_MAX, 0.319, 0.11, 0.42, 0.76, -6},
    {8, 2, (Nlm_FloatHi) INT2_MAX, 0.308, 0.090, 0.35, 0.89, -9},
    {7, 2, (Nlm_FloatHi) INT2_MAX, 0.293, 0.070, 0.27, 1.1, -14},
    {6, 2, (Nlm_FloatHi) INT2_MAX, 0.268, 0.045, 0.19, 1.4, -19},
    {11, 1, (Nlm_FloatHi) INT2_MAX, 0.314, 0.095, 0.35, 0.90, -9},
    {10, 1, (Nlm_FloatHi) INT2_MAX, 0.299, 0.071, 0.27, 1.1, -14},
    {9, 1, (Nlm_FloatHi) INT2_MAX, 0.279, 0.048, 0.20, 1.4, -19},
};

static Int4 blosum80_prefs[BLOSUM80_VALUES_MAX] = {
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_NOMINAL,
    BLAST_MATRIX_BEST,
    BLAST_MATRIX_NOMINAL
};

#define BLOSUM90_VALUES_MAX 8
static Nlm_FloatHi  blosum90_values[BLOSUM90_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.3346, 0.190, 0.7547, 0.4434, -1.4},
    {9, 2, (Nlm_FloatHi) INT2_MAX, 0.310, 0.12, 0.46, 0.67, -6},
    {8, 2, (Nlm_FloatHi) INT2_MAX, 0.300, 0.099, 0.39, 0.76, -7},
    {7, 2, (Nlm_FloatHi) INT2_MAX, 0.283, 0.072, 0.30, 0.93, -11},
    {6, 2, (Nlm_FloatHi) INT2_MAX, 0.259, 0.048, 0.22, 1.2, -16},
    {11, 1, (Nlm_FloatHi) INT2_MAX, 0.302, 0.093, 0.39, 0.78, -8},
    {10, 1, (Nlm_FloatHi) INT2_MAX, 0.290, 0.075, 0.28, 1.04, -15},
    {9, 1, (Nlm_FloatHi) INT2_MAX, 0.265, 0.044, 0.20, 1.3, -19},
};

static Int4 blosum90_prefs[BLOSUM90_VALUES_MAX] = {
	BLAST_MATRIX_NOMINAL,
	BLAST_MATRIX_NOMINAL,
	BLAST_MATRIX_NOMINAL,
	BLAST_MATRIX_NOMINAL,
	BLAST_MATRIX_NOMINAL,
	BLAST_MATRIX_NOMINAL,
	BLAST_MATRIX_BEST,
	BLAST_MATRIX_NOMINAL
};

#define PAM250_VALUES_MAX 16
static Nlm_FloatHi  pam250_values[PAM250_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.2252, 0.0868, 0.2223, 0.98, -5.0},
    {15, 3, (Nlm_FloatHi) INT2_MAX, 0.205, 0.049, 0.13, 1.6, -23},
    {14, 3, (Nlm_FloatHi) INT2_MAX, 0.200, 0.043, 0.12, 1.7, -26},
    {13, 3, (Nlm_FloatHi) INT2_MAX, 0.194, 0.036, 0.10, 1.9, -31},
    {12, 3, (Nlm_FloatHi) INT2_MAX, 0.186, 0.029, 0.085, 2.2, -41},
    {11, 3, (Nlm_FloatHi) INT2_MAX, 0.174, 0.020, 0.070, 2.5, -48},
    {17, 2, (Nlm_FloatHi) INT2_MAX, 0.204, 0.047, 0.12, 1.7, -28},
    {16, 2, (Nlm_FloatHi) INT2_MAX, 0.198, 0.038, 0.11, 1.8, -29},
    {15, 2, (Nlm_FloatHi) INT2_MAX, 0.191, 0.031, 0.087, 2.2, -44},
    {14, 2, (Nlm_FloatHi) INT2_MAX, 0.182, 0.024, 0.073, 2.5, -53},
    {13, 2, (Nlm_FloatHi) INT2_MAX, 0.171, 0.017, 0.059, 2.9, -64},
    {21, 1, (Nlm_FloatHi) INT2_MAX, 0.205, 0.045, 0.11, 1.8, -34},
    {20, 1, (Nlm_FloatHi) INT2_MAX, 0.199, 0.037, 0.10, 1.9, -35},
    {19, 1, (Nlm_FloatHi) INT2_MAX, 0.192, 0.029, 0.083, 2.3, -52},
    {18, 1, (Nlm_FloatHi) INT2_MAX, 0.183, 0.021, 0.070, 2.6, -60},
    {17, 1, (Nlm_FloatHi) INT2_MAX, 0.171, 0.014, 0.052, 3.3, -86},
};

static Int4 pam250_prefs[PAM250_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_BEST,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};

#define PAM30_VALUES_MAX 7
static Nlm_FloatHi  pam30_values[PAM30_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.3400, 0.283, 1.754, 0.1938, -0.3},
    {7, 2, (Nlm_FloatHi) INT2_MAX, 0.305, 0.15, 0.87, 0.35, -3},
    {6, 2, (Nlm_FloatHi) INT2_MAX, 0.287, 0.11, 0.68, 0.42, -4},
    {5, 2, (Nlm_FloatHi) INT2_MAX, 0.264, 0.079, 0.45, 0.59, -7},
    {10, 1, (Nlm_FloatHi) INT2_MAX, 0.309, 0.15, 0.88, 0.35, -3},
    {9, 1, (Nlm_FloatHi) INT2_MAX, 0.294, 0.11, 0.61, 0.48, -6},
    {8, 1, (Nlm_FloatHi) INT2_MAX, 0.270, 0.072, 0.40, 0.68, -10},
};

static Int4 pam30_prefs[PAM30_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_BEST,
BLAST_MATRIX_NOMINAL,
};


#define PAM70_VALUES_MAX 7
static Nlm_FloatHi  pam70_values[PAM70_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.3345, 0.229, 1.029, 0.3250,   -0.7},
    {8, 2, (Nlm_FloatHi) INT2_MAX, 0.301, 0.12, 0.54, 0.56, -5},
    {7, 2, (Nlm_FloatHi) INT2_MAX, 0.286, 0.093, 0.43, 0.67, -7},
    {6, 2, (Nlm_FloatHi) INT2_MAX, 0.264, 0.064, 0.29, 0.90, -12},
    {11, 1, (Nlm_FloatHi) INT2_MAX, 0.305, 0.12, 0.52, 0.59, -6},
    {10, 1, (Nlm_FloatHi) INT2_MAX, 0.291, 0.091, 0.41, 0.71, -9},
    {9, 1, (Nlm_FloatHi) INT2_MAX, 0.270, 0.060, 0.28, 0.97, -14},
};

static Int4 pam70_prefs[PAM70_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_BEST,
BLAST_MATRIX_NOMINAL
};



#define BLOSUM62_20_VALUES_MAX 65
static Nlm_FloatHi  blosum62_20_values[BLOSUM62_20_VALUES_MAX][8] = {
    {(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.03391, 0.125, 0.4544, 0.07462, -3.2},
    {100, 12, (Nlm_FloatHi) INT2_MAX, 0.0300, 0.056, 0.21, 0.14, -15},
    {95, 12, (Nlm_FloatHi) INT2_MAX, 0.0291, 0.047, 0.18, 0.16, -20},
    {90, 12, (Nlm_FloatHi) INT2_MAX, 0.0280, 0.038, 0.15, 0.19, -28},
    {85, 12, (Nlm_FloatHi) INT2_MAX, 0.0267, 0.030, 0.13, 0.21, -31},
    {80, 12, (Nlm_FloatHi) INT2_MAX, 0.0250, 0.021, 0.10, 0.25, -39},
    {105, 11, (Nlm_FloatHi) INT2_MAX, 0.0301, 0.056, 0.22, 0.14, -16},
    {100, 11, (Nlm_FloatHi) INT2_MAX, 0.0294, 0.049, 0.20, 0.15, -17},
    {95, 11, (Nlm_FloatHi) INT2_MAX, 0.0285, 0.042, 0.16, 0.18, -25},
    {90, 11, (Nlm_FloatHi) INT2_MAX, 0.0271, 0.031, 0.14, 0.20, -28},
    {85, 11, (Nlm_FloatHi) INT2_MAX, 0.0256, 0.023, 0.10, 0.26, -46},
    {115, 10, (Nlm_FloatHi) INT2_MAX, 0.0308, 0.062, 0.22, 0.14, -20},
    {110, 10, (Nlm_FloatHi) INT2_MAX, 0.0302, 0.056, 0.19, 0.16, -26},
    {105, 10, (Nlm_FloatHi) INT2_MAX, 0.0296, 0.050, 0.17, 0.17, -27},
    {100, 10, (Nlm_FloatHi) INT2_MAX, 0.0286, 0.041, 0.15, 0.19, -32},
    {95, 10, (Nlm_FloatHi) INT2_MAX, 0.0272, 0.030, 0.13, 0.21, -35},
    {90, 10, (Nlm_FloatHi) INT2_MAX, 0.0257, 0.022, 0.11, 0.24, -40},
    {85, 10, (Nlm_FloatHi) INT2_MAX, 0.0242, 0.017, 0.083, 0.29, -51},
    {115, 9, (Nlm_FloatHi) INT2_MAX, 0.0306, 0.061, 0.24, 0.13, -14},
    {110, 9, (Nlm_FloatHi) INT2_MAX, 0.0299, 0.053, 0.19, 0.16, -23},
    {105, 9, (Nlm_FloatHi) INT2_MAX, 0.0289, 0.043, 0.17, 0.17, -23},
    {100, 9, (Nlm_FloatHi) INT2_MAX, 0.0279, 0.036, 0.14, 0.20, -31},
    {95, 9, (Nlm_FloatHi) INT2_MAX, 0.0266, 0.028, 0.12, 0.23, -37},
    {120, 8, (Nlm_FloatHi) INT2_MAX, 0.0307, 0.062, 0.22, 0.14, -18},
    {115, 8, (Nlm_FloatHi) INT2_MAX, 0.0300, 0.053, 0.20, 0.15, -19},
    {110, 8, (Nlm_FloatHi) INT2_MAX, 0.0292, 0.046, 0.17, 0.17, -23},
    {105, 8, (Nlm_FloatHi) INT2_MAX, 0.0280, 0.035, 0.14, 0.20, -31},
    {100, 8, (Nlm_FloatHi) INT2_MAX, 0.0266, 0.026, 0.12, 0.23, -37},
    {125, 7, (Nlm_FloatHi) INT2_MAX, 0.0306, 0.058, 0.22, 0.14, -18},
    {120, 7, (Nlm_FloatHi) INT2_MAX, 0.0300, 0.052, 0.19, 0.16, -23},
    {115, 7, (Nlm_FloatHi) INT2_MAX, 0.0292, 0.044, 0.17, 0.17, -24},
    {110, 7, (Nlm_FloatHi) INT2_MAX, 0.0279, 0.032, 0.14, 0.20, -31},
    {105, 7, (Nlm_FloatHi) INT2_MAX, 0.0267, 0.026, 0.11, 0.24, -41},
    {120,10,5, 0.0298, 0.049, 0.19, 0.16, -21},
    {115,10,5, 0.0290, 0.042, 0.16, 0.18, -25},
    {110,10,5, 0.0279, 0.033, 0.13, 0.21, -32},
    {105,10,5, 0.0264, 0.024, 0.10, 0.26, -46},
    {100,10,5, 0.0250, 0.018, 0.081, 0.31, -56},
    {125,10,4, 0.0301, 0.053, 0.18, 0.17, -25},
    {120,10,4, 0.0292, 0.043, 0.15, 0.20, -33},
    {115,10,4, 0.0282, 0.035, 0.13, 0.22, -36},
    {110,10,4, 0.0270, 0.027, 0.11, 0.25, -41},
    {105,10,4, 0.0254, 0.020, 0.079, 0.32, -60},
    {130,10,3, 0.0300, 0.051, 0.17, 0.18, -27},
    {125,10,3, 0.0290, 0.040, 0.13, 0.22, -38},
    {120,10,3, 0.0278, 0.030, 0.11, 0.25, -44},
    {115,10,3, 0.0267, 0.025, 0.092, 0.29, -52},
    {110,10,3, 0.0252, 0.018, 0.070, 0.36, -70},
    {135,10,2, 0.0292, 0.040, 0.13, 0.22, -35},
    {130,10,2, 0.0283, 0.034, 0.10, 0.28, -51},
    {125,10,2, 0.0269, 0.024, 0.077, 0.35, -71},
    {120,10,2, 0.0253, 0.017, 0.059, 0.43, -90},
    {115,10,2, 0.0234, 0.011, 0.043, 0.55, -121},
    {100,14,3, 0.0258, 0.023, 0.087, 0.33, -59},
    {105,13,3, 0.0263, 0.024, 0.085, 0.31, -57},
    {110,12,3, 0.0271, 0.028, 0.093, 0.29, -54},
    {115,11,3, 0.0275, 0.030, 0.10, 0.27, -49},
    {125,9,3, 0.0283, 0.034, 0.12, 0.23, -38},
    {130,8,3, 0.0287, 0.037, 0.12, 0.23, -40},
    {125,7,3, 0.0287, 0.036, 0.12, 0.24, -44},
    {140,6,3, 0.0285, 0.033, 0.12, 0.23, -40},
    {105,14,3, 0.0270, 0.028, 0.10, 0.27, -46},
    {110,13,3, 0.0279, 0.034, 0.10, 0.27, -50},
    {115,12,3, 0.0282, 0.035, 0.12, 0.24, -42},
    {120,11,3, 0.0286, 0.037, 0.12, 0.24, -44},
};

static Int4 blosum62_20_prefs[BLOSUM62_20_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_BEST,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};

/*
	Allocates memory for the BLAST_ScoreBlkPtr.
*/

BLAST_ScoreBlkPtr LIBCALL
BLAST_ScoreBlkNew(Uint1 alphabet, Int2 number_of_contexts)

{
	BLAST_ScoreBlkPtr sbp;
	SeqCodeTablePtr sctp;

	if (alphabet != BLASTNA_SEQ_CODE)
	{
		sctp = SeqCodeTableFindObj(alphabet);
		if (sctp == NULL)
			return NULL;
	}

	sbp = (BLAST_ScoreBlkPtr) MemNew(sizeof(BLAST_ScoreBlk));

	if (sbp != NULL)
	{
		sbp->alphabet_code = alphabet;
		if (alphabet != BLASTNA_SEQ_CODE)
			sbp->alphabet_size = sctp->num;
		else
			sbp->alphabet_size = BLASTNA_SIZE;

/* set the alphabet type (protein or not). */
		switch (alphabet) {
			case Seq_code_iupacaa:
			case Seq_code_ncbi8aa:
			case Seq_code_ncbieaa:
			case Seq_code_ncbipaa:
			case Seq_code_iupacaa3:
				return NULL;

			case Seq_code_ncbistdaa:
				sbp->protein_alphabet = TRUE;
				break;
			case BLASTNA_SEQ_CODE:
				sbp->protein_alphabet = FALSE;
				break;
			default:
				break;
		}

		sbp->matrix_struct = BlastMatrixAllocate(sbp->alphabet_size);
		if (sbp->matrix_struct == NULL)
		{
			sbp = BLAST_ScoreBlkDestruct(sbp);
			return sbp;
		}
		sbp->matrix = sbp->matrix_struct->matrix;
		sbp->maxscore = (BLAST_ScorePtr) MemNew(BLAST_MATRIX_SIZE*sizeof(BLAST_Score));
		sbp->number_of_contexts = number_of_contexts;
		sbp->sfp = MemNew(sbp->number_of_contexts*sizeof(BLAST_ScoreFreqPtr));
		sbp->kbp_std = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
		sbp->kbp_gap_std = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
		sbp->kbp_psi = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
		sbp->kbp_gap_psi = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
	}

	return sbp;
}

BLAST_ScoreBlkPtr LIBCALL
BLAST_ScoreBlkDestruct(BLAST_ScoreBlkPtr sbp)

{
    Int4 index, rows;
    if (sbp == NULL)
        return NULL;
    
    for (index=0; index<sbp->number_of_contexts; index++) {
        if (sbp->sfp)
            sbp->sfp[index] = BlastScoreFreqDestruct(sbp->sfp[index]);
        if (sbp->kbp_std)
            sbp->kbp_std[index] = BlastKarlinBlkDestruct(sbp->kbp_std[index]);
        if (sbp->kbp_gap_std)
            sbp->kbp_gap_std[index] = BlastKarlinBlkDestruct(sbp->kbp_gap_std[index]);
        if (sbp->kbp_psi)
            sbp->kbp_psi[index] = BlastKarlinBlkDestruct(sbp->kbp_psi[index]);
        if (sbp->kbp_gap_psi)
            sbp->kbp_gap_psi[index] = BlastKarlinBlkDestruct(sbp->kbp_gap_psi[index]);
    }
    if (sbp->kbp_ideal)
        sbp->kbp_ideal = BlastKarlinBlkDestruct(sbp->kbp_ideal);
    sbp->sfp = MemFree(sbp->sfp);
    sbp->kbp_std = MemFree(sbp->kbp_std);
    sbp->kbp_psi = MemFree(sbp->kbp_psi);
    sbp->kbp_gap_std = MemFree(sbp->kbp_gap_std);
    sbp->kbp_gap_psi = MemFree(sbp->kbp_gap_psi);
    sbp->matrix_struct = BlastMatrixDestruct(sbp->matrix_struct);
    sbp->maxscore = MemFree(sbp->maxscore);
    sbp->comments = ValNodeFreeData(sbp->comments);
    sbp->name = MemFree(sbp->name);
    sbp->ambiguous_res = MemFree(sbp->ambiguous_res);
    
    /* Removing posMatrix and posFreqs if any */
    rows = sbp->query_length + 1;
    if(sbp->posMatrix != NULL) {
        for (index=0; index < rows; index++) {
            MemFree(sbp->posMatrix[index]);
        }
        MemFree(sbp->posMatrix);
    }
    
    if(sbp->posFreqs != NULL) {
        for (index = 0; index < rows; index++) {
            MemFree(sbp->posFreqs[index]);
        }
        MemFree(sbp->posFreqs);
    }
    
    sbp = MemFree(sbp);
    
    return sbp;
}

/* 
	Set the ambiguous residue (e.g, 'N', 'X') in the BLAST_ScoreBlkPtr.
	Convert from ncbieaa to sbp->alphabet_code (i.e., ncbistdaa) first.
*/
Int2 LIBCALL
BlastScoreSetAmbigRes(BLAST_ScoreBlkPtr sbp, Char ambiguous_res)

{
	Int2 index;
	SeqMapTablePtr smtp;
	Uint1Ptr ambig_buffer;

	if (sbp == NULL)
		return 1;

	if (sbp->ambig_occupy >= sbp->ambig_size)
	{
		sbp->ambig_size += 5;
		ambig_buffer = MemNew(sbp->ambig_size*sizeof(Uint1));
		for (index=0; index<sbp->ambig_occupy; index++)
		{
			ambig_buffer[index] = sbp->ambiguous_res[index];
		}
		sbp->ambiguous_res = MemFree(sbp->ambiguous_res);
		sbp->ambiguous_res = ambig_buffer;
	}

	if (sbp->alphabet_code == Seq_code_ncbistdaa)
	{
		smtp = SeqMapTableFind(sbp->alphabet_code, Seq_code_ncbieaa);
		sbp->ambiguous_res[sbp->ambig_occupy] = SeqMapTableConvert(smtp, ambiguous_res);
	}
	else if (sbp->alphabet_code == BLASTNA_SEQ_CODE)
	{
		smtp = SeqMapTableFind(Seq_code_ncbi4na, Seq_code_iupacna);
		sbp->ambiguous_res[sbp->ambig_occupy] = ncbi4na_to_blastna[SeqMapTableConvert(smtp, ambiguous_res)];
	}
	else if (sbp->alphabet_code == Seq_code_ncbi4na)
	{
		smtp = SeqMapTableFind(Seq_code_ncbi4na, Seq_code_iupacna);
		sbp->ambiguous_res[sbp->ambig_occupy] = SeqMapTableConvert(smtp, ambiguous_res);
	}

	(sbp->ambig_occupy)++;
	

	return 0;
}

/*
	Deallocates all data associated with the BLAST_MatrixPtr.
*/
BLAST_MatrixPtr LIBCALL
BLAST_MatrixDestruct(BLAST_MatrixPtr blast_matrix)

{
    Int4 index;

    if (blast_matrix == NULL)
        return NULL;
    
    /* We may have 2 different matrices in there */
    
    if(blast_matrix->original_matrix && 
       blast_matrix->original_matrix != blast_matrix->matrix) {
        for (index=0; index < 26; index++) {
            MemFree(blast_matrix->original_matrix[index]);
        }
        MemFree(blast_matrix->original_matrix);
    }
    
    blast_matrix->name = MemFree(blast_matrix->name);
    if (blast_matrix->matrix) {
	for (index=0; index<blast_matrix->rows; index++) {
	    MemFree(blast_matrix->matrix[index]);
	}
	MemFree(blast_matrix->matrix);
    }
    
    if(blast_matrix->posFreqs != NULL) {
        for (index = 0; index < blast_matrix->rows; index++) {
            MemFree(blast_matrix->posFreqs[index]);
        }
        MemFree(blast_matrix->posFreqs);
    }

    MemFree(blast_matrix);

    return NULL;
}


/*
	Allocates and fills the BLAST_MatrixPtr.  
	positionBased indicates that the posMatrix
	on BLAST_ScoreBlkPtr should be used rather
	than the 'normal' matrix.
*/
BLAST_MatrixPtr LIBCALL
BLAST_MatrixFill(BLAST_ScoreBlkPtr sbp, Boolean positionBased)

{
    BLAST_MatrixPtr blast_matrix;
    FloatHi karlinK = 0.0;
    Int4 index1, index2, dim1, dim2;
    Int4Ptr PNTR matrix, PNTR original_matrix;
    Nlm_FloatHi **posFreqs = NULL;

    if (sbp == NULL)
        return NULL;
    
    blast_matrix = (BLAST_MatrixPtr) MemNew(sizeof(BLAST_Matrix));
    
    dim1 = sbp->alphabet_size;
    dim2 = sbp->alphabet_size;
    original_matrix = sbp->matrix;

    if (sbp->kbp_gap_psi[0])
        karlinK = sbp->kbp_gap_psi[0]->K;
    
    matrix = MemNew(dim1*sizeof(Int4Ptr));
    for (index1=0; index1<dim1; index1++) {
        matrix[index1] = MemNew(dim2*sizeof(Int4));
        for (index2=0; index2<dim2; index2++) {
            matrix[index1][index2] = original_matrix[index1][index2];
        }
    }

    blast_matrix->original_matrix = matrix;
    
    /* For PSI BLAST blast_matrix->matrix will be position based */
    if (sbp->posMatrix) {
        dim1 = sbp->query_length + 1;
        dim2 = sbp->alphabet_size;
        original_matrix = sbp->posMatrix;
        matrix = MemNew(dim1*sizeof(Int4Ptr));
        for (index1=0; index1<dim1; index1++) {
           matrix[index1] = MemNew(dim2*sizeof(Int4));
           for (index2=0; index2<dim2; index2++) {
              matrix[index1][index2] = original_matrix[index1][index2];
           }
        }
    }
    blast_matrix->matrix = matrix;

    /* Copying posFreqs to the BLAST_Matrix */
    if ((sbp->posFreqs != NULL) && (sbp->posMatrix != NULL)) {
        posFreqs = MemNew(dim1*sizeof(Nlm_FloatHi *));
        for (index1 = 0; index1 < dim1; index1++) {
            posFreqs[index1] = MemNew(dim2*sizeof(Nlm_FloatHi));
            for (index2=0; index2 < dim2; index2++) {
                posFreqs[index1][index2] = sbp->posFreqs[index1][index2];
            }
        }
    }

    blast_matrix->is_prot = sbp->protein_alphabet;
    blast_matrix->name = StringSave(sbp->name);
    blast_matrix->rows = dim1;
    blast_matrix->columns = dim2;
    blast_matrix->matrix = matrix;
    blast_matrix->posFreqs = posFreqs;
    blast_matrix->karlinK = karlinK;
    
    return blast_matrix;
}

/* 
	Fill in the matrix for blastn using the penaly and rewards

	The query sequence alphabet is blastna, the subject sequence
	is ncbi2na.  The alphabet blastna is defined in blastkar.h
	and the first four elements of blastna are identical to ncbi2na.

	The query is in the first index, the subject is the second.
        if matrix==NULL, it is allocated and returned.
*/

BLAST_ScorePtr PNTR LIBCALL BlastScoreBlkMatCreateEx(BLAST_ScorePtr PNTR matrix,BLAST_Score penalty, BLAST_Score reward)
{

	Int2	index1, index2, degen;
	Int2 degeneracy[BLASTNA_SIZE+1];
	
        if(!matrix) {
            BLASTMatrixStructurePtr matrix_struct;
            matrix_struct =BlastMatrixAllocate((Int2) BLASTNA_SIZE);
            matrix = matrix_struct->matrix;
        }
	for (index1 = 0; index1<BLASTNA_SIZE; index1++) /* blastna */
		for (index2 = 0; index2<BLASTNA_SIZE; index2++) /* blastna */
			matrix[index1][index2] = 0;

	/* In blastna the 1st four bases are A, C, G, and T, exactly as it is ncbi2na. */
	/* ncbi4na gives them the value 1, 2, 4, and 8.  */

	/* Set the first four bases to degen. one */
	for (index1=0; index1<NUMBER_NON_AMBIG_BP; index1++)
		degeneracy[index1] = 1;

	for (index1=NUMBER_NON_AMBIG_BP; index1<BLASTNA_SIZE; index1++) /* blastna */
	{
		degen=0;
		for (index2=0; index2<NUMBER_NON_AMBIG_BP; index2++) /* ncbi2na */
		{
			if (blastna_to_ncbi4na[index1] & blastna_to_ncbi4na[index2])
				degen++;
		}
		degeneracy[index1] = degen;
	}


	for (index1=0; index1<BLASTNA_SIZE; index1++) /* blastna */
	{
		for (index2=index1; index2<BLASTNA_SIZE; index2++) /* blastna */
		{
			if (blastna_to_ncbi4na[index1] & blastna_to_ncbi4na[index2])
			{ /* round up for positive scores, down for negatives. */
				matrix[index1][index2] = Nlm_Nint( (double) ((degeneracy[index2]-1)*penalty + reward))/degeneracy[index2];
				if (index1 != index2)
				{
				      matrix[index2][index1] = matrix[index1][index2];
				}
			}
			else
			{
				matrix[index1][index2] = penalty;
				matrix[index2][index1] = penalty;
			}
		}
	}

        /* The value of 15 is a gap, which is a sentinel between strands in 
           the ungapped extension algorithm */
	for (index1=0; index1<BLASTNA_SIZE; index1++) /* blastna */
           matrix[BLASTNA_SIZE-1][index1] = INT4_MIN / 2;
	for (index1=0; index1<BLASTNA_SIZE; index1++) /* blastna */
           matrix[index1][BLASTNA_SIZE-1] = INT4_MIN / 2;

	return matrix;
}
/* 
	Fill in the matrix for blastn using the penaly and rewards on
	the BLAST_ScoreBlkPtr.

	The query sequence alphabet is blastna, the subject sequence
	is ncbi2na.  The alphabet blastna is defined in blastkar.h
	and the first four elements of blastna are identical to ncbi2na.

	The query is in the first index, the subject is the second.
*/


static Int2 BlastScoreBlkMatCreate(BLAST_ScoreBlkPtr sbp)
{
	Char buffer[25];

        sbp->matrix = BlastScoreBlkMatCreateEx(sbp->matrix,sbp->penalty, sbp->reward);
	sbp->mat_dim1 = BLASTNA_SIZE;
	sbp->mat_dim2 = BLASTNA_SIZE;

	sprintf(buffer, "blastn matrix:%ld %ld", (long) sbp->reward, (long) sbp->penalty);
	sbp->name = StringSave(buffer);

	return 0;
}





/*
	This function fills in the BLAST_ScoreBlk structure.  
	Tasks are:
		-read in the matrix
		-set maxscore
*/

Int2 LIBCALL
BlastScoreBlkMatFill(BLAST_ScoreBlkPtr sbp, CharPtr matrix)

{
    Char string[PATH_MAX] = "", alphabet_type[3] = "";
    CharPtr matrix_dir = NULL;
    Int2 status = 0;
    FILE *fp = NULL;
    
    if (sbp->read_in_matrix) {
        /* Convert matrix name to upper case. */
        matrix = Nlm_StrUpper(matrix);
        
        sbp->name = StringSave(matrix);	/* Save the name of the matrix. */

        /* 1. Try local directory */
        if(FileLength(matrix) > 0)
            fp = FileOpen(matrix, "r");
        
        /* 2. Try configuration file */
        if (fp == NULL) {
           if (sbp->protein_alphabet)
              Nlm_StringNCpy(alphabet_type, "aa", 2);
           else
              Nlm_StringNCpy(alphabet_type, "nt", 2);
           alphabet_type[2] = NULLB;

           if(FindPath("ncbi", "ncbi", "data", string, PATH_MAX)) {
               matrix_dir = StringSave(string);
               sprintf(string, "%s%s", matrix_dir, matrix);
               if(FileLength(string) > 0) {
                  fp = FileOpen(string, "r");
               } else {
                  sprintf(string, "%s%s%s%s", matrix_dir, 
                          alphabet_type, DIRDELIMSTR, matrix);
                  if(FileLength(string) > 0)
                     fp = FileOpen(string, "r");
               }
               matrix_dir = MemFree(matrix_dir);
            }
        }
        /* Trying to use local "data" directory */

        if(fp == NULL) {
            sprintf(string, "data%s%s", DIRDELIMSTR, matrix);
            if(FileLength(string) > 0)
                fp = FileOpen(string, "r");
        }
 
        /* Get the matrix locations from the environment for UNIX. */
        if (fp == NULL) {
            
            matrix_dir = getenv("BLASTMAT");
            if (matrix_dir != NULL) {
                sprintf(string, "%s%s%s%s%s", matrix_dir, DIRDELIMSTR,alphabet_type, DIRDELIMSTR, matrix); 
            }

            if(FileLength(string) > 0)
                fp = FileOpen(string, "r");
            
            /* Try again without "aa" or "nt" */
            if (fp == NULL) {
                if (matrix_dir != NULL) {
                    sprintf(string, "%s%s%s", matrix_dir, DIRDELIMSTR, matrix); 
                }

                if(FileLength(string) > 0)
                    fp = FileOpen(string, "r");
            }
        }

        if (fp == NULL) {
            ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", matrix);
            return 4;
        }
        
        if((status=BlastScoreBlkMatRead(sbp, fp)) != 0) {
            FileClose(fp);
            return status;
        }
        FileClose(fp);
    } else {
        if((status=BlastScoreBlkMatCreate(sbp)) != 0)
            return status;
    }
    
    if((status=BlastScoreBlkMaxScoreSet(sbp)) != 0)
        return status;
    
    return status;
}

/*
	Return the specified matrix.  Do this by setting up the ScoreBlkPtr
	and then fetching the matrix from disk.
*/

BLAST_MatrixPtr LIBCALL
BLAST_MatrixFetch(CharPtr matrix_name)

{
	BLAST_MatrixPtr matrix;
	BLAST_ScoreBlkPtr sbp;

	if (matrix_name == NULL)
		return NULL;

	sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa, 1);

	/* Read in for protein. */
	sbp->read_in_matrix = TRUE;

	BlastScoreBlkMatFill(sbp, matrix_name);

	matrix = BLAST_MatrixFill(sbp, FALSE);

	sbp = BLAST_ScoreBlkDestruct(sbp);	

	return matrix;
}


/*
  Calculate the Karlin parameters.  This function should be called once
  for each context, or frame translated.
  
  The rfp and stdrfp are calculated for each context, this should be
  fixed. 
*/

Int2 LIBCALL
BlastScoreBlkFill(BLAST_ScoreBlkPtr sbp, CharPtr query, Int4 query_length, Int2 context_number)
{
	BLAST_ResFreqPtr rfp, stdrfp;
	Int2 retval=0;

	rfp = BlastResFreqNew(sbp);
	stdrfp = BlastResFreqNew(sbp);
	BlastResFreqStdComp(sbp, stdrfp);
	BlastResFreqString(sbp, rfp, query, query_length);
	sbp->sfp[context_number] = BlastScoreFreqNew(sbp->loscore, sbp->hiscore);
	BlastScoreFreqCalc(sbp, sbp->sfp[context_number], rfp, stdrfp);
	sbp->kbp_std[context_number] = BlastKarlinBlkCreate();
	retval = BlastKarlinBlkCalc(sbp->kbp_std[context_number], sbp->sfp[context_number]);
	if (retval)
	{
		rfp = BlastResFreqDestruct(rfp);
		stdrfp = BlastResFreqDestruct(stdrfp);
		return retval;
	}
	sbp->kbp_psi[context_number] = BlastKarlinBlkCreate();
	retval = BlastKarlinBlkCalc(sbp->kbp_psi[context_number], sbp->sfp[context_number]);
	rfp = BlastResFreqDestruct(rfp);
	stdrfp = BlastResFreqDestruct(stdrfp);

	return retval;
}

/*
	Calculates the standard Karlin parameters.  This is used
	if the query is translated and the calculated (real) Karlin
	parameters are bad, as they're calculated for non-coding regions.
*/

BLAST_KarlinBlkPtr LIBCALL
BlastKarlinBlkStandardCalcEx(BLAST_ScoreBlkPtr sbp)

{
	BLAST_KarlinBlkPtr kbp_ideal;
	BLAST_ResFreqPtr stdrfp;
	BLAST_ScoreFreqPtr sfp;

	stdrfp = BlastResFreqNew(sbp);
	BlastResFreqStdComp(sbp, stdrfp);
	sfp = BlastScoreFreqNew(sbp->loscore, sbp->hiscore);
	BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
	kbp_ideal = BlastKarlinBlkCreate();
	BlastKarlinBlkCalc(kbp_ideal, sfp);
	stdrfp = BlastResFreqDestruct(stdrfp);

	sfp = BlastScoreFreqDestruct(sfp);

	return kbp_ideal;
}

Int2 LIBCALL
BlastKarlinBlkStandardCalc(BLAST_ScoreBlkPtr sbp, Int2 context_start, Int2 context_end)
{

	BLAST_KarlinBlkPtr kbp_ideal, kbp;
	Int2 index;

	kbp_ideal = BlastKarlinBlkStandardCalcEx(sbp);
/* Replace the calculated values with ideal ones for blastx, tblastx. */
	for (index=context_start; index<=context_end; index++)
	{
		kbp = sbp->kbp[index];	
		if (kbp->Lambda >= kbp_ideal->Lambda)
		{
			kbp->Lambda = kbp_ideal->Lambda;
			kbp->K = kbp_ideal->K;
			kbp->logK = kbp_ideal->logK;
			kbp->H = kbp_ideal->H;
		}
	}
	kbp_ideal = BlastKarlinBlkDestruct(kbp_ideal);

	return 0;
}

/*
	Creates the Karlin Block.
*/

BLAST_KarlinBlkPtr LIBCALL
BlastKarlinBlkCreate(void)

{
	BLAST_KarlinBlkPtr kbp;

	kbp = (BLAST_KarlinBlkPtr) MemNew(sizeof(BLAST_KarlinBlk));

	return kbp;
}

/*
	Deallocates the Karlin Block.
*/

BLAST_KarlinBlkPtr LIBCALL
BlastKarlinBlkDestruct(BLAST_KarlinBlkPtr kbp)

{
	kbp = MemFree(kbp);

	return kbp;
}


/* 
	Read in the matrix from the FILE *fp.

	This function ASSUMES that the matrices are in the ncbistdaa
	format.  BLAST should be able to use any alphabet, though it
	is expected that ncbistdaa will be used.
*/

static Char ASCII_TO_BLASTNA_CONVERT[128]={
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};

#define X_CODE 21
#define U_CODE 24

Int2 LIBCALL
BlastScoreBlkMatRead(BLAST_ScoreBlkPtr sbp, FILE *fp)
{
    Char	buf[512+3];
    Char	temp[512];
    CharPtr	cp, lp;
    Char		ch;
    BLAST_ScorePtr PNTR	matrix;
    BLAST_ScorePtr	m;
    BLAST_Score	score;
    Int2		a1cnt = 0, a2cnt = 0;
    Char    a1chars[BLAST_MAX_ALPHABET], a2chars[BLAST_MAX_ALPHABET];
    long	lineno = 0;
    Nlm_FloatHi	xscore;
    register int	index1, index2, total;
    SeqCodeTablePtr sctp;
    SeqMapTablePtr smtp=NULL;
    Int2 status;
    static TNlmMutex read_matrix_mutex;

    NlmMutexInit(&read_matrix_mutex);
    NlmMutexLock(read_matrix_mutex);
    
    matrix = sbp->matrix;	
    
    if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
        sctp = SeqCodeTableFindObj(sbp->alphabet_code);
        if(sctp == NULL) {
            NlmMutexUnlock(read_matrix_mutex); 
	    return 1;
        }

        total = sctp->start_at + sctp->num;
        for (index1 = sctp->start_at; index1 < total; index1++)
            for (index2 = sctp->start_at; index2 < total; index2++)
                matrix[index1][index2] = BLAST_SCORE_MIN;
        
        if (sbp->alphabet_code != Seq_code_ncbieaa) {
            smtp = SeqMapTableFind(sbp->alphabet_code, Seq_code_ncbieaa);
            if (smtp == NULL) {
                NlmMutexUnlock(read_matrix_mutex); 
                return 1;
            }
        }
    } else {
        /* Fill-in all the defaults ambiguity and normal codes */
        status=BlastScoreBlkMatCreate(sbp); 
        if(status != 0)
	{
        	NlmMutexUnlock(read_matrix_mutex); 
        	return status;
	}
    }
    
    /* Read the residue names for the second alphabet */
    while (Nlm_FileGets(buf, sizeof(buf), fp) != NULL) {
        ++lineno;
        if (Nlm_StrChr(buf, '\n') == NULL) {
            NlmMutexUnlock(read_matrix_mutex); 
            return 2;
        }

        if (buf[0] == COMMENT_CHR) {
            /* save the comment line in a linked list */
            *Nlm_StrChr(buf, '\n') = NULLB;
            ValNodeCopyStr(&sbp->comments, 0, buf+1);
            continue;
        }
        if ((cp = Nlm_StrChr(buf, COMMENT_CHR)) != NULL)
            *cp = NULLB;
        lp = (CharPtr)Nlm_StrTok(buf, TOKSTR);
        if (lp == NULL) /* skip blank lines */
            continue;
        while (lp != NULL) {
            
            if (smtp)
                ch = SeqMapTableConvert(smtp,  *lp);
            else {
                if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
                    ch = *lp;
                } else {
                    ch = ASCII_TO_BLASTNA_CONVERT[toupper(*lp)];
                }
            }
            a2chars[a2cnt++] = ch;
            lp = (CharPtr)Nlm_StrTok(NULL, TOKSTR);
        }
        
        break;	/* Exit loop after reading one line. */
    }
    
    if (a2cnt <= 1) { 
        NlmMutexUnlock(read_matrix_mutex); 
        return 2;
    }

    if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
        sbp->mat_dim2 = a2cnt;
    }
    while (Nlm_FileGets(buf, sizeof(buf), fp) != NULL)  {
        ++lineno;
        if ((cp = Nlm_StrChr(buf, '\n')) == NULL) {
            NlmMutexUnlock(read_matrix_mutex); 
            return 2;
        }
        if ((cp = Nlm_StrChr(buf, COMMENT_CHR)) != NULL)
            *cp = NULLB;
        if ((lp = (CharPtr)Nlm_StrTok(buf, TOKSTR)) == NULL)
            continue;
        ch = *lp;
        cp = (CharPtr) lp;
        if ((cp = Nlm_StrTok(NULL, TOKSTR)) == NULL) {
            NlmMutexUnlock(read_matrix_mutex); 
            return 2;
        }
        if (a1cnt >= DIM(a1chars)) {
            NlmMutexUnlock(read_matrix_mutex); 
            return 2;
        }

        if (smtp) {
            ch = SeqMapTableConvert(smtp,  ch);
            
        } else {
            if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
                ch = ASCII_TO_BLASTNA_CONVERT[toupper(ch)];
            }
        }
        a1chars[a1cnt++] = ch;
        m = &matrix[(int)ch][0];
        index2 = 0;
        while (cp != NULL) {
            if (index2 >= a2cnt) {
                NlmMutexUnlock(read_matrix_mutex); 
                return 2;
            }
            Nlm_StrCpy(temp, cp);
            
            if (Nlm_StrICmp(temp, "na") == 0)  {
                score = BLAST_SCORE_1MIN;
            } else  {
                if (sscanf(temp, "%lg", &xscore) != 1) {
                    NlmMutexUnlock(read_matrix_mutex); 
                    return 2;
                }
				/*xscore = MAX(xscore, BLAST_SCORE_1MIN);*/
                if (xscore > BLAST_SCORE_1MAX || xscore < BLAST_SCORE_1MIN) {
                    NlmMutexUnlock(read_matrix_mutex); 
                    return 2;
                }
                xscore += (xscore >= 0. ? 0.5 : -0.5);
                score = (BLAST_Score)xscore;
            }
            
            m[(int)a2chars[index2++]] = score;
            
            cp = Nlm_StrTok(NULL, TOKSTR);
        }
    }
    
    if (a1cnt <= 1) {
        NlmMutexUnlock(read_matrix_mutex); 
        return 2;
    }
    
    if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
        sbp->mat_dim1 = a1cnt;

        /* For protein matrices, copy the X scores to the
           U scores. Keeping the U scores as BLAST_SCORE_MIN
           means that U never aligns with another letter,
           even another U */
        for (index1 = 0; index1 < a2cnt; index1++) {
            matrix[U_CODE][index1] = matrix[X_CODE][index1];
            matrix[index1][U_CODE] = matrix[index1][X_CODE];
        }
    }
    
    NlmMutexUnlock(read_matrix_mutex); 
    return 0;
}

Int2 LIBCALL
BlastScoreBlkMaxScoreSet(BLAST_ScoreBlkPtr sbp)
{
	BLAST_Score score, maxscore;
	BLAST_ScorePtr PNTR  matrix; 
	Int2 index1, index2;

	sbp->loscore = BLAST_SCORE_1MAX;
        sbp->hiscore = BLAST_SCORE_1MIN;
	matrix = sbp->matrix;
	for (index1=0; index1<sbp->alphabet_size; index1++)
	{
		maxscore=BLAST_SCORE_MIN;
		for (index2=0; index2<sbp->alphabet_size; index2++)
		{
			score = matrix[index1][index2];
			if (score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX)
				continue;
			if (score > maxscore)
			{
				maxscore = score;
			}
			if (sbp->loscore > score)
				sbp->loscore = score;
			if (sbp->hiscore < score)
				sbp->hiscore = score;
		}
		sbp->maxscore[index1] = maxscore;
	}
/* If the lo/hi-scores are BLAST_SCORE_MIN/BLAST_SCORE_MAX, (i.e., for
gaps), then use other scores. */

	if (sbp->loscore < BLAST_SCORE_1MIN)
		sbp->loscore = BLAST_SCORE_1MIN;
	if (sbp->hiscore > BLAST_SCORE_1MAX)
		sbp->hiscore = BLAST_SCORE_1MAX;

	return 0;
}

/* maxscore for PSI Blast depends not on the residue, but on the position
   in the posMatrix, so maxscore array has size of the length of posMatrix */
/* SSH */
BLAST_ScorePtr BlastPSIMaxScoreGet(BLAST_ScorePtr PNTR posMatrix, 
                                   Int4 start, Int4 length)
{
    BLAST_Score score, maxscore;
    BLAST_ScorePtr maxscore_pos;
    Int4 index1, index2;
    
    if(posMatrix == NULL)
        return NULL;
    
    maxscore_pos = MemNew(length * sizeof(BLAST_Score));
    
    for (index1 = start; index1 < length; index1++) {
        maxscore = BLAST_SCORE_MIN;
        for (index2 = 0; index2 < PSI_ALPHABET_SIZE; index2++) {
            score = posMatrix[index1][index2];
            if (score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX)
                continue;
            if (score > maxscore) {
                maxscore = score;
            }
        }

        maxscore_pos[index1] = maxscore;
    }    
    return maxscore_pos;
}

static BLASTMatrixStructurePtr
BlastMatrixAllocate(Int2 alphabet_size)

{
	BLASTMatrixStructurePtr matrix_struct;
	Int2 index;

	if (alphabet_size <= 0 || alphabet_size >= BLAST_MATRIX_SIZE)
		return NULL;

	matrix_struct =	(BLASTMatrixStructurePtr) MemNew(sizeof(BLASTMatrixStructure));

	if (matrix_struct == NULL)
		return NULL;

	for (index=0; index<BLAST_MATRIX_SIZE-1; index++)
	{
		matrix_struct->matrix[index] = matrix_struct->long_matrix + index*BLAST_MATRIX_SIZE;
	}

	return matrix_struct;
}

static BLASTMatrixStructurePtr
BlastMatrixDestruct(BLASTMatrixStructurePtr matrix_struct)

{

	if (matrix_struct == NULL)
		return NULL;

	matrix_struct = MemFree(matrix_struct);

	return matrix_struct;
}

/* 
	Allocated the BLAST_ResCompPtr for a given alphabet.  Only the
	alphabets ncbistdaa and ncbi4na should be used by BLAST.
*/

BLAST_ResCompPtr LIBCALL
BlastResCompNew(BLAST_ScoreBlkPtr sbp)
{
	BLAST_ResCompPtr	rcp;

	rcp = (BLAST_ResCompPtr) MemNew(sizeof(BLAST_ResComp));
	if (rcp == NULL)
		return NULL;

	rcp->alphabet_code = sbp->alphabet_code;

/* comp0 has zero offset, comp starts at sctp->start_at, only one 
array is allocated.  */
	rcp->comp0 = (Int4Ptr) MemNew(BLAST_MATRIX_SIZE*sizeof(Int4));
	if (rcp->comp0 == NULL) 
	{
		rcp = BlastResCompDestruct(rcp);
		return rcp;
	}

	rcp->comp = rcp->comp0 - sbp->alphabet_start;

	return rcp;
}


BLAST_ResCompPtr LIBCALL
BlastResCompDestruct(BLAST_ResCompPtr rcp)
{
	if (rcp == NULL)
		return NULL;

	if (rcp->comp0 != NULL)
		rcp->comp0 = MemFree(rcp->comp0);

	return MemFree(rcp);
}

/*
	Store the composition of a (query) string.  
*/
Int2 LIBCALL
BlastResCompStr(BLAST_ScoreBlkPtr sbp, BLAST_ResCompPtr rcp, CharPtr str, Int4 length)
{
	CharPtr	lp, lpmax;
	Int2 index;
        Uint1 mask;

	if (sbp == NULL || rcp == NULL || str == NULL)
		return 1;

	if (rcp->alphabet_code != sbp->alphabet_code)
		return 1;

        /* For megablast, check only the first 4 bits of the sequence values */
        mask = (sbp->protein_alphabet ? 0xff : 0x0f);

/* comp0 starts at zero and extends for "num", comp is the same array, but 
"start_at" offset. */
	for (index=0; index<(sbp->alphabet_size); index++)
		rcp->comp0[index] = 0;

	for (lp = str, lpmax = lp+length; lp < lpmax; lp++)
	{
		++rcp->comp[(int)(*lp & mask)];
	}

	/* Don't count ambig. residues. */
	for (index=0; index<sbp->ambig_occupy; index++)
	{
		rcp->comp[sbp->ambiguous_res[index]] = 0;
	}

	return 0;
}

static Int2
BlastScoreChk(BLAST_Score lo, BLAST_Score hi)
{
	if (lo >= 0 || hi <= 0 ||
			lo < BLAST_SCORE_1MIN || hi > BLAST_SCORE_1MAX)
		return 1;

	if (hi - lo > BLAST_SCORE_RANGE_MAX)
		return 1;

	return 0;
}

BLAST_ScoreFreqPtr
BlastScoreFreqNew(BLAST_Score score_min, BLAST_Score score_max)
{
	BLAST_ScoreFreqPtr	sfp;
	BLAST_Score	range;

	if (BlastScoreChk(score_min, score_max) != 0)
		return NULL;

	sfp = (BLAST_ScoreFreqPtr) MemNew(sizeof(BLAST_ScoreFreq));
	if (sfp == NULL)
		return NULL;

	range = score_max - score_min + 1;
	sfp->sprob = (Nlm_FloatHi PNTR) MemNew(sizeof(Nlm_FloatHi) * range);
	if (sfp->sprob == NULL) 
	{
		BlastScoreFreqDestruct(sfp);
		return NULL;
	}

	sfp->sprob0 = sfp->sprob;
	sfp->sprob -= score_min;
	sfp->score_min = score_min;
	sfp->score_max = score_max;
	sfp->obs_min = sfp->obs_max = 0;
	sfp->score_avg = 0.0;
	return sfp;
}

BLAST_ScoreFreqPtr
BlastScoreFreqDestruct(BLAST_ScoreFreqPtr sfp)
{
	if (sfp == NULL)
		return NULL;

	if (sfp->sprob0 != NULL)
		sfp->sprob0 = MemFree(sfp->sprob0);
	sfp = MemFree(sfp);
	return sfp;
}

static Int2
BlastScoreFreqCalc(BLAST_ScoreBlkPtr sbp, BLAST_ScoreFreqPtr sfp, BLAST_ResFreqPtr rfp1, BLAST_ResFreqPtr rfp2)
{
	BLAST_ScorePtr PNTR	matrix;
	BLAST_Score	score, obs_min, obs_max;
	Nlm_FloatHi		score_sum, score_avg;
	Int2 		alphabet_start, alphabet_end, index1, index2;

	if (sbp == NULL || sfp == NULL)
		return 1;

	if (sbp->loscore < sfp->score_min || sbp->hiscore > sfp->score_max)
		return 1;

	for (score = sfp->score_min; score <= sfp->score_max; score++)
		sfp->sprob[score] = 0.0;

	matrix = sbp->matrix;

	alphabet_start = sbp->alphabet_start;
	alphabet_end = alphabet_start + sbp->alphabet_size;
	for (index1=alphabet_start; index1<alphabet_end; index1++)
	{
		for (index2=alphabet_start; index2<alphabet_end; index2++)
		{
			score = matrix[index1][index2];
			if (score >= sbp->loscore) 
			{
				sfp->sprob[score] += rfp1->prob[index1] * rfp2->prob[index2];
			}
		}
	}

	score_sum = 0.;
	obs_min = obs_max = BLAST_SCORE_MIN;
	for (score = sfp->score_min; score <= sfp->score_max; score++)
	{
		if (sfp->sprob[score] > 0.) 
		{
			score_sum += sfp->sprob[score];
			obs_max = score;
			if (obs_min == BLAST_SCORE_MIN)
				obs_min = score;
		}
	}
	sfp->obs_min = obs_min;
	sfp->obs_max = obs_max;

	score_avg = 0.0;
	if (score_sum > 0.0001 || score_sum < -0.0001)
	{
		for (score = obs_min; score <= obs_max; score++) 
		{
			sfp->sprob[score] /= score_sum;
			score_avg += score * sfp->sprob[score];
		}
	}
	sfp->score_avg = score_avg;

	return 0;
}

typedef struct {
		Char	ch;
		Nlm_FloatHi	p;
	} BLAST_LetterProb;

#define STD_AMINO_ACID_FREQS Robinson_prob

#if STD_AMINO_ACID_FREQS == Dayhoff_prob
/*  M. O. Dayhoff amino acid background frequencies   */
static BLAST_LetterProb	Dayhoff_prob[] = {
		{ 'A', 87.13 },
		{ 'C', 33.47 },
		{ 'D', 46.87 },
		{ 'E', 49.53 },
		{ 'F', 39.77 },
		{ 'G', 88.61 },
		{ 'H', 33.62 },
		{ 'I', 36.89 },
		{ 'K', 80.48 },
		{ 'L', 85.36 },
		{ 'M', 14.75 },
		{ 'N', 40.43 },
		{ 'P', 50.68 },
		{ 'Q', 38.26 },
		{ 'R', 40.90 },
		{ 'S', 69.58 },
		{ 'T', 58.54 },
		{ 'V', 64.72 },
		{ 'W', 10.49 },
		{ 'Y', 29.92 }
	};
#endif

#if STD_AMINO_ACID_FREQS == Altschul_prob
/* Stephen Altschul amino acid background frequencies */
static BLAST_LetterProb Altschul_prob[] = {
		{ 'A', 81.00 },
		{ 'C', 15.00 },
		{ 'D', 54.00 },
		{ 'E', 61.00 },
		{ 'F', 40.00 },
		{ 'G', 68.00 },
		{ 'H', 22.00 },
		{ 'I', 57.00 },
		{ 'K', 56.00 },
		{ 'L', 93.00 },
		{ 'M', 25.00 },
		{ 'N', 45.00 },
		{ 'P', 49.00 },
		{ 'Q', 39.00 },
		{ 'R', 57.00 },
		{ 'S', 68.00 },
		{ 'T', 58.00 },
		{ 'V', 67.00 },
		{ 'W', 13.00 },
		{ 'Y', 32.00 }
	};
#endif

#if STD_AMINO_ACID_FREQS == Robinson_prob
/* amino acid background frequencies from Robinson and Robinson */
static BLAST_LetterProb Robinson_prob[] = {
		{ 'A', 78.05 },
		{ 'C', 19.25 },
		{ 'D', 53.64 },
		{ 'E', 62.95 },
		{ 'F', 38.56 },
		{ 'G', 73.77 },
		{ 'H', 21.99 },
		{ 'I', 51.42 },
		{ 'K', 57.44 },
		{ 'L', 90.19 },
		{ 'M', 22.43 },
		{ 'N', 44.87 },
		{ 'P', 52.03 },
		{ 'Q', 42.64 },
		{ 'R', 51.29 },
		{ 'S', 71.20 },
		{ 'T', 58.41 },
		{ 'V', 64.41 },
		{ 'W', 13.30 },
		{ 'Y', 32.16 }
	};
#endif

static BLAST_LetterProb	nt_prob[] = {
		{ 'A', 25.00 },
		{ 'C', 25.00 },
		{ 'G', 25.00 },
		{ 'T', 25.00 }
	};

/*
	Allocates the BLAST_ResFreqPtr and then fills in the frequencies
	in the probabilities.
*/ 

BLAST_ResFreqPtr LIBCALL
BlastResFreqNew(BLAST_ScoreBlkPtr sbp)
{
	BLAST_ResFreqPtr	rfp;

	if (sbp == NULL)
	{
		return NULL;
	}

	rfp = (BLAST_ResFreqPtr) MemNew(sizeof(BLAST_ResFreq));
	if (rfp == NULL)
		return NULL;

	rfp->alphabet_code = sbp->alphabet_code;

	rfp->prob0 = (Nlm_FloatHi PNTR) MemNew(sizeof(Nlm_FloatHi) * sbp->alphabet_size);
	if (rfp->prob0 == NULL) 
	{
		rfp = BlastResFreqDestruct(rfp);
		return rfp;
	}
	rfp->prob = rfp->prob0 - sbp->alphabet_start;

	return rfp;
}

void LIBCALL BlastResFreqFree(BLAST_ResFreqPtr rfp)
{
    MemFree(rfp->prob0);
    MemFree(rfp);

    return;
}

/*
	Normalize the frequencies to "norm".
*/
Int2 LIBCALL
BlastResFreqNormalize(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp, Nlm_FloatHi norm)
{
	Int2	alphabet_stop, index;
	Nlm_FloatHi	sum = 0., p;

	if (norm == 0.)
		return 1;

	alphabet_stop = sbp->alphabet_start + sbp->alphabet_size;
	for (index=sbp->alphabet_start; index<alphabet_stop; index++)
	{
		p = rfp->prob[index];
		if (p < 0.)
			return 1;
		sum += p;
	}
	if (sum <= 0.)
		return 0;

	for (index=sbp->alphabet_start; index<alphabet_stop; index++)
	{
		rfp->prob[index] /= sum;
		rfp->prob[index] *= norm;
	}
	return 0;
}

/*
	Fills a buffer with the 'standard' alphabet (given by 
	STD_AMINO_ACID_FREQS[index].ch).

	Return value is the number of residues in alphabet.  
	Negative returns upon error.
*/

Int2 LIBCALL
BlastGetStdAlphabet (Uint1 alphabet_code, Uint1Ptr residues, Int4 residues_size)

{
	Int2 index;
	SeqMapTablePtr smtp=NULL;

	if (residues_size < DIM(STD_AMINO_ACID_FREQS))
		return -2;

	if (alphabet_code != Seq_code_ncbieaa)
	{
		smtp = SeqMapTableFind(alphabet_code, Seq_code_ncbieaa);
		if (smtp == NULL)
		{
			return -1;
		}
	}

	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
	{
		if (smtp)
		{
	 		residues[index] = SeqMapTableConvert(smtp,  STD_AMINO_ACID_FREQS[index].ch);
		}
		else
		{
			residues[index] = STD_AMINO_ACID_FREQS[index].ch;
		}
	}

	return index;
}

Int2 LIBCALL
BlastResFreqStdComp(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp)
{
	Int2 index, retval;
	Uint1Ptr residues;

	if (sbp->protein_alphabet == TRUE)
	{
		residues = (Uint1Ptr) MemNew(DIM(STD_AMINO_ACID_FREQS)*sizeof(Uint1));
		retval = BlastGetStdAlphabet(sbp->alphabet_code, residues, DIM(STD_AMINO_ACID_FREQS));
		if (retval < 1)
			return retval;

		for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
		{
			rfp->prob[residues[index]] = STD_AMINO_ACID_FREQS[index].p;
		}
		residues = MemFree(residues);
	}
	else
	{	/* beginning of blastna and ncbi2na are the same. */
		/* Only blastna used  for nucleotides. */
		for (index=0; index<DIM(nt_prob); index++) 
		{
			rfp->prob[index] = nt_prob[index].p;
		}
	}

	BlastResFreqNormalize(sbp, rfp, 1.0);

	return 0;
}

CharPtr  LIBCALL
BlastRepresentativeResidues(Int2 length)

{
	CharPtr buffer, ptr;
	Int2 index, total;
	Int4 number;
	total=0;

	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
	{
		total += (Int2) STD_AMINO_ACID_FREQS[index].p;
	}

	buffer = (CharPtr) MemNew((length+1)*sizeof(Char));

	ptr = buffer;
	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
	{
		number = Nint((STD_AMINO_ACID_FREQS[index].p)*((Nlm_FloatHi) length)/((Nlm_FloatHi) total));
		while (number > 0)
		{
			*ptr = STD_AMINO_ACID_FREQS[index].ch;
			ptr++;
			number--;
		}
	}

	return buffer;
}

Int2 LIBCALL
BlastResFreqString(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp, CharPtr string, Int4 length)
{
	BLAST_ResCompPtr rcp;
	
	rcp = BlastResCompNew(sbp);

	BlastResCompStr(sbp, rcp, string, length);

	BlastResFreqResComp(sbp, rfp, rcp);

	rcp = BlastResCompDestruct(rcp);

	return 0;
}

BLAST_ResFreqPtr LIBCALL
BlastResFreqDestruct(BLAST_ResFreqPtr rfp)
{
	if (rfp == NULL)
		return NULL;

	if (rfp->prob0 != NULL)
		MemFree(rfp->prob0);

	rfp = MemFree(rfp);

	return rfp;
}

/*
	Calculate the residue frequencies associated with the provided ResComp
*/
Int2 LIBCALL
BlastResFreqResComp(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp, BLAST_ResCompPtr rcp)
{
	Int2	alphabet_max, index;
	Nlm_FloatHi	sum = 0.;

	if (rfp == NULL || rcp == NULL)
		return 1;

	if (rfp->alphabet_code != rcp->alphabet_code)
		return 1;

	alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
	for (index=sbp->alphabet_start; index<alphabet_max; index++)
		sum += rcp->comp[index];

	if (sum == 0.) {
		BlastResFreqClr(sbp, rfp);
		return 0;
	}

	for (index=sbp->alphabet_start; index<alphabet_max; index++)
		rfp->prob[index] = rcp->comp[index] / sum;

	return 0;
}

Int2 LIBCALL
BlastResFreqClr(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp)
{
	Int2	alphabet_max, index;
 
	alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
	for (index=sbp->alphabet_start; index<alphabet_max; index++)
                rfp->prob[index] = 0.0;
 
        return 0;
}

/** Supported substitution and gap costs with corresponding quality values
 * for nucleotide sequence comparisons.
 * NB: the values 0 and 0 for the gap costs are treated as the defaults used for
 * the greedy gapped extension, i.e. 
 * gap opening = 0, 
 * gap extension = 1/2 match - mismatch.
 * 
 * The fields are:
 * 
 * 1. Gap opening cost,
 * 2. Gap extension cost,
 * 3. Lambda,
 * 4. K,
 * 5. H,
 * 6. Alpha,
 * 7. Beta,
 * 8. Theta
 */

/** Karlin-Altschul parameter values for substitution scores 1 and -5. */
static const array_of_8 blastn_values_1_5[] = {
    { 0, 0, 1.39, 0.747, 1.38, 1.00,  0, 100 },
    { 3, 3, 1.39, 0.747, 1.38, 1.00,  0, 100 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -4. */
static const array_of_8 blastn_values_1_4[] = {
    { 0, 0, 1.383, 0.738, 1.36, 1.02,  0, 100 },
    { 1, 2,  1.36,  0.67,  1.2,  1.1,  0,  98 }, 
    { 0, 2,  1.26,  0.43, 0.90,  1.4, -1,  91 },
    { 2, 1,  1.35,  0.61,  1.1,  1.2, -1,  98 },
    { 1, 1,  1.22,  0.35, 0.72,  1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -7. 
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_7[] = {
    { 0, 0,  0.69, 0.73, 1.34, 0.515,  0, 100 },
    { 2, 4,  0.68, 0.67,  1.2,  0.55,  0,  99 }, 
    { 0, 4,  0.63, 0.43, 0.90,   0.7, -1,  91 },
    { 4, 2, 0.675, 0.62,  1.1,   0.6, -1,  98 },
    { 2, 2,  0.61, 0.35, 0.72,   1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -3. */
static const array_of_8 blastn_values_1_3[] = {
    { 0, 0, 1.374, 0.711, 1.31, 1.05,  0, 100 },
    { 2, 2,  1.37,  0.70,  1.2,  1.1,  0,  99 },
    { 1, 2,  1.35,  0.64,  1.1,  1.2, -1,  98 },
    { 0, 2,  1.25,  0.42, 0.83,  1.5, -2,  91 },
    { 2, 1,  1.34,  0.60,  1.1,  1.2, -1,  97 },
    { 1, 1,  1.21,  0.34, 0.71,  1.7, -2,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -5.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_5[] = {
    { 0, 0, 0.675, 0.65,  1.1,  0.6, -1, 99 },
    { 2, 4,  0.67, 0.59,  1.1,  0.6, -1, 98 },
    { 0, 4,  0.62, 0.39, 0.78,  0.8, -2, 91 },
    { 4, 2,  0.67, 0.61,  1.0, 0.65, -2, 98 },
    { 2, 2,  0.56, 0.32, 0.59, 0.95, -4, 82 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -2. */
static const array_of_8 blastn_values_1_2[] = {
    { 0, 0, 1.28, 0.46, 0.85, 1.5, -2, 96 },
    { 2, 2, 1.33, 0.62,  1.1, 1.2,  0, 99 },
    { 1, 2, 1.30, 0.52, 0.93, 1.4, -2, 97 }, 
    { 0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89 },
    { 3, 1, 1.32, 0.57,  1.0, 1.3, -1, 99 }, 
    { 2, 1, 1.29, 0.49, 0.92, 1.4, -1, 96 }, 
    { 1, 1, 1.14, 0.26, 0.52, 2.2, -5, 85 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -3.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_3[] = {
    { 0, 0,  0.55, 0.21, 0.46,  1.2, -5, 87 },
    { 4, 4,  0.63, 0.42, 0.84, 0.75, -2, 99 },
    { 2, 4, 0.615, 0.37, 0.72, 0.85, -3, 97 },
    { 0, 4,  0.55, 0.21, 0.46,  1.2, -5, 87 },
    { 3, 3, 0.615, 0.37, 0.68,  0.9, -3, 97 },
    { 6, 2,  0.63, 0.42, 0.84, 0.75, -2, 99 },
    { 5, 2, 0.625, 0.41, 0.78,  0.8, -2, 99 },
    { 4, 2,  0.61, 0.35, 0.68,  0.9, -3, 96 },
    { 2, 2, 0.515, 0.14, 0.33, 1.55, -9, 81 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -4. */
static const array_of_8 blastn_values_3_4[] = {
    { 6, 3, 0.389, 0.25, 0.56, 0.7, -5, 95},
    { 5, 3, 0.375, 0.21, 0.47, 0.8, -6, 92},
    { 4, 3, 0.351, 0.14, 0.35, 1.0, -9, 86},
    { 6, 2, 0.362, 0.16, 0.45, 0.8, -4, 88},
    { 5, 2, 0.330, 0.092, 0.28, 1.2, -13, 81},
    { 4, 2, 0.281, 0.046, 0.16, 1.8, -23, 69}
};

/** Karlin-Altschul parameter values for substitution scores 4 and -5. */
static const array_of_8 blastn_values_4_5[] = {
    { 0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74 },
    { 6, 5, 0.28,  0.21, 0.47, 0.6 , -7, 93 },
    { 5, 5, 0.27,  0.17, 0.39, 0.7,  -9, 90 },
    { 4, 5, 0.25,  0.10, 0.31, 0.8, -10, 83 },
    { 3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -1. */
static const array_of_8 blastn_values_1_1[] = {
    { 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
    { 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 }, 
    { 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 }, 
    { 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
    { 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 }, 
    { 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 }, 
    { 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -2. */
static const array_of_8 blastn_values_3_2[] = {
    {  5,  5, 0.208, 0.030, 0.072, 2.9, -47, 77}
};

/** Karlin-Altschul parameter values for substitution scores 5 and -4. */
static const array_of_8 blastn_values_5_4[] = {
    { 10, 6, 0.163, 0.068, 0.16, 1.0, -19, 85 },
    {  8, 6, 0.146, 0.039, 0.11, 1.3, -29, 76 }
};

static Int2
s_SplitArrayOf8(const array_of_8* input, const array_of_8** normal, const array_of_8** non_affine, Boolean *split)
{

    if (input == NULL || normal == NULL || non_affine == NULL)
            return -1;

    *normal = NULL;
    *non_affine = NULL;

    if (input[0][0] == 0 && input[0][1] == 0)
    {
            *normal = input+1;
            *non_affine = input;
            *split = TRUE;
    }
    else
    {
            *normal = input;
            *split = FALSE;
    }
    return 0;

}

/** Returns the array of values corresponding to the given match/mismatch
 * scores, the number of supported gap cost combinations and thresholds for 
 * the gap costs, beyond which the ungapped statistics can be applied.
 * @param reward Match reward score [in]
 * @param penalty Mismatch penalty score [in]
 * @param array_size Number of supported combinations for this match/mismatch
 *                   pair [out]
 * @param normal the values for normal (e.g, "affine") gap costs [out]
 * @param non_affine specialized values used for megablast [out]
 * @param gap_open_max Gap opening cost threshold for infinite gap costs [out]
 * @param gap_extend_max Gap extension cost threshold for infinite gap costs [out]
 * @param round_down if set to TRUE only even scores should be used for calculation
 *    of expect value or bit scores [out]
 * @param error_return Pointer to error message [out]
 * @return zero on success, other values if error
 */
static Int2
s_GetNuclValuesArray(Int4 reward, Int4 penalty, Int4* array_size,
                     const array_of_8** normal, const array_of_8** non_affine,
                     Int4* gap_open_max, Int4* gap_extend_max, Boolean* round_down,
                     ValNodePtr* error_return)
{
    Int2 status = 0;
    const array_of_8 * kValues = NULL;
    const array_of_8 * kValues_non_affine = NULL;
    Boolean split = FALSE;

    *round_down = FALSE;

    *array_size = 0;

    if (reward == 1 && penalty == -5) {
        if ((status=s_SplitArrayOf8(blastn_values_1_5, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_1_5)/sizeof(array_of_8);
        *gap_open_max = 3;
        *gap_extend_max = 3;
    } else if (reward == 1 && penalty == -4) {
        if ((status=s_SplitArrayOf8(blastn_values_1_4, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *array_size = sizeof(blastn_values_1_4)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -7) {
        if ((status=s_SplitArrayOf8(blastn_values_2_7, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_7)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 4;
    } else if (reward == 1 && penalty == -3) { 
        if ((status=s_SplitArrayOf8(blastn_values_1_3, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *array_size = sizeof(blastn_values_1_3)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -5) {
        if ((status=s_SplitArrayOf8(blastn_values_2_5, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_5)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 4;
    } else if (reward == 1 && penalty == -2) {
        if ((status=s_SplitArrayOf8(blastn_values_1_2, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *array_size = sizeof(blastn_values_1_2)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -3) {
        if ((status=s_SplitArrayOf8(blastn_values_2_3, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_3)/sizeof(array_of_8);
        *gap_open_max = 6;
        *gap_extend_max = 4;
    } else if (reward == 3 && penalty == -4) {
        if ((status=s_SplitArrayOf8(blastn_values_3_4, &kValues, &kValues_non_affine, &split)))
           return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_3_4)/sizeof(array_of_8);
        *gap_open_max = 6;
        *gap_extend_max = 3;
    } else if (reward == 1 && penalty == -1) {
        if ((status=s_SplitArrayOf8(blastn_values_1_1, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *array_size = sizeof(blastn_values_1_1)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 2;
    } else if (reward == 3 && penalty == -2) {
        if ((status=s_SplitArrayOf8(blastn_values_3_2, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_3_2)/sizeof(array_of_8);
        *gap_open_max = 5;
        *gap_extend_max = 5;
    } else if (reward == 4 && penalty == -5) {
        if ((status=s_SplitArrayOf8(blastn_values_4_5, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *array_size = sizeof(blastn_values_4_5)/sizeof(array_of_8);
        *gap_open_max = 12;
        *gap_extend_max = 8;
    } else if (reward == 5 && penalty == -4) {
        if ((status=s_SplitArrayOf8(blastn_values_5_4, &kValues, &kValues_non_affine, &split)))
           return status;
        
        *array_size = sizeof(blastn_values_5_4)/sizeof(array_of_8);
        *gap_open_max = 25;
        *gap_extend_max = 10;
    } else  { /* Unsupported reward-penalty */
        status = -1;
        if (error_return) {
            char buffer[256];
            sprintf(buffer, "Substitution scores %d and %d are not supported", 
                reward, penalty);
            BlastConstructErrorMessage("s_GetNuclValuesArray", buffer, 2, error_return);
        }
    }
    if (split)
        (*array_size)--;

    *normal = kValues;
    *non_affine = kValues_non_affine;

    return status;
}


/** Returns the beta statistical parameter value, given the nucleotide 
 * substitution scores.
 * @param reward Match reward score [in]
 * @param penalty Mismatch penalty score [in]
 * @return The value of the beta parameter.
 */
static double s_GetUngappedBeta(Int4 reward, Int4 penalty)
{
    double beta = 0;
    if ((reward == 1 && penalty == -1) ||
        (reward == 2 && penalty == -3))
        beta = -2;

    return beta;
}

Int2 BlastKarlinGetNuclAlphaBeta(Int4 reward, Int4 penalty, Int4 gap_open, 
                            Int4 gap_extend, BLAST_KarlinBlkPtr kbp,
                            Boolean gapped_calculation,
                            double *alpha, double *beta)
{
    const int kGapOpenIndex = 0;
    const int kGapExtIndex = 1;
    const int kAlphaIndex = 5;
    const int kBetaIndex = 6;
    Int4 num_combinations = 0;
    Int4 gap_open_max = 0, gap_extend_max = 0;
    Int4 index = 0;
    const array_of_8* kNormal=NULL;
    const array_of_8* kNonAffine=NULL; 
    Boolean round_down = FALSE;
    Int2 status = s_GetNuclValuesArray(reward,
                                       penalty,
                                       &num_combinations, 
                                       &kNormal,
                                       &kNonAffine,
                                       &gap_open_max,
                                       &gap_extend_max,
                                       &round_down,
                                       NULL);
    
    if (status)
       return status;

    ASSERT(alpha && beta && kbp);

    /* For ungapped search return ungapped values of alpha and beta. */
    if (gapped_calculation && kNormal) {
        if (gap_open == 0 && gap_extend == 0 && kNonAffine)
        {
            *alpha = kNonAffine[0][kAlphaIndex];
             *beta = kNonAffine[0][kBetaIndex];
            return 0;
        }
        for (index = 0; index < num_combinations; ++index) {
            if (kNormal[index][kGapOpenIndex] == gap_open && 
                kNormal[index][kGapExtIndex] == gap_extend) {
                *alpha = kNormal[index][kAlphaIndex];
                *beta = kNormal[index][kBetaIndex];
                return 0;
            }
        }
    }

    /* If input values not found in tables, or if this is an ungapped search,
       return the ungapped values of alpha and beta. */
    *alpha = kbp->Lambda/kbp->H;
    *beta = s_GetUngappedBeta(reward, penalty);

    return 0;
}

/** Copies data from in to out.  Both in and out should be 
 * allocated.
 * @param in object to be copied [in]
 * @param out object to be copied to [in|out]
 */
static void
s_KarlinBlkCopy(BLAST_KarlinBlk* in, BLAST_KarlinBlk* out)
{
    ASSERT(in && out);

    out->Lambda = in->Lambda;
    out->K = in->K;
    out->logK = in->logK;
    out->H = in->H;
    out->paramC = in->paramC;

    return;
}

Int2
BlastKarlinBlkNuclGappedCalc(BLAST_KarlinBlk* kbp, Int4 gap_open, 
                              Int4 gap_extend, Int4 reward, Int4 penalty,
                              BLAST_KarlinBlk* kbp_ungap,
                              Boolean* round_down,
                              ValNodePtr* error_return)
{
    const int kGapOpenIndex = 0;
    const int kGapExtIndex = 1;
    const int kLambdaIndex = 2;
    const int kKIndex = 3;
    const int kHIndex = 4;
    int num_combinations = 0;
    int gap_open_max, gap_extend_max;
    const array_of_8* kNormal=NULL;
    const array_of_8* kNonAffine=NULL; 
    Int2 status = s_GetNuclValuesArray(reward,
                                       penalty,
                                       &num_combinations,
                                       &kNormal,
                                       &kNonAffine,
                                       &gap_open_max,
                                       &gap_extend_max,
                                       round_down,
                                       error_return);

    if (status)
       return status;

    ASSERT(kbp && kbp_ungap);


    /* Try to find the table entry corresponding to input gap costs. */
    if (gap_open == 0 && gap_extend == 0 && kNonAffine)
    {
        kbp->Lambda = kNonAffine[0][kLambdaIndex];
        kbp->K = kNonAffine[0][kKIndex];
        kbp->logK = log(kbp->K);
        kbp->H = kNonAffine[0][kHIndex];
    }
    else
    {
        int index=0;
        for (index = 0; index < num_combinations; ++index) {
            if (kNormal[index][kGapOpenIndex] == gap_open &&
                kNormal[index][kGapExtIndex] == gap_extend) {
                kbp->Lambda = kNormal[index][kLambdaIndex];
                kbp->K = kNormal[index][kKIndex];
                kbp->logK = log(kbp->K);
                kbp->H = kNormal[index][kHIndex];
                break;
            }
        }
    
        /* If gap costs are not found in the table, check if they belong to the
        infinite domain, where ungapped values of the parameters can be used. */
        if (index == num_combinations) {
        /* If gap costs are larger than maximal provided in tables, copy
           the values from the ungapped Karlin block. */
            if (gap_open >= gap_open_max && gap_extend >= gap_extend_max) {
                s_KarlinBlkCopy(kbp_ungap, kbp);
            } else if (error_return) {
                char buffer[8192];
                int i=0;
                int len=0;
                /* Unsupported gap costs combination. */
                sprintf(buffer, "Gap existence and extension values %ld and %ld "
                        "are not supported for substitution scores %ld and %ld\n", 
                        (long) gap_open, (long) gap_extend, (long) reward, (long) penalty);
                for (i = 0; i < num_combinations; ++i)
                {
                     len = strlen(buffer);
                     sprintf(buffer+len, "%ld and %ld are supported existence and extension values\n", 
                        (long) kNormal[i][kGapOpenIndex],  (long) kNormal[i][kGapExtIndex]);
                }
                len = strlen(buffer);
                sprintf(buffer+len, "%ld and %ld are supported existence and extension values\n", 
                     (long) gap_open_max, (long) gap_extend_max);
                len = strlen(buffer);
                sprintf(buffer+len, "Any values more stringent than %ld and %ld are supported\n", 
                     (long) gap_open_max, (long) gap_extend_max);
                BlastConstructErrorMessage("BlastKarlinBlkNuclGappedCalc", buffer, 2, error_return);
                return 1;
            }
        }
    }

    return 0;
}


/*
	Deallocates MatrixInfoPtr
*/

static MatrixInfoPtr
MatrixInfoDestruct(MatrixInfoPtr matrix_info)

{
	if (matrix_info == NULL)
		return NULL;

	MemFree(matrix_info->name);
	return MemFree(matrix_info);
}

/*
	Makes New MatrixInfoPtr
*/

static MatrixInfoPtr
MatrixInfoNew(CharPtr name, array_of_8 *values, Int4Ptr prefs, Int4 max_number)

{
	MatrixInfoPtr matrix_info;

	matrix_info = (MatrixInfoPtr) MemNew(sizeof(MatrixInfo));
	matrix_info->name = StringSave(name);
	matrix_info->values = values;
	matrix_info->prefs = prefs;
	matrix_info->max_number_values = max_number;

	return matrix_info;
}

static ValNodePtr
BlastMatrixValuesDestruct(ValNodePtr vnp)

{
	ValNodePtr head;

	head = vnp;
	while (vnp)
	{
		MatrixInfoDestruct((MatrixInfoPtr) vnp->data.ptrvalue);
		vnp = vnp->next;
	}

	head = ValNodeFree(head);

	return head;
}

/*
	ValNodePtr BlastLoadMatrixValues (void)
	
	Loads all the matrix values, returns a ValNodePtr chain that contains
	MatrixInfoPtr's.

*/
static ValNodePtr 
BlastLoadMatrixValues (void)

{
	MatrixInfoPtr matrix_info;
	ValNodePtr retval=NULL;

	matrix_info = MatrixInfoNew("BLOSUM80", blosum80_values, blosum80_prefs, BLOSUM80_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("BLOSUM62", blosum62_values, blosum62_prefs, BLOSUM62_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("BLOSUM50", blosum50_values, blosum50_prefs, BLOSUM50_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("BLOSUM45", blosum45_values, blosum45_prefs, BLOSUM45_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("PAM250", pam250_values, pam250_prefs, PAM250_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("BLOSUM62_20", blosum62_20_values, blosum62_20_prefs, BLOSUM62_20_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("BLOSUM90", blosum90_values, blosum90_prefs, BLOSUM90_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("PAM30", pam30_values, pam30_prefs, PAM30_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	matrix_info = MatrixInfoNew("PAM70", pam70_values, pam70_prefs, PAM70_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	return retval;
}
/*
Int2 LIBCALL
BlastKarlinGetMatrixValues(CharPtr matrix, Int4Ptr open, Int4Ptr extension, FloatHiPtr lambda, FloatHiPtr K, FloatHiPtr H)
	
Obtains arrays of the allowed opening and extension penalties for gapped BLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  Any of these fields that
are not required should be set to NULL.  The Int2 return value is the length of the
arrays.
*/

Int2 LIBCALL
BlastKarlinGetMatrixValues(CharPtr matrix, Int4Ptr PNTR open, Int4Ptr PNTR extension, FloatHiPtr PNTR lambda, FloatHiPtr PNTR K, FloatHiPtr PNTR H, Int4Ptr PNTR pref_flags)

{
	return BlastKarlinGetMatrixValuesEx2(matrix, open, extension, NULL, lambda, K, H, NULL, NULL, pref_flags);

}

/*
Int2 LIBCALL
BlastKarlinGetMatrixValuesEx(CharPtr matrix, Int4Ptr open, Int4Ptr extension, FloatHiPtr lambda, FloatHiPtr K, FloatHiPtr H)
	
Obtains arrays of the allowed opening and extension penalties for gapped BLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  Any of these fields that
are not required should be set to NULL.  The Int2 return value is the length of the
arrays.
*/

Int2 LIBCALL
BlastKarlinGetMatrixValuesEx(CharPtr matrix, Int4Ptr PNTR open, Int4Ptr PNTR extension, Int4Ptr PNTR decline_align, FloatHiPtr PNTR lambda, FloatHiPtr PNTR K, FloatHiPtr PNTR H, Int4Ptr PNTR pref_flags)

{
	return BlastKarlinGetMatrixValuesEx2(matrix, open, extension, decline_align, lambda, K, H, NULL, NULL, pref_flags);

}

/*
Int2 LIBCALL
BlastKarlinGetMatrixValuesEx2(CharPtr matrix, Int4Ptr open, Int4Ptr extension, Int4Ptr decline_align, FloatHiPtr lambda, FloatHiPtr K, FloatHiPtr H)
	
Obtains arrays of the allowed opening and extension penalties for gapped BLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  Any of these fields that
are not required should be set to NULL.  The Int2 return value is the length of the
arrays.
*/

Int2 LIBCALL
BlastKarlinGetMatrixValuesEx2(CharPtr matrix, Int4Ptr PNTR open, Int4Ptr PNTR extension, Int4Ptr PNTR decline_align, FloatHiPtr PNTR lambda, FloatHiPtr PNTR K, FloatHiPtr PNTR H, FloatHiPtr PNTR alpha, FloatHiPtr PNTR beta, Int4Ptr PNTR pref_flags)

{
	array_of_8 *values;
	Boolean found_matrix=FALSE;
	Int4 index, max_number_values=0;
	Int4Ptr open_array=NULL, extension_array=NULL, decline_align_array=NULL, pref_flags_array=NULL, prefs;
	Nlm_FloatHiPtr lambda_array=NULL, K_array=NULL, H_array=NULL, alpha_array=NULL, beta_array=NULL;
	MatrixInfoPtr matrix_info;
	ValNodePtr vnp, head;

	if (matrix == NULL)
		return 0;

	vnp = head = BlastLoadMatrixValues();

	while (vnp)
	{
		matrix_info = vnp->data.ptrvalue;
		if (StringICmp(matrix_info->name, matrix) == 0)
		{
			values = matrix_info->values;
			max_number_values = matrix_info->max_number_values;
			prefs = matrix_info->prefs;
			found_matrix = TRUE;
			break;
		}
		vnp = vnp->next;
	}

	if (found_matrix)
	{
		if (open)
			*open = open_array = MemNew(max_number_values*sizeof(Int4));
		if (extension)
			*extension = extension_array = MemNew(max_number_values*sizeof(Int4));
		if (decline_align)
			*decline_align = decline_align_array = MemNew(max_number_values*sizeof(Int4));
		if (lambda)
			*lambda = lambda_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (K)
			*K = K_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (H)
			*H = H_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (alpha)
			*alpha = alpha_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (beta)
			*beta = beta_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (pref_flags)
			*pref_flags = pref_flags_array = MemNew(max_number_values*sizeof(Int4));

		for (index=0; index<max_number_values; index++)
		{
			if (open)
				open_array[index] = values[index][0];
			if (extension)
				extension_array[index] = values[index][1];
			if (decline_align)
				decline_align_array[index] = values[index][2];
			if (lambda)
				lambda_array[index] = values[index][3];
			if (K)
				K_array[index] = values[index][4];
			if (H)
				H_array[index] = values[index][5];
			if (alpha)
				alpha_array[index] = values[index][6];
			if (beta)
				beta_array[index] = values[index][7];
			if (pref_flags)
				pref_flags_array[index] = prefs[index];
		}
	}

	BlastMatrixValuesDestruct(head);

	return max_number_values;
}

/*Extract the alpha and beta settings for this matrixName, and these
  gap open and gap extension costs*/
void LIBCALL getAlphaBeta(CharPtr matrixName, Nlm_FloatHi *alpha,
Nlm_FloatHi *beta, Boolean gapped, Int4 gap_open, Int4 gap_extend)
{
   Int4Ptr gapOpen_arr, gapExtend_arr, pref_flags;
   FloatHiPtr alpha_arr, beta_arr;
   Int2 num_values;
   Int4 i; /*loop index*/

   num_values = BlastKarlinGetMatrixValuesEx2(matrixName, &gapOpen_arr, 
     &gapExtend_arr, NULL, NULL, NULL, NULL,  &alpha_arr, &beta_arr, 
     &pref_flags);

   if (gapped) {
     if ((0 == gap_open) && (0 == gap_extend)) {
       for(i = 1; i < num_values; i++) {
	 if(pref_flags[i]==BLAST_MATRIX_BEST) {
	   (*alpha) = alpha_arr[i];
	   (*beta) = beta_arr[i];
	   break;
	 }
       }
     }
     else {
       for(i = 1; i < num_values; i++) {
	 if ((gapOpen_arr[i] == gap_open) &&
	     (gapExtend_arr[i] == gap_extend)) {
	   (*alpha) = alpha_arr[i];
	   (*beta) = beta_arr[i];
	   break;
	 }
       }
     }
   }
   else if (num_values > 0) {
     (*alpha) = alpha_arr[0];
     (*beta) = beta_arr[0];
   }

   MemFree(gapOpen_arr);
   MemFree(gapExtend_arr);
   MemFree(pref_flags);
   MemFree(alpha_arr);
   MemFree(beta_arr);
}

/*
  Conveniently return default/best Karling-Altschul parameters for a given matrix.
  
 */

Int2 LIBCALL
BlastKarlinGetDefaultMatrixValues(CharPtr matrix, Int4Ptr open, Int4Ptr extension, FloatHiPtr lambda, FloatHiPtr K, FloatHiPtr H) {
    Int4Ptr gapOpen_arr, gapExtend_arr, pref_flags;
    FloatHiPtr Lambda_arr, Kappa_arr, H_arr;
    Int4 i,n;
    if(matrix==NULL)
        matrix = "BLOSUM62";
    n = BlastKarlinGetMatrixValues(matrix, &gapOpen_arr, &gapExtend_arr, &Lambda_arr, &Kappa_arr, &H_arr,&pref_flags);
    if(n) {
        *open = gapOpen_arr[0];
        *extension = gapExtend_arr[0];
        *K = Kappa_arr[0];
        *lambda = Lambda_arr[0];
        *H = H_arr[0];
        for(i=0;i<n;i++) {
            if(pref_flags[i]==BLAST_MATRIX_PREFERRED) {
                *open = gapOpen_arr[i];
                *extension = gapExtend_arr[i];
                *K = Kappa_arr[i];
                *lambda = Lambda_arr[i];
                *H = H_arr[i];
            } else if(pref_flags[i]==BLAST_MATRIX_BEST) {
                *open = gapOpen_arr[i];
                *extension = gapExtend_arr[i];
                *K = Kappa_arr[i];
                *lambda = Lambda_arr[i];
                *H = H_arr[i];
                i+=n;
            }
        }
        MemFree(gapOpen_arr);
        MemFree(gapExtend_arr);
        MemFree(Kappa_arr);
        MemFree(Lambda_arr);
        MemFree(H_arr);
        MemFree(pref_flags);
        return 1;
    } else
        return 0;
}
/*
	Supplies lambda, H, and K values, as calcualted by Stephen Altschul 
	in Methods in Enzy. (vol 266, page 474).

	if kbp is NULL, then a validation is perfomed.
*/

Int2 LIBCALL
BlastKarlinBlkGappedCalc(BLAST_KarlinBlkPtr kbp, Int4 gap_open, Int4 gap_extend, CharPtr matrix_name, ValNodePtr PNTR error_return)

{

	return BlastKarlinBlkGappedCalcEx(kbp, gap_open, gap_extend, (FloatHi) INT2_MAX, matrix_name, error_return);

}
	
/*
	Supplies lambda, H, and K values, as calcualted by Stephen Altschul 
	in Methods in Enzy. (vol 266, page 474).

	if kbp is NULL, then a validation is perfomed.
*/

Int2 LIBCALL
BlastKarlinBlkGappedCalcEx(BLAST_KarlinBlkPtr kbp, Int4 gap_open, Int4 gap_extend, Int4 decline_align, CharPtr matrix_name, ValNodePtr PNTR error_return)

{
	Char buffer[256];
	Int2 status = BlastKarlinkGapBlkFill(kbp, gap_open, gap_extend, decline_align, matrix_name);

	if (status && error_return)
	{
		if (status == 1)
		{
			MatrixInfoPtr matrix_info;
			ValNodePtr vnp, head; 		

			vnp = head = BlastLoadMatrixValues();

			sprintf(buffer, "%s is not a supported matrix", matrix_name);
			BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);

			while (vnp)
			{
				matrix_info = vnp->data.ptrvalue;
				sprintf(buffer, "%s is a supported matrix", matrix_info->name);
				BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
				vnp = vnp->next;
			}

			BlastMatrixValuesDestruct(head);
		}
		else if (status == 2)
		{
			if (decline_align == INT2_MAX)
				sprintf(buffer, "Gap existence and extension values of %ld and %ld not supported for %s", (long) gap_open, (long) gap_extend, matrix_name);
			else
				sprintf(buffer, "Gap existence, extension and decline-to-align values of %ld, %ld and %ld not supported for %s", (long) gap_open, (long) gap_extend, (long) decline_align, matrix_name);
			BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
			BlastKarlinReportAllowedValues(matrix_name, error_return);
		}
	}

	return status;
}

/*
	Attempts to fill KarlinBlk for given gap opening, extensions etc.  
	Will return non-zero status if that fails.

	return values:	-1 if matrix_name is NULL;
			1 if matrix not found
			2 if matrix found, but open, extend etc. values not supported.
*/
Int2 LIBCALL
BlastKarlinkGapBlkFill(BLAST_KarlinBlkPtr kbp, Int4 gap_open, Int4 gap_extend, Int4 decline_align, CharPtr matrix_name)
{
	Boolean found_matrix=FALSE, found_values=FALSE;
	array_of_8 *values;
	Int2 status=0;
	Int4 index, max_number_values=0;
	MatrixInfoPtr matrix_info;
	ValNodePtr vnp, head;
	
	if (matrix_name == NULL)
		return -1;

	values = NULL;

	vnp = head = BlastLoadMatrixValues();
	while (vnp)
	{
		matrix_info = vnp->data.ptrvalue;
		if (StringICmp(matrix_info->name, matrix_name) == 0)
		{
			values = matrix_info->values;
			max_number_values = matrix_info->max_number_values;
			found_matrix = TRUE;
			break;
		}
		vnp = vnp->next;
	}


	if (found_matrix)
	{
		for (index=0; index<max_number_values; index++)
		{
			if (Nint(values[index][0]) == gap_open &&
				Nint(values[index][1]) == gap_extend &&
				(Nint(values[index][2]) == INT2_MAX || Nint(values[index][2]) == decline_align))
			{
				if (kbp)
				{
					kbp->Lambda = values[index][3];
					kbp->K = values[index][4];
					kbp->logK = log(kbp->K);
					kbp->H = values[index][5];
				}
				found_values = TRUE;
				break;
			}
		}

		if (found_values == TRUE)
		{
			status = 0;
		}
		else
		{
			status = 2;
		}
	}
	else
	{
		status = 1;
	}

	BlastMatrixValuesDestruct(head);

	return status;
}

CharPtr
PrintMatrixMessage(const Char *matrix_name)
{
	CharPtr buffer=MemNew(1024*sizeof(Char));
	CharPtr ptr;
	MatrixInfoPtr matrix_info;
        ValNodePtr vnp, head;

	ptr = buffer;
        sprintf(ptr, "%s is not a supported matrix, supported matrices are:\n", matrix_name);

	ptr += StringLen(ptr);

        vnp = head = BlastLoadMatrixValues();

        while (vnp)
        {
        	matrix_info = vnp->data.ptrvalue;
        	sprintf(ptr, "%s \n", matrix_info->name);
		ptr += StringLen(ptr);
		vnp = vnp->next;
        }
        BlastMatrixValuesDestruct(head);

	return buffer;
}

CharPtr
PrintAllowedValuesMessage(const Char *matrix_name, Int4 gap_open, Int4 gap_extend, Int4 decline_align) 
{
	array_of_8 *values;
	Boolean found_matrix=FALSE;
	CharPtr buffer, ptr;
	Int4 index, max_number_values=0;
	MatrixInfoPtr matrix_info;
	ValNodePtr vnp, head;

	ptr = buffer = MemNew(2048*sizeof(Char));

        if (decline_align == INT2_MAX)
              sprintf(ptr, "Gap existence and extension values of %ld and %ld not supported for %s\nsupported values are:\n", 
		(long) gap_open, (long) gap_extend, matrix_name);
        else
               sprintf(ptr, "Gap existence, extension and decline-to-align values of %ld, %ld and %ld not supported for %s\n supported values are:\n", 
		(long) gap_open, (long) gap_extend, (long) decline_align, matrix_name);

	ptr += StringLen(ptr);

	vnp = head = BlastLoadMatrixValues();
	while (vnp)
	{
		matrix_info = vnp->data.ptrvalue;
		if (StringICmp(matrix_info->name, matrix_name) == 0)
		{
			values = matrix_info->values;
			max_number_values = matrix_info->max_number_values;
			found_matrix = TRUE;
			break;
		}
		vnp = vnp->next;
	}

	if (found_matrix)
	{
		for (index=0; index<max_number_values; index++)
		{
			if (Nint(values[index][2]) == INT2_MAX)
				sprintf(ptr, "%ld, %ld\n", (long) Nint(values[index][0]), (long) Nint(values[index][1]));
			else
				sprintf(ptr, "%ld, %ld, %ld\n", (long) Nint(values[index][0]), (long) Nint(values[index][1]), (long) Nint(values[index][2]));
			ptr += StringLen(ptr);
		}
	}

	BlastMatrixValuesDestruct(head);

	return buffer;
}
	

Int2 LIBCALL
BlastKarlinReportAllowedValues(const Char *matrix_name, ValNodePtr PNTR error_return)
{
	array_of_8 *values;
	Boolean found_matrix=FALSE;
	Char buffer[256];
	Int4 index, max_number_values=0;
	MatrixInfoPtr matrix_info;
	ValNodePtr vnp, head;

	vnp = head = BlastLoadMatrixValues();
	while (vnp)
	{
		matrix_info = vnp->data.ptrvalue;
		if (StringICmp(matrix_info->name, matrix_name) == 0)
		{
			values = matrix_info->values;
			max_number_values = matrix_info->max_number_values;
			found_matrix = TRUE;
			break;
		}
		vnp = vnp->next;
	}

	if (found_matrix)
	{
		for (index=0; index<max_number_values; index++)
		{
			if (Nint(values[index][2]) == INT2_MAX)
				sprintf(buffer, "Gap existence and extension values of %ld and %ld are supported", (long) Nint(values[index][0]), (long) Nint(values[index][1]));
			else
				sprintf(buffer, "Gap existence, extension and decline-to-align values of %ld, %ld and %ld are supported", (long) Nint(values[index][0]), (long) Nint(values[index][1]), (long) Nint(values[index][2]));
			BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
		}
	}

	BlastMatrixValuesDestruct(head);

	return 0;
}

/*
	BlastKarlinLtoH

	Calculate H, the relative entropy of the p's and q's
*/
static Nlm_FloatHi LIBCALL
BlastKarlinLtoH(BLAST_ScoreFreqPtr sfp, Nlm_FloatHi	lambda)
{
	BLAST_Score	score;
	Nlm_FloatHi	H, etonlam, sum, scale;

	Nlm_FloatHi PNTR probs = sfp->sprob;
	BLAST_Score low   = sfp->obs_min,  high  = sfp->obs_max;

	if (lambda < 0.) {
		return -1.;
	}
	if (BlastScoreChk(low, high) != 0) return -1.;

	etonlam = exp( - lambda );
  sum = low * probs[low];
  for( score = low + 1; score <= high; score++ ) {
    sum = score * probs[score] + etonlam * sum;
  }

  scale = Nlm_Powi( etonlam, high );
  if( scale > 0.0 ) {
    H = lambda * sum/scale;
  } else { /* Underflow of exp( -lambda * high ) */
    H = lambda * exp( lambda * high + log(sum) );
  }
	return H;
}

/* 
        Everything below here was (more or less) copied from the old 
        karlin.c and could work separately from the stuff above. 
*/ 

/**************** Statistical Significance Parameter Subroutine ****************

    Version 1.0     February 2, 1990
    Version 1.2     July 6,     1990

    Program by:     Stephen Altschul

    Address:        National Center for Biotechnology Information
                    National Library of Medicine
                    National Institutes of Health
                    Bethesda, MD  20894

    Internet:       altschul@ncbi.nlm.nih.gov

See:  Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
    Significance of Molecular Sequence Features by Using General Scoring
    Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

    Computes the parameters lambda and K for use in calculating the
    statistical significance of high-scoring segments or subalignments.

    The scoring scheme must be integer valued.  A positive score must be
    possible, but the expected (mean) score must be negative.

    A program that calls this routine must provide the value of the lowest
    possible score, the value of the greatest possible score, and a pointer
    to an array of probabilities for the occurrence of all scores between
    these two extreme scores.  For example, if score -2 occurs with
    probability 0.7, score 0 occurs with probability 0.1, and score 3
    occurs with probability 0.2, then the subroutine must be called with
    low = -2, high = 3, and pr pointing to the array of values
    { 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
    pointers to lambda and K; the subroutine will then calculate the values
    of these two parameters.  In this example, lambda=0.330 and K=0.154.

    The parameters lambda and K can be used as follows.  Suppose we are
    given a length N random sequence of independent letters.  Associated
    with each letter is a score, and the probabilities of the letters
    determine the probability for each score.  Let S be the aggregate score
    of the highest scoring contiguous segment of this sequence.  Then if N
    is sufficiently large (greater than 100), the following bound on the
    probability that S is greater than or equal to x applies:

            P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].

    In other words, the p-value for this segment can be written as
    1-exp[-KN*exp(-lambda*S)].

    This formula can be applied to pairwise sequence comparison by assigning
    scores to pairs of letters (e.g. amino acids), and by replacing N in the
    formula with N*M, where N and M are the lengths of the two sequences
    being compared.

    In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
    distinct segments all with score >= S is given by:

                           2             m-1           -y
            1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

    Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
    as the previous formula.

*******************************************************************************/

Int2 LIBCALL
BlastKarlinBlkCalc(BLAST_KarlinBlkPtr kbp, BLAST_ScoreFreqPtr sfp)
{
	

	if (kbp == NULL || sfp == NULL)
		return 1;

	/* Calculate the parameter Lambda */

	kbp->Lambda = BlastKarlinLambdaNR(sfp);
	if (kbp->Lambda < 0.)
		goto ErrExit;


	/* Calculate H */

	kbp->H = BlastKarlinLtoH(sfp, kbp->Lambda);
	if (kbp->H < 0.)
		goto ErrExit;


	/* Calculate K and log(K) */

	kbp->K = BlastKarlinLHtoK(sfp, kbp->Lambda, kbp->H);
	if (kbp->K < 0.)
		goto ErrExit;
	kbp->logK = log(kbp->K);

	/* Normal return */
	return 0;

ErrExit:
	kbp->Lambda = kbp->H = kbp->K = -1.;
#ifdef BLASTKAR_HUGE_VAL
	kbp->logK = BLASTKAR_HUGE_VAL;
#else
	kbp->logK = 1.e30;
#endif
	return 1;
}
#define DIMOFP0	(iterlimit*range + 1)
#define DIMOFP0_MAX (BLAST_KARLIN_K_ITER_MAX*BLAST_SCORE_RANGE_MAX+1)


#define smallLambdaThreshold 20 /*defines special case in K computation*/
                                /*threshold is on exp(-Lambda)*/

/*The following procedure computes K. The input includes Lambda, H,
 *  and an array of probabilities for each score.
 *  There are distinct closed form for three cases:
 *  1. high score is 1 low score is -1
 *  2. high score is 1 low score is not -1
 *  3. low score is -1, high score is not 1
 *
 * Otherwise, in most cases the value is computed as:
 * -exp(-2.0*outerSum) / ((H/lambda)*(exp(-lambda) - 1)
 * The last term (exp(-lambda) - 1) can be computed in two different
 * ways depending on whether lambda is small or not.
 * outerSum is a sum of the terms
 * innerSum/j, where j is denoted by iterCounter in the code.
 * The sum is truncated when the new term innersum/j i sufficiently small.
 * innerSum is a weighted sum of the probabilities of
 * of achieving a total score i in a gapless alignment,
 * which we denote by P(i,j).
 * of exactly j characters. innerSum(j) has two parts
 * Sum over i < 0  P(i,j)exp(-i * lambda) +
 * Sum over i >=0  P(i,j)
 * The terms P(i,j) are computed by dynamic programming.
 * An earlier version was flawed in that ignored the special case 1
 * and tried to replace the tail of the computation of outerSum
 * by a geometric series, but the base of the geometric series
 * was not accurately estimated in some cases.
 */

Nlm_FloatHi
BlastKarlinLHtoK(BLAST_ScoreFreqPtr sfp, Nlm_FloatHi    lambda, Nlm_FloatHi H)
{
    /*The next array stores the probabilities of getting each possible
      score in an alignment of fixed length; the array is shifted
      during part of the computation, so that
      entry 0 is for score 0.  */
    Nlm_FloatHi     PNTR alignmentScoreProbabilities = NULL;
    BLAST_Score     low;    /* Lowest score (must be negative) */
    BLAST_Score     high;   /* Highest score (must be positive) */
    BLAST_Score     range;  /* range of scores, computed as high - low*/
    Nlm_FloatHi     K;      /* local copy of K  to return*/
    Int4            i;   /*loop index*/
    Int4            iterCounter; /*counter on iterations*/
    BLAST_Score     divisor; /*candidate divisor of all scores with
                               non-zero probabilities*/
    /*highest and lowest possible alignment scores for current length*/
    BLAST_Score   lowAlignmentScore, highAlignmentScore;
    BLAST_Score   first, last; /*loop indices for dynamic program*/
    register Nlm_FloatHi    innerSum;
    Nlm_FloatHi    oldsum, oldsum2;  /* values of innerSum on previous
                                        iterations*/
    Nlm_FloatHi     outerSum;        /* holds sum over j of (innerSum
                                        for iteration j/j)*/

    Nlm_FloatHi    score_avg; /*average score*/
    /*first term to use in the closed form for the case where
      high == 1 or low == -1, but not both*/
    Nlm_FloatHi     firstTermClosedForm;  /*usually store H/lambda*/
    Int4            iterlimit; /*upper limit on iterations*/
    Nlm_FloatHi     sumlimit; /*lower limit on contributions
                                to sum over scores*/

    /*array of score probabilities reindexed so that low is at index 0*/
    Nlm_FloatHi     PNTR probArrayStartLow;

    /*pointers used in dynamic program*/
    Nlm_FloatHi     PNTR ptrP, PNTR ptr1, PNTR ptr2, PNTR ptr1e;
    Nlm_FloatHi     expMinusLambda; /*e^^(-Lambda) */

    if (lambda <= 0. || H <= 0.) {
        /* Theory dictates that H and lambda must be positive, so
         * return -1 to indicate an error */
        return -1.;
    }

    /*Karlin-Altschul theory works only if the expected score
      is negative*/
    if (sfp->score_avg >= 0.0) {
        return -1.;
    }

    low   = sfp->obs_min;
    high  = sfp->obs_max;
    range = high - low;

    probArrayStartLow = &sfp->sprob[low];
    /* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
       Karlin&Altschul (1990) */
    for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) {
        if (probArrayStartLow[i] != 0.0)
            divisor = Nlm_Gcd(divisor, i);
    }

    high   /= divisor;
    low    /= divisor;
    lambda *= divisor;

    range = high - low;

    firstTermClosedForm = H/lambda;
    expMinusLambda      = exp((Nlm_FloatHi) -lambda);

    if (low == -1 && high == 1) {
        K = (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) *
            (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) / 
            sfp->sprob[low*divisor];
        return(K);
    }

    if (low == -1 || high == 1) {
        if (high != 1) {
            score_avg = sfp->score_avg / divisor;
            firstTermClosedForm
                = (score_avg * score_avg) / firstTermClosedForm;
        }
        return firstTermClosedForm * (1.0 - expMinusLambda);
    }

    sumlimit  = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    iterlimit = BLAST_KARLIN_K_ITER_MAX;

    if (DIMOFP0 > DIMOFP0_MAX) {
        return -1.;
    }
    alignmentScoreProbabilities =
        (Nlm_FloatHi PNTR)MemNew(DIMOFP0 *
                                 sizeof(*alignmentScoreProbabilities));
    if (alignmentScoreProbabilities == NULL)
        return -1.;

    outerSum = 0.;
    lowAlignmentScore = highAlignmentScore = 0;
    alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

    for (iterCounter = 0;
         ((iterCounter < iterlimit) && (innerSum > sumlimit));
         outerSum += innerSum /= ++iterCounter) {
        first = last = range;
        lowAlignmentScore  += low;
        highAlignmentScore += high;
        /*dynamic program to compute P(i,j)*/
        for (ptrP = alignmentScoreProbabilities +
                 (highAlignmentScore-lowAlignmentScore);
             ptrP >= alignmentScoreProbabilities;
             *ptrP-- =innerSum) {
            ptr1  = ptrP - first;
            ptr1e = ptrP - last;
            ptr2  = probArrayStartLow + first;
            for (innerSum = 0.; ptr1 >= ptr1e; )
                innerSum += *ptr1--  *  *ptr2++;
            if (first)
                --first;
            if (ptrP - alignmentScoreProbabilities <= range)
                --last;
        }
        /* Horner's rule */
        innerSum = *++ptrP;
        for( i = lowAlignmentScore + 1; i < 0; i++ ) {
            innerSum = *++ptrP + innerSum * expMinusLambda;
        }
        innerSum *= expMinusLambda;

        for (; i <= highAlignmentScore; ++i)
            innerSum += *++ptrP;
        oldsum2 = oldsum;
        oldsum  = innerSum;
    }

#ifdef ADD_GEOMETRIC_TERMS_TO_K
    /*old code assumed that the later terms in sum were
      asymptotically comparable to those of a geometric
      progression, and tried to speed up convergence by
      guessing the estimated ratio between sucessive terms
      and using the explicit terms of a geometric progression
      to speed up convergence. However, the assumption does not
      always hold, and convergenece of the above code is fast
      enough in practice*/
    /* Terms of geometric progression added for correction */
    {
        Nlm_FloatHi     ratio;  /* fraction used to generate the
                                   geometric progression */
        
        ratio = oldsum / oldsum2;
        if (ratio >= (1.0 - sumlimit*0.001)) {
            K = -1.;
            if (alignmentScoreProbabilities != NULL)
                MemFree(alignmentScoreProbabilities);
            return K;
        }
        sumlimit *= 0.01;
        while (innerSum > sumlimit) {
            oldsum   *= ratio;
            outerSum += innerSum = oldsum / ++iterCounter;
        }
    }
#endif

    if (expMinusLambda <  smallLambdaThreshold ) {
        K = -exp((double)-2.0*outerSum) /
            (firstTermClosedForm*(expMinusLambda - 1.0));
    } else {
        K = -exp((double)-2.0*outerSum) /
            (firstTermClosedForm*Expm1(-(double)lambda));
    }

    if (alignmentScoreProbabilities != NULL)
        MemFree(alignmentScoreProbabilities);

    return K;
}


/**
 * Find positive solution to 
 *
 *     sum_{i=low}^{high} exp(i lambda) * probs[i] = 1.
 * 
 * Note that this solution does not exist unless the average score is
 * negative and the largest score that occurs with nonzero probability
 * is positive.
 * 
 * @param probs         probabilities of a score occurring 
 * @param d             the gcd of the possible scores. This equals 1 if
 *                      the scores are not a lattice
 * @param low           the lowest possible score that occurs with
 *                      nonzero probability
 * @param high          the highest possible score that occurs with
 *                      nonzero probability.
 * @param lambda0       an initial guess for lambda
 * @param tolx          the tolerance to which lambda must be computed
 * @param itmax         the maximum number of times the function may be
 *                      evaluated
 * @param maxNewton     the maximum permissible number of Newton
 *                      iterations; after that the computation will proceed
 *                      by bisection.
 * @param *itn          the number of iterations needed to compute Lambda,
 *                      or itmax if Lambda could not be computed.
 *
 * Let phi(lambda) =  sum_{i=low}^{high} exp(i lambda) - 1. Then
 * phi(lambda) may be written
 *
 *     phi(lamdba) = exp(u lambda) f( exp(-lambda) )
 *
 * where f(x) is a polynomial that has exactly two zeros, one at x = 1
 * and one at x = exp(-lamdba).  It is simpler to solve this problem
 * in x = exp(-lambda) than it is to solve it in lambda, because we
 * know that for x, a solution lies in [0,1], and because Newton's
 * method is generally more stable and efficient for polynomials than
 * it is for exponentials.
 * 
 * For the most part, this function is a standard safeguarded Newton
 * iteration: define an interval of uncertainty [a,b] with f(a) > 0
 * and f(b) < 0 (except for the initial value b = 1, where f(b) = 0);
 * evaluate the function and use the sign of that value to shrink the
 * interval of uncertainty; compute a Newton step; and if the Newton
 * step suggests a point outside the interval of uncertainty or fails
 * to decrease the function sufficiently, then bisect.  There are
 * three further details needed to understand the algorithm:
 *
 * 1)  If y the unique solution in [0,1], then f is positive to the left of
 *     y, and negative to the right.  Therefore, we may determine whether
 *     the Newton step -f(x)/f'(x) is moving toward, or away from, y by
 *     examining the sign of f'(x).  If f'(x) >= 0, we bisect instead
 *     of taking the Newton step.
 * 2)  There is a neighborhood around x = 1 for which f'(x) >= 0, so
 *     (1) prevents convergence to x = 1 (and for a similar reason
 *     prevents convergence to x = 0, if the function is incorrectly
 *     called with probs[high] == 0).
 * 3)  Conditions like  fabs(p) < lambda_tolerance * x * (1-x) are used in
 *     convergence criteria because these values translate to a bound
 *     on the relative error in lambda.  This is proved in the
 *     "Blast Scoring Parameters" document that accompanies the BLAST
 *     code.
 *
 * The iteration on f(x) is robust and doesn't overflow; defining a
 * robust safeguarded Newton iteration on phi(lambda) that cannot
 * converge to lambda = 0 and that is protected against overflow is
 * more difficult.  So (despite the length of this comment) the Newton
 * iteration on f(x) is the simpler solution.
 */

static Nlm_FloatHi 
NlmKarlinLambdaNR(Nlm_FloatHi PNTR probs, BLAST_Score d, BLAST_Score low,
                  BLAST_Score high, Nlm_FloatHi lambda0, Nlm_FloatHi tolx,
                  int itmax, int maxNewton, int * itn) 
{
  int k;
  Nlm_FloatHi x0, x, a = 0, b = 1;
  Nlm_FloatHi f = 4;  /* Larger than any possible value of the poly in [0,1] */
  int isNewton = 0; /* we haven't yet taken a Newton step. */

  assert( d > 0 );

	x0 = exp( -lambda0 );
  x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;
  
  for( k = 0; k < itmax; k++ ) { /* all iteration indices k */
    int i;
    Nlm_FloatHi g, fold = f;
    int wasNewton = isNewton; /* If true, then the previous step was a */
                              /* Newton step */
    isNewton  = 0;            /* Assume that this step is not */
    
    /* Horner's rule for evaluating a polynomial and its derivative */
    g = 0;
    f = probs[low];
    for( i = low + d; i < 0; i += d ) {
      g = x * g + f;
      f = f * x + probs[i];
    }
    g = x * g + f;
    f = f * x + probs[0] - 1;
    for( i = d; i <= high; i += d ) {
      g = x * g + f;
      f = f * x + probs[i];
    }
    /* End Horner's rule */

    if( f > 0 ) {
      a = x; /* move the left endpoint */
    } else if( f < 0 ) { 
      b = x; /* move the right endpoint */
    } else { /* f == 0 */
      break; /* x is an exact solution */
    }
    if( b - a < 2 * a * ( 1 - b ) * tolx ) {
      /* The midpoint of the interval converged */
      x = (a + b) / 2; break;
    }

    if( k >= maxNewton ||
        /* If convergence of Newton's method appears to be failing; or */
				( wasNewton && fabs( f ) > .9 * fabs(fold) ) ||  
        /* if the previous iteration was a Newton step but didn't decrease 
         * f sufficiently; or */
        g >= 0 
        /* if a Newton step will move us away from the desired solution */
        ) { /* then */
      /* bisect */
      x = (a + b)/2;
    } else {
      /* try a Newton step */
      double p = - f/g;
      double y = x + p;
      if( y <= a || y >= b ) { /* The proposed iterate is not in (a,b) */
        x = (a + b)/2;
      } else { /* The proposed iterate is in (a,b). Accept it. */
        isNewton = 1;
        x = y;
        if( fabs( p ) < tolx * x * (1-x) ) break; /* Converged */
      } /* else the proposed iterate is in (a,b) */
    } /* else try a Newton step. */ 
  } /* end for all iteration indices k */
	*itn = k; 
  return -log(x)/d;
}


Nlm_FloatHi
BlastKarlinLambdaNR(BLAST_ScoreFreqPtr sfp)
{
	BLAST_Score	low;			/* Lowest score (must be negative)  */
	BLAST_Score	high;			/* Highest score (must be positive) */
	int		itn;
	BLAST_Score	i, d;
	Nlm_FloatHi PNTR	sprob;
	Nlm_FloatHi	returnValue;

	low = sfp->obs_min;
	high = sfp->obs_max;
	if (sfp->score_avg >= 0.) {	/* Expected score must be negative */
		return -1.0;
	}
	if (BlastScoreChk(low, high) != 0) return -1.;
	
	sprob = sfp->sprob;
	/* Find greatest common divisor of all scores */
	for (i = 1, d = -low; i <= high-low && d > 1; ++i) {
		if (sprob[i+low] != 0.0) {
			d = Nlm_Gcd(d, i);
		}
	}
	returnValue =
		NlmKarlinLambdaNR( sprob, d, low, high,
											 BLAST_KARLIN_LAMBDA0_DEFAULT,
											 BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
											 20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn );


	return returnValue;
}

Nlm_FloatHi LIBCALL
impalaKarlinLambdaNR(BLAST_ScoreFreqPtr sfp, Nlm_FloatHi initialLambda)
{
	Nlm_FloatHi returnValue;
	int itn;
	Nlm_FloatHi PNTR	sprob = sfp->sprob;

	if (sfp->score_avg >= 0.) {	/* Expected score must be negative */
		return -1.0;
	}
	{
		Boolean foundPositive = FALSE;
		BLAST_Score	j;
		for(j = 1; j <=sfp->obs_max; j++) {
			if (sprob[j] > 0.0) {
				foundPositive = TRUE;
				break;
			}
		}
		if (!foundPositive) return(-1);
	}

	returnValue =
		NlmKarlinLambdaNR( sprob, 1, sfp->obs_min, sfp->obs_max,
			                 initialLambda, BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
											 20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn );

	return returnValue;
}



/* Compute a divisor used to weight the evalue of a collection of
 * "nsegs" distinct alignments.  These divisors are used to compensate
 * for the effect of choosing the best among multiple collections of
 * alignments.  See
 *
 * Stephen F. Altschul. Evaluating the statitical significance of
 * multiple distinct local alignments. In Suhai, editior, Theoretical
 * and Computational Methods in Genome Research, pages 1-14. Plenum
 * Press, New York, 1997.
 *
 * The "decayrate" parameter of this routine is a value in the
 * interval (0,1). Typical values of decayrate are .1 and .5. */

Nlm_FloatHi LIBCALL
BlastGapDecayDivisor(Nlm_FloatHi decayrate, unsigned nsegs )
{
    return (1. - decayrate) * Nlm_Powi(decayrate, nsegs - 1);
}

/*
	BlastCutoffs

	Calculate the cutoff score, S, and the highest expected score.

	WRG (later modified by TLM).
*/
Int2 LIBCALL
BlastCutoffs(BLAST_ScorePtr S, /* cutoff score */
	Nlm_FloatHi PNTR E, /* expected no. of HSPs scoring at or above S */
	BLAST_KarlinBlkPtr kbp,
	Nlm_FloatHi searchsp, /* size of search space. */
	Nlm_Boolean dodecay,  /* TRUE ==> use gapdecay feature */
  Nlm_FloatHi gap_decay_rate )
{
	BLAST_Score	s = *S, es;
	Nlm_FloatHi	e = *E, esave;
	Boolean		s_changed = FALSE;

	if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.)
		return 1;

	/*
	Calculate a cutoff score, S, from the Expected
	(or desired) number of reported HSPs, E.
	*/
	es = 1;
	esave = e;
	if (e > 0.) 
	{
        if (dodecay) {
            /* Invert the adjustment to the e-value that will be applied
             * to compensate for the effect of choosing the best among
             * multiple alignments */
            if( gap_decay_rate > 0 && gap_decay_rate < 1 ) {
                e *= BlastGapDecayDivisor(gap_decay_rate, 1);
            }
        }
		es = BlastKarlinEtoS_simple(e, kbp, searchsp);
	}
	/*
	Pick the larger cutoff score between the user's choice
	and that calculated from the value of E.
	*/
	if (es > s) {
		s_changed = TRUE;
		*S = s = es;
	}

	/*
		Re-calculate E from the cutoff score, if E going in was too high
	*/
	if (esave <= 0. || !s_changed) 
	{
		e = BlastKarlinStoE_simple(s, kbp, searchsp);
        if (dodecay) {
            /* Weight the e-value to compensate for the effect of
             * choosing the best of more than one collection of
             * distinct alignments */
            if( gap_decay_rate > 0 && gap_decay_rate < 1 ) {
                e /= BlastGapDecayDivisor(gap_decay_rate, 1);
            }
        }
		*E = e;
	}

	return 0;
}

/*
	BlastKarlinEtoS() -- given an Expect value, return the associated cutoff score

	Error return value is BLAST_SCORE_MIN
*/
BLAST_Score LIBCALL
BlastKarlinEtoS(Nlm_FloatHi	E,	/* Expect value */
	BLAST_KarlinBlkPtr	kbp,
	Nlm_FloatHi	qlen,	/* length of query sequence */
	Nlm_FloatHi	dblen)	/* length of database */
{
	return BlastKarlinEtoS_simple(E, kbp, qlen*dblen);
}


/* Smallest float that might not cause a floating point exception in
	S = (BLAST_Score) (ceil( log((Nlm_FloatHi)(K * searchsp / E)) / Lambda ));
below.
*/
#define BLASTKAR_SMALL_FLOAT 1.0e-297
BLAST_Score LIBCALL
BlastKarlinEtoS_simple(Nlm_FloatHi	E,	/* Expect value */
	BLAST_KarlinBlkPtr	kbp,
	Nlm_FloatHi	searchsp)	/* size of search space */
{

	Nlm_FloatHi	Lambda, K, H; /* parameters for Karlin statistics */
	BLAST_Score	S;

	Lambda = kbp->Lambda;
	K = kbp->K;
	H = kbp->H;
	if (Lambda < 0. || K < 0. || H < 0.) 
	{
		return BLAST_SCORE_MIN;
	}

	E = MAX(E, BLASTKAR_SMALL_FLOAT);

	S = (BLAST_Score) (ceil( log((Nlm_FloatHi)(K * searchsp / E)) / Lambda ));
	return S;
}

/*
	BlastKarlinStoP
	
	Calculate the probability (as opposed to expectation)
	of achieving a particular score.

	On error, return value is -1. (same as BlastKarlinStoE()).
*/
Nlm_FloatHi LIBCALL
BlastKarlinStoP(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	qlen,	/* length of query sequence */
		Nlm_FloatHi	dblen)	/* length of database */
{
	return BlastKarlinStoP_simple(S, kbp, qlen*dblen);
}

Nlm_FloatHi LIBCALL
BlastKarlinStoP_simple(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	searchsp)	/* size of search space. */
{
	Nlm_FloatHi	x, p;

	x = BlastKarlinStoE_simple(S, kbp, searchsp);
	if (x == -1.)
		return x;
	p = BlastKarlinEtoP(x);

	return p;
}

/*
	BlastKarlinStoE() -- given a score, return the associated Expect value

	Error return value is -1.
*/
Nlm_FloatHi LIBCALL
BlastKarlinStoE(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	qlen,	/* length of query sequence */
		Nlm_FloatHi	dblen)	/* length of database */
{
	return BlastKarlinStoE_simple(S, kbp, qlen*dblen);
}

Nlm_FloatHi LIBCALL
BlastKarlinStoE_simple(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	searchsp)	/* size of search space. */
{
	Nlm_FloatHi	Lambda, K, H; /* parameters for Karlin statistics */

	Lambda = kbp->Lambda;
	K = kbp->K;
	H = kbp->H;
	if (Lambda < 0. || K < 0. || H < 0.) {
		return -1.;
	}

	return searchsp * exp((Nlm_FloatHi)(-Lambda * S) + kbp->logK);
}

/*
BlastKarlinPtoE -- convert a P-value to an Expect value
When using BlastKarlinPtoE in the context of a database search,
the returned E-value should be multiplied by the effective
length of the database and divided by the effective lnegth of
the subject.
*/
Nlm_FloatHi LIBCALL
BlastKarlinPtoE(Nlm_FloatHi p)
{
        if (p < 0. || p > 1.0) 
	{
                return INT4_MIN;
        }

	if (p == 1)
		return INT4_MAX;

        return -Nlm_Log1p(-p);
}

/*
BlastKarlinEtoP -- convert an Expect value to a P-value
When using BlastKarlinEtoP in the context of a database search,
the input paramter  E-value should be divided by the effective
length of the database and multiplied by the effective lnegth of
the subject, before BlastKarlinEtoP is called.
*/
Nlm_FloatHi LIBCALL
BlastKarlinEtoP(Nlm_FloatHi x)
{
	return -Nlm_Expm1((Nlm_FloatHi)-x);
}

/*
	BlastKarlinStoLen()

	Given a score, return the length expected for an HSP of that score
*/
Nlm_FloatHi LIBCALL
BlastKarlinStoLen(BLAST_KarlinBlkPtr kbp, BLAST_Score S)
{
	return kbp->Lambda * S / kbp->H;
}

static Nlm_FloatHi	tab2[] = { /* table for r == 2 */
0.01669,  0.0249,   0.03683,  0.05390,  0.07794,  0.1111,   0.1559,   0.2146,   
0.2890,   0.3794,   0.4836,   0.5965,   0.7092,   0.8114,   0.8931,   0.9490,   
0.9806,   0.9944,   0.9989
		};

static Nlm_FloatHi	tab3[] = { /* table for r == 3 */
0.9806,   0.9944,   0.9989,   0.0001682,0.0002542,0.0003829,0.0005745,0.0008587,
0.001278, 0.001893, 0.002789, 0.004088, 0.005958, 0.008627, 0.01240,  0.01770,  
0.02505,  0.03514,  0.04880,  0.06704,  0.09103,  0.1220,   0.1612,   0.2097,   
0.2682,   0.3368,   0.4145,   0.4994,   0.5881,   0.6765,   0.7596,   0.8326,   
0.8922,   0.9367,   0.9667,   0.9846,   0.9939,   0.9980
		};

static Nlm_FloatHi	tab4[] = { /* table for r == 4 */
2.658e-07,4.064e-07,6.203e-07,9.450e-07,1.437e-06,2.181e-06,3.302e-06,4.990e-06,
7.524e-06,1.132e-05,1.698e-05,2.541e-05,3.791e-05,5.641e-05,8.368e-05,0.0001237,
0.0001823,0.0002677,0.0003915,0.0005704,0.0008275,0.001195, 0.001718, 0.002457,
0.003494, 0.004942, 0.006948, 0.009702, 0.01346,  0.01853,  0.02532,  0.03431,
0.04607,  0.06128,  0.08068,  0.1051,   0.1352,   0.1719,   0.2157,   0.2669,
0.3254,   0.3906,   0.4612,   0.5355,   0.6110,   0.6849,   0.7544,   0.8168,
0.8699,   0.9127,   0.9451,   0.9679,   0.9827,   0.9915,   0.9963
		};

static Nlm_FloatHi PNTR table[] = { tab2, tab3, tab4 };
static short tabsize[] = { DIM(tab2)-1, DIM(tab3)-1, DIM(tab4)-1 };

static Nlm_FloatHi LIBCALL f PROTO((Nlm_FloatHi,Nlm_VoidPtr));
static Nlm_FloatHi LIBCALL g PROTO((Nlm_FloatHi,Nlm_VoidPtr));


/*
    Estimate the Sum P-value by calculation or interpolation, as appropriate.
	Approx. 2-1/2 digits accuracy minimum throughout the range of r, s.
	r = number of segments
	s = total score (in nats), adjusted by -r*log(KN)
*/
Nlm_FloatHi LIBCALL
BlastSumP(Int4 r, Nlm_FloatHi s)
{
	Int4		i, r1, r2;
	Nlm_FloatHi	a;

	if (r == 1)
		return -Nlm_Expm1(-exp(-s));

	if (r <= 4) {
		if (r < 1)
			return 0.;
		r1 = r - 1;
		if (s >= r*r+r1) {
			a = Nlm_LnGammaInt(r+1);
			return r * exp(r1*log(s)-s-a-a);
		}
		if (s > -2*r) {
			/* interpolate */
			i = (Int4) (a = s+s+(4*r));
			a -= i;
			i = tabsize[r2 = r - 2] - i;
			return a*table[r2][i-1] + (1.-a)*table[r2][i];
		}
		return 1.;
	}

	return BlastSumPCalc(r, s);
}

/*
    BlastSumPCalc

    Evaluate the following Nlm_FloatHi integral, where r = number of segments
    and s = the adjusted score in nats:

                    (r-2)         oo           oo
     Prob(r,s) =   r              -            -   (r-2)
                 -------------   |   exp(-y)  |   x   exp(-exp(x - y/r)) dx dy
                 (r-1)! (r-2)!  U            U
                                s            0
*/
static Nlm_FloatHi
BlastSumPCalc(int r, Nlm_FloatHi s)
{
	int		r1, itmin;
	Nlm_FloatHi	t, d, epsilon;
	Nlm_FloatHi	est_mean, mean, stddev, stddev4;
	Nlm_FloatHi	xr, xr1, xr2, logr;
	Nlm_FloatHi	args[6];

	epsilon = BLAST_SUMP_EPSILON_DEFAULT; /* accuracy for SumP calcs. */

	if (r == 1) {
		if (s > 8.)
			return exp(-s);
		return -Nlm_Expm1(-exp(-s));
	}
	if (r < 1)
		return 0.;

	xr = r;

	if (r < 8) {
		if (s <= -2.3*xr)
			return 1.;
	}
	else if (r < 15) {
			if (s <= -2.5*xr)
				return 1.;
	}
	else if (r < 27) {
			if (s <= -3.0*xr)
				return 1.;
	}
	else if (r < 51) {
			if (s <= -3.4*xr)
				return 1.;
	}
	else if (r < 101) {
			if (s <= -4.0*xr)
				return 1.;
	}

	/* stddev in the limit of infinite r, but quite good for even small r */
	stddev = sqrt(xr);
	stddev4 = 4.*stddev;
	xr1 = r1 = r - 1;

	if (r > 100) {
		/* Calculate lower bound on the mean using inequality log(r) <= r */
		est_mean = -xr * xr1;
		if (s <= est_mean - stddev4)
			return 1.;
	}

	/* mean is rather close to the mode, and the mean is readily calculated */
	/* mean in the limit of infinite r, but quite good for even small r */
	logr = log(xr);
	mean = xr * (1. - logr) - 0.5;
	if (s <= mean - stddev4)
		return 1.;

	if (s >= mean) {
		t = s + 6.*stddev;
		itmin = 1;
	}
	else {
		t = mean + 6.*stddev;
		itmin = 2;
	}

#define ARG_R args[0]
#define ARG_R2 args[1]
#define ARG_ADJ1 args[2]
#define ARG_ADJ2 args[3]
#define ARG_SDIVR args[4]
#define ARG_EPS args[5]

	ARG_R = xr;
	ARG_R2 = xr2 = r - 2;
	ARG_ADJ1 = xr2*logr - Nlm_LnGammaInt(r1) - Nlm_LnGammaInt(r);
	ARG_EPS = epsilon;

	do {
		d = Nlm_RombergIntegrate(g, args, s, t, epsilon, 0, itmin);
#ifdef BLASTKAR_HUGE_VAL
		if (d == BLASTKAR_HUGE_VAL)
			return d;
#endif
	} while (s < mean && d < 0.4 && itmin++ < 4);

	return (d < 1. ? d : 1.);
}

static Nlm_FloatHi LIBCALL
g(Nlm_FloatHi	s, Nlm_VoidPtr	vp)
{
	register Nlm_FloatHi PNTR	args = vp;
	Nlm_FloatHi	mx;
	
	ARG_ADJ2 = ARG_ADJ1 - s;
	ARG_SDIVR = s / ARG_R;	/* = s / r */
	mx = (s > 0. ? ARG_SDIVR + 3. : 3.);
	return Nlm_RombergIntegrate(f, vp, 0., mx, ARG_EPS, 0, 1);
}

static Nlm_FloatHi LIBCALL
f(Nlm_FloatHi	x, Nlm_VoidPtr	vp)
{
	register Nlm_FloatHi PNTR	args = vp;
	register Nlm_FloatHi	y;

	y = exp(x - ARG_SDIVR);
#ifdef BLASTKAR_HUGE_VAL
	if (y == BLASTKAR_HUGE_VAL)
		return 0.;
#endif
	if (ARG_R2 == 0.)
		return exp(ARG_ADJ2 - y);
	if (x == 0.)
		return 0.;
	return exp(ARG_R2*log(x) + ARG_ADJ2 - y);
}

/*
    Calculates the e-value for alignments with "small" gaps (typically
    under fifty residues/basepairs) following ideas of Stephen Altschul's.
*/

Nlm_FloatHi LIBCALL
BlastSmallGapSumE(
    Int4 starting_points,       /* the number of starting points
                                 * permitted between adjacent
                                 * alignments;
                                 * max_overlap + max_gap + 1 */
    Int2 num,                   /* the number of distinct alignments in this
                                 * collection */
    Nlm_FloatHi xsum,           /* the sum of the scores of these alignments
                                 * each weighted by an appropriate value of
                                 * Lambda and logK */
    Int4 query_length,          /* the effective len of the query seq */
    Int4 subject_length,        /* the effective len of the database seq */
    Int8 dblen_eff,             /* the effective len of the database */
    Nlm_FloatHi weight_divisor) /* a divisor used to weight the e-value
                                 * when multiple collections of alignments
                                 * are being considered by the calling
                                 * routine */
{
    Nlm_FloatHi search_space;   /* The effective size of the search space */
    Nlm_FloatHi sum_e;          /* The e-value of this set of alignments */

    if(num == 1) {
        search_space = (Nlm_FloatHi) query_length * (Nlm_FloatHi)dblen_eff;

        sum_e = search_space * exp(-xsum);
    } else {
        Nlm_FloatHi sum_p;      /* The p-value of this set of alignments */

        search_space = (Nlm_FloatHi)subject_length * (Nlm_FloatHi)query_length;

        xsum -=
          log(search_space) + 2 * (num-1)*log((Nlm_FloatHi)starting_points);
        xsum -= LnFactorial((Nlm_FloatHi) num);

        sum_p = BlastSumP(num, xsum);
        sum_e = BlastKarlinPtoE(sum_p) *
            ((Nlm_FloatHi) dblen_eff / (Nlm_FloatHi) subject_length);
    }
    if( weight_divisor == 0 || (sum_e /= weight_divisor) > INT4_MAX ) {
        sum_e = INT4_MAX;
    }

    return sum_e;
}


/*
    Calculates the e-value of a collection multiple distinct
    alignments with asymmetric gaps between the alignments. The gaps
    in one (protein) sequence are typically small (like in
    BlastSmallGapSumE) gap an the gaps in the other (translated DNA)
    sequence are possibly large (up to 4000 bp.)  This routine is used
    for linking HSPs representing exons in the DNA sequence that are
    separated by introns.
*/

Nlm_FloatHi LIBCALL
BlastUnevenGapSumE(
    Int4 query_start_points,    /* the number of starting points in
                                 * the query sequence permitted
                                 * between adjacent alignments;
                                 * max_overlap + max_gap + 1 */
    Int4 subject_start_points,  /* the number of starting points in
                                 * the subject sequence permitted
                                 * between adjacent alignments */
    Int2 num,                   /* the number of distinct alignments in this
                                 * collection */
    Nlm_FloatHi xsum,           /* the sum of the scores of these alignments
                                 * each weighted by an appropriate value of
                                 * Lambda and logK */
    Int4 query_length,          /* the effective len of the query seq */
    Int4 subject_length,        /* the effective len of the database seq */
    Int8 dblen_eff,             /* the effective len of the database */
    Nlm_FloatHi weight_divisor) /* a divisor used to weight the e-value
                                 * when multiple collections of alignments
                                 * are being considered by the calling
                                 * routine */
{
    Nlm_FloatHi search_space;   /* The effective size of the search space */
    Nlm_FloatHi sum_e;          /* The e-value of this set of alignments */

    if( num == 1 ) {
        search_space = (Nlm_FloatHi)query_length * (Nlm_FloatHi)dblen_eff;

        sum_e = search_space * exp(-xsum);
    } else {
        Nlm_FloatHi sum_p;        /* The p-value of this set of alignments */

        search_space = (Nlm_FloatHi)subject_length * (Nlm_FloatHi)query_length;

        xsum -= log(search_space) +
            (num-1)*(log((Nlm_FloatHi) query_start_points) +
                     log((Nlm_FloatHi) subject_start_points));
        xsum -= LnFactorial((Nlm_FloatHi) num);

        sum_p = BlastSumP(num, xsum);

        sum_e = BlastKarlinPtoE(sum_p) *
            ((Nlm_FloatHi) dblen_eff / (Nlm_FloatHi) subject_length);
    }
    if( weight_divisor == 0 || (sum_e /= weight_divisor) > INT4_MAX ) {
        sum_e = INT4_MAX;
    }

    return sum_e;
}

/*
    Calculates the e-value if a collection of distinct alignments with
    arbitrarily large gaps between the alignments
*/

Nlm_FloatHi LIBCALL
BlastLargeGapSumE(
    Int2 num,                   /* the number of distinct alignments in this
                                 * collection */
    Nlm_FloatHi xsum,           /* the sum of the scores of these alignments
                                 * each weighted by an appropriate value of
                                 * Lambda and logK */
    Int4 query_length,          /* the effective len of the query seq */
    Int4 subject_length,        /* the effective len of the database seq */
    Int8 dblen_eff,             /* the effective len of the database */
    Nlm_FloatHi weight_divisor) /* a divisor used to weight the e-value
                                 * when multiple collections of alignments
                                 * are being considered by the calling
                                 * routine */
{
    Nlm_FloatHi sum_p;          /* The p-value of this set of alignments */
    Nlm_FloatHi sum_e;          /* The e-value of this set of alignments */

/* The next two variables are for compatability with Warren's code. */
    Nlm_FloatHi lcl_subject_length;     /* query_length as a float */
    Nlm_FloatHi lcl_query_length;       /* subject_length as a float */

    lcl_query_length = (Nlm_FloatHi) query_length;
    lcl_subject_length = (Nlm_FloatHi) subject_length;

    if( num == 1 ) {
        sum_e = lcl_query_length * (Nlm_FloatHi) dblen_eff * exp(-xsum);
    } else {
        xsum -= num*log(lcl_subject_length*lcl_query_length)
            - LnFactorial((Nlm_FloatHi) num);

        sum_p = BlastSumP(num, xsum);

        sum_e = BlastKarlinPtoE(sum_p) *
            ((Nlm_FloatHi) dblen_eff / (Nlm_FloatHi) subject_length);
    }
    if( weight_divisor == 0 || (sum_e /= weight_divisor) > INT4_MAX ) {
        sum_e = INT4_MAX;
    }

    return sum_e;
}


/********************************************************************
*
*	The following function, from Stephen Altschul, calculates a 
*	pseudoscore from a p-vlaue and n, the product of the database
*	and query lengths.
*	double	p;		 p-value	
*	double	n;		 search space 
*********************************************************************/
/* Move the following constant into blast.h??, or only the last one. */
#define		PSCALE 20.0
#define 	PSEUDO_SCORE_MAX 32767
#define 	SYBASE_MIN 1.0e-300

Int2 
ConvertPtoPseudoS(Nlm_FloatHi p, Nlm_FloatHi n)
{
	Int2	s;
	Nlm_FloatHi	E;

/* If p is 1.0, then E is very large and E/n is about one. */
	if (p > 0.99)
		return 0.5;
/* If p is very small, the highest score should be returned. */
	else if (p < SYBASE_MIN)
		return PSEUDO_SCORE_MAX;

/* E = -ln(1-p); the following calculates it. */
        E = -Nlm_Log1p(-p);

	s= ConvertEtoPseudoS (E, n);

	return s;
}

/*******************************************************************
*
*	This function calculates a pseudo-score from an E value.
*	The searchsp is the product of the query and database
*	lengths.  As the E value is directly related to the search
*	space, this effectively scales the database size out of
*	the calculation of the pseudo-score.
*******************************************************************/
Int2 
ConvertEtoPseudoS(Nlm_FloatHi E, Nlm_FloatHi searchsp)
{
	Int2	s;

/* If E is very small, the highest score should be returned. */
	if (E < SYBASE_MIN)
		return PSEUDO_SCORE_MAX;

	if (E > searchsp)
		return 0;

/* 0.5 is added to make sure this is rounded up, not down. */
	s= 0.5-PSCALE*log(E/searchsp);

	return s;
}


 
/*
Given a pseudoscore, a subroutine for calculating an E-value for a comparison
of size n (the product of the sequence length) is:
*/

Nlm_FloatHi 
ConvertPseudoStoE(Int2 s, Nlm_FloatHi n)
{
	return n*exp(-s/PSCALE);
}

/****************************************************************************
*
*	This function attaches a new ScorePtr to the "*old" one passed
*	in.  If *old is NULL, then *old is set equal to the new ptr.
*	The new pointer (NOT the head of the chain) is returned.
*
*	If the value is an int, then set prob to zero and pass the value in
*	as score; if the value is a Nlm_FloatHi, pass it in as prob.
*
*	The type of score is stored in an Object-id.str
****************************************************************************/

ScorePtr
MakeBlastScore (ScorePtr PNTR old, CharPtr scoretype, Nlm_FloatHi prob, Int4 score)

{
	ScorePtr	scrp, scrp1;

	scrp = ScoreNew();

	if (score) 
	{
		scrp->choice = 1;
		scrp->value.intvalue = score;
	}
	else 
	{
		scrp->choice = 2;
		scrp->value.realvalue = (Nlm_FloatHi) prob;
	}

	scrp->id = ObjectIdNew();
	scrp->id->str = StringSave(scoretype);


	if (*old == NULL)
		*old = scrp;
	else
	{
		for (scrp1=*old; scrp1->next; scrp1=scrp1->next)
			;
		scrp1->next = scrp;	
	}
		

	return scrp;
}

#ifndef TX_MATRIX_SIZE
#define TX_MATRIX_SIZE 128
#endif

Int4Ptr PNTR LIBCALL BlastMatrixToTxMatrix(BLAST_MatrixPtr blast_matrix)
{
   Uint1 i, j, index1, index2;
   Int4Ptr PNTR matrix = blast_matrix->original_matrix;
   Int4Ptr PNTR txmatrix;
   SeqMapTablePtr smtp;
   SeqCodeTablePtr sctp;

   if (!blast_matrix->is_prot || matrix == NULL) 
      return NULL;

   sctp = SeqCodeTableFindObj(Seq_code_ncbistdaa);
   smtp = SeqMapTableFind(Seq_code_ncbieaa, Seq_code_ncbistdaa);

   txmatrix = Malloc(TX_MATRIX_SIZE*sizeof(Int4Ptr));

   for (i=0; i<TX_MATRIX_SIZE; i++) {
      txmatrix[i] = Malloc(TX_MATRIX_SIZE*sizeof(Int4));
      for (j=0; j<TX_MATRIX_SIZE; j++) 
         txmatrix[i][j] = BLAST_SCORE_MIN;
   }
   
   for (i=sctp->start_at; i < sctp->start_at + sctp->num; i++) {
      for (j=sctp->start_at; j < sctp->start_at + sctp->num; j++) { 
         index1 = SeqMapTableConvert(smtp, i);
         index2 = SeqMapTableConvert(smtp, j);
         txmatrix[index1][index2] = matrix[i][j];
      }
   }

   return txmatrix;
}

Int4Ptr PNTR LIBCALL TxMatrixDestruct(Int4Ptr PNTR txmatrix) 
{
   Int2 i;

   if (txmatrix == NULL)
      return NULL;

   for (i=0; i<TX_MATRIX_SIZE; i++)
      MemFree(txmatrix[i]);

   MemFree(txmatrix);
   return NULL;
}


/** 
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues. 
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta + 
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 * 
 * The value of the length adjustment computed by this routine, A, 
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that 
 *
 *    K * (m - A) * (n - N * A) > min(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge. 
 *
 * @param K      the statistical parameter K
 * @param logK   the natural logarithm of K
 * @param alpha_d_lambda    the ratio of the statistical parameters 
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used)
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0)
 * @param query_length      the length of the query sequence
 * @param db_length         the length of the database
 * @param db_num_seq        the number of sequences in the database
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
 */
Int4
BlastComputeLengthAdjustment(Nlm_FloatHi K,
                             Nlm_FloatHi logK,
                             Nlm_FloatHi alpha_d_lambda,
                             Nlm_FloatHi beta,
                             Int4 query_length,
                             Int8 db_length,
                             Int4 db_num_seqs,
                             Int4 * length_adjustment)
{
    Int4 i;                     /* iteration index */
    const Int4 maxits = 20;     /* maximum allowed iterations */
    Nlm_FloatHi m = query_length, n = db_length, N = db_num_seqs;

    Nlm_FloatHi ell;            /* A float value of the length adjustment */
    Nlm_FloatHi ss;             /* effective size of the search space */
    Nlm_FloatHi ell_min = 0, ell_max;   /* At each iteration i,
                                         * ell_min <= ell <= ell_max. */
    Boolean converged    = FALSE;       /* True if the iteration converged */
    Nlm_FloatHi ell_next = 0;   /* Value the variable ell takes at iteration
                                 * i + 1 */
    /* Choose ell_max to be the largest nonnegative value that satisfies
     *
     *    K * (m - ell) * (n - N * ell) > max(m,n)
     *
     * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */
    { /* scope of a, mb, and c, the coefficients in the quadratic formula
       * (the variable mb is -b) */
        Nlm_FloatHi a  = N;
        Nlm_FloatHi mb = m * N + n;
        Nlm_FloatHi c  = n * m - MAX(m, n) / K;

        if(c < 0) {
            *length_adjustment = 0;
            return 1;
        } else {
            ell_max = 2 * c / (mb + sqrt(mb * mb - 4 * a * c));
        }
    } /* end scope of a, mb and c */

    for(i = 1; i <= maxits; i++) {      /* for all iteration indices */
        Nlm_FloatHi ell_bar;    /* proposed next value of ell */
        ell      = ell_next;
        ss       = (m - ell) * (n - N * ell);
        ell_bar  = alpha_d_lambda * (logK + log(ss)) + beta;
        if(ell_bar >= ell) { /* ell is no bigger than the true fixed point */
            ell_min = ell;
            if(ell_bar - ell_min <= 1.0) {
                converged = TRUE;
                break;
            }
            if(ell_min == ell_max) { /* There are no more points to check */
                break;
            }
        } else { /* else ell is greater than the true fixed point */
            ell_max = ell;
        }
        if(ell_min <= ell_bar && ell_bar <= ell_max) { 
          /* ell_bar is in range. Accept it */
            ell_next = ell_bar;
        } else { /* else ell_bar is not in range. Reject it */
            ell_next = (i == 1) ? ell_max : (ell_min + ell_max) / 2;
        }
    } /* end for all iteration indices */
    if(converged) { /* the iteration converged */
        /* If ell_fixed is the (unknown) true fixed point, then we
         * wish to set (*length_adjustment) to floor(ell_fixed).  We
         * assume that floor(ell_min) = floor(ell_fixed) */
        *length_adjustment = (Int4) ell_min;
        /* But verify that ceil(ell_min) != floor(ell_fixed) */
        ell = ceil(ell_min);
        if( ell <= ell_max ) {
          ss = (m - ell) * (n - N * ell);
          if(alpha_d_lambda * (logK + log(ss)) + beta >= ell) {
            /* ceil(ell_min) == floor(ell_fixed) */
            *length_adjustment = (Int4) ell;
          }
        }
    } else { /* else the iteration did not converge. */
        /* Use the best value seen so far */
        *length_adjustment = (Int4) ell_min;
    }

    return converged ? 0 : 1;
}
