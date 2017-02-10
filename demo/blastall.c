static char const rcsid[] = "$Id: blastall.c,v 6.177 2006/04/26 12:47:48 madden Exp $";

/* $Id: blastall.c,v 6.177 2006/04/26 12:47:48 madden Exp $
**************************************************************************
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
*                                                                         *
************************************************************************** 
 * 
 * $Log: blastall.c,v $
 * Revision 6.177  2006/04/26 12:47:48  madden
 * Use SBlastMessage in place of Blast_Message
 *
 * Revision 6.176  2006/04/25 18:00:19  papadopo
 * change signature of BlastTabularFormatDataNew
 *
 * Revision 6.175  2006/04/21 14:34:50  madden
 * BLAST_GetQuerySeqLoc prototype change
 *
 * Revision 6.174  2006/04/20 15:32:36  papadopo
 * if query IDs are actually used, verify that there are no duplicate IDs
 *
 * Revision 6.173  2006/04/04 13:13:48  madden
 * 1.) Add range check to ARG_FORMAT argument.
 * 2.) Rework error reporting so as to not truncate results if there is a warning.
 *
 * Revision 6.172  2006/03/08 15:41:26  coulouri
 * tune query concatenation limits
 *
 * Revision 6.171  2006/01/24 18:33:44  papadopo
 * from Mike Gertz: Use enumerated values, rather than #define'd constants, to specify the composition adjustment method
 *
 * Revision 6.170  2006/01/23 16:44:05  papadopo
 * change signature of FillHitSavingOptions
 *
 * Revision 6.169  2006/01/13 15:39:17  madden
 * - Enabled the use of composition-based statistics for tblastn (and
 *   only tblastn) using the new engine.
 * - Disabled setting evalue to a larger number when mode > 1
 *   composition-based statistics is used by setting EVALUE_EXPAND to 1.
 * - In Main_new, don't create a BlastSeqSrc.  It is apparently never
 *   used and is a memory leak.  (from Mike Gertz).
 *
 * Revision 6.168  2006/01/10 20:44:10  madden
 * Use SBlastSeqalignArray
 *
 * Revision 6.167  2005/12/22 14:22:19  papadopo
 * change signature of BLAST_FillLookupTableOptions
 *
 * Revision 6.166  2005/12/16 18:30:08  coulouri
 * disable new engine for smith-waterman and composition-based statistics until they are implemented
 *
 * Revision 6.165  2005/12/14 14:43:12  coulouri
 * enable new engine by default
 *
 * Revision 6.164  2005/12/12 13:42:59  madden
 * SBlastOptionsSetRewardPenaltyAndGapCosts now has new greedy Boolean, BLAST_FillScoringOptions no longer called
 *
 * Revision 6.163  2005/10/31 14:15:10  madden
 * Call SBlastOptionsSetRewardPenaltyAndGapCosts
 *
 * Revision 6.162  2005/10/17 14:07:13  madden
 * Use -1 rather than zero for unset gap parameters
 *
 * Revision 6.161  2005/09/29 17:39:29  coulouri
 * from mike gertz:
 *    - Removed the unused static routing GetLambdaFast.
 *    - Removed unused variables from s_FillOptions.
 *    - Removed unused variables from Main_new.
 *    - For tblastn, enabled query concatenation when composition-based
 *      statistics or Smith-Waterman are used.
 *    - Free seq_annot_arr only if query concatenation is being used.
 *    - In Nlm_Main, add preprocessor directives around the declaration
 *      of use_new_engine to suppress a compiler warning.
 *
 * Revision 6.159  2005/09/26 15:02:58  morgulis
 * Fixing some memort leaks when using query concatenation in blastn and tblastn.
 *
 * Revision 6.158  2005/09/16 14:10:03  madden
 * Print out more informative message when Blast_DatabaseSearch has non-zero return if available, add call to SBlastOptionsSetRewardPenaltyAndGapCosts
 *
 * Revision 6.157  2005/09/13 17:39:05  kans
 * include repeats_filter.h
 *
 * Revision 6.156  2005/09/08 14:02:06  coulouri
 * From Mike Gertz:
 *   - Introduced the new options -C for using composition-based
 *     statistics with tblastn; and -s for using Smith-Waterman alignments
 *     with tblastn.
 *   - Forbid the use of the -B option when -C or -s is present; we
 *     expect to remove this restriction.
 *
 * Revision 6.155  2005/09/01 12:28:52  madden
 * 1.) add new function Main_new and put old Main in Main_old, which one is called
 * depends upon the -V option as well as the other params.
 * 2.) Main_new now runs searches with the new engine.
 * 3.) Add headers to allow new engine.
 * 4.) all of the above can be turned off at compile time with a BLASTALL_TOOLS_ONLY define
 *
 * Revision 6.154  2005/08/17 12:42:31  madden
 * Set TXALIGN_SHOW_NO_OF_SEGS for tblastx
 *
 * Revision 6.153  2005/08/08 15:47:41  dondosha
 * Added call to TransTableFreeAll, fixing a memory leak
 *
 * Revision 6.152  2005/06/15 21:37:23  dondosha
 * Do not trigger on-the-fly output with -m8 option for megablast
 *
 * Revision 6.151  2005/05/05 14:41:32  coulouri
 * plug object manager entity id leak - rt ticket 15084082
 *
 * Revision 6.150  2005/02/07 15:30:39  dondosha
 * Removed restriction on the value of longest intron option
 *
 * Revision 6.149  2005/01/10 18:52:28  coulouri
 * fixes from morgulis to allow concatenation of >255 queries in [t]blastn
 *
 * Revision 6.148  2004/09/28 16:06:38  papadopo
 * From Michael Gertz:
 * 1. Disabled ungapped psitblastn.
 * 2. The longest_intron parameter no longer has a minimum value of 4000.
 * 3. Changed the command line help for the longest_intron parameter.
 *
 * Revision 6.147  2004/08/17 17:22:33  madden
 * Add BlastArguments enum for command-line arguments
 *
 * Revision 6.146  2004/07/29 00:05:57  coulouri
 * fix blastcl3 umr
 *
 * Revision 6.145  2004/07/28 18:49:56  coulouri
 * fix printf specifier
 *
 * Revision 6.144  2004/06/30 12:33:30  madden
 * Add include for blfmtutl.h
 *
 * Revision 6.143  2004/05/13 18:42:44  coulouri
 * disable -B for blastcl3
 *
 * Revision 6.142  2004/04/29 19:56:00  dondosha
 * Mask filtered locations in query sequence lines in XML output
 *
 * Revision 6.141  2004/04/20 14:55:47  morgulis
 * 1. Fixed query offsets in results when -B option is used.
 * 2. Fixes for lower case masking handling with -B option.
 *
 * Revision 6.140  2004/03/26 21:42:19  coulouri
 * remove unused variables
 *
 * Revision 6.139  2004/03/18 15:14:21  coulouri
 * do not dereference null seqalignptr
 *
 * Revision 6.138  2004/02/27 14:22:47  coulouri
 * Correct typo
 *
 * Revision 6.137  2004/02/10 18:49:06  coulouri
 * do not allow 1-hit blastn searches
 *
 * Revision 6.136  2003/11/05 22:28:06  dondosha
 * No need to shift subsequence coordinates in tabular output, since they are already shifted in the seqalign
 *
 * Revision 6.135  2003/08/21 15:37:54  dondosha
 * Corrections for out-of-frame tabular output and megablast XML output
 *
 * Revision 6.134  2003/05/30 17:31:09  coulouri
 * add rcsid
 *
 * Revision 6.133  2003/05/09 18:44:49  coulouri
 * make ErrPostEx(SEV_FATAL, ...) exit with nonzero status
 *
 * Revision 6.132  2003/05/06 18:57:46  dondosha
 * Do not set cutoff_s for megablast, it is not needed
 *
 * Revision 6.131  2003/04/08 17:33:42  dondosha
 * Scale the default values of gap costs if match reward is > 1
 *
 * Revision 6.130  2003/04/07 14:46:25  madden
 * Disallow query concatenation if XML, tabular, or ASN.1
 *
 * Revision 6.129  2003/04/01 22:40:09  dondosha
 * Check lower case masking option if megablast option is on
 *
 * Revision 6.128  2003/03/25 15:28:08  dondosha
 * Print tabular output header before checking if seqalign is NULL
 *
 * Revision 6.127  2003/03/24 21:17:08  madden
 * XML fix, remove random printf statements
 *
 * Revision 6.126  2003/03/24 19:43:05  madden
 * Changes to support query concatenation for blastn and tblastn
 *
 * Revision 6.125  2003/03/20 13:44:23  madden
 * Fix -m 10/11 output to make them SeqAnnots
 *
 * Revision 6.124  2002/12/31 22:47:16  boemker
 * Added support for printing output as ASN (text, with -m 10, or binary, with
 * -m 11).
 *
 * Revision 6.123  2002/09/18 20:34:30  camacho
 * Restored -P option
 *
 * Revision 6.122  2002/08/23 16:45:36  madden
 * Issue WARNING for out-of-frame alignments
 *
 * Revision 6.121  2002/08/14 15:09:59  camacho
 * Only change default window size if its command-line value is non-zero
 *
 * Revision 6.120  2002/08/09 19:41:25  camacho
 * 1) Added blast version number to command-line options
 * 2) Added explanations for some default parameters
 *
 * Revision 6.119  2002/06/19 22:50:17  dondosha
 * Added all queries information for tabular output with multiple queries
 *
 * Revision 6.118  2002/05/09 15:37:52  dondosha
 * Call BLASTOptionNewEx instead of BLASTOptionNew, so megablast defaults are set in a central place
 *
 * Revision 6.117  2002/05/04 13:04:43  madden
 * Unsuppress options
 *
 * Revision 6.116  2002/04/29 19:55:26  madden
 * Use ARG_FLOAT for db length
 *
 * Revision 6.115  2002/04/25 21:57:45  madden
 * Strip options for release
 *
 * Revision 6.114  2002/04/25 21:49:28  madden
 * Reset mask_loc_start to NULL for every query
 *
 * Revision 6.113  2002/04/24 19:55:13  madden
 * Rolled back last change
 *
 * Revision 6.112  2002/04/23 20:58:52  madden
 * Suppress options for release
 *
 * Revision 6.111  2002/04/18 20:18:22  dondosha
 * Separate mask locations when formatting results for multiple queries
 *
 * Revision 6.110  2002/04/16 21:10:58  madden
 * Change placement of ReadDBBioseqFetchEnable so db open only once (for HPUX)
 *
 * Revision 6.109  2002/04/16 14:06:00  madden
 * Do not print headers for XML or tabular output
 *
 * Revision 6.108  2002/03/19 23:29:38  dondosha
 * Do not increment options->wordsize by 4 for megablast any more
 *
 * Revision 6.107  2002/02/19 23:21:45  dondosha
 * Fix for XML output if megablast option is used
 *
 * Revision 6.106  2001/12/20 21:51:06  madden
 * Uncomment DO_NOT_SUPPRESS_BLAST_OP
 *
 * Revision 6.105  2001/12/17 20:23:44  madden
 * comment out DO_NOT_SUPPRESS_BLAST_OP
 *
 * Revision 6.104  2001/09/06 20:24:34  dondosha
 * Removed threshold_first
 *
 * Revision 6.103  2001/08/28 17:34:34  madden
 * Add -m 9 as tabular output with comments
 *
 * Revision 6.102  2001/08/28 16:23:12  madden
 * Do not suppress args
 *
 * Revision 6.101  2001/07/27 21:47:35  dondosha
 * Fixed dummy variable declaration for call to StringToInt8
 *
 * Revision 6.100  2001/07/26 18:21:04  dondosha
 * Dummy variable type correction
 *
 * Revision 6.99  2001/07/20 13:31:23  dondosha
 * Undeclared variable correction
 *
 * Revision 6.98  2001/07/19 22:05:47  dondosha
 * Made db_length option a string, to convert to Int8 value
 *
 * Revision 6.97  2001/07/05 15:40:33  madden
 * Comment out DO_NOT_SUPPRESS_BLAST_OP for release
 *
 * Revision 6.96  2001/07/03 20:50:33  madden
 * Commented out call to PrintTabularOutputHeader
 *
 * Revision 6.95  2001/06/21 21:49:55  dondosha
 * No need to declare extra variable vnp
 *
 * Revision 6.94  2001/06/21 21:29:08  dondosha
 * Fixed memory leaks: destroy all error returns, free private_slp
 *
 * Revision 6.93  2001/06/15 21:20:19  dondosha
 * Moved -m9 option to -m8; added header for tabular output
 *
 * Revision 6.92  2001/06/07 19:30:03  dondosha
 * Pass believe query argument to BlastPrintTabulatedResults
 *
 * Revision 6.91  2001/06/06 21:22:44  dondosha
 * Added (query) Bioseq and SeqLoc arguments to function BlastPrintTabulatedResults
 *
 * Revision 6.90  2001/05/25 19:26:36  vakatov
 * Nested comment typo fixed
 *
 * Revision 6.89  2001/05/23 22:38:47  dondosha
 * Added option -m 9 to print post-search tabulated output
 *
 * Revision 6.88  2001/04/10 19:20:52  madden
 * Unsuppress some options suppressed for the release
 *
 * Revision 6.87  2001/04/02 13:52:15  madden
 * Fix for last checkin, properly suppress some options
 *
 * Revision 6.85  2001/03/19 22:39:24  dondosha
 * Allow location on the first query sequence for megablast
 *
 * Revision 6.84  2001/03/13 21:58:23  madden
 * add support for multiple hits blastn, add option for window size
 *
 * Revision 6.83  2001/02/22 20:26:03  dondosha
 * If location stop is -1, make it end of sequence
 *
 * Revision 6.82  2001/02/22 20:11:58  dondosha
 * Previous change reversed; added option to set location on query sequence
 *
 * Revision 6.81  2001/02/22 16:16:43  shavirin
 * Added options for required start and required stop of the query to be
 * used in the Blast search.
 *
 * Revision 6.80  2001/02/22 15:38:48  dondosha
 * Corrected the argument number for longest intron length
 *
 * Revision 6.79  2001/02/09 22:22:36  madden
 * Do not use BlastPruneHitsFromSeqAlign for printing DefLines
 *
 * Revision 6.78  2001/02/08 20:41:17  dondosha
 * Implemented tabulated output for all translated programs
 *
 * Revision 6.77  2001/02/07 21:17:22  dondosha
 * Added support to produce tabulated output (-m 8 option)
 *
 * Revision 6.76  2001/01/19 20:03:47  dondosha
 * Uninitialized variable seqannot caused core dump with XML output
 *
 * Revision 6.75  2000/12/19 18:40:47  madden
 * Add calls to BlastSetUserErrorString and BlastDeleteUserErrorString
 *
 * Revision 6.74  2000/12/15 21:32:12  dondosha
 * Appended getargs explanation of new tblastn (-t) option
 *
 * Revision 6.73  2000/11/21 15:47:21  dondosha
 * Corrected default wordsize for megablast option
 *
 * Revision 6.72  2000/11/17 21:56:26  dondosha
 * Do not free query_lcase_mask in client-server case - already freed
 *
 * Revision 6.71  2000/11/17 20:56:50  dondosha
 * Returned Mega BLAST option which existed in blastcl3 and was removed
 *
 * Revision 6.70  2000/11/17 17:54:50  dondosha
 * Added argument to allow greedy (a la Mega BLAST) extension in blastn
 *
 * Revision 6.69  2000/11/15 15:10:27  shavirin
 * This revision is result of merge between blastall.c and blastcl3.c
 * programs. Using define BLAST_CS_API - client/server version may be
 * created.
 *
 * Revision 6.68  2000/11/09 15:01:00  dondosha
 * Set longest intron length in options in nucleotide coordinates
 *
 * Revision 6.67  2000/11/08 22:24:07  dondosha
 * Enabled new tblastn by adding longest intron option
 *
 * Revision 6.66  2000/11/01 16:26:50  madden
 * Changes from Futamura for psitblastn
 *
 * Revision 6.65  2000/10/27 19:14:40  madden
 * Change description of -b option
 *
 * Revision 6.64  2000/10/23 22:14:04  shavirin
 * Added possibility to pass valid error message into XML output in case
 * of failure or no hits.
 *
 * Revision 6.63  2000/10/23 19:58:22  dondosha
 * Open and close AsnIo outside of call(s) to BXMLPrintOutput
 *
 * Revision 6.62  2000/10/17 19:37:41  shavirin
 * Fixed compilation problems detected on Mac.
 *
 * Revision 6.61  2000/10/17 17:19:49  shavirin
 * Temporary - for toolkit release - commented OOF shift penalty parameter.
 *
 * Revision 6.60  2000/10/06 17:54:28  shavirin
 * Added usage of correct matrix in case of OOF alignment.
 *
 * Revision 6.59  2000/09/26 15:48:15  dondosha
 * Put back printing of header before results of every search when multiple queries are submitted
 *
 * Revision 6.58  2000/09/13 22:26:23  dondosha
 * Removed extra </PRE> that is now printed in PrintDefLinesFromSeqAlign
 *
 * Revision 6.57  2000/09/13 21:39:31  dondosha
 * Corrected html output when input contains multiple queries
 *
 * Revision 6.56  2000/09/12 16:08:43  dondosha
 * Create txalign style matrix from search matrix
 *
 * Revision 6.55  2000/09/12 16:02:13  madden
 * do not allow -P with blastn, fix typo
 *
 * Revision 6.54  2000/09/07 20:25:59  madden
 * Remove L option, turn off K (culling) by default, add -P option
 *
 * Revision 6.53  2000/09/07 16:27:07  shavirin
 * Added option for OOF gap alignment for blastx.
 *
 * Revision 6.52  2000/08/24 14:13:23  shavirin
 * Added return 1 if database do not exists on any path.
 *
 * Revision 6.51  2000/08/11 18:03:58  shavirin
 * Added possibility to make blastx and tblastx with XML output.
 *
 * Revision 6.50  2000/08/11 17:54:08  shavirin
 * Added possibility to print XML output (with -m 7 option)
 *
 * Revision 6.49  2000/08/01 16:35:34  madden
 * Append Seq-annot, do not overwrite
 *
 * Revision 6.48  2000/06/27 15:25:18  madden
 * Changed master-slave to query-anchored
 *
 * Revision 6.47  2000/06/13 19:38:46  shavirin
 * Added ability to print XML Blast output.
 *
 * Revision 6.46  2000/06/05 19:31:31  madden
 * Free query->lcase_mask between searches
 *
 * Revision 6.45  2000/05/26 19:28:44  shavirin
 * Added adjustment of dropoff_1st_pass if dropoff_1st_pass > dropoff_2nd_pass
 *
 * Revision 6.44  2000/05/26 18:48:23  shavirin
 * Added two new parameters; '-y' and '-Z'
 *
 * Revision 6.43  2000/05/09 15:57:26  shavirin
 * Added call to the function ReadDBBioseqSetDbGeneticCode().
 *
 * Revision 6.42  2000/04/25 20:50:45  dondosha
 * Removed unavailable option to use greedy algorithm
 *
 * Revision 6.41  2000/04/13 13:34:19  shavirin
 * Added call to ObjMgrFreeCache() back after fixes in API.
 *
 * Revision 6.40  2000/04/04 18:29:13  shavirin
 * Added some missing HTML tags.
 *
 * Revision 6.39  2000/03/31 19:13:33  dondosha
 * Changed some names related to MegaBlast
 *
 * Revision 6.38  2000/03/24 21:49:30  madden
 * Comment out ObjMgrFreeCache
 *
 * Revision 6.37  2000/03/02 21:06:09  shavirin
 * Added -U option, that allows to consider low characters in FASTA files
 * as filtered regions (for blastn, blastp and tblastn).
 *
 * Revision 6.36  2000/02/01 20:05:31  dondosha
 * Added option -B: use greedy basic alignment search if set to T
 *
 * Revision 6.35  2000/01/28 16:46:54  madden
 * Added function BlastGetMaskingLoc
 *
 * Revision 6.34  1999/12/17 20:48:53  egorov
 * Fix 'gcc -Wall' warnings and remove old stuff.
 *
 * Revision 6.33  1999/10/12 19:35:26  madden
 * Deallocate Mask information
 *
 * Revision 6.32  1999/08/26 14:58:06  madden
 * Use float for db length
 *
 * Revision 6.31  1999/05/26 13:12:56  madden
 * Initialized matrix to NULL
 *
 * Revision 6.30  1999/03/31 16:58:04  madden
 * Removed static FindProt and FindNuc
 *
 * Revision 6.29  1999/02/10 21:12:26  madden
 * Added HTML and GI list option, fixed filtering
 *
 * Revision 6.28  1999/01/22 17:24:51  madden
 * added line breaks for alignment views
 *
 * Revision 6.27  1998/12/31 18:18:27  madden
 * Added strand option
 *
 * Revision 6.26  1998/12/29 20:03:14  kans
 * calls UseLocalAsnloadDataAndErrMsg at startup
 *
 * Revision 6.25  1998/11/19 14:04:34  madden
 * Changed message level to SEV_WARNING
 *
 * Revision 6.24  1998/11/16 16:29:19  madden
 * Added ErrSetMessageLevel(SEV_INFO)
 *
 * Revision 6.23  1998/07/17 15:41:36  madden
 * Added effective search space flag
 *
 * Revision 6.22  1998/06/29 13:02:01  madden
 * Deallocate matrix
 *
 * Revision 6.21  1998/06/10 13:33:14  madden
 * Change -K from zero to 100
 *
 * Revision 6.20  1998/06/05 21:48:42  madden
 * Added -K and -L options
 *
 * Revision 6.19  1998/05/18 18:01:04  madden
 * Changed args to allow filter options to be changed
 *
 * Revision 6.18  1998/05/01 18:31:02  egorov
 * Add new parametes to BLASTOptionSetGapParam()
 *
 * Revision 6.17  1998/04/30 14:32:32  madden
 * init_buff_ex arg changed to 90 for reference
 *
 * Revision 6.16  1998/04/29 14:29:30  madden
 * Made reference line longer
 *
 * Revision 6.15  1998/04/01 22:49:12  madden
 * Print No hits found message
 *
 * Revision 6.14  1998/02/25 20:50:48  madden
 * Added arg for db length
 *
 * Revision 6.13  1998/02/24 22:48:34  madden
 * Removed options for culling
 *
 * Revision 6.12  1998/01/31 21:35:17  madden
 * zeroed out values between searches
 *
 * Revision 6.11  1997/12/31 17:48:52  madden
 * Added wordsize option
 *
 * Revision 6.10  1997/12/23 21:09:47  madden
 * Added -K and -L for range-dependent blast
 *
 * Revision 6.9  1997/11/19 14:26:43  madden
 * Removed extra break statement
 *
 * Revision 6.8  1997/11/18 22:24:22  madden
 * Added call to BLASTOptionSetGapParams
 *
 * Revision 6.7  1997/10/27 22:26:52  madden
 * Added call to ObjMgrFreeCache(0)
 *
 * Revision 6.6  1997/10/23 20:26:12  madden
 * Use of init_buff_ex rather than init_buff
 *
 * Revision 6.5  1997/10/22 21:56:04  madden
 * Added matrix option
 *
 * Revision 6.3  1997/10/07 21:33:38  madden
 * Added BLUNT option
 *
 * Revision 6.2  1997/09/23 22:13:19  madden
 * enabled descriptions and alignment options
 *
 * Revision 6.1  1997/09/16 16:34:32  madden
 * Dbinfo printing changed for multiple db searches
 *
 * Revision 6.0  1997/08/25 18:19:14  madden
 * Revision changed to 6.0
 *
 * Revision 1.16  1997/07/29 19:33:02  madden
 * Added TXALIGN_SHOW_QS flag
 *
 * Revision 1.15  1997/07/28 17:01:23  madden
 * Added include for simutil.h
 *
 * Revision 1.14  1997/07/28 14:31:09  madden
 * Changes for masking alignments.
 *
 * Revision 1.13  1997/07/22 19:06:35  madden
 * Option changes, Printing of verison info
 *
 * Revision 1.12  1997/07/18 20:09:22  madden
 * Conversion from blast2 output to new output
 *
 * Revision 1.3  1997/02/24  22:08:38  madden
 * Added reward and penalty for match and mismatch.
 *
 * Revision 1.2  1997/02/23  16:48:52  madden
 * Call to AcknowledgeBlastQuery added.
 *
 * Revision 1.1  1997/02/19  21:44:28  madden
 * Initial revision
 *
 *
*/


#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <xmlblast.h>
#include <mblast.h>
#include <blfmtutl.h>
#include <algo/blast/composition_adjustment/composition_constants.h>
#ifdef BLAST_CS_API
#include <objblst3.h>
#include <netblap3.h>
#endif
#ifndef BLASTALL_TOOLS_ONLY
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/hspstream_collector.h>
#include <algo/blast/core/blast_stat.h>
#include <algo/blast/api/hspstream_queue.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/api/blast_input.h>
#include <algo/blast/api/blast_format.h>
#include <algo/blast/api/blast_seqalign.h>
#include <algo/blast/api/seqsrc_readdb.h>
#include <algo/blast/api/blast_tabular.h>
#include <algo/blast/api/blast_mtlock.h>
#include <algo/blast/api/blast_prelim.h>
#include <algo/blast/api/blast_api.h>
#include <algo/blast/api/repeats_filter.h>
#endif   /* BLASTALL_TOOLS_ONLY */

#define DEFLINE_BUF 255


/* Used by the callback function. */
FILE *global_fp=NULL;
/*
	Callback to print out ticks, in UNIX only due to file systems
	portability issues.
*/

#ifdef BLAST_CS_API
static  Boolean LIBCALLBACK
tick_callback (BlastResponsePtr brp, Boolean PNTR cancel)
{
    
#if 0
    fprintf(global_fp, ".");
    fflush(global_fp);
#endif

    return TRUE;    
}

#else
static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{
#ifdef OS_UNIX
    /* #ifndef BLAST_CS_API */
    fprintf(global_fp, "%s", ".");
    fflush(global_fp);
    /* #endif */
#endif
    return 0;
}
#endif

static Int2
BlastGetMaskingLoc(FILE *infp, FILE *outfp, CharPtr instructions)
{
	BioseqPtr bsp;
	Char buffer[50];
	SeqEntryPtr sep;
	SeqLocPtr slp, slp_start, tmp_slp;

	if (infp == NULL || outfp == NULL || instructions == NULL)
		return 1;

	while ((sep=FastaToSeqEntryEx(infp, TRUE, NULL, TRUE)) != NULL) 
	{
		bsp = NULL;
		SeqEntryExplore(sep, &bsp, FindNuc);

		if (bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
	   		return 2;
		}
		SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, 50);
		fprintf(outfp, ">%s\n", buffer);
		slp_start = slp = BlastBioseqFilter(bsp, instructions);
        	while (slp)
        	{
               		tmp_slp=NULL;
               		while((tmp_slp = SeqLocFindNext(slp, tmp_slp))!=NULL)
               	 	{
				fprintf(outfp, "%ld %ld\n", (long) (1+SeqLocStart(tmp_slp)), (long) (1+SeqLocStop(tmp_slp)));
                 	}
                	slp = slp->next;
        	}

/* used for debugging. */
#if 0
{{
	BioseqPtr bsp_tmp;
	ByteStorePtr byte_sp;
	Int4 index;
	SeqLocPtr tmp_slp_1, tmp_filter_slp;
	SeqPortPtr spp;
	Uint1Ptr tmp_query_seq, tmp_query_seq_start;
	Uint1 residue;
	FILE *tmp_fp;

		spp = SeqPortNew(bsp, 0, -1, 0, Seq_code_iupacna);
                SeqPortSet_do_virtual(spp, TRUE);
		tmp_query_seq_start = (Uint1Ptr) MemNew(((BioseqGetLen(bsp))+2)*sizeof(Uint1));
		tmp_query_seq_start[0] = NULLB;
		tmp_query_seq = tmp_query_seq_start+1;
		index=0;
		while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{

			if (IS_residue(residue))
			{
				tmp_query_seq[index] = residue;
				index++;
			}
		}
		BlastMaskTheResidues(tmp_query_seq, BioseqGetLen(bsp), 78, slp_start, FALSE, 0);
		bsp_tmp = BioseqNew();
		bsp_tmp->length = BioseqGetLen(bsp);
		byte_sp = BSNew(1);
		BSWrite(byte_sp, tmp_query_seq, bsp->length);
		bsp_tmp->seq_data = byte_sp;
		bsp_tmp->repr = Seq_repr_raw;
		bsp_tmp->seq_data_type = Seq_code_iupacna;
		bsp_tmp->mol = 1;

		bsp_tmp->id = bsp->id;
		bsp_tmp->descr = bsp->descr;

		tmp_fp = FileOpen("masked.fsa", "w");
		BioseqRawToFastaExtra(bsp_tmp, tmp_fp, 50);

		bsp_tmp->id = NULL;
		bsp_tmp->descr = NULL;

		spp = SeqPortFree(spp);
		bsp_tmp = BioseqFree(bsp_tmp);
		tmp_query_seq_start = MemFree(tmp_query_seq_start);
		FileClose(tmp_fp);

		tmp_filter_slp = slp_start;
		tmp_fp = FileOpen("locations.msk", "w");
        	while (tmp_filter_slp)
        	{
               	 tmp_slp_1=NULL;
               	 while((tmp_slp_1 = SeqLocFindNext(tmp_filter_slp, tmp_slp_1))!=NULL)
               	 {
			fprintf(tmp_fp, "%ld %ld\n", (long) (1+SeqLocStart(tmp_slp_1)), (long) (1+SeqLocStop(tmp_slp_1)));

                 }
                	tmp_filter_slp = tmp_filter_slp->next;
        	}


		FileClose(tmp_fp);
}}
#endif
		slp_start = SeqLocSetFree(slp_start);
		sep = SeqEntryFree(sep);
	}

	return 0;
}

/* Breaks up a location like "2000 3000" into two integers 
   that are returned.

   If location is NULL then the integers are set to 0.
*/

/* FIXME: better name, move to API directory?? */
static Boolean
sGetLoc(char* location, Int4* start, Int4* end)
{
        CharPtr delimiters = " ,;";

        if (start == NULL || end == NULL)
           return FALSE;

        *start = 0;
        *end = 0;

        if (location == NULL)
           return TRUE;

        *start =  atoi(StringTokMT(location, delimiters, &location));
        *end = atoi(location);

        return TRUE;
}

typedef enum {
ARG_PROGRAM = 0,
ARG_DB,
ARG_QUERY,
ARG_EVALUE,
ARG_FORMAT,
ARG_OUT,
ARG_FILTER,
ARG_GAPOPEN,
ARG_GAPEXT,
ARG_XDROP,
ARG_SHOWGIS,
ARG_MISMATCH,
ARG_MATCH,
ARG_DESCRIPTIONS,
ARG_ALIGNMENTS,
ARG_THRESHOLD,
ARG_GAPPED,
ARG_QGENETIC_CODE,
ARG_DBGENCODE,
ARG_THREADS, 
ARG_ASNOUT,
ARG_BELIEVEQUERY,
ARG_MATRIX,
ARG_WORDSIZE,
ARG_DBSIZE,
ARG_BESTHITS,
ARG_MULTIPLEHITS,
ARG_SEARCHSP,
ARG_STRAND,
ARG_HTML,
#ifdef BLAST_CS_API
ARG_ENTREZQ,
#else
ARG_GILIST,
#endif
ARG_LCASE,
ARG_XDROP_UNGAPPED,
ARG_XDROP_FINAL,
#ifdef BLAST_CS_API
ARG_RPSBLAST,
#else
ARG_PSITCHKPNT,
#endif
ARG_USEMEGABLAST,
ARG_QUERYLOC,
ARG_WINDOW,
ARG_FRAMESHIFT,
ARG_INTRON,
#ifndef BLAST_CS_API
ARG_NUMQUERIES,
#ifndef BLASTALL_TOOLS_ONLY
ARG_FORCE_OLD,
#endif
#endif
ARG_COMP_BASED_STATS,
ARG_SMITH_WATERMAN
} BlastArguments;

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs[] = {
    { "Program Name",           
      NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},    /* ARG_PROGRAM */
    { "Database",               
      "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},    /* ARG_DB */
    { "Query File",            
      "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL}, /* ARG_QUERY */
    { "Expectation value (E)",  
      "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},    /* ARG_EVALUE */
    { "alignment view options:\n0 = pairwise,\n1 = query-anchored showing identities,\n2 = query-anchored no identities,\n3 = flat query-anchored, show identities,\n4 = flat query-anchored, no identities,\n5 = query-anchored no identities and blunt ends,\n6 = flat query-anchored, no identities and blunt ends,\n7 = XML Blast output,\n8 = tabular, \n9 tabular with comment lines\n10 ASN, text\n11 ASN, binary", /* 4 */
      "0", "0", "11", FALSE, 'm', ARG_INT, 0.0, 0, NULL},         /* ARG_FORMAT */
    { "BLAST report Output File", 
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL}, /* ARG_OUT */
    { "Filter query sequence (DUST with blastn, SEG with others)", 
      "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},       /* ARG_FILTER */
    { "Cost to open a gap (-1 invokes default behavior)", 
      "-1", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},          /* ARG_GAPOPEN */
    { "Cost to extend a gap (-1 invokes default behavior)", 
      "-1", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},          /* ARG_GAPEXT */
    { "X dropoff value for gapped alignment (in bits) (zero invokes default "
      "behavior)\n      blastn 30, megablast 20, tblastx 0, all others 15", 
      "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},          /* ARG_XDROP */
    { "Show GI's in deflines",  /* 10 */
      "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},      /* ARG_SHOWGIS */
    { "Penalty for a nucleotide mismatch (blastn only)", 
      "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},         /* ARG_MISMATCH */
    { "Reward for a nucleotide match (blastn only)", 
      "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},          /* ARG_MATCH */
    { "Number of database sequences to show one-line descriptions for (V)", 
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},         /*  ARG_DESCRIPTIONS */
    { "Number of database sequence to show alignments for (B)", 
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},        /* ARG_ALIGNMENTS */
    { "Threshold for extending hits, default if zero\n" 
      "      blastp 11, blastn 0, blastx 12, tblastn 13\n"
      "      tblastx 13, megablast 0",
      "0", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},           /* ARG_THRESHOLD */
    { "Perform gapped alignment (not available with tblastx)", 
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},     /* ARG_GAPPED */
    { "Query Genetic code to use", /* 17 */
      "1", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},           /* ARG_QGENETIC_CODE */
    { "DB Genetic code (for tblast[nx] only)", /* 18 */
      "1", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},           /* ARG_DBGENCODE */
    { "Number of processors to use", /* 19 */
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},           /* ARG_THREADS */
    { "SeqAlign file",          /* 20 */
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},      /* ARG_ASNOUT */
    { "Believe the query defline", /* 21 */
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},        /* ARG_BELIEVEQUERY */
    { "Matrix",                 /* 22 */
      "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},  /* ARG_MATRIX */
    { "Word size, default if zero (blastn 11, megablast 28, "
        "all others 3)", /* 23 */
      "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},            /* ARG_WORDSIZE */
    { "Effective length of the database (use zero for the real size)", 
      "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},          /* ARG_DBSIZE */
    { "Number of best hits from a region to keep (off by default, if used a value of 100 is recommended)", 
      "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},            /* ARG_BESTHITS */
    { "0 for multiple hit, 1 for single hit (does not apply to blastn)",
       "0",  NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},           /* ARG_MULTIPLEHITS */
    { "Effective length of the search space (use zero for the real size)", 
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},           /* ARG_SEARCHSP */
    { "Query strands to search against database (for blast[nx], and tblastx)\n"
      "       3 is both, 1 is top, 2 is bottom", 
      "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},             /* ARG_STRAND */
    { "Produce HTML output",    /* 29 */
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},         /* ARG_HTML */
#ifdef BLAST_CS_API
    { "Restrict search of database to results of Entrez2 lookup", 
      NULL, NULL, NULL, TRUE, 'u', ARG_STRING, 0.0, 0, NULL},          /* ARG_ENTREZQ */
#else
    { "Restrict search of database to list of GI's",             
      NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},          /* ARG_GILIST */
#endif
    {"Use lower case filtering of FASTA sequence", 
     NULL, NULL, NULL, TRUE, 'U', ARG_BOOLEAN, 0.0, 0, NULL},          /* ARG_LCASE */
    { "X dropoff value for ungapped extensions in bits (0.0 invokes default "
      "behavior)\n      blastn 20, megablast 10, all others 7", 
      "0.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},         /* ARG_XDROP_UNGAPPED */       
    { "X dropoff value for final gapped alignment in bits " 
      "(0.0 invokes default behavior)\n"
      "      blastn/megablast 50, tblastx 0, all others 25",
      "0", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},             /* ARG_XDROP_FINAL */
#ifdef BLAST_CS_API
    { "RPS Blast search",            /* 34 */
      "F", NULL, NULL, FALSE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},          /* ARG_RPSBLAST */
#else
    { "PSI-TBLASTN checkpoint file", /* 34 */
      NULL, NULL, NULL, TRUE, 'R', ARG_FILE_IN, 0.0, 0, NULL},         /* ARG_PSITCHKPNT */
#endif
    { "MegaBlast search",       /* 35 */
      "F", NULL, NULL, FALSE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},         /* ARG_USEMEGABLAST */
    { "Location on query sequence",/* 36 */
      NULL, NULL, NULL, TRUE, 'L', ARG_STRING, 0.0, 0, NULL},          /* ARG_QUERYLOC */
    { "Multiple Hits window size, default if zero (blastn/megablast 0, "
        "all others 40", /* 37 */
      "0", NULL, NULL, FALSE, 'A', ARG_INT, 0.0, 0, NULL},             /* ARG_WINDOW */
    { "Frame shift penalty (OOF algorithm for blastx)", 
      "0", NULL, NULL, FALSE, 'w', ARG_INT, 0.0, 0, NULL},             /* ARG_FRAMESHIFT */
    { "Length of the largest intron allowed in a translated nucleotide "
      "sequence when "
      "linking multiple distinct alignments. (0 invokes default behavior; a "
      "negative value disables linking.)", 
      "0", NULL, NULL, FALSE, 't', ARG_INT, 0.0, 0, NULL},             /* ARG_INTRON */
/*--KM
   seems ok to add another param b/c NUMARG is defined based on 
    sizeof(myargs) itself
   made optional=TRUE but this may change?
*/
#ifndef BLAST_CS_API
    { "Number of concatenated queries, for blastn and tblastn", 
      "0", NULL, NULL, TRUE, 'B', ARG_INT, 0.0, 0, NULL},               /* ARG_NUMQUERIES */
#ifndef BLASTALL_TOOLS_ONLY
    { "Force use of the legacy BLAST engine", 
      "F", NULL, NULL, TRUE, 'V', ARG_BOOLEAN, 0.0, 0, NULL},              /* ARG_FORCE_OLD */
#endif  /* BLASTALL_TOOLS_ONLY */
#endif
    { "Use composition-based statistics for tblastn:\n"                /* ARG_COMP_BASED_STATS */
      "      D or d: default (equivalent to F)\n"
      "      0 or F or f: no composition-based statistics\n"
      "      1 or T or t: Composition-based statistics as in "
      "NAR 29:2994-3005, 2001\n"
      "      2: Composition-based score adjustment as in "
      "Bioinformatics 21:902-911,\n"
      "          2005, conditioned on sequence properties\n"
      "      3: Composition-based score adjustment as in "
      "Bioinformatics 21:902-911,\n"
      "          2005, unconditionally\n"
      "      For programs other than tblastn, must either be absent "
      "or be D, F or 0.\n     ",
      "D", NULL, NULL, FALSE, 'C', ARG_STRING, 0.0, 0, NULL},
    { "Compute locally optimal Smith-Waterman alignments "
        "(This option is only\n"
      "      available for gapped tblastn.)",                          /* ARG_SMITH_WATERMAN */
      "F", NULL, NULL, FALSE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
};

#ifdef BLAST_CS_API
static BlastNet3Hptr BNETInitializeBlast(CharPtr database, CharPtr program, 
                                  FILE *outfp, Boolean db_is_na,
                                  Boolean is_rps_blast, Boolean html, Boolean header)
{
    BlastNet3Hptr    bl3hp;
    BlastResponsePtr response = NULL;
    BlastVersionPtr	blast_version;

    if (! BlastInit("blastcl3", &bl3hp, &response)) {
        ErrPostEx(SEV_FATAL, 1, 0, "Unable to initialize BLAST service");
        return NULL;
    }
    if (response && response->choice == BlastResponse_init) {
        blast_version = response->data.ptrvalue;
    } else {
        ErrPostEx(SEV_FATAL, 1, 0, "Unable to connect to the BLAST service");
        return NULL;
    }
    
    BlastNetBioseqFetchEnable(bl3hp, database, db_is_na, TRUE);
    
    if(is_rps_blast == TRUE && header)
    {
        BlastPrintVersionInfoEx("RPS-BLAST", html, blast_version->version, 
                                blast_version->date, outfp);
    }
    else if (header) 
    {
	init_buff_ex(90);
        BlastPrintVersionInfoEx(program, html, blast_version->version, 
                                blast_version->date, outfp);
        fprintf(outfp, "\n");
        BlastPrintReference(html, 80, outfp);
	free_buff();
    }

    BlastResponseFree(response);

    return bl3hp;
}
#endif

/* Needed for Mega BLAST only */
#define MAX_NUM_QUERIES 16383 /* == 1/2 INT2_MAX */

#ifndef BLASTALL_TOOLS_ONLY

/** Fills all the options structures with user defined values. Uses the 
 * myargs global structure obtained from GetArgs.
 * @param lookup_options Lookup table options [in]
 * @param query_setup_options Query options [in]
 * @param word_options Initial word processing options [in]
 * @param ext_options Extension options [in]
 * @param hit_options Hit saving options [out]
 * @param score_options Scoring options [out]
 * @param eff_len_options Effective length options [out]
 * @param psi_options Protein BLAST options [out]
 * @param db_options BLAST database options [out]
 */
static Int2 
s_FillOptions(SBlastOptions* options)
{
   LookupTableOptions* lookup_options = options->lookup_options;
   QuerySetUpOptions* query_setup_options = options->query_options; 
   BlastInitialWordOptions* word_options = options->word_options;
   BlastExtensionOptions* ext_options = options->ext_options;
   BlastHitSavingOptions* hit_options = options->hit_options ;
   BlastScoringOptions* score_options = options->score_options;
   BlastEffectiveLengthsOptions* eff_len_options = options->eff_len_options;

   Boolean mb_lookup = FALSE;
   Boolean greedy = FALSE;
   Boolean is_gapped = FALSE;
   EBlastProgramType program_number = options->program;

   if (myargs[ARG_USEMEGABLAST].intvalue != 0)
   {
       greedy = TRUE;
       mb_lookup = TRUE;
   }

   BLAST_FillLookupTableOptions(lookup_options, program_number, mb_lookup,
      myargs[ARG_THRESHOLD].intvalue, (Int2)myargs[ARG_WORDSIZE].intvalue);

   BLAST_FillQuerySetUpOptions(query_setup_options, program_number, 
      myargs[ARG_FILTER].strvalue, (Uint1)myargs[ARG_STRAND].intvalue);

   if (myargs[ARG_QGENETIC_CODE].intvalue &&
       (program_number == eBlastTypeBlastx || 
        program_number == eBlastTypeTblastx))
      query_setup_options->genetic_code = myargs[ARG_QGENETIC_CODE].intvalue;

   BLAST_FillInitialWordOptions(word_options, program_number, 
      greedy, myargs[ARG_WINDOW].intvalue, myargs[ARG_XDROP_UNGAPPED].intvalue);

   BLAST_FillExtensionOptions(ext_options, program_number, greedy, 
      myargs[ARG_XDROP].intvalue, myargs[ARG_XDROP_FINAL].intvalue);

   /* if both gap_open and gap_extend are zero then they are set to suggested values */
   SBlastOptionsSetMatrixAndGapCosts(options, myargs[ARG_MATRIX].strvalue,
        myargs[ARG_GAPOPEN].intvalue, myargs[ARG_GAPEXT].intvalue);

   SBlastOptionsSetRewardPenaltyAndGapCosts(options,
        myargs[ARG_MATCH].intvalue,
        myargs[ARG_MISMATCH].intvalue,
        myargs[ARG_GAPOPEN].intvalue,
        myargs[ARG_GAPEXT].intvalue,
        FALSE);

   if (myargs[ARG_WINDOW].intvalue < 0)
       word_options->window_size = 0;
   else
       SBlastOptionsSetWindowSize(options, myargs[ARG_WINDOW].intvalue);

   SBlastOptionsSetThreshold(options, myargs[ARG_THRESHOLD].intvalue);

   if (program_number != eBlastTypeTblastx)
      is_gapped = myargs[ARG_GAPPED].intvalue;
   else
      is_gapped = FALSE;

   score_options->gapped_calculation = is_gapped;
   if (myargs[ARG_FRAMESHIFT].intvalue) {
      score_options->shift_pen = myargs[ARG_FRAMESHIFT].intvalue;
      score_options->is_ooframe = TRUE;
   }

   BLAST_FillHitSavingOptions(hit_options, 
      myargs[ARG_EVALUE].floatvalue, 
      MAX(myargs[ARG_DESCRIPTIONS].intvalue, 
          myargs[ARG_ALIGNMENTS].intvalue),
          is_gapped, 
      0,                /* culling limit */
      0);               /* min diag separation */
 
   hit_options->longest_intron = MIN(myargs[ARG_INTRON].intvalue, MAX_INTRON_LENGTH);

   if (myargs[ARG_SEARCHSP].floatvalue != 0)
      eff_len_options->searchsp_eff = (Int8) myargs[ARG_SEARCHSP].floatvalue; 

   if (program_number == eBlastTypeTblastn ||
       program_number == eBlastTypeRpsTblastn ||
       program_number == eBlastTypeTblastx) {
       SBlastOptionsSetDbGeneticCode(options, myargs[ARG_DBGENCODE].intvalue);
   }
   if (program_number == eBlastTypeTblastn && is_gapped) {
       /* Set options specific to gapped tblastn */
       switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
       case 'D': case 'd':
       case '0': case 'F': case 'f':
           ext_options->compositionBasedStats = eNoCompositionBasedStats;
           break;
       case '1': case 'T': case 't':
           ext_options->compositionBasedStats = eCompositionBasedStats;
           break;
       case '2':
           ErrPostEx(SEV_WARNING, 1, 0, "the -C 2 argument "
                     "is currently experimental\n");
           ext_options->compositionBasedStats = eCompositionMatrixAdjust;
           break;
       case '3':
           ErrPostEx(SEV_WARNING, 1, 0, "the -C 3 argument "
                     "is currently experimental\n");
           ext_options->compositionBasedStats = eCompoForceFullMatrixAdjust;
           break;
       default:
           ErrPostEx(SEV_FATAL, 1, 0, "invalid argument for composition-"
                     "based statistics; see -C options\n");
           break;
       }
       if (myargs[ARG_SMITH_WATERMAN].intvalue) {
           ext_options->eTbackExt = eSmithWatermanTbck;
       }
   } else {
       /* Make sure tblastn and blastp parameters were not set for
        * other programs */
       
       switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
       case '0': case 'D': case 'd': case 'F': case 'f':
           break;
       default:
           ErrPostEx(SEV_FATAL, 1, 0,
                     "Invalid option -C: only gapped tblastn may use"
                     " composition based statistics.");
           break;
       }
       if(myargs[ARG_SMITH_WATERMAN].intvalue) {
           ErrPostEx(SEV_FATAL, 1, 0,
                     "Invalid option -s: Smith-Waterman alignments are only "
                     "available for gapped tblastn and blastp.");
       }
   }

   return 0;
}

#ifndef TX_MATRIX_SIZE
#define TX_MATRIX_SIZE 128
#endif

Int4** LIBCALL BlastMatrixConvert(Int4** old)
{
   Int4 i, j, index1, index2;
   Int4** new;
   SeqMapTablePtr smtp;
   SeqCodeTablePtr sctp;

   if (!old)
      return NULL;

   sctp = SeqCodeTableFindObj(Seq_code_ncbistdaa);
   smtp = SeqMapTableFind(Seq_code_ncbieaa, Seq_code_ncbistdaa);

   new = malloc(TX_MATRIX_SIZE*sizeof(Int4Ptr));

   for (i=0; i<TX_MATRIX_SIZE; i++) {
      new[i] = malloc(TX_MATRIX_SIZE*sizeof(Int4));
      for (j=0; j<TX_MATRIX_SIZE; j++)
         new[i][j] = BLAST_SCORE_MIN;
   }

   for (i=sctp->start_at; i < sctp->start_at + sctp->num; i++) {
      for (j=sctp->start_at; j < sctp->start_at + sctp->num; j++) {
         index1 = SeqMapTableConvert(smtp, i);
         index2 = SeqMapTableConvert(smtp, j);
         new[index1][index2] = old[i][j];
      }
   }

   return new;
}

Int2 Main_new (void)

{
   Boolean query_is_na;
   Boolean db_is_na;
   Boolean believe_query = FALSE;
   EBlastProgramType program_number;
   Int2 status = 0;
   Int4 start=0, end=0;   /* start and end of sequence to be searched as specified by ARG_QUERYLOC */
   FILE *infp=NULL, *outfp=NULL;
   SBlastOptions* options = NULL;
   BlastFormattingInfo* format_info = NULL;
   Int4 ctr = 1;
   Boolean tabular_output = FALSE;
   Blast_SummaryReturn* sum_returns = Blast_SummaryReturnNew();
   Blast_SummaryReturn* full_sum_returns = NULL;
   char* blast_program = myargs[ARG_PROGRAM].strvalue;
   char* dbname = myargs[ARG_DB].strvalue;
   Int4 maxquery = 0; /* maximum number of bases/residues to concatenate per
                         database pass */

   status = SBlastOptionsNew(blast_program, &options, sum_returns);

   if (status) {
       if (sum_returns->error) {
           SBlastMessageErrPost(sum_returns->error);
           sum_returns = Blast_SummaryReturnFree(sum_returns);
       }
       return -1;
   }

   s_FillOptions(options);
   program_number = options->program;

   switch(program_number) {
       case eBlastTypeBlastn:
           maxquery = 40000;
           break;
       case eBlastTypeTblastn:
           maxquery = 20000;
           break;
       case eBlastTypeBlastp:
       case eBlastTypeBlastx:
       case eBlastTypeTblastx:
       default:
           maxquery = 10000;
   }

   BlastGetTypes(myargs[ARG_PROGRAM].strvalue, &query_is_na, &db_is_na);

   if (myargs[ARG_BELIEVEQUERY].intvalue != 0)
        believe_query = TRUE;

   SBlastOptionsSetBelieveQuery(options, believe_query);

   if (myargs[ARG_FORMAT].intvalue == 8 && myargs[ARG_USEMEGABLAST].intvalue)
        tabular_output = TRUE;

   if (!tabular_output) {
       Int2 finfo_status = BlastFormattingInfoNew(myargs[ARG_FORMAT].intvalue, options,
                              blast_program, dbname,
                              myargs[ARG_OUT].strvalue, &format_info);
       if (finfo_status != 0)
       {
           ErrPostEx(SEV_FATAL, 1, 0, "BlastFormattingInfoNew returned non-zero status");
       }

       /* Pass TRUE for the "is megablast" argument. Since megablast is always
          gapped, pass FALSE for the "is ungapped" argument. */
       BlastFormattingInfoSetUpOptions(format_info,
                                       myargs[ARG_DESCRIPTIONS].intvalue,
                                       myargs[ARG_ALIGNMENTS].intvalue,
                                       (Boolean) myargs[ARG_HTML].intvalue,
                                       FALSE,
                                       (Boolean) myargs[ARG_SHOWGIS].intvalue,
                                       believe_query);

      if (myargs[ARG_DB].strvalue) {
         BLAST_PrintOutputHeader(format_info); 
      }
   }
   else
   { /* tabular output requires raw FILE*. */
       if ((outfp = FileOpen(myargs[ARG_OUT].strvalue, "w")) == NULL) {
            ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", 
                myargs[ARG_OUT].strvalue);
            return (1);
       }
       believe_query = TRUE;
       /* FetchEnable/Disable called in blast_format.c for non-tabular output. */
       ReadDBBioseqFetchEnable ("blastall", myargs[ARG_DB].strvalue, db_is_na, TRUE);
   }


   if ((infp = FileOpen(myargs[ARG_QUERY].strvalue, "r")) == NULL) {
      ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open input file %s\n", 
                myargs[ARG_QUERY].strvalue);
      return (1);
   }

   sGetLoc(myargs[ARG_QUERYLOC].strvalue, &start, &end);

   /* Get the query (queries), loop if necessary. */
   while (1) {
      SBlastSeqalignArray* seqalign_arr=NULL;
      BlastTabularFormatData* tf_data = NULL;
      SeqLoc* lcase_mask = NULL;
      SeqLoc* repeat_mask = NULL; /* Repeat mask locations */
      SeqLoc* query_slp = NULL;
      SeqLoc* filter_loc=NULL;	/* All masking locations */
      Int4 num_queries; /* Number of queries read this time. */
      Int4  letters_read;  /* number of letters (bases/residues) read. */

      if ((Boolean)myargs[ARG_LCASE].intvalue) {
         letters_read = BLAST_GetQuerySeqLoc(infp, query_is_na, 
                   myargs[ARG_STRAND].intvalue, maxquery, start, end,
                   &lcase_mask, &query_slp, &ctr, &num_queries, believe_query,
                   myargs[ARG_QGENETIC_CODE].intvalue);
      } else {
         letters_read = BLAST_GetQuerySeqLoc(infp, query_is_na,
                   myargs[ARG_STRAND].intvalue, maxquery, start, end, 
                   NULL, &query_slp, &ctr, &num_queries, believe_query,
                   myargs[ARG_QGENETIC_CODE].intvalue);
      }

      if (letters_read == 0)
          break;

      if (letters_read < 0)
      {
	   ErrPostEx(SEV_FATAL, 1, 0, "BLAST_GetQuerySeqLoc returned an error\n");
           return -1;
      }

      if (believe_query && BlastSeqlocsHaveDuplicateIDs(query_slp)) {
         ErrPostEx(SEV_FATAL, 1, 0, 
                 "Duplicate IDs detected; please ensure that "
                 "all query sequence identifiers are unique");
      }

      if (tabular_output) {
           EBlastTabularFormatOptions tab_option = eBlastTabularDefault;
           if (tabular_output == 2) {
               if (program_number == eBlastTypeBlastn) {
                   tab_option = eBlastTabularAddSequences;
               } else {
                   fprintf(stderr, 
                           "WARNING: Sequences printout in tabular output"
                           " allowed only for blastn\n");
               }
           } 
           
           /* Print the header of tabular output. */
           PrintTabularOutputHeader(dbname, NULL, query_slp, 
                                    blast_program, 0, FALSE, outfp);
           
           tf_data = BlastTabularFormatDataNew(outfp, query_slp, 
                                               tab_option, believe_query);
           tf_data->show_gi = (Boolean) myargs[ARG_SHOWGIS].intvalue;
           tf_data->show_accession = TRUE;
      }

      options->num_cpus = myargs[ARG_THREADS].intvalue;

      /* Find repeat mask, if necessary */
      if ((status = Blast_FindRepeatFilterSeqLoc(query_slp, myargs[ARG_FILTER].strvalue,
                                &repeat_mask, &sum_returns->error)) != 0)
      {
            if (sum_returns && sum_returns->error)
            {
                   ErrSev max_sev = SBlastMessageErrPost(sum_returns->error);
                   if (max_sev >= SEV_ERROR)
                         return status;
            }
      }

      /* Combine repeat mask with lower case mask */
      if (repeat_mask)
          lcase_mask = ValNodeLink(&lcase_mask, repeat_mask);


      if ((status = Blast_DatabaseSearch(query_slp, dbname, lcase_mask, options, 
                                    tf_data, &seqalign_arr, &filter_loc, 
                                    sum_returns)) != 0)
      {
            /* Jump out if fatal error or unknown reason for exit. */
            if (sum_returns && sum_returns->error)
            {
                ErrSev max_severity = SBlastMessageErrPost(sum_returns->error);
                if (max_severity >= SEV_ERROR)
                   return status;
            }
            else if (!sum_returns || !sum_returns->error)
            {
                   ErrPostEx(SEV_ERROR, 1, 0, "Non-zero return from Blast_DatabaseSearch\n");
                   return status;
            }
      }

       /* Deallocate the formatting thread data structure. */
       if (tabular_output)
           BlastTabularFormatDataFree(tf_data);

       /* Free the lower case mask in SeqLoc form. */
       lcase_mask = Blast_ValNodeMaskListFree(lcase_mask);

      /* If masking was done for lookup table only, free the masking locations,
          because they will not be used for formatting. */
       if (SBlastOptionsGetMaskAtHash(options))
           filter_loc = Blast_ValNodeMaskListFree(filter_loc);

       /* Post warning or error messages, no matter what the search status was. */
       SBlastMessageErrPost(sum_returns->error);

       if (!status && !tabular_output) {
/*   FIXME:
           Int4** ascii_matrix = BlastMatrixConvert(sbp->matrix);
*/
           if (myargs[ARG_ASNOUT].strvalue) {
                   /* This just prints out the ASN.1 to a secondary file. */
                   BlastFormattingInfo* asn_format_info = NULL;
                   BlastFormattingInfoNew(eAlignViewAsnText, options,
                              blast_program, dbname,
                              myargs[ARG_ASNOUT].strvalue, &asn_format_info);

                   BlastFormattingInfoSetUpOptions(asn_format_info,
                                       myargs[ARG_DESCRIPTIONS].intvalue,
                                       myargs[ARG_ALIGNMENTS].intvalue,
                                       (Boolean) myargs[ARG_HTML].intvalue,
                                       FALSE,
                                       (Boolean) myargs[ARG_SHOWGIS].intvalue,
                                       believe_query);
                   status = 
                       BLAST_FormatResults(seqalign_arr, num_queries, query_slp, 
                                   NULL, asn_format_info, sum_returns);
                   asn_format_info = BlastFormattingInfoFree(asn_format_info);
           }
           
           /* Format the results */
           status = 
               BLAST_FormatResults(seqalign_arr, num_queries, query_slp, 
                                   filter_loc, format_info, sum_returns);
       }

       seqalign_arr = SBlastSeqalignArrayFree(seqalign_arr);
       /* Update the cumulative summary returns structure and clean the returns
          substructures for the current search iteration. */
       Blast_SummaryReturnUpdate(sum_returns, &full_sum_returns);
       Blast_SummaryReturnClean(sum_returns);
       filter_loc = Blast_ValNodeMaskListFree(filter_loc);
       query_slp = SeqLocSetFree(query_slp);
   } /* End loop on sets of queries */

   Blast_PrintOutputFooter(format_info, full_sum_returns);

   options = SBlastOptionsFree(options);
   sum_returns = Blast_SummaryReturnFree(sum_returns);
   full_sum_returns = Blast_SummaryReturnFree(full_sum_returns);

   if (!tabular_output)
      format_info = BlastFormattingInfoFree(format_info);
   else
   {
      FileClose(outfp);
      /* FetchEnable/Disable called in blast_format.c for non-tabular output. */
      ReadDBBioseqFetchDisable();
   }

   if (infp)
      FileClose(infp);
   
   return status;
}

#endif /* BLASTALL_TOOLS_ONLY */


/* Amount to relax the evalue threshold for preliminary alignments
 * when compositionally adjusted score matrices are used. */
#define EVALUE_EXPAND 1


Int2 Main_old (void)
 
{
    AsnIoPtr aip, xml_aip;
    BioseqPtr fake_bsp = NULL, query_bsp, bsp;
    BioSourcePtr source;
    BLAST_MatrixPtr matrix;
    Int4Ptr PNTR txmatrix;
    BLAST_OptionsBlkPtr options;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
    BlastPruneSapStructPtr prune;
    Boolean db_is_na, query_is_na, show_gi, believe_query=FALSE;
    Boolean html = FALSE;
    CharPtr params_buffer=NULL;
    Int4 number_of_descriptions, number_of_alignments;
    SeqAlignPtr  seqalign;
    SeqAnnotPtr seqannot = NULL;
    SeqEntryPtr sep;
    TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
    Uint1 align_type, align_view, err_ticket;
    Uint4 align_options, print_options;
    ValNodePtr mask_loc, mask_loc_start = NULL, vnp, next_mask_loc = NULL;
    ValNodePtr other_returns, error_returns;
    CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
    FILE *infp, *outfp;
    /* Mega BLAST related variables */
    SeqAlignPtr sap, next_seqalign, PNTR seqalignp;
    Int4 num_bsps, index;
    SeqLocPtr last_mask, mask_slp, slp = NULL, tmp_slp;
    Int2 ctr = 1;
    Char prefix[2];
    Boolean done = TRUE;
    int (LIBCALLBACK *handle_results)(VoidPtr srch);       
    Int4 from = 0, to = -1;
    Uint4 num_queries;		/*--KM for concatenated queries in blastn, tblastn */
    Uint4 num_iters;
    Uint4 sap_iter;
    SeqAlignPtr curr_seqalign;
    SeqAlignPtrArray sap_array;		/*--KM for separating seqaligns to test concat printing, temporary?*/
    SeqAnnotPtr curr_seqannot;
    SeqAnnotPtrArray seq_annot_arr;
    Uint4 bsp_iter;
    BspArray fake_bsp_arr;	/*--KM the array of fake_bsps for indiv. queries */ 
    SeqLocPtr PNTR lcase_mask_arr = NULL;	/* AM: information about lower case masked parts of queries */
    Boolean concat_done, nuc_concat;
    QueriesPtr mult_queries = NULL;	/*--KM, AM: stores information related to 
                                                    query multipolexing, to put in search */
    BioseqPtr curr_bsp;

    /* AM: Support for query multiplexing. */
    Uint4 num_spacers;
    ValNodePtr orig_mask_loc = NULL;

#ifdef BLAST_CS_API
    BlastNet3Hptr    bl3hp;
    Boolean status;
#endif
    
    blast_program = myargs[ARG_PROGRAM].strvalue;

#ifdef BLAST_CS_API
    /* For RPS Blast - anything not "blastp" - is "tblastn" */    
    if(myargs[ARG_RPSBLAST].intvalue) {
        if(StringICmp(blast_program, "blastp")) {
            StringCpy(blast_program, "blastx");
        }
    }
#endif

    blast_database = myargs[ARG_DB].strvalue;
    blast_inputfile = myargs[ARG_QUERY].strvalue;
    blast_outputfile = myargs[ARG_OUT].strvalue;

    if (myargs[ARG_HTML].intvalue)
        html = TRUE;
    
    if ((infp = FileOpen(blast_inputfile, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open input file %s\n", blast_inputfile);
        return (1);
    }

    align_view = (Int1) myargs[ARG_FORMAT].intvalue;
    outfp = NULL;
    if (align_view != 7 && align_view != 10 && align_view != 11 && blast_outputfile != NULL) {
        if ((outfp = FileOpen(blast_outputfile, "w")) == NULL) {
            ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", blast_outputfile);
            return (1);
        }
    }
    
    if (StringCmp("filter", blast_program) == 0) {
        BlastGetMaskingLoc(infp, outfp, myargs[ARG_FILTER].strvalue);
        FileClose(outfp);
        FileClose(infp);	
        return 0;
    }
    
    align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);

    if(align_view < 7) {
        if (StringICmp("blastx", blast_program) == 0) {
            if (align_view != 0) {
                ErrPostEx(SEV_FATAL, 1, 0, "This option is not available with blastx");
                return 1;
            }
        } else if (StringICmp("tblastx", blast_program) == 0) {
            if (align_view != 0) {
                ErrPostEx(SEV_FATAL, 1, 0, "This option is not available with tblastx");
                return 1;
            }
        }
    }
    
    believe_query = FALSE;
    if (myargs[ARG_BELIEVEQUERY].intvalue != 0)
        believe_query = TRUE;
    
    if (believe_query == FALSE && (myargs[ARG_ASNOUT].strvalue || align_view == 10 || align_view ==11)) {
        ErrPostEx(SEV_FATAL, 1, 0, "-J option must be TRUE to produce a SeqAlign file");
    }
    
    options = BLASTOptionNewEx(blast_program, (Boolean) myargs[ARG_GAPPED].intvalue, (Boolean) myargs[ARG_USEMEGABLAST].intvalue);
    if (options == NULL)
        return 3;

#ifdef BLAST_CS_API
    if(myargs[ARG_RPSBLAST].intvalue) 
        options->is_rps_blast = TRUE;
#endif
    
    handle_results = NULL;

    BLASTOptionSetGapParams(options, myargs[ARG_MATRIX].strvalue, 0, 0); 
    options->kappa_expect_value =
        options->expect_value  = (Nlm_FloatHi) myargs[ARG_EVALUE].floatvalue;
    number_of_descriptions = myargs[ARG_DESCRIPTIONS].intvalue;	
    number_of_alignments = myargs[ARG_ALIGNMENTS].intvalue;	
    options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);

    if (StringICmp("blastn", blast_program) == 0) {
        options->penalty = myargs[ARG_MISMATCH].intvalue;
        options->reward = myargs[ARG_MATCH].intvalue;
        if (options->reward > 1) {
           /* Scale the default values for gap costs; will be overridden
              later, if command line values are non-zero */
           options->gap_open *= options->reward;
           options->gap_extend *= options->reward;
        }
    } else {
        if (myargs[ARG_THRESHOLD].intvalue != 0) {
            options->threshold_second = myargs[ARG_THRESHOLD].intvalue;
        }
    }
    
    if (myargs[ARG_GAPOPEN].intvalue >= 0)
        options->gap_open = myargs[ARG_GAPOPEN].intvalue;
    if (myargs[ARG_GAPEXT].intvalue >= 0)
        options->gap_extend = myargs[ARG_GAPEXT].intvalue;
    if (myargs[ARG_XDROP].intvalue != 0)
        options->gap_x_dropoff = myargs[ARG_XDROP].intvalue;

    /* use one-hit if specified or it's a blastn search */
    if ( (myargs[ARG_MULTIPLEHITS].intvalue == 1) || (StringICmp("blastn", blast_program) == 0 ) )
      {
        options->two_pass_method  = FALSE;
        options->multiple_hits_only  = FALSE;
      }
    /* otherwise, use two-hit */
    else
      { 
        /* all other inputs, including the default 0 use 2-hit method */
        options->two_pass_method  = FALSE;
        options->multiple_hits_only  = TRUE;
      }
    
    if(myargs[ARG_XDROP_FINAL].intvalue != 0) 
        options->gap_x_dropoff_final = myargs[ARG_XDROP_FINAL].intvalue;

    if (StringICmp(myargs[ARG_FILTER].strvalue, "T") == 0) {
        if (StringICmp("blastn", blast_program) == 0)
            options->filter_string = StringSave("D");
        else
            options->filter_string = StringSave("S");
    } else {
        options->filter_string = StringSave(myargs[ARG_FILTER].strvalue);
    }
    
    show_gi = (Boolean) myargs[ARG_SHOWGIS].intvalue;

    options->genetic_code = myargs[ARG_QGENETIC_CODE].intvalue;
    options->db_genetic_code = myargs[ARG_DBGENCODE].intvalue;
    options->number_of_cpus = myargs[ARG_THREADS].intvalue;
    if (myargs[ARG_WORDSIZE].intvalue != 0) {
        options->wordsize = myargs[ARG_WORDSIZE].intvalue;
    }
    
    if (options->is_megablast_search) {
       options->cutoff_s2 = options->wordsize*options->reward;
    }

    options->db_length = (Int8) myargs[ARG_DBSIZE].floatvalue;
    
    options->hsp_range_max  = myargs[ARG_BESTHITS].intvalue;
    if (options->hsp_range_max != 0)
        options->perform_culling = TRUE;
    if (myargs[ARG_SEARCHSP].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[ARG_SEARCHSP].floatvalue;

    if (0 != StringICmp("tblastn", blast_program) ||
       !options->gapped_calculation) {
        /* Set some gapped tblastn-specific options to the correct
         * defaults for non-tblastn or non-gapped modes of operation.
         */
        options->tweak_parameters = eNoCompositionBasedStats;
        options->smith_waterman = 0;
        
        switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
        case '0': case 'D': case 'd': case 'F': case 'f':
            options->tweak_parameters = eNoCompositionBasedStats;
            break;
        default:
            ErrPostEx(SEV_FATAL, 1, 0,
               "Invalid option -C: only gapped tblastn may use"
               " composition based statistics.");
            break;
        }
        if(myargs[ARG_SMITH_WATERMAN].intvalue) {
            ErrPostEx(SEV_FATAL, 1, 0,
               "Invalid option -s: Smith-Waterman alignments are only "
               "available for gapped tblastn.");
        }
    } else {
        /* Set options specific to gapped tblastn */
        switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
        case 'D': case 'd':
        case '0': case 'F': case 'f':
            options->tweak_parameters = eNoCompositionBasedStats;
            break;
        case '1': case 'T': case 't':
            options->tweak_parameters = eCompositionBasedStats;
            break;
        case '2':
            ErrPostEx(SEV_WARNING, 1, 0, "the -C 2 argument "
                      "is currently experimental\n");
            options->tweak_parameters = eCompositionMatrixAdjust;
            break;
        case '3':
            ErrPostEx(SEV_WARNING, 1, 0, "the -C 3 argument "
                      "is currently experimental\n");
            options->tweak_parameters = eCompoForceFullMatrixAdjust;
        break;
        default:
            ErrPostEx(SEV_FATAL, 1, 0, "invalid argument for composition-"
                      "based statistics; see -C options\n");
            break;
        }
        options->smith_waterman =
            (Boolean) myargs[ARG_SMITH_WATERMAN].intvalue;
    }
    if (options->tweak_parameters > 1) {
        /* Compositionally adjusted score matrices are being used, and
         * these can improve evalue, so relax the evalue cutoff for
         * the preliminary alignments.  (Note that traditional
         * composition based statistics can only make evalues larger.)
         */
        options->expect_value *= EVALUE_EXPAND;
    }

    options->strand_option = myargs[ARG_STRAND].intvalue;

    if(myargs[ARG_XDROP_UNGAPPED].floatvalue != 0.0) {
        options->dropoff_2nd_pass  = myargs[ARG_XDROP_UNGAPPED].floatvalue;
        if(options->dropoff_1st_pass > options->dropoff_2nd_pass)
            options->dropoff_1st_pass = options->dropoff_2nd_pass;
    }

    if (myargs[ARG_WINDOW].intvalue != 0)
        options->window_size = myargs[ARG_WINDOW].intvalue;

    print_options = 0;
    align_options = 0;
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;
    if (StringICmp("blastx", blast_program) == 0) {
        align_options += TXALIGN_BLASTX_SPECIAL;
    }
    if (show_gi) {
        align_options += TXALIGN_SHOW_GI;
        print_options += TXALIGN_SHOW_GI;
    }
    if (myargs[ARG_GAPPED].intvalue == 0 || StringICmp("tblastx", blast_program) == 0)
        print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
    if (align_view) {
        align_options += TXALIGN_MASTER;
        if (align_view == 1 || align_view == 3)
            align_options += TXALIGN_MISMATCH;
        if (align_view == 3 || align_view == 4 || align_view == 6)
            align_options += TXALIGN_FLAT_INS;
        if (align_view == 5 || align_view == 6)
            align_options += TXALIGN_BLUNT_END;
    } else {
        align_options += TXALIGN_MATRIX_VAL;
        align_options += TXALIGN_SHOW_QS;
    }
    
    if (html) {
        align_options += TXALIGN_HTML;
        print_options += TXALIGN_HTML;
    }

#ifdef BLAST_CS_API
    if(myargs[ARG_ENTREZQ].strvalue)
        options->entrez_query = StringSave(myargs[ARG_ENTREZQ].strvalue);
#else    
    if (myargs[ARG_GILIST].strvalue) {
        options->gifile = StringSave(myargs[ARG_GILIST].strvalue);
    }
#endif
    
    /* 
       Out-of-frame option is valid only for blastx, tblastn and 
       psitblastnsearches
    */

    if(myargs[ARG_FRAMESHIFT].intvalue > 0) {
        if (!StringICmp("blastx", blast_program) || 
            !StringICmp("tblastn", blast_program)||
	    !StringICmp("psitblastn", blast_program)) {
           if (!StringICmp("blastx", blast_program)) {
              options->is_ooframe = TRUE;
              options->shift_pen = myargs[ARG_FRAMESHIFT].intvalue;
           }
        }
    }
        
    /* Input longest intron length is in nucleotide scale; in the lower level
       code it will be used in protein scale */
    options->longest_intron = myargs[ARG_INTRON].intvalue;

    aip = NULL;
    if (myargs[ARG_ASNOUT].strvalue != NULL) {
        if ((aip = AsnIoOpen (myargs[ARG_ASNOUT].strvalue,"w")) == NULL) {
                ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", myargs[ARG_ASNOUT].strvalue);
                return 1;
        }
    }
    else if (align_view == 10 || align_view == 11) 
    {
        const char* mode = (align_view == 10) ? "w" : "wb";
        if ((aip = AsnIoOpen (blast_outputfile, (char*) mode)) == NULL) {
                ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", myargs[ARG_ASNOUT].strvalue);
                return 1;
        }
    }

    if(align_view < 7) {
       if (html) {
          fprintf(outfp, "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
          fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                  "VLINK=\"#660099\" ALINK=\"#660099\">\n");
          fprintf(outfp, "<PRE>\n");
       }
    } else if (align_view == 7 ) {
        xml_aip = AsnIoOpen(blast_outputfile, "wx");
    }

#ifndef BLAST_CS_API
    if(align_view >= 7 && myargs[ARG_NUMQUERIES].intvalue > 1)
    {
      ErrPostEx(SEV_FATAL, 1, 0, 
                 "blast: Query concatenation is currently not supported with -m > 7");
      return 1;
    }
#endif


                  /* Futamura: Setting up the psitblastn options */
#ifndef BLAST_CS_API
    if (NULL != myargs[ARG_PSITCHKPNT].strvalue) {
          options->recoverCheckpoint = TRUE;
          options->freqCheckpoint = TRUE;
    }
    options->CheckpointFileName=myargs[ARG_PSITCHKPNT].strvalue;
#endif

#ifdef BLAST_CS_API
    if (align_view < 7)
    	bl3hp = BNETInitializeBlast(blast_database, blast_program, outfp, 
                                db_is_na, options->is_rps_blast, html, TRUE);
    else
    	bl3hp = BNETInitializeBlast(blast_database, blast_program, outfp, 
                                db_is_na, options->is_rps_blast, html, FALSE);
#endif

    /*--KM get number of queries for concatenated blastn/tblastn queries */

#ifndef BLAST_CS_API
    options->NumQueries=myargs[ARG_NUMQUERIES].intvalue;  
#endif

    num_queries = options->NumQueries;
    if (num_queries>0 && 
	!( (StringICmp("blastn",  blast_program) == 0) || 
	   (StringICmp("tblastn", blast_program) == 0)   ) ) {

	ErrPostEx(SEV_FATAL, 1, 0, "blast: Can't concat with program %s\n", myargs[ARG_PROGRAM].strvalue);
       return 1;
    }
    
    /* AM: Query concatenation is not consistent with ungapped search */
    if( num_queries > 0 && !myargs[ARG_GAPPED].intvalue )
    {
      ErrPostEx(SEV_FATAL, 1, 0, 
                 "blast: Query concatenation is inconsistent with ungapped search\n" );
      return 1;
    }
    if( !myargs[ARG_GAPPED].intvalue &&
        0 == StringCmp("psitblastn", blast_program ) ) {
      ErrPostEx(SEV_FATAL, 1, 0,"blast: Ungapped alignment is not appropriate "
                "for PSI-tBLASTn.\n" );
    }

    /* --KM set bool value if DNA and concat needed, need for Fasta->seq functions */
    if (num_queries>0 && query_is_na == TRUE) {
        nuc_concat = TRUE;
    } else {
        nuc_concat = FALSE;
    }
 
    /* --- Main loop over all FASTA entries in the input file ---- */

    concat_done = FALSE;	/*--KM */

    sGetLoc(myargs[ARG_QUERYLOC].strvalue, &from, &to);

    while (TRUE) {
       if (options->is_megablast_search) {
          StrCpy(prefix, "");
          slp = NULL;
	  num_bsps = 0;
          done = TRUE;
	  SeqMgrHoldIndexing(TRUE);
	  mask_slp = last_mask = NULL;
	  while ((sep=FastaToSeqEntryForDb(infp, query_is_na, NULL,
					   believe_query, prefix, &ctr, 
					   &mask_slp)) != NULL) {
	     if ((Boolean)myargs[ARG_LCASE].intvalue) {
                if (mask_slp) {
                   if (!last_mask)
                      options->query_lcase_mask = last_mask = mask_slp;
                   else {
                      last_mask->next = mask_slp;
                      last_mask = last_mask->next;
                   }
                   mask_slp = NULL;
                }
             } else {
                mask_slp = SeqLocSetFree(mask_slp);
             }
	     query_bsp = NULL;
	     if (query_is_na) 
		SeqEntryExplore(sep, &query_bsp, FindNuc);
	     else
		SeqEntryExplore(sep, &query_bsp, FindProt);
	     
	     if (query_bsp == NULL) {
		ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
		return 2;
	     }
	     
             /* Only for the first query */
             if (num_bsps == 0) {
                 to = MIN(to, query_bsp->length - 1);
                 
                 /* -1 means end of sequence */
                 if (to < 0)
                     to = query_bsp->length - 1;
                 if (from >= query_bsp->length || to < 0) {
                     ErrPostEx(SEV_FATAL, 1, 0, 
                               "Location outside of the query sequence range\n");
                     return 3;
                 }
                 slp = SeqLocIntNew(from, to, options->strand_option, 
                                    SeqIdFindBest(query_bsp->id, SEQID_GI));
             } else 
                 ValNodeAddPointer(&slp, SEQLOC_WHOLE,
                                   SeqIdDup(SeqIdFindBest(query_bsp->id,
                                                          SEQID_GI)));
	     num_bsps++;
	     if (num_bsps >= MAX_NUM_QUERIES) {
                done = FALSE;
		break;
	     }
	     /*sep = MemFree(sep);*/ /* Do not free the underlying Bioseq */
	  }
	  SeqMgrHoldIndexing(FALSE);
          if (num_bsps == 0) 
             break;
       } else {
          /* not megablast */

          /*--KM make array of fake_bsp's if concat. query */
          if (concat_done)
             break;
          if (num_queries > 0)  {
             fake_bsp_arr = (BspArray) MemNew(sizeof(BioseqPtr)*num_queries); 

	     if( myargs[ARG_LCASE].intvalue )
	       lcase_mask_arr = (SeqLocPtr PNTR)MemNew( sizeof( SeqLocPtr )*num_queries );
          }
          num_iters = (num_queries>0) ? num_queries : 1; 
          for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {

             if(myargs[ARG_LCASE].intvalue) {
	        /* AM: query multiplexing */
		if( !num_queries )
                  sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, believe_query, NULL, NULL, &options->query_lcase_mask);
                else
		  sep = FastaToSeqEntryInternalEx( infp, FASTA_FILE_IO, NULL, query_is_na, NULL, believe_query,
		                                   NULL, NULL, NULL, lcase_mask_arr + bsp_iter );
                
             } else {
                sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query);
             }
          
             /* if concat and num_queries has not been reached and sep is NULL, crap out */
             if (sep == NULL && bsp_iter < num_queries) {   /* implies num_queries>0 */
                ErrPostEx(SEV_FATAL, 1, 0, "blast: Only %d queries found!\n", bsp_iter); 
                return (1);
             }
               
             if(sep == NULL)
                break;  /* no more queries, can go to finish with next break */
          
             query_bsp = NULL;
             if (query_is_na) {
                SeqEntryExplore(sep, &query_bsp, FindNuc);
             } else {
                SeqEntryExplore(sep, &query_bsp, FindProt);
             }
          
             if (query_bsp == NULL) {
                ErrPostEx(SEV_FATAL, 1, 0, "Unable to obtain bioseq\n");
                return 2;
             }

             if (num_queries>0) {
                *(fake_bsp_arr + bsp_iter) = query_bsp;
             }
          }
          if ( (sep == NULL && num_queries ==0) || (num_queries>0 && concat_done) )
             break;  /* go to finish */

          /* --KM */
          
          if (num_queries>0) {
             concat_done = TRUE;   /* --KM to prevent futher looping */

             /* AM: Determine the number of query separators. */
	     num_spacers = GetNumSpacers( options, believe_query, fake_bsp_arr ); 

	     if( num_spacers%2 ) ++num_spacers;

             /* --KM make the concatenated fake_bsp */
	     /* AM: Added num_spacers. */
	     if( query_is_na )
               fake_bsp = (BioseqPtr) 
                          BlastMakeFakeBspConcat(fake_bsp_arr, num_queries, query_is_na, num_spacers); 
             else
               fake_bsp = (BioseqPtr) 
                          BlastMakeFakeBspConcat(fake_bsp_arr, num_queries, query_is_na, num_spacers); 
             
             /* construct the MultQueries struct here*/
             mult_queries = (QueriesPtr) BlastMakeMultQueries(fake_bsp_arr, num_queries, query_is_na, num_spacers,
	                                                      lcase_mask_arr);
          } else {
             if(believe_query)
                fake_bsp = query_bsp;
             else 
                fake_bsp = BlastMakeFakeBioseq(query_bsp, NULL);
          }

	  err_ticket = BlastSetUserErrorString(NULL, query_bsp->id, believe_query);
        
          /* If fake_bsp created mask should be updated to use it's id */
	  /* AM: query multiplexing */
	  if( !mult_queries )
            BLASTUpdateSeqIdInSeqInt(options->query_lcase_mask, fake_bsp->id);
          else for( bsp_iter = 0; bsp_iter < num_iters; ++bsp_iter )
	         if( mult_queries->LCaseMasks )
	           BLASTUpdateSeqIdInSeqInt( mult_queries->LCaseMasks[bsp_iter],
		                             mult_queries->FakeBsps[bsp_iter]->id );
        
          source = BioSourceNew();
          source->org = OrgRefNew();
          source->org->orgname = OrgNameNew();
          source->org->orgname->gcode = options->genetic_code;
          ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);

       /* free sep later when done. --KM remember to free all if array*/
       }

       global_fp = outfp;
          
       if(align_view < 7) {
#ifndef BLAST_CS_API
           init_buff_ex(90);
           BlastPrintVersionInfo(blast_program, html, outfp);
           fprintf(outfp, "\n");
           BlastPrintReference(html, 90, outfp);
           fprintf(outfp, "\n");
#else
           fprintf(outfp, "\n");
#endif            
           if (!options->is_megablast_search) {
              /* KM added loop here for concat case */
              num_iters = (num_queries>0) ? num_queries : 1;
              for (bsp_iter=0; bsp_iter<num_iters; bsp_iter++) {
                 curr_bsp = (num_queries>0) ? *(fake_bsp_arr + bsp_iter) : query_bsp; 
                 AcknowledgeBlastQuery(curr_bsp, 70, outfp, believe_query, html);
              }
           }

            /* Here we first check, that database do no exists */

#ifndef BLAST_CS_API
           if(!PrintDbInformation(blast_database, !db_is_na, 70, outfp, html))
                return 1;
#else

            {{
                BlastDbinfoPtr dbinfo;
                static Boolean not_first_time;

                /* For CS version we will print database info ones to 
                   decrease network traffic */                

                if(!not_first_time) {
                    dbinfo = BlastRequestDbInfo(bl3hp, blast_database, !db_is_na);
                    if (dbinfo)
                        PrintDbInformationBasic(blast_database, !db_is_na, 70, dbinfo->definition, dbinfo->number_seqs, dbinfo->total_length, outfp, html);
                    dbinfo = BlastDbinfoFree(dbinfo);
                    not_first_time = TRUE;
                }
            }}
#endif    /* BLAST_CS_API */        
            free_buff();
		if (options->is_ooframe)
        		ErrPostEx(SEV_WARNING, 0, 0, "Out-of-frame option selected, Expect values are only approximate and calculated not assuming out-of-frame alignments");
        }
#ifdef OS_UNIX
        if(align_view < 7) { /*--KM why not fold into previous if statement? */
#ifdef BLAST_CS_API
            fprintf(global_fp, "%s", "Searching... please wait.. ");
#else
            fprintf(global_fp, "%s", "Searching");
#endif
        }
#endif
        other_returns = NULL;
        error_returns = NULL;

        if (options->is_megablast_search) {
#ifdef BLAST_CS_API
           seqalign = MegaBlastSeqLocNetCore(bl3hp, slp, blast_program, 
                                   blast_database, options, 
                                   &other_returns, &error_returns,
                                   align_view < 7 ? tick_callback : NULL,
                                   &status);
#else
           seqalignp = BioseqMegaBlastEngineByLoc(slp, blast_program,
			           blast_database, options, &other_returns, 
                                   &error_returns, 
                                   align_view < 7 ? tick_callback : NULL,
                                   NULL, NULL, 0, handle_results);
           seqalign = NULL;
	   for (index=0; index<num_bsps; index++) { 
	      if (seqalignp && seqalignp[index]) {
                 if (seqalign == NULL) 
                    sap = seqalign = seqalignp[index];
                 else
                    sap->next = seqalignp[index];
                 while (sap->next != NULL)
                    sap = sap->next;
              }
           }
           seqalignp = MemFree(seqalignp);
#endif
        } else if (!myargs[ARG_QUERYLOC].strvalue) {       
#ifdef BLAST_CS_API
           seqalign = BlastBioseqNetCore(bl3hp, fake_bsp, blast_program, 
                                      blast_database, options,
				      &other_returns, &error_returns,
                                      align_view < 7 ? tick_callback : NULL,
				      NULL, &status);
#else
           /* KM added mult_queries param */
           seqalign = BioseqBlastEngineWithCallbackMult(fake_bsp, blast_program, blast_database, options, &other_returns, &error_returns, align_view < 7 ? tick_callback : NULL, handle_results, mult_queries);
#endif
        } else { /* Location on query provided */
           to = MIN(to, fake_bsp->length - 1);
           
           /* -1 means end of sequence */
           if (to < 0)
              to = fake_bsp->length - 1;
           if (from >= fake_bsp->length || to < 0) {
              ErrPostEx(SEV_FATAL, 1, 0, 
                        "Location outside of the query sequence range\n");
              return 3;
           }
           slp = SeqLocIntNew(from, to, options->strand_option, 
                              fake_bsp->id);
           
#ifdef BLAST_CS_API
           seqalign = BlastSeqLocNetCore(bl3hp, slp, blast_program, 
                                         blast_database, options,
                                         &other_returns, &error_returns,
                                         align_view < 7 ? tick_callback : NULL,
                                         NULL, &status);
#else
           seqalign = BioseqBlastEngineByLocWithCallbackMult(slp, blast_program, blast_database, options, &other_returns, &error_returns, align_view < 7 ? tick_callback : NULL, NULL, NULL, 0, handle_results, mult_queries);
#endif
           
        }
#if 0
        seqalign = BLASTFilterOverlapRegions(seqalign, 0, !db_is_na, 
                                             options->is_ooframe, FALSE);
#endif
        
        BlastErrorPrint(error_returns);

        dbinfo = NULL;
        ka_params = NULL;
        ka_params_gap = NULL;
        params_buffer = NULL;
        mask_loc = NULL;
        matrix = NULL;
        txmatrix = NULL;
        for (vnp=other_returns; vnp; vnp = vnp->next) {
            switch (vnp->choice) {
            case TXDBINFO:
                dbinfo = vnp->data.ptrvalue;
                break;
            case TXKABLK_NOGAP:
                ka_params = vnp->data.ptrvalue;
                break;
            case TXKABLK_GAP:
                ka_params_gap = vnp->data.ptrvalue;
                break;
            case TXPARAMETERS:
                params_buffer = vnp->data.ptrvalue;
                break;
            case TXMATRIX:
                matrix = vnp->data.ptrvalue;
                if (matrix)
                   txmatrix = BlastMatrixToTxMatrix(matrix);
                break;
            case SEQLOC_MASKING_NOTSET:
            case SEQLOC_MASKING_PLUS1:
            case SEQLOC_MASKING_PLUS2:
            case SEQLOC_MASKING_PLUS3:
            case SEQLOC_MASKING_MINUS1:
            case SEQLOC_MASKING_MINUS2:
            case SEQLOC_MASKING_MINUS3:
                ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
                break;
            default:
                break;
            }
        }	

#ifdef OS_UNIX
        fflush(global_fp);
#endif
        
#ifdef OS_UNIX
        if(align_view < 7) {
            fprintf(global_fp, "%s", "done");
        }
#endif
        
#ifndef BLAST_CS_API
    ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
#endif
        ReadDBBioseqSetDbGeneticCode(options->db_genetic_code);

        tmp_slp = slp;
        if (slp)
           query_bsp = NULL;

        if (getenv("POST_BLAST_CLUSTER_HITS") != NULL)
           BlastClusterHitsFromSeqAlign(seqalign, blast_program, blast_database, 
                                        options, 0.9, 1.6, 0.5, TRUE);

        if (mask_loc) {
           mask_loc_start = mask_loc;
        }
	else
	{	/* Could have become non-NUll for last query. */
           mask_loc_start = NULL;
	}
        /* Print header in any case */
        if (align_view == 9) {
           PrintTabularOutputHeader(blast_database, query_bsp, slp, 
              blast_program, 0, believe_query, global_fp);
        }

        if (seqalign) {
	   if (num_queries > 0) { /* AM: Support for query multiplexing. */
	      sap_array = mult_queries->sap_array_data->sap_array;
	   }   
	
	   if (align_view == 8 || align_view == 9) {
/* --KM need to put a loop around this. seqaligns already broken up
   note the method for looping if num_aligns > 0 - reuse this method everywhere */
	      num_iters = (num_queries>0) ? num_queries : 1;
	      for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
	         curr_seqalign = (num_queries>0) ? *(sap_array + sap_iter) : seqalign;
	         BlastPrintTabularResults(curr_seqalign, query_bsp, slp, 
               number_of_alignments, blast_program, 
               !options->gapped_calculation, options->is_ooframe,
               believe_query, 0, 0, global_fp, NULL, (align_view == 9));

	         SeqAlignSetFree(curr_seqalign);
	      }
           } else {
           while (seqalign) {
	    		
              if (!options->is_megablast_search){
                 next_seqalign = NULL;
	      } else {
                 SeqIdPtr sip, next_sip = NULL;
                 
                 sap = seqalign;
                 sip = TxGetQueryIdFromSeqAlign(seqalign);
                 while (sap != NULL) { 
                    if (sap->next != NULL) {
                       next_sip = TxGetQueryIdFromSeqAlign(sap->next);

                       if (SeqIdComp(sip, next_sip) != SIC_YES) {
                          next_seqalign = sap->next;
                          sap->next = NULL;
                       }
                    } else
                       next_seqalign = NULL;
                    sap = sap->next;
                 }
                 
                 while (tmp_slp && SeqIdComp(sip, SeqLocId(tmp_slp)) != SIC_YES)
                    tmp_slp = tmp_slp->next;
                 if (tmp_slp == NULL) /* Should never happen */
                    break;
                 /* Separate the mask locations list for this query */
                 if (!mask_loc && next_mask_loc) {
                    mask_loc = next_mask_loc;
                    next_mask_loc = NULL;
                 }
                 if (mask_loc) {
                    if (next_mask_loc) {
                       mask_loc->next = next_mask_loc;
                       mask_loc = next_mask_loc;
                    }
                    mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
                    next_mask_loc = mask_loc;
                    while (SeqIdComp(SeqLocId(mask_slp), sip) != SIC_YES) {
                       mask_loc = mask_loc->next;
                       if (!mask_loc)
                          break;
                       mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
                    }
                    if (mask_loc) {
                       next_mask_loc = mask_loc->next;
                       mask_loc->next = NULL;
                    }
                 }
                 if (align_view < 7) {
                     bsp = BioseqLockById(SeqLocId(tmp_slp));
                     init_buff_ex(85);
                     fprintf(outfp, "\n");
                     AcknowledgeBlastQuery(bsp, 70, outfp, believe_query, 
                                           html);
                     free_buff();
                     BioseqUnlock(bsp);
                 }
              }
              if((align_view == 7) && !options->is_ooframe) {
                 if (options->is_megablast_search) {
                    bsp = BioseqLockById(SeqLocId(tmp_slp));
                    BXMLPrintOutput(xml_aip, seqalign, 
                                    options, blast_program, blast_database, 
                                    bsp, other_returns, 0, NULL, mask_loc);
                    BioseqUnlock(bsp);
                    AsnIoReset(xml_aip);
                    SeqAlignSetFree(seqalign);
                 } else {
	            num_iters = (num_queries>0) ? num_queries : 1;
                    for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                       curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
                       BXMLPrintOutput(xml_aip, curr_seqalign, 
                                    options, blast_program, blast_database, 
                                    fake_bsp, other_returns, 0, NULL, mask_loc);
                       AsnIoReset(xml_aip);
                       SeqAlignSetFree(curr_seqalign);
                    } /* for loop over sap-array (concat) */
                 } /* not MBlast case */
              } else {
	         /* create the array of SeqAnnotPtrs, if necessary */

	         num_iters = (num_queries > 0) ? num_queries : 1; 
                 for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                    curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
                    if ( (num_queries > 0) && (sap_iter == 0) ) {
                       seq_annot_arr = (SeqAnnotPtrArray) MemNew(sizeof(SeqAnnotPtr)*num_queries);
                    }
                    seqannot = SeqAnnotNew();
                    seqannot->type = 2;
                    AddAlignInfoToSeqAnnot(seqannot, align_type);
                    seqannot->data = curr_seqalign;
                    if (aip) {
                       SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
                       AsnIoReset(aip);
                    }
                    if (num_queries > 0) {
                       *(seq_annot_arr + sap_iter) = seqannot;
                    }
                 } /* make seqannots over the sap_iters from concat, or the single seqalign */
                    
                 if (outfp) { /* Uncacheing causes problems with ordinal nos. vs. gi's. */
                    ObjMgrSetHold();
                    /* print deflines */
                    for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
                       curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;

                       init_buff_ex(85);

                       PrintDefLinesFromSeqAlignEx2(curr_seqalign, 80, outfp, 
					print_options, FIRST_PASS, NULL, 
					number_of_descriptions, NULL, NULL);
                       free_buff();
                    } /* print deflines, looped if concat */

                    for (sap_iter=0; sap_iter < num_iters; sap_iter++) {
		       /* AM: Query concatenation. */
		       if( mult_queries && mask_loc )
		       {
		         orig_mask_loc = mask_loc;
			 
			 if( !mask_loc->data.ptrvalue ) mask_loc = NULL;
		       }

                       curr_seqalign = (num_queries > 0) ? *(sap_array + sap_iter) : seqalign;
                       curr_seqannot = (num_queries > 0) ? *(seq_annot_arr + sap_iter) : seqannot;

                       prune = BlastPruneHitsFromSeqAlign(curr_seqalign,
					number_of_alignments, NULL);
                       curr_seqannot->data = prune->sap;

                       if(options->is_ooframe) {
                          OOFShowBlastAlignment(curr_seqalign, /*mask*/ NULL,
                       			outfp, align_options, txmatrix);
                       } else {
                          if (align_view != 0)
                             ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL, 
                       			align_options, txmatrix, mask_loc, NULL);
                          else
                             ShowTextAlignFromAnnot(curr_seqannot, 60, outfp, NULL, NULL,
                        		align_options, txmatrix, mask_loc, 
                       			FormatScoreFunc);
                       }
                    
                       curr_seqannot->data = curr_seqalign;
                       prune = BlastPruneSapStructDestruct(prune);

                       /* AM: Query concatenation. */
		       if( mult_queries && orig_mask_loc ) 
		       {
		         mask_loc = orig_mask_loc;
		         mask_loc = mask_loc->next;
                       }
                    } /* show text align, loop over seqalign/seqannots for concat */
                    ObjMgrClearHold();
                 } /* if outfp */
                 for (sap_iter=0; sap_iter < num_queries; sap_iter++) {
                    /* upper bound is num_queries, take care not to do this unless concat */
                    *(seq_annot_arr + sap_iter) = SeqAnnotFree(*(seq_annot_arr + sap_iter)); 
                 }
                 if (mult_queries) 
                     seq_annot_arr = MemFree(seq_annot_arr);
	/*--KM free seqalign array and all seqaligns?? */

              } /* end of else (not XML Printing) */
              if (options->is_megablast_search)
                 tmp_slp = tmp_slp->next;
        /* --KM watch for memory leaks */
              if (seqannot && num_queries == 0)   
                 seqannot = SeqAnnotFree(seqannot);
              seqalign = next_seqalign;
           } /* End of loop on all seqaligns */
           if (mask_loc && next_mask_loc)
              mask_loc->next = next_mask_loc;

           } /* end of align_view not tabular case */
        } else {         /* seqalign is NULL */
           if((align_view == 7) && !options->is_ooframe) {
              BlastErrorMsgPtr error_msg;
              CharPtr message;
              
              if (error_returns == NULL) {
                 message = "No hits found";
              } else {
                 error_msg = error_returns->data.ptrvalue;
                 message = error_msg->msg;
              }
              if (options->is_megablast_search) {
                 bsp = BioseqLockById(SeqLocId(tmp_slp));
                 BXMLPrintOutput(xml_aip, seqalign, 
                                 options, blast_program, blast_database, 
                                 bsp, other_returns, 0, NULL, mask_loc);
                 BioseqUnlock(bsp);
              } else {
                 BXMLPrintOutput(xml_aip, NULL, options, blast_program, 
                                 blast_database, fake_bsp, other_returns, 0, 
                                 message, mask_loc);
              }
              AsnIoReset(xml_aip);
           } else if (align_view < 8) {
              fprintf(outfp, "\n\n ***** No hits found ******\n\n");
           }
           if (error_returns != NULL) {
              for (vnp = error_returns; vnp; vnp = vnp->next) {
                 BlastDestroyErrorMessage((BlastErrorMsgPtr)vnp->data.ptrvalue);
              }
              ValNodeFree(error_returns);
           }
        }
        
        slp = SeqLocSetFree(slp);
        matrix = BLAST_MatrixDestruct(matrix);
        if (txmatrix)
           txmatrix = TxMatrixDestruct(txmatrix);
        
        if(html) {
           fprintf(outfp, "<PRE>\n");
        }
        
        init_buff_ex(85);
        dbinfo_head = dbinfo;
        
        if(align_view < 7 && done) {
           while (dbinfo) {
              PrintDbReport(dbinfo, 70, outfp);
              dbinfo = dbinfo->next;
           }
        }
        dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
        
        if (ka_params) {
           if(align_view < 7 && done) {
              PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
           }
           MemFree(ka_params);
        }
        
        if (ka_params_gap) {
           if(align_view < 7 && done) {
              PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
           }
           MemFree(ka_params_gap);
        }
        
        if(align_view < 7 && done) {
           PrintTildeSepLines(params_buffer, 70, outfp);
        }
        
        MemFree(params_buffer);
        free_buff();
        mask_loc = mask_loc_start;
        while (mask_loc) {
           SeqLocSetFree(mask_loc->data.ptrvalue);
           mask_loc = mask_loc->next;
        }
        ValNodeFree(mask_loc_start);
        
        if(num_queries > 0) { /* AM: query concatenation */
            BSFree(fake_bsp->seq_data);
            fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
        } else if(!believe_query ) {
            fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
        }
        other_returns = ValNodeFree(other_returns);
        if (done) 
           sep = SeqEntryFree(sep);
#ifndef BLAST_CS_API
        /* This is freed earlier in client-server case */
        options->query_lcase_mask = SeqLocSetFree(options->query_lcase_mask);
        /* Free the database translation tables, if applicable. */
        TransTableFreeAll();
        ReadDBBioseqFetchDisable();
#endif
        if (html)
           fprintf(outfp, "</PRE>\n<P><HR><BR>\n<PRE>");
        
        if (!options->is_megablast_search) 
           BlastDeleteUserErrorString(err_ticket);
                    
        ObjMgrFreeCache(0);
    } /* while(TRUE)  - main loop of the program over all FASTA entries */
    
#ifdef BLAST_CS_API
    BlastNetBioseqFetchDisable(bl3hp, blast_database, db_is_na);
    BlastFini(bl3hp);
#endif
    
    aip = AsnIoClose(aip);
    
    if(align_view < 7) {
        if (html) {
            fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");
        }
    } else if (align_view == 7)
        xml_aip = AsnIoClose(xml_aip);
    
    /* AM: query concatenation. */
    mult_queries = BlastMultQueriesDestruct( mult_queries );

    options = BLASTOptionDelete(options);
    FileClose(infp);
    return 0;
}

/*
        This function decides whether the new blast code can handle this database or not.
        Currently it should return FALSE for any database that uses a gilist.
        This implementation only works for nucleotide databases.

        If it is not possible to initialize the database or some error condition exists then FALSE
        will also be returned and the old engine should deal with this.
*/
static Boolean
readdb_use_new_blast(char* dbname)
{
      Boolean db_is_na, query_is_na;
      Boolean retval=TRUE;
      ReadDBFILEPtr rdfp=NULL;
      ReadDBFILEPtr rdfp_var=NULL;

      if (!dbname)
           return FALSE;

      BlastGetTypes(myargs[ARG_PROGRAM].strvalue, &query_is_na, &db_is_na);
      rdfp = readdb_new(dbname, !db_is_na);
      if (!rdfp)
           return FALSE;

      rdfp_var = rdfp;
      while (rdfp_var)
      {
            if (rdfp_var->gilist != NULL)
            {
                   retval = FALSE;
                   break;  /* Break out and free rdfp. */
            }
            rdfp_var = rdfp_var->next;
      }
      rdfp = readdb_destruct(rdfp);
      return retval;
}

Int2 Nlm_Main(void)
{
#ifndef BLASTALL_TOOLS_ONLY
    Boolean use_new_engine=FALSE;
#endif
    char buf[256] = { '\0' };

#ifdef BLAST_CS_API
    StringCpy(buf, "blastcl3 ");
    StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf)-1);
    if (! GetArgs (buf, NUMARG, myargs)) {
        return (1);
    }
#else
    StringCpy(buf, "blastall ");
    StringNCat(buf, BlastGetVersionNumber(), sizeof(buf)-StringLen(buf));
    if (! GetArgs (buf, NUMARG, myargs)) {
        return (1);
    }
#endif

    UseLocalAsnloadDataAndErrMsg ();

    if (! SeqEntryLoad())
                return 1;

    ErrSetMessageLevel(SEV_WARNING);

#ifdef BLAST_CS_API
    return Main_old();
#else
#ifndef BLASTALL_TOOLS_ONLY
    if (myargs[ARG_FORCE_OLD].intvalue == 0 &&
                myargs[ARG_PSITCHKPNT].strvalue == NULL &&
                      myargs[ARG_GILIST].strvalue == NULL)
          use_new_engine = readdb_use_new_blast(myargs[ARG_DB].strvalue);

    if (use_new_engine)
        return Main_new();
    else
#endif /* BLASTALL_TOOLS_ONLY */
        return Main_old();
#endif
}
