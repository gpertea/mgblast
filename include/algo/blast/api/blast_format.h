/* $Id: blast_format.h,v 1.44 2006/04/25 17:56:54 papadopo Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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

/** @file blast_format.h
 * Functions needed for formatting of BLAST results
 */

#ifndef __BLAST_FORMAT__
#define __BLAST_FORMAT__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <ncbi.h>
#include <asn.h>
#include <bxmlobj.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_diagnostics.h>   
#include <algo/blast/api/twoseq_api.h>
#include <algo/blast/api/blast_options_api.h>
#include <algo/blast/api/blast_seqalign.h>
#include <xmlblast.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

/** Enums for selecting view type (query-anchored, pairwise, etc.)
*/
typedef enum EAlignView {
    /** Pairwise. */
    eAlignViewPairwise                   = 0,
    /** Query anchored with identities. */
    eAlignViewQueryAnchoredIdent         = 1,
    /** Query anchored without identities. */
    eAlignViewQueryAnchoredNoIdent       = 2,
    /** Flat query anchored with identities. */
    eAlignViewFlatQueryAnchoredIdent     = 3,
    /** Flat query anchored without identities. */
    eAlignViewFlatQueryAnchoredNoIdent   = 4,
    /** Query anchored with blunt ends. */
    eAlignViewQueryAnchoredBluntEnds     = 5,
    /** Flat query anchored with blunt ends. */ 
    eAlignViewFlatQueryAnchoredBluntEnds = 6,
    /** XML. */
    eAlignViewXml                        = 7,
    /** Tabular without comments. */
    eAlignViewTabular                    = 8,
    /** Tabular with per-query comments. */
    eAlignViewTabularWithComments        = 9,
    /** ASN.1 in text form. */
    eAlignViewAsnText                    = 10,
    /** ASN.1 in binary form. */
    eAlignViewAsnBinary                  = 11,
    /** Sentinel value, binding the allowed range. */
    eAlignViewMax
} EAlignView;


/** Options for formatting BLAST results 
 */
typedef struct BlastFormattingOptions {
   EAlignView align_view;       /**< How to show alignments? */
   Uint4 align_options;         /**< Options for showing alignments. */
   Uint4 print_options;         /**< Printing options. */
   Boolean believe_query;       /**< Should query defline be trusted? */
   Boolean html;                /**< Create an HTML output? */
   Int4 number_of_descriptions; /**< Number of descriptions to show. */
   Int4 number_of_alignments;   /**< Number of alignments to show. */
   Boolean is_megablast;        /**< Is this a megablast search? Needed for
                                   determination of which reference to use. */
} BlastFormattingOptions;

/** Structure containing all information necessary for producing the
 * formatted BLAST output. 
 */
typedef struct BlastFormattingInfo { 
    const SBlastOptions* search_options; /**< Search options. */
    char* program_name;          /**< Program name as a sting. */
    char* db_name;               /**< Database name. Null for a 2 sequences 
                                    search. */
    int num_formatted;           /**< Number of queries formatted so far. */
    BlastFormattingOptions* format_options; /**< Formatting options. */
    FILE* outfp;                 /**< Output stream for text output. */
    AsnIo* aip;                  /**< Output stream for ASN.1 output */
    MBXml* xmlp;                 /**< Output stream for XML output */
    Boolean is_seqalign_null;    /**< flag indicating absence of seqalign */
} BlastFormattingInfo;

/** Allocates and initializes the formatting information structure.
 * @param align_view Alignment view option [in]
 * @param search_options BLAST search options [in]
 * @param program_name BLAST program name [in]
 * @param db_name BLAST database name (NULL for non-database search) [in]
 * @param outfile_name Name of the output file. [in]
 * @param info_ptr Initialized structure. [out]
 */
Int2 BlastFormattingInfoNew(EAlignView align_view,
                            const SBlastOptions* search_options,
                            const char* program_name, const char* db_name, 
                            const char* outfile_name, 
                            BlastFormattingInfo* *info_ptr);

/** Allocate and initialize the formatting information structure, given the 
 * basic two sequences search options.
 * @param align_view What kind of formatted output to show? [in]
 * @param basic_options Basic two sequences search options [in]
 * @param query_seqloc Query Seq-loc [in]
 * @param outfile_name Name of the output file [in]
 * @param advanced_options Advanced BLAST search options structure. Must be 
 *                         returned to avoid memory leak. [out]
 * @param info_ptr The initialized formatting info structure [out]
*/
Int2 BlastFormattingInfoNewBasic(EAlignView align_view, 
                                 const BLAST_SummaryOptions* basic_options,
                                 SeqLoc* query_seqloc,
                                 const char* outfile_name, 
                                 SBlastOptions* *advanced_options,
                                 BlastFormattingInfo* *info_ptr);

/** Fill the formatting options in the formatting info structure.
 * @param info The structure to work with [in]
 * @param num_descriptions Number of descriptions to show in the formatted
 *                         output. [in]
 * @param num_alignments Number of alignments to show in the formatted
 *                       output. [in]
 * @param html Format with HTML tags or in plain text? [in]
 * @param is_megablast Is this a megablast search, needed for correct 
 *                     reference? [in]
 * @param show_gi Should gis be shown in database sequence descriptions? [in]
 * @param believe_defline Should query def-lines be parsed? [in]
 */ 
void 
BlastFormattingInfoSetUpOptions(BlastFormattingInfo* info, Int4 num_descriptions,
                                Int4 num_alignments, Boolean html, 
                                Boolean is_megablast, Boolean show_gi,
                                Boolean believe_defline);

/** Deallocates the formatting information structure. In particular, closes the
 * output stream(s).
 * @param info The structure to free. [in]
 */
BlastFormattingInfo* BlastFormattingInfoFree(BlastFormattingInfo* info);

/** Print formatted output.
 * @param seqalign_arr object with arrays of seqaligns [in]
 * @param num_queries Number of query sequences [in]
 * @param query_slp Linked list of query SeqLocs [in]
 * @param mask_loc Masking locations for all queries [in]
 * @param format_info Formatting options and other information [in]
 * @param sum_returns Summary data returned from the search. [in]
 */
Int2 BLAST_FormatResults(SBlastSeqalignArray* seqalign_arr, Int4 num_queries, 
                         SeqLocPtr query_slp, SeqLoc* mask_loc, 
                         BlastFormattingInfo* format_info,
                         Blast_SummaryReturn* sum_returns);

/** Print the summary at the end of the BLAST report.
 * @param format_info Formatting options and other information. [in]
 * @param sum_returns infor from inside blast engine [in]
 */
Int2 Blast_PrintOutputFooter(const BlastFormattingInfo* format_info, 
                             const Blast_SummaryReturn* sum_returns);

/** Prints the top part of the traditional BLAST output, including version, 
 * reference(s) and database information.
 * @param format_info Formatting options and other information [in]
 */
Int2 BLAST_PrintOutputHeader(const BlastFormattingInfo* format_info);

/** Prints the "log info" at the bottom of the megablast output, looks something
 * like:
 * "Mega BLAST run finished, processed 2 queries"
 *
 * @param outfp output file pointer [in]
 * @param num_processed number of records proceessed [in]
 *                    stream pointer. [in]
 * @return 0 on success 
 */
Int2 BlastPrintLogReport(FILE* outfp, Int4 num_processed);

/** Given a Seq-id structure, returns a buffer in requested form.
 * @param sip Seq-id to get information from. [in]
 * @param buffer_ptr Buffer to fill. [out]
 * @param ncbi_gi Should the NCBI gi be included in the id buffer? [in]
 * @param accession_only Should only accession be returned (only gi if ncbi_gi
 *                       is TRUE)? [in]
 * @param search_for_id If TRUE, attempt to find a sequence identifier.
 *                  If FALSE, the first token in the defline is treated
 *                  as an identifier, even if it is not an actual ID [in]
 */
void 
Blast_SeqIdGetDefLine(SeqIdPtr sip, char** buffer_ptr, Boolean ncbi_gi, 
                      Boolean accession_only, Boolean search_for_id);

/** Posts an error message from the summary returns, if any, to a provided
 * output stream. If no stream provided, nothing is written. 
 * @param sum_return Summary returns structure. [in]
 */
Int2 
Blast_SummaryReturnsPostError(Blast_SummaryReturn* sum_return);

/** Format the PHI BLAST report.
 * @param phivnps List of ValNode's, containing PHI BLAST Seq-aligns 
 *                corresponding to different pattern occurrences in query. [in]
 * @param query_slp Query Seq-loc [in]
 * @param format_info Formatting information structure [in]
 * @param sum_returns Search summary returns [in]
 */
Int2 PHIBlastFormatResults(ValNode* phivnps, SeqLocPtr query_slp,
                           const BlastFormattingInfo* format_info, 
                           Blast_SummaryReturn* sum_returns);

/** Frees a chain of ValNode's containing PHI BLAST results.
 * @param phivnps Head of chain to be freed. [in]
 * @return NULL.
 */
ValNode* PHIBlastResultsFree(ValNode* phivnps);

/* @} */

#ifdef __cplusplus
}
#endif

#endif /* !__BLAST_FORMAT__ */

