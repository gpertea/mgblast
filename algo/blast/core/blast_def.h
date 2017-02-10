/* $Id: blast_def.h,v 1.70 2006/04/11 15:47:19 camacho Exp $
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

/** @file blast_def.h
 * Definitions of major structures used throughout BLAST
 */

#ifndef __BLAST_DEF__
#define __BLAST_DEF__

#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_export.h>

#ifdef __cplusplus
extern "C" {
#endif

/****************************** Constants *********************************/

extern const int kDustLevel;  /**< Level parameter used by dust. */
extern const int kDustWindow; /**< Window parameter used by dust. */
extern const int kDustLinker; /**< Parameter used by dust to link together close low-complexity segments. */

extern const int kSegWindow;  /**< Window that SEG examines at once. */
extern const double kSegLocut;   /**< Locut parameter for SEG. */
extern const double kSegHicut;   /**< Hicut parameter for SEG. */

/** Maximum number of HPSs to be saved in an ungapped search.
 * Value defined in blast_options.c
 */
extern const int kUngappedHSPNumMax; 

/******************** Preprocessor definitions ******************************/

/** Codons are always of length 3 */
#ifndef CODON_LENGTH
#define CODON_LENGTH 3
#endif

/** For translated gapped searches, this is the default value in
 * nucleotides of longest_intron (for ungapped translated searches,
 * the default value of longest_intron is zero, which causes a legacy
 * method of HSP linking that does not use longest_intron to be
 * invoked).
 *
 * The value 122 corresponds to 40 amino acids: 40 codons * 3
 * nucleotides per codon + up to 2 frame shifts.  40 amino acids is
 * the maximum gap size in the untranslated sequence, so
 * DEFAULT_LONGEST_INTRON makes these two gap sizes equal.
 */ 
#ifndef DEFAULT_LONGEST_INTRON
#define DEFAULT_LONGEST_INTRON 122
#endif

/** Compression ratio of nucleotide bases (4 bases in 1 byte) */
#ifndef COMPRESSION_RATIO
#define COMPRESSION_RATIO 4
#endif

/** Number of frames to which we translate in translating searches */
#ifndef NUM_FRAMES
#define NUM_FRAMES 6
#endif

/** Number of frames in a nucleotide sequence */
#ifndef NUM_STRANDS
#define NUM_STRANDS 2
#endif

/**
 * A macro expression that returns 1, 0, -1 if a is greater than,
 * equal to or less than b, respectively.  This macro evaluates its
 * arguments more than once.
 */
#ifndef BLAST_CMP
#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))
#endif

/** Safe free a pointer: belongs to a higher level header. */
#ifndef sfree
#define sfree(x) __sfree((void**)&(x))
#endif

/** Implemented in lookup_util.c. @sa sfree */
NCBI_XBLAST_EXPORT
void __sfree(void** x);

/********************* Structure definitions ********************************/

/** A structure containing two integers, used e.g. for locations for the 
 * lookup table.
 */
typedef struct SSeqRange {
   Int4 left;  /**< left endpoint of range (zero based) */
   Int4 right;  /**< right endpoint of range (zero based) */
} SSeqRange;

/** Used to hold a set of positions, mostly used for filtering. 
 * oid holds the index of the query sequence.
*/
typedef struct BlastSeqLoc {
        struct BlastSeqLoc *next;  /**< next in linked list */
        SSeqRange *ssr;            /**< location data on the sequence. */
} BlastSeqLoc;

/** Structure for keeping the query masking information */
typedef struct BlastMaskLoc {
   /** Total size of the BlastSeqLoc array below. This is always the number 
     of queries times the number of contexts. Note that in the case of 
     translated query searches, these locations must be provided in protein 
     coordinates to BLAST_MainSetUp.
     @sa BLAST_GetNumberOfContexts 
     @sa BlastMaskLocDNAToProtein
    */
   Int4 total_size; 

   /** Array of masked locations. 
     Every query is allocated the number of contexts associated with the 
     program. In the case of nucleotide searches, the strand(s) to search 
     dictatate which elements of the array for a given query are filled. For 
     translated searches, this should also be the same (by design) but the 
     C toolkit API does NOT implement this, it rather fills all elements 
     for all queries with masked locations in protein coordinates (if any). 
     The C++ API does follow the convention which populates each element, only
     if so dictated by the strand(s) to search for each query.
     @sa BLAST_GetNumberOfContexts
    */
   BlastSeqLoc** seqloc_array; 
} BlastMaskLoc;

/** Structure to hold a sequence. */
typedef struct BLAST_SequenceBlk {
   Uint1* sequence; /**< Sequence used for search (could be translation). */
   Uint1* sequence_start; /**< Start of sequence, usually one byte before 
                               sequence as that byte is a NULL sentinel byte.*/
   Int4     length;         /**< Length of sequence. */
   Int4 context; /**< Context of the query, needed for multi-query searches */
   Int2 frame; /**< Frame of the query, needed for translated searches */
   Int4 oid; /**< The ordinal id of the current sequence */
   Boolean sequence_allocated; /**< TRUE if memory has been allocated for 
                                  sequence */
   Boolean sequence_start_allocated; /**< TRUE if memory has been allocated 
                                        for sequence_start */
   Uint1* oof_sequence; /**< Mixed-frame protein representation of a
                             nucleotide sequence for out-of-frame alignment */
   Boolean oof_sequence_allocated; /**< TRUE if memory has been allocated 
                                        for oof_sequence */
   BlastMaskLoc* lcase_mask; /**< Locations to be masked from operations on 
                                this sequence: lookup table for query; 
                                scanning for subject. */
   Boolean lcase_mask_allocated; /**< TRUE if memory has been allocated for 
                                    lcase_mask */
} BLAST_SequenceBlk;

/** Information about a single pattern occurence in the query. */
typedef struct SPHIPatternInfo {
    Int4 offset;  /**< Starting offset of this pattern occurrence. */
    Int4 length;  /**< Length of this pattern occurrence. */
} SPHIPatternInfo;

/** In PHI BLAST, structure containing information about all pattern 
 * occurrences in query.
 */
typedef struct SPHIQueryInfo {
    Int4 num_patterns;  /**< Number of pattern occurrences in query. */
    SPHIPatternInfo *occurrences; /**< Array of pattern occurrence information
                                        structures. */
    Int4 allocated_size; /**< Allocated size of the occurrences array. */
    double probability; /**< Probability of the pattern */
} SPHIQueryInfo;

/************************* Progress monitoring/interruptible API *************/

/** Enumeration for the stages in the BLAST search */
typedef enum EBlastStage {
    ePrelimSearch,
    eTracebackSearch
} EBlastStage;

/** Progress monitoring structure. This is updated by the engine to provided to
 * the user as an argument to the user-supplied callback function 
 * (TInterruptFnPtr). This function then can assess whether the search 
 * should proceed or exit prematurely.
 * @sa TInterruptFnPtr
 */
typedef struct SBlastProgress {
    EBlastStage stage;      /**< Stage of the BLAST search currently in
                              progress */
    void* user_data;        /**< Pointer to user-provided data */
} SBlastProgress;

/** Prototype for function pointer to determine whether the BLAST search
 * should proceed or be interrupted. If this function returns true, all 
 * processing must stop and the search must discard all interim results 
 * @note In order to avoid undue overhead, this function should not perform any
 * time consuming operations and should always return (i.e.: it should never 
 * block)
 */
typedef Boolean (*TInterruptFnPtr) (SBlastProgress* progress_info);

/** Allocates and initializes a new SBlastProgress structure.
 * @param user_data user-provided data (not owned by the resulting structure)
 * [in]
 * Implemented in blast_util.c 
 */
SBlastProgress* SBlastProgressNew(void* user_data);

/** Deallocates a SBlastProgress structure.
 * Implemented in blast_util.c */
SBlastProgress* SBlastProgressFree(SBlastProgress* progress_info);

/** Resets the progress structure to its original state (as if newly allocated)
 * for a fresh start without touching the user_data field */
void SBlastProgressReset(SBlastProgress* progress_info);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_DEF__ */
