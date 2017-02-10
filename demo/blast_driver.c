/* $Id: blast_driver.c,v 1.110 2006/04/26 12:48:06 madden Exp $
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
* ===========================================================================*/

/*****************************************************************************

File name: blast_driver.c

Author: Ilya Dondoshansky

Contents: Main function for running BLAST

******************************************************************************
 * $Revision: 1.110 $
 * */

static char const rcsid[] = "$Id: blast_driver.c,v 1.110 2006/04/26 12:48:06 madden Exp $";

#include <ncbi.h>
#include <sqnutils.h>
#include <readdb.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/api/blast_input.h>
#include <algo/blast/api/blast_format.h>
#include <algo/blast/api/blast_seqalign.h>
#include <algo/blast/api/blast_tabular.h>
#include <algo/blast/api/blast_api.h>
#include <algo/blast/api/repeats_filter.h>

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

typedef enum {
   ARG_PROGRAM = 0,
   ARG_DB,
   ARG_QUERY,
   ARG_QUERY_LOC,
   ARG_SUBJECT,
   ARG_SUBJECT_LOC,
   ARG_STRAND,
   ARG_GENCODE,
   ARG_DBGENCODE,
   ARG_FILTER,
   ARG_LCASE,
   ARG_LOOKUP,
   ARG_MATRIX,
   ARG_MISMATCH,
   ARG_MATCH,
   ARG_WORDSIZE,
   ARG_TEMPL_LEN,
   ARG_TEMPL_TYPE,
   ARG_EVERYBASE,
   ARG_PHI,
   ARG_THRESHOLD,
   ARG_WINDOW,
   ARG_XDROP_UNGAPPED,
   ARG_UNGAPPED,
   ARG_GREEDY,
   ARG_GAPOPEN,
   ARG_GAPEXT,
   ARG_FRAMESHIFT,
   ARG_XDROP,
   ARG_XDROP_FINAL,
   ARG_EVALUE,
   ARG_SEARCHSP,
   ARG_PERC_IDENT,
   ARG_INTRON,
   ARG_DESCRIPTIONS,
   ARG_CULLING,
   ARG_ALIGNMENTS,
   ARG_OUT,
   ARG_FORMAT,
   ARG_HTML,
   ARG_TABULAR,
   ARG_THREADS,
   ARG_SHOWGI,
   ARG_ACCESSION,
   ARG_COMP_BASED_STATS,
   ARG_SMITH_WATERMAN
} BlastArguments;

static Args myargs[] = {
   { "Program Name",           /* ARG_PROGRAM */
      NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
   { "Database name (if not set, second sequence FASTA must be provided)",
     NULL, NULL, NULL, TRUE, 'd', ARG_STRING, 0.0, 0, NULL}, /* ARG_DB */
   { "Query File",               /* ARG_QUERY */
     "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
   { "Query location offsets; format: \"start stop\";\n"
     "Applies only if query file contains 1 sequence", /* ARG_QUERY_LOC */
     NULL, NULL, NULL, TRUE, 'I', ARG_STRING, 0.0, 0, NULL},
   { "Subject File (for sequence sets comparison)", /* ARG_SUBJECT */
     "stdin", NULL, NULL, FALSE, 'j', ARG_FILE_IN, 0.0, 0, NULL},
   { "Subject location offsets; format: \"start stop\";\n"/* ARG_SUBJECT_LOC */
     "Applies only if subject file (-j) contains 1 sequence",
     NULL, NULL, NULL, TRUE, 'J', ARG_STRING, 0.0, 0, NULL},
   { "Query strands to search against database: 0 or 3 is both, 1 is top, "
     "2 is bottom", /* ARG_STRAND */
     "0", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
   { "Genetic code for translation of the query sequence", /* ARG_GENCODE */
     "0", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
   { "Genetic code for translation of the database", /* ARG_DBGENCODE */
     "0", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},
   { "Filter query sequence (DUST with blastn, SEG with others)", /* ARG_FILTER */
     "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
   { "Mask lower case",
     "F", NULL, NULL, FALSE, 'c', ARG_BOOLEAN, 0.0, 0, NULL}, /* ARG_LCASE */
   { "Use (classical Mega BLAST) lookup table with width 12", 
     "F", NULL, NULL, FALSE, 'L', ARG_BOOLEAN, 0.0, 0, NULL},/* ARG_LOOKUP */
   { "Matrix",                 /* ARG_MATRIX */
     "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
   { "Penalty for a nucleotide mismatch (blastn only)", /* ARG_MISMATCH */
     "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},
   { "Reward for a nucleotide match (blastn only)", /* ARG_MATCH */
     "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},
   { "Word size, default if 0 (blastn 11, others 3) ", /* ARG_WORDSIZE */  
     "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL}, 
   { "Length of a discontiguous word template (contiguous word if 0)",
     "0", NULL, NULL, FALSE, 't', ARG_INT, 0.0, 0, NULL}, /* ARG_TEMPL_LEN */
   { "Type of a discontiguous word template (0 - coding, 1 - optimal, "
     "2 - two simultaneous", /* ARG_TEMPL_TYPE */
     "0", NULL, NULL, FALSE, 'T', ARG_INT, 0.0, 0, NULL},
   {"Generate words for every base of the database (default is every 4th base; may only be used with discontiguous words)",
        "F", NULL, NULL, TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},    /* ARG_EVERYBASE */
   { "Pattern for PHI BLAST",
     NULL, NULL, NULL, TRUE, 'k', ARG_STRING, 0.0, 0, NULL}, /* ARG_PHI */
   { "Threshold for extending hits, default if zero\n" /* ARG_THRESHOLD */
     "      blastp 11, blastn 0, blastx 12, tblastn 13\n"
     "      tblastx 13, megablast 0",
     "0", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
   { "Window size (max. allowed distance between a pair of initial hits;\n"
     "      0 causes default behavior, -1 turns off multiple hits)", 
     "0", NULL, NULL, FALSE, 'w', ARG_INT, 0.0, 0, NULL}, /* ARG_WINDOW */
   { "X dropoff value for ungapped extensions in bits (0 invokes default "
     "behavior)\n      blastn 20, others 7",/*ARG_XDROP_UNGAPPED*/
      "0", NULL, NULL, FALSE, 'y', ARG_INT, 0.0, 0, NULL},
   { "Do only ungapped alignment (always TRUE for tblastx)",/*ARG_UNGAPPED*/
     "F", NULL, NULL, FALSE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
   { "Use greedy algorithm for gapped extensions:\n      0 no, 1 one-step, "
     "2 two-step, 3 two-step with ungapped", /* ARG_GREEDY */
     "0", NULL, NULL, FALSE, 'g', ARG_INT, 0.0, 0, NULL},
   { "Gap open penalty (-1 means default: non-affine if greedy; 5 if dyn. prog.)", 
     "-1", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL}, /* ARG_GAPOPEN */
   { "Gap extension penalty (-1 means default: non-affine if greedy; 2 otherwise)",
     "-1", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL}, /* ARG_GAPEXT */
   { "Frame shift penalty for out-of-frame gapping (blastx, tblastn only)",
     "0", NULL, NULL, FALSE, 'h', ARG_INT, 0.0, 0, NULL}, /* ARG_FRAMESHIFT */
   { "X dropoff value for gapped alignment (in bits) (zero invokes default "
     "behavior)\n      blastn 30, tblastx 0, others 15", /* ARG_XDROP */
     "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
   { "X dropoff value for final gapped alignment in bits "
     "(0 invokes default behavior)\n"
     "      blastn 50, tblastx 0, others 25",  /* ARG_XDROP_FINAL */
     "0", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
   { "Expected value",                         /* ARG_EVALUE */
     "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
   { "Effective length of the search space (use zero for the real size)", 
     "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL}, /* ARG_SEARCHSP */
   {"Identity percentage cut-off",  /* ARG_PERC_IDENT */
    "0", NULL, NULL, FALSE, 'P', ARG_FLOAT, 0.0, 0, NULL},
   { "Longest intron length for uneven gap HSP linking (tblastn only)",
     "0", NULL, NULL, FALSE, 'z', ARG_INT, 0.0, 0, NULL}, /* ARG_INTRON */
   { "Number of database sequences to show one-line descriptions for (V)",
     "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL}, /* ARG_DESCRIPTIONS */
   { "Remove hits whose query range is contained in at least "
      "this many higher-scoring hits (ignored if zero)", /* ARG_CULLING */
     "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},
   { "Number of database sequence to show alignments for (B)", /* ARG_ALIGNMENTS */
     "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
   { "Final output file name",             /* ARG_OUT */
     "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL}, 
   { "alignment view options:\n0 = pairwise,\n1 = query-anchored showing "
     "identities,\n2 = query-anchored no identities,\n3 = flat "
     "query-anchored, show identities,\n4 = flat query-anchored, no "
     "identities,\n5 = query-anchored no identities and blunt ends,\n6 = "
     "flat query-anchored, no identities and blunt ends,\n7 = XML Blast "
     "output,\n8 = tabular, \n9 tabular with comment lines\n10 ASN, text\n"
     "11 ASN, binary",                         /* ARG_FORMAT */
     "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
   { "Produce HTML output",                    /* ARG_HTML */
     "F", NULL, NULL, FALSE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
   { "Produce on-the-fly tabular output; 1 - just offsets and quality values;\n"
     "2 - add sequence data.",
     "0", NULL, NULL, FALSE, 'B', ARG_INT, 0.0, 0, NULL}, /* ARG_TABULAR */
   { "Number of threads to use in preliminary search stage",
     "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL}, /* ARG_THREADS */
   { "Show gis in sequence ids",
     "F", NULL, NULL, FALSE, 'n', ARG_BOOLEAN, 0.0, 0, NULL}, /* ARG_SHOWGI */
   { "Show only accessions for sequence ids in tabular output",
     "F", NULL, NULL, FALSE, 'N', ARG_BOOLEAN, 0.0, 0, NULL},/* ARG_ACCESSION */
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
     "F", NULL, NULL, FALSE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
};

extern void PrintTabularOutputHeader PROTO((char* blast_database, 
               BioseqPtr query_bsp, SeqLocPtr query_slp, char* blast_program, 
               Int4 iteration, Boolean believe_query, FILE *outfp));

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
   Int4 greedy_extension = 0;
   Int4 diag_separation = 0;
   Boolean greedy_with_ungapped = FALSE;
   Boolean is_gapped = FALSE;
   EBlastProgramType program_number = options->program;

   /* The following options are for blastn only */
   if (program_number == eBlastTypeBlastn) {
      if (myargs[ARG_TEMPL_LEN].intvalue == 0) {
         mb_lookup = (Boolean) myargs[ARG_LOOKUP].intvalue;
      } else {
         /* Discontiguous words */
         mb_lookup = TRUE;
      }
      greedy_extension = MIN(myargs[ARG_GREEDY].intvalue, 2);
      greedy_with_ungapped = (myargs[ARG_GREEDY].intvalue == 3);
   }

   BLAST_FillLookupTableOptions(lookup_options, program_number, mb_lookup,
      myargs[ARG_THRESHOLD].intvalue, (Int2)myargs[ARG_WORDSIZE].intvalue);
   /* Fill the rest of the lookup table options */
   lookup_options->mb_template_length = 
      (Uint1) myargs[ARG_TEMPL_LEN].intvalue;
   lookup_options->mb_template_type = 
      (Uint1) myargs[ARG_TEMPL_TYPE].intvalue;
   if (myargs[ARG_EVERYBASE].intvalue)
      lookup_options->full_byte_scan = FALSE;

   if (myargs[ARG_PHI].strvalue) {
      lookup_options->phi_pattern = strdup(myargs[ARG_PHI].strvalue);
      /* Set the lookup table type, and also change program type to 
         indicate a PHI BLAST search. */
      if (program_number == eBlastTypeBlastn) {
          lookup_options->lut_type = PHI_NA_LOOKUP;
          options->program = eBlastTypePhiBlastn;
      } else {
          lookup_options->lut_type = PHI_AA_LOOKUP;
          options->program = eBlastTypePhiBlastp;
      }
   }

   BLAST_FillQuerySetUpOptions(query_setup_options, program_number, 
      myargs[ARG_FILTER].strvalue, (Uint1)myargs[ARG_STRAND].intvalue);

   if (myargs[ARG_GENCODE].intvalue &&
       (program_number == eBlastTypeBlastx || 
        program_number == eBlastTypeTblastx))
      query_setup_options->genetic_code = myargs[ARG_GENCODE].intvalue;

   BLAST_FillInitialWordOptions(word_options, program_number, 
      (Boolean)(greedy_extension && !greedy_with_ungapped), 
      myargs[ARG_WINDOW].intvalue, myargs[ARG_XDROP_UNGAPPED].intvalue);

   if (!greedy_extension)
      word_options->ungapped_extension = TRUE;

   if (myargs[ARG_WINDOW].intvalue < 0)
       word_options->window_size = 0;

   BLAST_FillExtensionOptions(ext_options, program_number, greedy_extension, 
      myargs[ARG_XDROP].intvalue, myargs[ARG_XDROP_FINAL].intvalue);

   BLAST_FillScoringOptions(score_options, program_number, 
       (Boolean)greedy_extension, myargs[ARG_MISMATCH].intvalue, 
        myargs[ARG_MATCH].intvalue, myargs[ARG_MATRIX].strvalue, 
        myargs[ARG_GAPOPEN].intvalue, myargs[ARG_GAPEXT].intvalue);

   if (program_number != eBlastTypeTblastx)
      is_gapped = !myargs[ARG_UNGAPPED].intvalue;

   score_options->gapped_calculation = is_gapped;
   if (myargs[ARG_FRAMESHIFT].intvalue) {
      score_options->shift_pen = myargs[ARG_FRAMESHIFT].intvalue;
      score_options->is_ooframe = TRUE;
   }

   if (mb_lookup)
       diag_separation = 6;

   BLAST_FillHitSavingOptions(hit_options, 
      myargs[ARG_EVALUE].floatvalue, 
      MAX(myargs[ARG_DESCRIPTIONS].intvalue, 
          myargs[ARG_ALIGNMENTS].intvalue),
      is_gapped,
      myargs[ARG_CULLING].intvalue,
      diag_separation);
 
   hit_options->percent_identity = myargs[ARG_PERC_IDENT].floatvalue;
   hit_options->longest_intron = myargs[ARG_INTRON].intvalue;

   if (myargs[ARG_SEARCHSP].floatvalue != 0) {
      eff_len_options->searchsp_eff = (Int8) myargs[ARG_SEARCHSP].floatvalue; 
   }
   if ((program_number == eBlastTypeTblastn ||
        program_number == eBlastTypeBlastp) && is_gapped) {
       /* Set options specific to gapped tblastn */
       switch (myargs[ARG_COMP_BASED_STATS].strvalue[0]) {
       case 'D': case 'd':
       case '0': case 'F': case 'f':
           ext_options->compositionBasedStats = 0;
           break;
        case '1': case 'T': case 't':
            ext_options->compositionBasedStats = 1;
            break;
       case '2':
           ErrPostEx(SEV_WARNING, 1, 0, "the -C 2 argument "
                     "is currently experimental\n");
           ext_options->compositionBasedStats = 2;
           break;
       case '3':
           ErrPostEx(SEV_WARNING, 1, 0, "the -C 3 argument "
                     "is currently experimental\n");
           ext_options->compositionBasedStats = 3;
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
   if (program_number == eBlastTypeTblastn ||
       program_number == eBlastTypeRpsTblastn ||
       program_number == eBlastTypeTblastx) {
       SBlastOptionsSetDbGeneticCode(options, myargs[ARG_DBGENCODE].intvalue);
   }

   return 0;
}

/** Parses an argument specifying an interval on a sequence to search.
 * Argument has form "Start[ ,;]Stop"
 * @param arg Argument string (no restriction if NULL) [in]
 * @param from_ptr Start of the location [out]
 * @param to_ptr End of the location [out]
 * @return -1 if invalid location, otherwise 0.
 */
static Int2
s_ParseIntervalLocationArgument(char* arg, Int4* from_ptr, Int4* to_ptr)
{
   const char* delimiters = " ,;";
   Int4 from = 0, to = 0;

   if (!arg)
      return 0;
   
   from = atoi(StringTokMT(arg, delimiters, &arg));
   to = atoi(arg);
      
   from = MAX(from, 0);
   to = MAX(to, 0);

   *from_ptr = from;
   *to_ptr = to;

   if (to > 0 && from > to)
      return -1;
   else
      return 0;
}

Int2 Nlm_Main(void)
{
   SeqLoc* subject_slp = NULL; /* SeqLoc for the subject sequence in two
                                    sequences case */
   Boolean query_is_na, db_is_na;
   char buf[256] = { '\0' };
   char* blast_program;
   EBlastProgramType program_number;
   char* dbname = NULL;
   Int2 status = 0;
   SeqLoc* lcase_mask = NULL;
   SeqLoc* query_slp = NULL;
   FILE *infp, *outfp = NULL;
   SBlastOptions* options = NULL;
   BlastFormattingInfo* format_info = NULL;
   Int4 ctr = 1;
   Int4 num_queries=0;
   int tabular_output = FALSE;
   BlastTabularFormatData* tf_data = NULL;
   Boolean believe_query = FALSE;
   Int4 q_from = 0, q_to = 0;
   Blast_SummaryReturn* sum_returns = NULL;
   Boolean phi_blast = FALSE;
   ValNode* phivnps = NULL;
   Blast_SummaryReturn* full_sum_returns = NULL;

   if (! GetArgs (buf, NUMARG, myargs))
      return (1);
   
   UseLocalAsnloadDataAndErrMsg ();
   
   if (! SeqEntryLoad())
      return 1;
   
   ErrSetMessageLevel(SEV_WARNING);
   
   blast_program = myargs[ARG_PROGRAM].strvalue;

   sum_returns = Blast_SummaryReturnNew();
   status = SBlastOptionsNew(blast_program, &options, sum_returns);

   if (status) {
       if (sum_returns->error) {
           SBlastMessageErrPost(sum_returns->error);
           sum_returns = Blast_SummaryReturnFree(sum_returns);
       }
       return -1;
   }

   program_number = options->program;

   db_is_na = (program_number == eBlastTypeBlastn || 
               program_number == eBlastTypeTblastn || 
               program_number == eBlastTypeTblastx ||
               program_number == eBlastTypePhiBlastn);
   query_is_na = (program_number == eBlastTypeBlastn || 
                  program_number == eBlastTypeBlastx || 
                  program_number == eBlastTypeRpsTblastn || 
                  program_number == eBlastTypeTblastx ||
                  program_number == eBlastTypePhiBlastn);

   phi_blast = Blast_ProgramIsPhiBlast(program_number);

   if ((dbname = myargs[ARG_DB].strvalue) == NULL) {
      Int4 letters_read;
      FILE *infp2;
      Int4 s_from = 0, s_to = 0;
      char *subject_file = strdup(myargs[ARG_SUBJECT].strvalue);
      if ((infp2 = FileOpen(subject_file, "r")) == NULL) {
         ErrPostEx(SEV_FATAL, 1, 0, 
                   "blast: Unable to open second input file %s\n", 
                   subject_file);
         return (1);
      }
      sfree(subject_file);

      if (s_ParseIntervalLocationArgument(myargs[ARG_SUBJECT_LOC].strvalue,
                                          &s_from, &s_to)) {
         ErrPostEx(SEV_FATAL, 1, 0, "Invalid subject sequence location\n");
         return -1;
      }

      letters_read = 
         BLAST_GetQuerySeqLoc(infp2, db_is_na, 0, 0, s_from, s_to, NULL, 
                              &subject_slp, &ctr, NULL, FALSE,
                              myargs[ARG_DBGENCODE].intvalue);

      if (letters_read <= 0)
      {
           ErrPostEx(SEV_FATAL, 1, 0, "Unable to read subject sequence\n");
           return -1;
      }
      FileClose(infp2);
   }

   s_FillOptions(options);

   if ((infp = FileOpen(myargs[ARG_QUERY].strvalue, "r")) == NULL) {
      ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open input file %s\n", 
                myargs[ARG_QUERY].strvalue);
      return (1);
   }

   if (s_ParseIntervalLocationArgument(myargs[ARG_QUERY_LOC].strvalue,
                                       &q_from, &q_to)) {
         ErrPostEx(SEV_FATAL, 1, 0, "Invalid query sequence location\n");
         return -1;
   }

   tabular_output = myargs[ARG_TABULAR].intvalue;
   
   /* For on-the-fly tabular and for ASN.1 output, the believe query defline 
      option must be set to true, otherwise it should be false. */
   believe_query = 
       (myargs[ARG_FORMAT].intvalue == eAlignViewAsnText ||
        myargs[ARG_FORMAT].intvalue == eAlignViewAsnBinary);


   if (tabular_output) {
      if ((outfp = FileOpen(myargs[ARG_OUT].strvalue, "w")) == NULL) {
         ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", 
                   myargs[ARG_OUT].strvalue);
        return (1);
      }
   } else {
       BlastFormattingInfoNew(myargs[ARG_FORMAT].intvalue, options, 
                              blast_program, dbname, 
                              myargs[ARG_OUT].strvalue, &format_info);

       BlastFormattingInfoSetUpOptions(format_info, 
                                       myargs[ARG_DESCRIPTIONS].intvalue, 
                                       myargs[ARG_ALIGNMENTS].intvalue,
                                       (Boolean) myargs[ARG_HTML].intvalue,
                                       (Boolean) myargs[ARG_GREEDY].intvalue,
                                       (Boolean) myargs[ARG_SHOWGI].intvalue,
                                       believe_query);

       if (dbname)
           BLAST_PrintOutputHeader(format_info);
   }

   /* Get the query (queries), loop if necessary. */
   while (1) {
       SBlastSeqalignArray* seqalign_arr = NULL;
       SeqLoc* filter_loc=NULL;	/* All masking locations */
       SeqLoc* repeat_mask = NULL; /* Repeat mask locations */
       Int4 letters_read;

       if ((Boolean)myargs[ARG_LCASE].intvalue) {
           letters_read = 
               BLAST_GetQuerySeqLoc(infp, query_is_na, 
                   (Uint1)myargs[ARG_STRAND].intvalue, 10000, q_from, q_to, &lcase_mask,
                   &query_slp, &ctr, &num_queries, believe_query,
                   myargs[ARG_GENCODE].intvalue);
       } else {
           letters_read = 
               BLAST_GetQuerySeqLoc(infp, query_is_na,
                   (Uint1)myargs[ARG_STRAND].intvalue, 0, q_from, q_to, NULL, 
                   &query_slp, &ctr, &num_queries, believe_query,
                   myargs[ARG_GENCODE].intvalue);
       }

       /* If there is no sequence data left in the input file, break out of 
          the loop. */
       if (letters_read == 0)
           break;

       if (letters_read < 0) {
           ErrPostEx(SEV_FATAL, 1, 0, "Unable to read query sequence(s)\n");
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
           tf_data->show_gi = (Boolean) myargs[ARG_SHOWGI].intvalue;
           tf_data->show_accession = (Boolean) myargs[ARG_ACCESSION].intvalue;
       }

       options->num_cpus = myargs[ARG_THREADS].intvalue;

       /* Find repeat mask, if necessary */
       if ((status = Blast_FindRepeatFilterSeqLoc(query_slp, myargs[ARG_FILTER].strvalue,
                                &repeat_mask, &sum_returns->error)) != 0)
       {
            if (sum_returns && sum_returns->error && sum_returns->error->sev >= SEV_ERROR)
            {
                   ErrPostEx(SEV_ERROR, 1, 0, sum_returns->error->message);
                   return status;
            }
       }

       /* Combine repeat mask with lower case mask */
       if (repeat_mask)
           lcase_mask = ValNodeLink(&lcase_mask, repeat_mask);

       /* Do the main search */
       if (phi_blast) {
           status = PHIBlastRunSearch(query_slp, dbname, lcase_mask, options, 
                                      &phivnps, &filter_loc, sum_returns);
       } else if (dbname) {
           status = 
               Blast_DatabaseSearch(query_slp, dbname, lcase_mask, options, 
                                    tf_data, &seqalign_arr, &filter_loc, 
                                    sum_returns);
       } else {
           status = 
               Blast_TwoSeqLocSetsAdvanced(query_slp, subject_slp, lcase_mask, 
                                           options, tf_data, &seqalign_arr, 
                                           &filter_loc, sum_returns);
       }

       /* Deallocate the data structure used for tabular formatting. */
       BlastTabularFormatDataFree(tf_data);

       /* Free the lower case mask in SeqLoc form. */
       lcase_mask = Blast_ValNodeMaskListFree(lcase_mask);

       /* If masking was done for lookup table only, free the masking locations,
          because they will not be used for formatting. */
       if (SBlastOptionsGetMaskAtHash(options))
           filter_loc = Blast_ValNodeMaskListFree(filter_loc);

       /* Post warning or error messages, no matter what the search status was. */
       SBlastMessageErrPost(sum_returns->error);

       if (!status) {
           if (phi_blast) {
               status =
                   PHIBlastFormatResults(phivnps, query_slp, format_info, 
                                         sum_returns);
               phivnps = PHIBlastResultsFree(phivnps);
           } else if (!tabular_output) {
               
               /* Format the results */
               status = 
                   BLAST_FormatResults(seqalign_arr, num_queries, query_slp, 
                                       filter_loc, format_info, sum_returns);
               
               seqalign_arr = SBlastSeqalignArrayFree(seqalign_arr);
           }

       }
       /* Update the cumulative summary returns structure and clean the returns
          substructures for the current search iteration. */
       Blast_SummaryReturnUpdate(sum_returns, &full_sum_returns);
       Blast_SummaryReturnClean(sum_returns);
       filter_loc = Blast_ValNodeMaskListFree(filter_loc);
       query_slp = SeqLocSetFree(query_slp);
   } /* End loop on sets of queries */
   
   if (infp)
      FileClose(infp);
   
   subject_slp = SeqLocSetFree(subject_slp);

   /* Print the footer with cumulative summary information. */
   Blast_PrintOutputFooter(format_info, full_sum_returns);

   sum_returns = Blast_SummaryReturnFree(sum_returns);
   full_sum_returns = Blast_SummaryReturnFree(full_sum_returns);

   if (!tabular_output) 
       format_info = BlastFormattingInfoFree(format_info);       
   else
       FileClose(outfp);

   options = SBlastOptionsFree(options);

   return status;
}




