/* $Id: blast_api.c,v 1.33 2006/04/26 12:44:37 madden Exp $
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
*
* ===========================================================================
*
* Author:  Ilya Dondoshansky
*
* ===========================================================================
*/

/** @file blast_api.c
 * Functions for C toolkit applications to perform a BLAST search 
 * against a BLAST database, using the rewritten blast engine.
 */

#include <algo/blast/api/blast_api.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/blast_traceback.h>
#include <algo/blast/core/hspstream_collector.h>
#include <algo/blast/core/phi_lookup.h>
#include <algo/blast/api/hspstream_queue.h>
#include <algo/blast/api/blast_mtlock.h>
#include <algo/blast/api/blast_prelim.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/api/seqsrc_readdb.h>
#include <algo/blast/api/seqsrc_multiseq.h>
#include <algo/blast/api/blast_seqalign.h>
#include <algo/blast/api/dust_filter.h>
#include <algo/blast/api/blast_message_api.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

/** Initializes the auxiliary structure with RPS BLAST database information.
 * @param ppinfo Resulting structure. [out]
 * @param rps_mmap Memory mapped lookup table [out]
 * @param rps_pssm_mmap Memory mapped PSSM [out]
 * @param dbname Name of the database [in]
 */
static Int2 
s_BlastRPSInfoInit(BlastRPSInfo **ppinfo, Nlm_MemMap **rps_mmap,
                   Nlm_MemMap **rps_pssm_mmap, const char* dbname)
{
   char filename[PATH_MAX];
   char pathname[PATH_MAX];
   BlastRPSInfo *info;
   FILE *auxfile;
   Int4 i;
   Int4 seq_size;
   Int4 num_db_seqs;
   Nlm_MemMapPtr lut_mmap;
   Nlm_MemMapPtr pssm_mmap;
   char buffer[PATH_MAX];
   ReadDBFILEPtr rdfp;
   char *tmp_dbname;

   info = (BlastRPSInfo *)malloc(sizeof(BlastRPSInfo));
   if (info == NULL)
      ErrPostEx(SEV_FATAL, 1, 0, "Memory allocation failed");

   /* find the path to the RPS database */
   tmp_dbname = strdup(dbname);
   rdfp = readdb_new_ex2(tmp_dbname, READDB_DB_IS_PROT, 
                         READDB_NEW_DO_REPORT, NULL, NULL);
   sfree(tmp_dbname);
   if (rdfp == NULL)
      ErrPostEx(SEV_FATAL, 1, 0, "Cannot map RPS BLAST database");
   sprintf(pathname, "%s", rdfp->full_filename);
   rdfp = readdb_destruct(rdfp);

   sprintf(filename, "%s.loo", (char *)pathname);
   lut_mmap = Nlm_MemMapInit(filename);
   if (lut_mmap == NULL)
      ErrPostEx(SEV_FATAL, 1, 0, "Cannot map RPS BLAST lookup file");

   info->lookup_header = (BlastRPSLookupFileHeader *)lut_mmap->mmp_begin;
   if (info->lookup_header->magic_number != RPS_MAGIC_NUM) {
       if (Nlm_SwitchUint4(info->lookup_header->magic_number)==RPS_MAGIC_NUM) {
           ErrPostEx(SEV_FATAL, 1, 0, "RPS BLAST lookup file was created "
                           "on an incompatible platform");
       }
       else {
           ErrPostEx(SEV_FATAL, 1, 0, "RPS BLAST lookup file is corrupt");
       }
   }

   sprintf(filename, "%s.rps", (char *)pathname);
   pssm_mmap = Nlm_MemMapInit(filename);
   if (pssm_mmap == NULL)
      ErrPostEx(SEV_FATAL, 1, 0, "Cannot map RPS BLAST profile file");

   info->profile_header = (BlastRPSProfileHeader *)pssm_mmap->mmp_begin;
   if (info->profile_header->magic_number != RPS_MAGIC_NUM) {
       if (Nlm_SwitchUint4(info->profile_header->magic_number)==RPS_MAGIC_NUM) {
           ErrPostEx(SEV_FATAL, 1, 0, "RPS BLAST profile file was created "
                           "on an incompatible platform");
       }
       else {
           ErrPostEx(SEV_FATAL, 1, 0, "RPS BLAST profile file is corrupt");
       }
   }

   num_db_seqs = info->profile_header->num_profiles;

   sprintf(filename, "%s.aux", (char *)pathname);
   auxfile = FileOpen(filename, "r");
   if (auxfile == NULL)
      ErrPostEx(SEV_FATAL, 1, 0,"Cannot open RPS BLAST parameters file");

   fscanf(auxfile, "%s", buffer);
   info->aux_info.orig_score_matrix = strdup(buffer);
   fscanf(auxfile, "%d", &info->aux_info.gap_open_penalty);
   fscanf(auxfile, "%d", &info->aux_info.gap_extend_penalty);
   fscanf(auxfile, "%le", &info->aux_info.ungapped_k);
   fscanf(auxfile, "%le", &info->aux_info.ungapped_h);
   fscanf(auxfile, "%d", &info->aux_info.max_db_seq_length);
   fscanf(auxfile, "%d", &info->aux_info.db_length);
   fscanf(auxfile, "%lf", &info->aux_info.scale_factor);

   info->aux_info.karlin_k = (double *)malloc(num_db_seqs * sizeof(double));
   for (i = 0; i < num_db_seqs && !feof(auxfile); i++) {
      fscanf(auxfile, "%d", &seq_size); /* not used */
      fscanf(auxfile, "%le", &info->aux_info.karlin_k[i]);
   }

   if (i < num_db_seqs)
      ErrPostEx(SEV_FATAL, 1, 0, "Missing Karlin parameters");

   FileClose(auxfile);
   *ppinfo = info;
   *rps_mmap = lut_mmap;
   *rps_pssm_mmap = pssm_mmap;
   return 0;
}

/** Initializes and populates the RPS BLAST specific structures. 
 * @param seq_src Database sequences source [in]
 * @param options All search options [in]
 * @param rps_options Copy of options, with RPS-specific modifications for 
 *                    scoring and hit saving options. All other options pointers
 *                    are left the same as in input "options". [out]
 * @param rps_info_out Auxiliary structure with RPS-specific information [out]
 * @param rps_mmap Memory mapped lookup table [out]
 * @param rps_pssm_mmap Memory mapped PSSM [out]
 * @param scale_factor Scaling factor for RPS matrix. [out]
 * @param extra_returns Structure containing error information [in] [out]
 * @return Status.
 */
static Int2 
s_RPSExtraStructsSetUp(const BlastSeqSrc* seq_src, const SBlastOptions* options,
                       SBlastOptions* *rps_options, BlastRPSInfo* *rps_info_out,
                       Nlm_MemMapPtr *rps_mmap, Nlm_MemMapPtr *rps_pssm_mmap,
                       double *scale_factor, Blast_SummaryReturn* extra_returns)
{
    const char* kDbName;
    BlastRPSInfo* rps_info = NULL;
    BlastScoringOptions* rps_score_options;
    BlastHitSavingOptions* rps_hit_options;
    Int2 status = 0;

    /* The caller has already checked these, so we are just asserting it here. */
    ASSERT(seq_src && options && extra_returns);

    kDbName = BlastSeqSrcGetName(seq_src);

    if (kDbName == NULL ||
        (status = s_BlastRPSInfoInit(&rps_info, rps_mmap, 
                                     rps_pssm_mmap, kDbName)) != 0) {
        SBlastMessageWrite(&extra_returns->error, SEV_WARNING, "RPS BLAST database setup failed", NULL, FALSE);
        return status;
    }
    *rps_info_out = rps_info;
    *scale_factor = rps_info->aux_info.scale_factor;
    rps_score_options = (BlastScoringOptions*)
        BlastMemDup(options->score_options, sizeof(BlastScoringOptions));
    rps_hit_options = (BlastHitSavingOptions*)
        BlastMemDup(options->hit_options, sizeof(BlastHitSavingOptions));
    rps_score_options->gap_open = 
        rps_info->aux_info.gap_open_penalty;
    rps_score_options->gap_extend = 
        rps_info->aux_info.gap_extend_penalty;
    rps_score_options->matrix = 
        strdup(rps_info->aux_info.orig_score_matrix);

    *rps_options = (SBlastOptions*) BlastMemDup(options, sizeof(SBlastOptions));
    (*rps_options)->score_options = rps_score_options;
    (*rps_options)->hit_options = rps_hit_options;

    return 0;
}

/** Frees the RPS BLAST specific extra structures. 
 * @param rps_info Auxiliary structure with RPS-specific information [in]
 * @param rps_mmap Memory mapped lookup table [in]
 * @param rps_pssm_mmap Memory mapped PSSM [in]
 * @param options Copy of the options wrapper structure, containing scoring 
 *                and hit saving options, specially modified for RPS search.
 *                All other options are the same as in the original structure,
 *                so they should not be freed here. [in]
 */
static void
s_RPSExtraStructsFree(BlastRPSInfo* rps_info, Nlm_MemMapPtr rps_mmap,
                      Nlm_MemMapPtr rps_pssm_mmap, SBlastOptions* options)
{
    Nlm_MemMapFini(rps_mmap);
    Nlm_MemMapFini(rps_pssm_mmap);

    if (rps_info) {
        sfree(rps_info->aux_info.karlin_k);
        sfree(rps_info->aux_info.orig_score_matrix);
        sfree(rps_info);
    }
    if (options) {
        if (options->score_options) {
            sfree(options->score_options->matrix);
            sfree(options->score_options);
        }
        sfree(options->hit_options);
        sfree(options);
    }
}

/** Sets up the HSP stream, depending on whether the search is single or
 * multithreaded, and whether on-the-fly tabular output option is set. 
 */
static Int2
s_BlastHSPStreamSetUp(BLAST_SequenceBlk* query, BlastQueryInfo* query_info, 
                      const BlastSeqSrc* seq_src, const SBlastOptions* options, 
                      BlastScoreBlk* sbp, BlastTabularFormatData* tf_data, 
                      BlastHSPStream* *hsp_stream,
                      Blast_SummaryReturn* extra_returns)
{
    Int2 status = 0;

    /* If any of the required inputs were NULL, the caller would have exited 
       before getting to this point. ASSERT this here. */
    ASSERT(query && query_info && seq_src && options && sbp);

    if (!tf_data) {
        const Int4 kNumResults = query_info->num_queries;
        /* Results in the collector stream should be sorted only for a
           database search. The latter is true if and only if the sequence
           source has non-zero database length. */
        const Boolean kSortOnRead = (BlastSeqSrcGetTotLen(seq_src) != 0);
        SBlastHitsParameters* blasthit_params=NULL;
        MT_LOCK lock = NULL;
        if (options->num_cpus > 1)
            lock = Blast_MT_LOCKInit();
        
        SBlastHitsParametersNew(options->hit_options, options->ext_options,
                                options->score_options, &blasthit_params);
        *hsp_stream =
            Blast_HSPListCollectorInitMT(options->program, blasthit_params,
                                         kNumResults, kSortOnRead, lock);
    } else {
        /* Initialize the queue HSP stream for tabular formatting. */
        *hsp_stream = Blast_HSPListQueueInit();
        if ((status = Blast_TabularFormatDataSetUp(tf_data, options->program,
                          *hsp_stream, seq_src, query, query_info,
                          options->score_options, sbp, options->eff_len_options,
                          options->ext_options, options->hit_options, 
                          options->db_options)) != 0) {
            SBlastMessageWrite(&extra_returns->error, SEV_WARNING,
                "Failed to set up tabular formatting data structure", NULL, FALSE);
            return status;
        }
    }
    return status;
}

/** Starts and joins all threads performing a multi-threaded search, with or 
 * without on-the-fly output, or performs a single-threaded search.
 */
static Int2
s_BlastThreadManager(BLAST_SequenceBlk* query, BlastQueryInfo* query_info, 
                     const BlastSeqSrc* seq_src, const SBlastOptions* options, 
                     LookupTableWrap* lookup_wrap, BlastScoreBlk* sbp, 
                     BlastHSPStream* hsp_stream, BlastRPSInfo* rps_info, 
                     BlastTabularFormatData* tf_data, BlastHSPResults **results,
                     Blast_SummaryReturn* extra_returns)
{
    Int2 status = 0;
    /* The options input cannot be NULL here. The program would have exited 
       before entering this function if it was. */
    const BlastInitialWordOptions* word_options = options->word_options;
    const BlastScoringOptions* score_options = options->score_options;
    const BlastExtensionOptions* ext_options = options->ext_options;
    const BlastHitSavingOptions* hit_options = options->hit_options;
    const BlastEffectiveLengthsOptions* eff_len_options = 
        options->eff_len_options;
    const PSIBlastOptions* psi_options = options->psi_options;
    const BlastDatabaseOptions* db_options = options->db_options;
    TNlmThread format_thread = NULL;
    BlastDiagnostics* diagnostics = NULL;
    const EBlastProgramType kProgram = options->program;
    const int kNumCpus = options->num_cpus;

    /* Assert that all required inputs are not NULL. They must be - otherwise 
       the program should have exited before entering this function. */
    ASSERT(query && query_info && seq_src && lookup_wrap && sbp && 
           hsp_stream && extra_returns);

    /* Start the formatting thread */
    if(tf_data && NlmThreadsAvailable() &&
       (format_thread =
        NlmThreadCreate(Blast_TabularFormatThread, (void*) tf_data))
       == NULL_thread) {
        SBlastMessageWrite(&extra_returns->error, SEV_WARNING,
                           "Cannot create thread for formatting tabular output\n", NULL, options->believe_query);
        return -1;
    }

    if (NlmThreadsAvailable() && kNumCpus > 1) {
        TNlmThread* thread_array =
            (TNlmThread*) calloc(kNumCpus, sizeof(TNlmThread));
        BlastPrelimSearchThreadData* search_data = NULL;
        void* join_status = NULL;
        int index;
        
        diagnostics = Blast_DiagnosticsInitMT(Blast_MT_LOCKInit());

        for (index = 0; index < kNumCpus; index++) {
            search_data = 
                BlastPrelimSearchThreadDataInit(kProgram, query, 
                    query_info, seq_src, lookup_wrap, score_options, 
                    word_options, ext_options, hit_options, eff_len_options, 
                    psi_options, db_options, sbp, diagnostics, hsp_stream);

            thread_array[index] =
               NlmThreadCreate(Blast_PrelimSearchThreadRun, 
                               (void*) search_data);
        }
        for (index = 0; index < kNumCpus; index++)
            NlmThreadJoin(thread_array[index], &join_status);
  
        MemFree(thread_array);
        
        if (!tf_data) {
            SPHIPatternSearchBlk* pattern_blk = NULL;
            if (Blast_ProgramIsPhiBlast(kProgram)) {
                PHIPatternSpaceCalc(query_info, diagnostics);
                pattern_blk = (SPHIPatternSearchBlk*) lookup_wrap->lut;
            }

            if ((status = Blast_RunTracebackSearch(kProgram, query, 
                             query_info, seq_src, score_options, 
                             ext_options, hit_options, eff_len_options, 
                             db_options, psi_options, sbp, hsp_stream, 
                             rps_info, pattern_blk, results)) != 0) {
                SBlastMessageWrite(&extra_returns->error, SEV_ERROR,
                                   "Traceback engine failed\n", NULL, options->believe_query);
            }
        }
    } else {
        diagnostics = Blast_DiagnosticsInit();

        if (tf_data) { /* Single thread, tabular */
            if ((status = 
                 Blast_RunPreliminarySearch(kProgram, query, query_info, 
                     seq_src, score_options, sbp, lookup_wrap, word_options, 
                     ext_options, hit_options, eff_len_options, psi_options, 
                     db_options, hsp_stream, diagnostics)) != 0) {
                SBlastMessageWrite(&extra_returns->error, SEV_ERROR,
                                   "Preliminary search engine failed\n", NULL, options->believe_query);
            }
        } else { /* Single thread, non-tabular */
            if ((status=Blast_RunFullSearch(kProgram, query, query_info, 
                            seq_src, sbp, score_options, lookup_wrap, 
                            word_options, ext_options, hit_options, 
                            eff_len_options, psi_options, db_options, hsp_stream,
                            rps_info, diagnostics, results, 0, 0)) != 0) {
                SBlastMessageWrite(&extra_returns->error, SEV_ERROR,  
                                    "Blast_RunFullSearch failed\n", NULL, options->believe_query);
            }
        }
    }

    if (tf_data) {
        void* join_status = NULL;
        BlastHSPStreamClose(hsp_stream);
        NlmThreadJoin(format_thread, &join_status);
        /* Free the internally allocated structures used for tabular
           formatting. */
        BlastTabularFormatDataClean(tf_data);
    }

    hsp_stream = BlastHSPStreamFree(hsp_stream);
    Blast_SummaryReturnFill(kProgram, score_options, sbp, options->lookup_options, 
                            word_options, ext_options, hit_options,
                            eff_len_options, options->query_options, query_info, 
                            seq_src, &diagnostics, extra_returns);

    return status;
}

Int2 
Blast_RunSearch(SeqLoc* query_seqloc, 
                 const BlastSeqSrc* seq_src,
                 SeqLoc* masking_locs,
                 const SBlastOptions* options,
                 BlastTabularFormatData* tf_data,
                 BlastHSPResults **results,
                 SeqLoc** filter_out,
                 Blast_SummaryReturn* extra_returns)
{
    Int2 status = 0;
    BLAST_SequenceBlk *query = NULL;
    BlastQueryInfo* query_info = NULL;
    double scale_factor = 1.0;
    BlastSeqLoc* lookup_segments = NULL;
    BlastScoreBlk* sbp = NULL;
    LookupTableWrap* lookup_wrap = NULL;
    BlastMaskLoc* mask_loc = NULL;
    BlastHSPStream* hsp_stream = NULL;
    const EBlastProgramType kProgram = options->program;
    const Boolean kRpsBlast = 
        (kProgram == eBlastTypeRpsBlast ||
         kProgram == eBlastTypeRpsTblastn);
    BlastRPSInfo* rps_info = NULL;
    Nlm_MemMapPtr rps_mmap = NULL;
    Nlm_MemMapPtr rps_pssm_mmap = NULL;
    const QuerySetUpOptions* query_options = options->query_options;
    const LookupTableOptions* lookup_options = options->lookup_options;
    const BlastScoringOptions* score_options = options->score_options;
    const BlastHitSavingOptions* hit_options = options->hit_options;
    SBlastOptions* rps_options = NULL;
    const Boolean kPhiBlast = Blast_ProgramIsPhiBlast(kProgram);
    const Uint1 kDeallocateMe = 253;
    Blast_Message *core_msg = NULL;

    if (!query_seqloc || !seq_src || !options || !extra_returns) 
        return -1;

    if ((status = 
         BLAST_ValidateOptions(kProgram, options->ext_options, score_options, 
                               lookup_options, options->word_options, hit_options, 
                               &core_msg)) != 0) {
         extra_returns->error = Blast_MessageToSBlastMessage(core_msg, NULL, NULL, options->believe_query);
         core_msg = Blast_MessageFree(core_msg);
         
        return status;
    }

    if (options->program == eBlastTypeBlastn)
    {
         SeqLoc* dust_mask = NULL; /* Dust mask locations */
         Blast_FindDustSeqLoc(query_seqloc, options, &dust_mask);
         /* Combine dust mask with lower case mask 
            The dust mask will be deallocated by the end of this function
            though as it's copied in BLAST_MainSetUp 
            Not deallocating it will result in a memory leak if masking_locs
            was NULL at the start of this function */
         if (dust_mask)
         {
            SeqLoc* dust_mask_var = dust_mask;
            while (dust_mask_var)
            {
               dust_mask_var->choice = kDeallocateMe;
               dust_mask_var = dust_mask_var->next;
            }
            ValNodeLink(&masking_locs, dust_mask);
         }
    }

    if (kRpsBlast) {
        if ((status = 
             s_RPSExtraStructsSetUp(seq_src, options, &rps_options, &rps_info, 
                                    &rps_mmap, &rps_pssm_mmap, &scale_factor, 
                                    extra_returns)))
            return status;
        score_options = rps_options->score_options;
        hit_options = rps_options->hit_options;
        options = rps_options; /* This will not change the caller's pointer. */
    }

    if ((status = BLAST_SetUpQuery(kProgram, query_seqloc, query_options, 
                                   masking_locs, &query_info, &query))) {
        SBlastMessageWrite(&extra_returns->error, SEV_ERROR,  
                "BLAST_SetUpQuery returned non-zero status\n", NULL, FALSE);
        return status;
    }

    status = 
        BLAST_MainSetUp(kProgram, query_options, score_options, query, 
                        query_info, scale_factor, &lookup_segments, &mask_loc,
                        &sbp, &core_msg);
    if (core_msg)
    {
       extra_returns->error = Blast_MessageToSBlastMessage(core_msg, query_seqloc, query_info, options->believe_query);
       core_msg = Blast_MessageFree(core_msg);
    }

    if (status)
        return status;

    if (filter_out) {
        *filter_out = 
            BlastMaskLocToSeqLoc(kProgram, mask_loc, query_seqloc);
    }

    /* Mask locations in BlastMaskLoc form are no longer needed. */
    BlastMaskLocFree(mask_loc);

    if (masking_locs)
    {
          SeqLocPtr slp_var = masking_locs;
          SeqLocPtr last = NULL;
          while (slp_var)
          {
              if (slp_var->choice == kDeallocateMe)
              {
                  if (last == NULL)
                  {
                     masking_locs = slp_var->next;
                     slp_var->next = NULL;
                     Blast_ValNodeMaskListFree(slp_var);
                     slp_var = masking_locs;
                  }
                  else
                  {
                     last->next = slp_var->next;
                     slp_var->next = NULL;
                     Blast_ValNodeMaskListFree(slp_var);
                     slp_var = last->next;
                  }
              } 
              else
              {
                  last = slp_var;
                  slp_var = slp_var->next;
              }
          }
    }

    if ((status = LookupTableWrapInit(query, lookup_options, 
                        lookup_segments, sbp, &lookup_wrap, rps_info)))
        return status;

    /* For PHI BLAST, save information about pattern occurrences in
       query in the BlastQueryInfo structure. */
    if (kPhiBlast) {
        SPHIPatternSearchBlk* pattern_blk = 
            (SPHIPatternSearchBlk*) lookup_wrap->lut;
        Blast_SetPHIPatternInfo(kProgram, pattern_blk, query, lookup_segments, 
                                query_info);
    }
    /* Only need for the setup of lookup table. */
    lookup_segments = BlastSeqLocFree(lookup_segments);

    if ((status = s_BlastHSPStreamSetUp(query, query_info, seq_src, options, sbp,
                                        tf_data, &hsp_stream, extra_returns)))
        return status;

    if ((status = s_BlastThreadManager(query, query_info, seq_src, options,
                                       lookup_wrap, sbp, hsp_stream, rps_info, 
                                       tf_data, results, extra_returns)))
        return status;
    
    lookup_wrap = LookupTableWrapFree(lookup_wrap);
    
    query = BlastSequenceBlkFree(query);
    query_info = BlastQueryInfoFree(query_info);
    BlastScoreBlkFree(sbp);
    
    if (kRpsBlast)
        s_RPSExtraStructsFree(rps_info, rps_mmap, rps_pssm_mmap, rps_options);
    
    return status;
}

Int2 
Blast_DatabaseSearch(SeqLoc* query_seqloc, char* db_name,
                     SeqLoc* masking_locs,
                     const SBlastOptions* options,
                     BlastTabularFormatData* tf_data,
                     SBlastSeqalignArray* *seqalign_arr,
                     SeqLoc** filter_out,
                     Blast_SummaryReturn* extra_returns)
{
    BlastSeqSrc *seq_src = NULL;
    Boolean db_is_prot;
    Int2 status = 0;
    BlastHSPResults* results = NULL;
    ReadDBFILE* rdfp = NULL;

    if (!options || !query_seqloc || !db_name || !extra_returns)
        return -1;

    db_is_prot = 
        (options->program == eBlastTypeBlastp   ||
         options->program == eBlastTypeBlastx   ||
         options->program == eBlastTypeRpsBlast ||
         options->program == eBlastTypeRpsTblastn);

    rdfp = readdb_new(db_name, db_is_prot);

    seq_src = ReaddbBlastSeqSrcAttach(rdfp);

    if (seq_src == NULL) {
        SBlastMessageWrite(&extra_returns->error, SEV_WARNING,
                           "Initialization of subject sequences source failed", NULL, options->believe_query);
    } else {
        char* error_str = BlastSeqSrcGetInitError(seq_src);
        if (error_str)
            SBlastMessageWrite(&extra_returns->error, SEV_WARNING, error_str, NULL, options->believe_query); 
    }

    /* If there was an error initializing the sequence source, return without 
       doing the search. */
    if (extra_returns->error)
        return -1;

    status = 
        Blast_RunSearch(query_seqloc, seq_src, masking_locs, options, tf_data,
                         &results, filter_out, extra_returns);

    /* The ReadDBFILE structure will not be destroyed here, because the 
       initialising function used readdb_attach */
    BlastSeqSrcFree(seq_src);

    if (!status && !tf_data) {
        status = 
            BLAST_ResultsToSeqAlign(options->program, &results, 
                                    query_seqloc, rdfp, NULL, 
                                    options->score_options->gapped_calculation,
                                    options->score_options->is_ooframe, 
                                    seqalign_arr);
    }

    readdb_destruct(rdfp);

    if (status)
        return status;

    return status;
}

/** Splits the PHI BLAST results corresponding to different pattern occurrences
 * in query, converts them to Seq-aligns and puts in a list of ValNodes.
 * @param results All results from different pattern occurrences 
 *                mixed together. On return points to NULL. [in]
 * @param pattern_info Query pattern occurrences information [in]
 * @param program Program type (phiblastp or phiblastn) [in]
 * @param query_seqloc List of query locations [in]
 * @param rdfp blast db object [in]
 * @param phivnps List of ValNodes containing Seq-aligns. [out]
 * @return Status, 0 on success, -1 on failure.
 */
static Int2
s_PHIResultsToSeqAlign(const BlastHSPResults* results, 
                       const SPHIQueryInfo* pattern_info,
                       EBlastProgramType program, SeqLoc* query_seqloc, 
                       ReadDBFILE* rdfp, ValNode* *phivnps)
{
    Int2 status = 0;
    /* Split results into an array of BlastHSPResults structures corresponding
       to different pattern occurrences. */
    BlastHSPResults* *phi_results = 
        PHIBlast_HSPResultsSplit(results, pattern_info);

    if (phi_results) {
        int pattern_index; /* Index over pattern occurrences. */

        for (pattern_index = 0; pattern_index < pattern_info->num_patterns;
             ++pattern_index) {
            SBlastSeqalignArray* seqalign_arr = NULL;
            SeqAlign* seqalign = NULL;
            BlastHSPResults* one_phi_results = phi_results[pattern_index];

            if (one_phi_results) {
                /* PHI BLAST is always gapped, and never out-of-frame, hence
                 * TRUE and FALSE values for the respective booleans in the next
                 * call.
                 */
                status =
                    BLAST_ResultsToSeqAlign(program, &one_phi_results, 
                                            query_seqloc, rdfp, NULL, TRUE, 
                                            FALSE, &seqalign_arr);
                if (seqalign_arr)
                {
                    seqalign = seqalign_arr->array[0];
                    seqalign_arr->array[0] = NULL;
                    SBlastSeqalignArrayFree(seqalign_arr);
                }
                ValNodeAddPointer(phivnps, pattern_index, seqalign);
            }
        }
        sfree(phi_results);
    }
    return status;
}

Int2
PHIBlastRunSearch(SeqLoc* query_seqloc, char* db_name, SeqLoc* masking_locs,
                  const SBlastOptions* options, ValNode* *phivnps,
                  SeqLoc** filter_out, Blast_SummaryReturn* extra_returns)
{
    BlastSeqSrc *seq_src = NULL;
    Boolean is_prot;
    Int2 status = 0;
    BlastHSPResults* results = NULL;
    ReadDBFILE* rdfp = NULL;

    if (!options || !query_seqloc || !db_name || !extra_returns || !phivnps)
        return -1;

    ASSERT(Blast_ProgramIsPhiBlast(options->program));

    is_prot = (options->program == eBlastTypePhiBlastp);

    rdfp = readdb_new(db_name, is_prot);

    seq_src = ReaddbBlastSeqSrcAttach(rdfp);

    if (seq_src == NULL) {
        SBlastMessageWrite(&extra_returns->error, SEV_WARNING,
                           "Initialization of subject sequences source failed", NULL, options->believe_query);
    } else {
        char* error_str = BlastSeqSrcGetInitError(seq_src);
        if (error_str)
            SBlastMessageWrite(&extra_returns->error, SEV_WARNING, error_str, NULL, options->believe_query);
    }

    /* If there was an error initializing the sequence source, return without
       doing the search. */
    if (extra_returns->error)
        return -1;

    /* Masking at hash and on-the-fly tabular output are not applicable for 
       PHI BLAST, so pass NULL in corresponding arguments. */
    status =
        Blast_RunSearch(query_seqloc, seq_src, masking_locs, options, NULL,
                        &results, filter_out, extra_returns);

    /* The ReadDBFILE structure will not be destroyed here, because the
       initialising function used readdb_attach */
    BlastSeqSrcFree(seq_src);

    *phivnps = NULL;

    if (!status) {
        status = 
            s_PHIResultsToSeqAlign(results, extra_returns->pattern_info, 
                                   options->program, query_seqloc, rdfp,
                                   phivnps);
    }

    results = Blast_HSPResultsFree(results);

    readdb_destruct(rdfp);
    return status;
}

Int2 
Blast_TwoSeqLocSetsAdvanced(SeqLoc* query_seqloc, 
                            SeqLoc* subject_seqloc,
                            SeqLoc* masking_locs,
                            const SBlastOptions* options,
                            BlastTabularFormatData* tf_data,
                            SBlastSeqalignArray* *seqalign_arr,
                            SeqLoc** filter_out,
                            Blast_SummaryReturn* extra_returns)
{
    BlastSeqSrc *seq_src = NULL;
    Int2 status = 0;
    BlastHSPResults* results = NULL;


    if (!options || !query_seqloc || !subject_seqloc || !extra_returns)
        return -1;


    seq_src = MultiSeqBlastSeqSrcInit(subject_seqloc, options->program);

    if (seq_src == NULL) {
        SBlastMessageWrite(&extra_returns->error, SEV_WARNING,
                           "Initialization of subject sequences source failed", NULL, options->believe_query);
    } else {
        char* error_str = BlastSeqSrcGetInitError(seq_src);
        if (error_str)
            SBlastMessageWrite(&extra_returns->error, SEV_WARNING, error_str, NULL, options->believe_query); 
    }

    /* If there was an error initializing the sequence source, return without 
       doing the search. */
    if (extra_returns->error)
        return -1;

    status = 
        Blast_RunSearch(query_seqloc, seq_src, masking_locs, options, tf_data, 
                         &results, filter_out, extra_returns);

    /* The ReadDBFILE structure will not be destroyed here, because the 
       initialising function used readdb_attach */
    BlastSeqSrcFree(seq_src);

    if (!status) {
        status = 
            BLAST_ResultsToSeqAlign(options->program, &results, query_seqloc, 
                                    NULL, subject_seqloc, 
                                    options->score_options->gapped_calculation,
                                    options->score_options->is_ooframe, 
                                    seqalign_arr);
    }

    if (status)
        return status;

    return status;
}

/* @} */



