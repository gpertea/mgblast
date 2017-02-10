static char const rcsid[] = "$Id: seg.c,v 6.19 2003/05/30 17:25:38 coulouri Exp $";

/* 
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


/*--------------------------------------------------------------------------*/

#include "seg.h"
#include "lnfac.h"

static char _this_module[] = "seg";
#undef THIS_MODULE
#define THIS_MODULE _this_module

/*------------------------------------------------------(local functions)---*/

SequencePtr seq_phase(SequencePtr seq, Int4 phase, Int4 period);
SegPtr per_mergesegs(SequencePtr seq, SegPtr *persegs);
FloatHiPtr seqent(SequencePtr seq, Int4 window, Int4 maxbogus);

SequencePtr openwin(SequencePtr seq, Int4 start, Int4 length);
Int4 shiftwin1(SequencePtr win);
void closewin(SequencePtr win);
void compon(SequencePtr win);
void stateon(SequencePtr win);
void enton(SequencePtr win);
FloatHi entropy(Int4Ptr sv);
static int LIBCALLBACK state_cmp(VoidPtr s1, VoidPtr s2);
Int4 findlo(Int4 i, Int4 limit, FloatHi hicut, FloatHiPtr H);
Int4 findhi(Int4 i, Int4 limit, FloatHi hicut, FloatHiPtr H);
void trim(SequencePtr seq, Int4Ptr leftend, Int4Ptr rightend, 
          SegParamsPtr sparamsp);
FloatHi getprob(Int4Ptr sv, Int4 total, AlphaPtr palpha);
FloatHi lnperm(Int4Ptr sv, Int4 tot);
FloatHi lnass(Int4Ptr sv, Int4 alphasize);
void mergesegs(SequencePtr seq, SegPtr segs, Boolean overlap);
void decrementsv(Int4Ptr sv, Int4 class);
void incrementsv(Int4Ptr sv, Int4 class);
static void appendseg(SegPtr segs, SegPtr seg);
static Boolean hasdash(SequencePtr win);
AlphaPtr AA20alpha (void);
AlphaPtr NA4alpha (void);
void AlphaFree(AlphaPtr palpha);
AlphaPtr AlphaCopy(AlphaPtr palpha);

SeqLocPtr SegsToSeqLoc(BioseqPtr bsp, SegPtr segs);
SeqLocPtr SeqlocSegsToSeqLoc(SeqLocPtr slp, SegPtr segs);

/*------------------------------------------------------------(BioseqSeg)---*/

SeqLocPtr BioseqSeg (BioseqPtr bsp, SegParamsPtr sparamsp)

  {
   SeqLocPtr	slp = NULL;

/* error msg stuff */

   ErrSetOptFlags (EO_MSG_CODES);

/* bail on null bioseq */

   if (!bsp)
     {
       ErrPostEx (SEV_ERROR, 1, 1, "no bioseq");
       ErrShow ();
       return (slp);
     }

   if (ISA_na (bsp->mol))
     {
      slp = BioseqSegNa (bsp, sparamsp);
     }

   if (ISA_aa (bsp->mol))
     {
      slp = BioseqSegAa (bsp, sparamsp);
     }

/* clean up & return */
   return (slp);
  }

/*----------------------------------------------------------(BioseqSegNa)---*/

SeqLocPtr BioseqSegNa (BioseqPtr bsp, SegParamsPtr sparamsp)

  {
   SeqPortPtr   spp=NULL;
   SeqLocPtr	slp = NULL;
   SequencePtr seqwin;
   SegPtr segs;
   Int4 len, index;
   CharPtr seq;
   Uint1       residue;

/* error msg stuff */

   ErrSetOptFlags (EO_MSG_CODES);
   SegParamsCheck (sparamsp);

/* bail on null bioseq */

   if (!bsp)
     {
       ErrPostEx (SEV_ERROR, 2, 1, "no bioseq");
       ErrShow ();
       return (slp);
     }

   if (!ISA_na (bsp->mol))
     {
      ErrPostEx (SEV_WARNING, 2, 2, "not nucleic acid sequence");
      ErrShow ();
      return (slp);
     }

   if (!sparamsp)
     {
      sparamsp = SegParamsNewNa();
      if (!sparamsp)
        {
         ErrPostEx (SEV_WARNING, 0, 0, "null parameters object");
         ErrShow();
         return(slp);
        }
     }
   SegParamsCheck (sparamsp);

/* make an old-style genwin sequence window object */

   len = bsp->length;
   seq = (CharPtr) MemNew((len+1)*sizeof(Char));
   spp = SeqPortNew(bsp, 0, -1, bsp->strand, Seq_code_ncbi8na);
   if (spp == NULL) {
      ErrPostEx (SEV_ERROR, 0, 0, "SeqPortNew failure");
      ErrShow();
   }

   if (!seq)
     {
      ErrPostEx (SEV_ERROR, 0, 0, "memory allocation failure");
      ErrShow();
      return(slp);
     }

   index = 0;
   while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
     {
      if (IS_residue(residue))
        {
         seq[index] = residue;
         index++;
        }
      else if (residue == SEQPORT_EOS)
        {
         continue; /*[Segment boundary]*/
        }
      else if (residue == SEQPORT_VIRT)
        {
         continue; /*[Virtual Sequence]*/
        } 
     }

   seq[index] = NULLB;

   seqwin = SeqNew();
   seqwin->seq = (CharPtr) seq;
   seqwin->length = bsp->length;
   seqwin->palpha = AlphaCopy(sparamsp->palpha);
   
/* seg the sequence */

   segs = (SegPtr) NULL;
   SegSeq (seqwin, sparamsp, &segs, 0);

/* merge the segment if desired. */

   if (sparamsp->overlaps)
	mergesegs(seqwin, segs, sparamsp->overlaps);

/* convert segs to seqlocs */

   slp = SegsToSeqLoc(bsp, segs);   

/* clean up & return */

   SeqFree (seqwin);
   SegFree (segs);
   SeqPortFree (spp);
   return (slp);
  }

/*----------------------------------------------------------(BioseqSegAa)---*/

SeqLocPtr BioseqSegAa (BioseqPtr bsp, SegParamsPtr sparamsp)

  {
   SeqPortPtr	spp=NULL;
   SeqLocPtr	slp = NULL;
   SequencePtr seqwin;
   SegPtr segs;
   Int4 index, len;
   CharPtr seq;
   Uint1       residue;
   Boolean params_allocated = FALSE;

/* error msg stuff */

   ErrSetOptFlags (EO_MSG_CODES);
   SegParamsCheck (sparamsp);

/* bail on null bioseq */

   if (!bsp)
     {
       ErrPostEx (SEV_ERROR, 0, 0, "no bioseq");
       ErrShow ();
       return (slp);
     }

   if (!ISA_aa (bsp->mol))
     {
      ErrPostEx (SEV_WARNING, 0, 0, "not protein sequence");
      ErrShow ();
      return (slp);
     }

   if (!sparamsp)
     {
      sparamsp = SegParamsNewAa();
      SegParamsCheck (sparamsp);
      params_allocated = TRUE;
      if (!sparamsp)
        {
         ErrPostEx (SEV_WARNING, 0, 0, "null parameters object");
         ErrShow();
         return(slp);
        }
     }

/* make an old-style genwin sequence window object */

   len = bsp->length;
   seq = (CharPtr) MemNew((len+1)*sizeof(Char));
   spp = SeqPortNew(bsp, 0, -1, bsp->strand, Seq_code_ncbieaa);
   if (spp == NULL) {
      ErrPostEx (SEV_ERROR, 0, 0, "SeqPortNew failure");
      ErrShow();
   }

   if (!seq)
   {
      ErrPostEx (SEV_ERROR, 0, 0, "memory allocation failure");
      ErrShow();
	  SeqPortFree(spp);
      return(slp);
   }

   index = 0;
   while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
     {
      if (IS_residue(residue))
        {
         seq[index] = residue;
         index++;
        }
      else if (residue == SEQPORT_EOS)
        {
         continue; /*[Segment boundary]*/
        }
      else if (residue == SEQPORT_VIRT)
        {
         continue; /*[Virtual Sequence]*/
        } 
     }

   seq[index] = NULLB;


   seqwin = SeqNew();
   seqwin->seq = (CharPtr) seq;
   seqwin->length = bsp->length;
   seqwin->palpha = AlphaCopy(sparamsp->palpha);
   
/* seg the sequence */

   segs = (SegPtr) NULL;
   SegSeq (seqwin, sparamsp, &segs, 0);

/* merge the segment if desired. */

   if (sparamsp->overlaps)
	mergesegs(seqwin, segs, sparamsp->overlaps);

/* convert segs to seqlocs */

   slp = SegsToSeqLoc(bsp, segs);   

/* clean up & return */

   SeqFree (seqwin);
   SegFree (segs);
   SeqPortFree(spp);
   if(params_allocated)
       SegParamsFree(sparamsp);
   
   return (slp);
  }

/*----------------------------------------------------------(SeqlocSegAa)---*/

SeqLocPtr SeqlocSegAa (SeqLocPtr slpin, SegParamsPtr sparamsp)

{
    SeqPortPtr	spp=NULL;
    SeqLocPtr	slp = NULL;
    SequencePtr seqwin;
    SegPtr segs;
    Int4 index, len;
    Int4 start, stop, temp;
    CharPtr seq;
    Uint1       residue;
    Boolean params_allocated = FALSE;
    
    /* error msg stuff */
    
    ErrSetOptFlags (EO_MSG_CODES);
    SegParamsCheck (sparamsp);
    
    /* bail on null bioseq */
    
    if (!slpin) {
        ErrPostEx (SEV_ERROR, 0, 0, "no seqloc");
        ErrShow ();
        return (slp);
    }
    
    /* only simple intervals */
    
    if (slpin->choice != SEQLOC_INT &&
        slpin->choice != SEQLOC_WHOLE) return(slp);
    
    /* get coordinate range */
    
    start = SeqLocStart(slpin);
    stop = SeqLocStop(slpin);
    if (stop<start) {
        temp = start;
        start = stop;
        stop = temp;
    }
    
    len = stop - start + 1;
    
    /* check seg parameters */
    
    if (!sparamsp) {
        params_allocated = TRUE;
        sparamsp = SegParamsNewAa();
        SegParamsCheck (sparamsp);
        if (!sparamsp) {
            ErrPostEx (SEV_WARNING, 0, 0, "null parameters object");
            ErrShow();
            return(slp);
        }
    }
    
    /* make an old-style genwin sequence window object */
    
    seq = (CharPtr) MemNew((len+1)*sizeof(Char));
    spp = SeqPortNewByLoc(slpin, Seq_code_ncbieaa);
    if (spp == NULL) {
        ErrPostEx (SEV_ERROR, 0, 0, "SeqPortNew failure");
        ErrShow();
    }
    
    if (!seq) {
        ErrPostEx (SEV_ERROR, 0, 0, "memory allocation failure");
        ErrShow();
        SeqPortFree(spp);
        return(slp);
    }
    
   index = 0;
   while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
       if (IS_residue(residue)) {
           seq[index] = residue;
           index++;
       } else if (residue == SEQPORT_EOS) {
           continue; /*[Segment boundary]*/
       } else if (residue == SEQPORT_VIRT) {
           continue; /*[Virtual Sequence]*/
       } 
   }
   
   seq[index] = NULLB;
   
   seqwin = SeqNew();
   seqwin->seq = (CharPtr) seq;
   seqwin->length = len;
   seqwin->palpha = AlphaCopy(sparamsp->palpha);
   
   /* seg the sequence */
   
   segs = (SegPtr) NULL;
   SegSeq (seqwin, sparamsp, &segs, 0);

/* merge the segment if desired. */

   if (sparamsp->overlaps)
	mergesegs(seqwin, segs, sparamsp->overlaps);

   /* convert segs to seqlocs */
   
   slp = SeqlocSegsToSeqLoc(slpin, segs);   
   
   /* clean up & return */
   
   SeqFree (seqwin);
   SegFree (segs);
   SeqPortFree(spp);

   if(params_allocated)
       SegParamsFree(sparamsp);
   
   return (slp);
}

/*---------------------------------------------------------(SegsToSeqLoc)---*/

SeqLocPtr SegsToSeqLoc(BioseqPtr bsp, SegPtr segs)

  {
   SeqLocPtr slp, slp_last, slp_next, slp_packed;
   SeqIntPtr sip;

   if (bsp==NULL || segs==NULL) return ((SeqLocPtr) NULL);

   slp = (SeqLocPtr) ValNodeNew(NULL);
   slp->choice = SEQLOC_INT;

   sip = SeqIntNew();
   sip->from = segs->begin;
   sip->to = segs->end;
   sip->strand = bsp->strand;
   sip->id = SeqIdDup(bsp->id);

   slp->data.ptrvalue = sip;

/*---                             SEQLOC_INT for a single segment  ---*/

   if (segs->next==NULL)
     {
      slp->next = NULL;
      return(slp);
     }

/*---                     SEQLOC_PACKED_INT for multiple segments  ---*/

   slp_last = slp;
   for (segs=segs->next; segs!=NULL; segs=segs->next)
     {
      slp_next = (SeqLocPtr) ValNodeNew(NULL);
      slp_next->choice = SEQLOC_INT;

      sip = SeqIntNew();
      sip->from = segs->begin;
      sip->to = segs->end;
      sip->strand = bsp->strand;
      sip->id = SeqIdDup(bsp->id);

      slp_next->data.ptrvalue = sip;
      slp_last->next = slp_next;
      slp_last = slp_next;
     }

   slp_last->next = NULL;

   slp_packed = (SeqLocPtr) ValNodeNew(NULL);
   slp_packed->choice = SEQLOC_PACKED_INT;
   slp_packed->data.ptrvalue = slp;
   slp_packed->next = NULL;

   return(slp_packed);
  }

/*---------------------------------------------------(SeqlocSegsToSeqLoc)---*/

SeqLocPtr SeqlocSegsToSeqLoc(SeqLocPtr slpin, SegPtr segs)

  {
   SeqLocPtr slp, slp_last, slp_next, slp_packed;
   SeqIntPtr sip;
   Int4 start, stop;
   Uint1 strand;
   SeqIdPtr id;

   if (slpin==NULL || segs==NULL) return ((SeqLocPtr) NULL);
   if (slpin->choice != SEQLOC_INT &&
       slpin->choice != SEQLOC_WHOLE) return ((SeqLocPtr) NULL);

   start = SeqLocStart(slpin);
   stop = SeqLocStop(slpin);
   if (stop < start) return ((SeqLocPtr) NULL);
   strand = SeqLocStrand(slpin);
   id = SeqLocId(slpin);

   slp = (SeqLocPtr) ValNodeNew(NULL);
   slp->choice = SEQLOC_INT;

   sip = SeqIntNew();
   sip->from = segs->begin + start;
   sip->to = segs->end + start;
   sip->strand = strand;
   sip->id = SeqIdDup(id);

   slp->data.ptrvalue = sip;

/*---                             SEQLOC_INT for a single segment  ---*/

   if (segs->next==NULL)
     {
      slp->next = NULL;
      return(slp);
     }

/*---                     SEQLOC_PACKED_INT for multiple segments  ---*/

   slp_last = slp;
   for (segs=segs->next; segs!=NULL; segs=segs->next)
     {
      slp_next = (SeqLocPtr) ValNodeNew(NULL);
      slp_next->choice = SEQLOC_INT;

      sip = SeqIntNew();
      sip->from = segs->begin + start;
      sip->to = segs->end + start;
      sip->strand = strand;
      sip->id = SeqIdDup(id);

      slp_next->data.ptrvalue = sip;
      slp_last->next = slp_next;
      slp_last = slp_next;
     }

   slp_last->next = NULL;

   slp_packed = (SeqLocPtr) ValNodeNew(NULL);
   slp_packed->choice = SEQLOC_PACKED_INT;
   slp_packed->data.ptrvalue = slp;
   slp_packed->next = NULL;

   return(slp_packed);
  }

/*---------------------------------------------------------------(SeqNew)---*/

SequencePtr SeqNew(void)

  {
   SequencePtr seq;

   seq = (SequencePtr) MemNew(sizeof(Sequence));
   if (seq==NULL)
     {
      /* raise error flag and etc. */
      return(seq);
     }

   seq->parent = (SequencePtr) NULL;
   seq->seq = (CharPtr) NULL;
   seq->palpha = (AlphaPtr) NULL;
   seq->start = seq->length = 0;
   seq->bogus = seq->punctuation = FALSE;
   seq->composition = seq->state = (Int4Ptr) NULL;
   seq->entropy = (FloatHi) 0.0;

   return(seq);
  }

/*--------------------------------------------------------------(SeqFree)---*/

void SeqFree(SequencePtr seq)

  {
   if (seq==NULL) return;

   MemFree(seq->seq);
   AlphaFree(seq->palpha);
   MemFree(seq->composition);
   MemFree(seq->state);
   MemFree(seq);
   return;
  }

/*--------------------------------------------------------------(SegFree)---*/

void SegFree(SegPtr seg)

  {
   SegPtr nextseg;

   while (seg)
     {
      nextseg = seg->next;
      MemFree(seg);
      seg = nextseg;
     }

   return;
  }

/*---------------------------------------------------------------(SegSeq)---*/

void SegSeq(SequencePtr seq, SegParamsPtr sparamsp, SegPtr *segs,
            Int4 offset)

  {
   SegPtr seg, leftsegs;
   SequencePtr leftseq;
   Int4 window;
   FloatHi locut, hicut;
   Int4 maxbogus;
   Int4 downset, upset;
   Int4 first, last, lowlim;
   Int4 loi, hii, i;
   Int4 leftend, rightend, lend, rend;
   FloatHiPtr H;

   if (sparamsp->window<=0) return;
   if (sparamsp->locut<=0.) sparamsp->locut = 0.;
   if (sparamsp->hicut<=0.) sparamsp->hicut = 0.;

   window = sparamsp->window;
   locut = sparamsp->locut;
   hicut = sparamsp->hicut;
   downset = (window+1)/2 - 1;
   upset = window - downset;
   maxbogus = sparamsp->maxbogus;
      
   H = seqent(seq, window, maxbogus);
   if (H==NULL) return;

   first = downset;
   last = seq->length - upset;
   lowlim = first;

   for (i=first; i<=last; i++)
     {
      if (H[i]<=locut && H[i]!=-1)
        {
         loi = findlo(i, lowlim, hicut, H);
         hii = findhi(i, last, hicut, H);

         leftend = loi - downset;
         rightend = hii + upset - 1;

         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend,
              sparamsp);

         if (i+upset-1<leftend)   /* check for trigger window in left trim */
           {
            lend = loi - downset;
            rend = leftend - 1;

            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (SegPtr) NULL;
            SegSeq(leftseq, sparamsp, &leftsegs, offset+lend);
            if (leftsegs!=NULL)
              {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
              }
            closewin(leftseq);

/*          trim(openwin(seq, lend, rend-lend+1), &lend, &rend);
            seg = (SegPtr) MemNew(sizeof(Seg));
            seg->begin = lend;
            seg->end = rend;
            seg->next = (SegPtr) NULL;
            if (segs==NULL) segs = seg;
            else appendseg(segs, seg);  */
           }

         seg = (SegPtr) MemNew(sizeof(Seg));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (SegPtr) NULL;

         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);

         i = MIN(hii, rightend+downset);
         lowlim = i + 1;
/*       i = hii;     this ignores the trimmed residues... */
        }
     }

   MemFree(H);
   return;
  }

/*------------------------------------------------------------(seq_phase)---*/

SequencePtr seq_phase (SequencePtr seq, Int4 phase, Int4 period)

  {
   SequencePtr perseq;
   Int4 len, i, j;

   perseq = (SequencePtr) MemNew(sizeof(Sequence));

   len = ((seq->length)/period) + 1;
   perseq->seq = (CharPtr) MemNew((len+1)*sizeof(Char));

   perseq->length = 0;
   for (i=0, j=phase; j<seq->length; i++, j+=period)
     {
      perseq->seq[i] = seq->seq[j];
      perseq->length++;
     }
   perseq->seq[i] = '\0';

   return(perseq);
  }

/*---------------------------------------------------------------(seqent)---*/

FloatHiPtr seqent(SequencePtr seq, Int4 window, Int4 maxbogus)

  {
   SequencePtr win;
   FloatHiPtr H;
   Int4 i, first, last, downset, upset;

   downset = (window+1)/2 - 1;
   upset = window - downset;

   if (window>seq->length)
     {
      return((FloatHiPtr) NULL);
     }

   H = (FloatHiPtr) MemNew(seq->length*sizeof(FloatHi));

   for (i=0; i<seq->length; i++)
     {
      H[i] = -1.;
     }

   win = openwin(seq, 0, window);
   enton(win);

   first = downset;
   last = seq->length - upset;

   for (i=first; i<=last; i++)
     {
      if (seq->punctuation && hasdash(win))
        {
         H[i] = -1.;
         shiftwin1(win);
         continue;
        }
      if (win->bogus > maxbogus)
        {
         H[i] = -1.;
         shiftwin1(win);
         continue;
        }
      H[i] = win->entropy;
      shiftwin1(win);
     }

   closewin(win);
   return(H);
  }

/*--------------------------------------------------------------(hasdash)---*/

static Boolean hasdash(SequencePtr win)


{
	register char	*seq, *seqmax;

	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		if (*seq++ == '-')
			return TRUE;
	}
	return FALSE;
}

/*---------------------------------------------------------------(findlo)---*/

Int4 findlo(Int4 i, Int4 limit, FloatHi hicut, FloatHiPtr H)

  {
   Int4 j;

   for (j=i; j>=limit; j--)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j+1);
  }

/*---------------------------------------------------------------(findhi)---*/

Int4 findhi(Int4 i, Int4 limit, FloatHi hicut, FloatHiPtr H)

  {
   Int4 j;

   for (j=i; j<=limit; j++)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j-1);
  }

/*-----------------------------------------------------------------(trim)---*/

void trim(SequencePtr seq, Int4Ptr leftend, Int4Ptr rightend,
          SegParamsPtr sparamsp)

  {
   SequencePtr win;
   FloatHi prob, minprob;
   Int4 shift, len, i;
   Int4 lend, rend;
   Int4 minlen;
   Int4 maxtrim;

/* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   maxtrim = sparamsp->maxtrim;
   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

   minprob = 1.;
/*   fprintf(stderr, "\n");                                           -*/
   for (len=seq->length; len>minlen; len--)
     {
/*      fprintf(stderr, "%5d ", len);                                 -*/
      win = openwin(seq, 0, len);
      i = 0;

      shift = TRUE;
      while (shift)
        {
         prob = getprob(win->state, len, sparamsp->palpha);
/*         realprob = exp(prob);                 (for tracing the trim)-*/
/*         fprintf(stderr, "%2e ", realprob);                          -*/
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
           }
         shift = shiftwin1(win);
         i++;
        }
      closewin(win);
/*      fprintf(stderr, "\n");                                       -*/
     }

/*   fprintf(stderr, "%d-%d ", *leftend, *rightend);                 -*/

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/*   fprintf(stderr, "%d-%d\n", *leftend, *rightend);                -*/

   closewin(seq);
   return;
  }

/*--------------------------------------------------------------(getprob)---*/

FloatHi getprob(Int4Ptr sv, Int4 total, AlphaPtr palpha)

  {
   FloatHi ans, ans1, ans2, totseq;

   /* #define LN20	2.9957322735539909 */
   /* #define LN4	1.3862943611198906 */

   totseq = ((double) total) * palpha->lnalphasize;

/*
   ans = lnass(sv, palpha->alphasize) + lnperm(sv, total) - totseq;
*/
   ans1 = lnass(sv, palpha->alphasize);
   if (ans1 > -100000.0 && sv[0] != INT4_MIN)
   {
	ans2 = lnperm(sv, total);
   }
   else
   {
        ErrPostEx (SEV_ERROR, 0, 0, "Illegal value returned by lnass");
   }
   ans = ans1 + ans2 - totseq;
/*
fprintf(stderr, "%lf %lf %lf\n", ans, ans1, ans2);
*/

   return ans;
  }

/*---------------------------------------------------------------(lnperm)---*/

FloatHi lnperm(Int4Ptr sv, Int4 tot)

  {
   FloatHi ans;
   Int4 i;

   ans = lnfac[tot];

   for (i=0; sv[i]!=0; i++) 
     {
      ans -= lnfac[sv[i]];
     }

   return(ans);
  }

/*----------------------------------------------------------------(lnass)---*/

FloatHi lnass(Int4Ptr sv, Int4 alphasize)

{
	double	ans;
	int	svi, svim1;
	int	class, total;
	int    i;

	ans = lnfac[alphasize];
	if (sv[0] == 0)
		return ans;

	total = alphasize;
	class = 1;
	svi = *sv;
	svim1 = sv[0];
	for (i=0;; svim1 = svi) {
	        if (++i==alphasize) {
		        ans -= lnfac[class];
			break;
		      }
		else if ((svi = *++sv) == svim1) {
			class++;
			continue;
		}
		else {
			total -= class;
			ans -= lnfac[class];
			if (svi == 0) {
				ans -= lnfac[total];
				break;
			}
			else {
				class = 1;
				continue;
			}
		}
	}

	return ans;
}

/*------------------------------------------------------------(mergesegs)---*/
/* merge together overlapping segments, 
	hilenmin also does something, but we need to ask Scott Federhen what?
*/

void mergesegs(SequencePtr seq, SegPtr segs, Boolean overlaps)

  {SegPtr seg, nextseg;
   Int4 hilenmin;		/* hilenmin yet unset */
   Int4 len;

   hilenmin = 0;               /* hilenmin - temporary default */

   if (overlaps == FALSE) 
	return;
	
   if (segs==NULL) return;

   if (segs->begin<hilenmin) segs->begin = 0;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL)
     {
      if (seg->end>=nextseg->begin && seg->end>=nextseg->end)
        {
         seg->next = nextseg->next;
         MemFree(nextseg);
         nextseg = seg->next;
         continue;
        }
      if (seg->end>=nextseg->begin)               /* overlapping segments */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         MemFree(nextseg);
         nextseg = seg->next;
         continue;
        }
      len = nextseg->begin - seg->end - 1;
      if (len<hilenmin)                            /* short hient segment */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         MemFree(nextseg);
         nextseg = seg->next;
         continue;
        }
      seg = nextseg;
      nextseg = seg->next;
     }

   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;

   return;
  }

/*--------------------------------------------------------(per_mergesegs)---*/

SegPtr per_mergesegs(SequencePtr seq, SegPtr *persegs)

  {SegPtr *localsegs;
   SegPtr firstseg, segs, seg;
   Int4 period;                           /* period yet unset */
   Int4 first;
   Int4 phase, savephase;

   period = 1;                            /* period set temporarily */

   segs = (SegPtr) NULL;
   if (persegs==NULL) return(segs);

   localsegs = (SegPtr *) MemNew(period*sizeof(FloatHi));
   for (phase=0; phase<period; phase++) localsegs[phase] = persegs[phase];

   while (TRUE)
     {
      firstseg = (SegPtr) NULL;
      first = -1;
      savephase = -1;

      for (phase=0; phase<period; phase++)
        {
         if (localsegs[phase]==NULL) continue;
         if (first==-1 || localsegs[phase]->begin<first)
           {
            savephase = phase;
            firstseg = localsegs[phase];
            first = firstseg->begin;
           }
        }

      if (firstseg==NULL) break;

      seg = (SegPtr) MemNew(sizeof(Seg));
      seg->begin = ((firstseg->begin)*period)+savephase;
      seg->end = ((firstseg->end)*period)+savephase;
      seg->next = (SegPtr) NULL;
      if (segs==NULL) segs = seg; else appendseg(segs, seg);

      localsegs[savephase] = localsegs[savephase]->next;
     }

   MemFree(localsegs);
   mergesegs(seq, segs, FALSE);
   return(segs);
  }



/*------------------------------------------------------------(appendseg)---*/

static void
appendseg(SegPtr segs, SegPtr seg)

  {SegPtr temp;

   temp = segs;
   while (TRUE)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }

/*-------------------------------------------------------(SegParamsNewAa)---*/

SegParamsPtr SegParamsNewAa (void)

  {
   SegParamsPtr sparamsp;

   sparamsp = (SegParamsPtr) MemNew (sizeof(SegParams));

   sparamsp->window = 12;
   sparamsp->locut = 2.2;
   sparamsp->hicut = 2.5;
   sparamsp->period = 1;
   sparamsp->hilenmin = 0;
   sparamsp->overlaps = FALSE;
   sparamsp->maxtrim = 50;
   sparamsp->maxbogus = 2;
   sparamsp->palpha = AA20alpha();

   return (sparamsp);
  }

/*-------------------------------------------------------(SegParamsNewNa)---*/

SegParamsPtr SegParamsNewNa (void)

  {
   SegParamsPtr sparamsp;

   sparamsp = (SegParamsPtr) MemNew (sizeof(SegParams));

   sparamsp->window = 21;
   sparamsp->locut = 1.4;
   sparamsp->hicut = 1.6;
   sparamsp->period = 1;
   sparamsp->hilenmin = 0;
   sparamsp->overlaps = FALSE;
   sparamsp->maxtrim = 100;
   sparamsp->maxbogus = 3;
   sparamsp->palpha = NA4alpha();

   return (sparamsp);
  }

/*-------------------------------------------------------(SegParamsCheck)---*/

void SegParamsCheck (SegParamsPtr sparamsp)

  {
   if (!sparamsp) return;

   if (sparamsp->window <= 0) sparamsp->window = 12;

   if (sparamsp->locut < 0.0) sparamsp->locut = 0.0;
/* if (sparamsp->locut > sparamsp->palpha->lnalphasize)
       sparamsp->locut = sparamsp->palpha->lnalphasize; */

   if (sparamsp->hicut < 0.0) sparamsp->hicut = 0.0;
/* if (sparamsp->hicut > sparamsp->palpha->lnalphasize)
       sparamsp->hicut = sparamsp->palpha->lnalphasize; */

   if (sparamsp->locut > sparamsp->hicut)
       sparamsp->hicut = sparamsp->locut;

   if (sparamsp->maxbogus < 0)
       sparamsp->maxbogus = 0;
   if (sparamsp->maxbogus > sparamsp->window)
       sparamsp->maxbogus = sparamsp->window;

   if (sparamsp->period <= 0) sparamsp->period = 1;
   if (sparamsp->maxtrim < 0) sparamsp->maxtrim = 0;

   return;
  }

/*------------------------------------------------------------(AA20alpha)---*/

AlphaPtr AA20alpha (void)

  {
   AlphaPtr palpha;
   Int4Ptr alphaindex;
   BoolPtr alphaflag;
   CharPtr alphachar;
   Int4 i;
   Char c;

   palpha = (AlphaPtr) MemNew (sizeof(Alpha));

   palpha->alphabet = AA20;
   palpha->alphasize = 20;
   palpha->lnalphasize = LN20;

   alphaindex = (Int4Ptr) MemNew (CHAR_SET * sizeof(Int4));
   alphaflag = (BoolPtr) MemNew (CHAR_SET * sizeof(Boolean));
   alphachar = (CharPtr) MemNew (palpha->alphasize * sizeof(Char));

   for (i=0; i<128; i++)
     {
      c = (Char) i;
      if (c=='a' || c=='A')      
        {alphaflag[i] = FALSE; alphaindex[i] = 0; alphachar[0] = c;}
      else if (c=='c' || c=='C') 
        {alphaflag[i] = FALSE; alphaindex[i] = 1; alphachar[1] = c;}
      else if (c=='d' || c=='D')
        {alphaflag[i] = FALSE; alphaindex[i] = 2; alphachar[2] = c;} 
      else if (c=='e' || c=='E')
        {alphaflag[i] = FALSE; alphaindex[i] = 3; alphachar[3] = c;} 
      else if (c=='f' || c=='F')
        {alphaflag[i] = FALSE; alphaindex[i] = 4; alphachar[4] = c;} 
      else if (c=='g' || c=='G')
        {alphaflag[i] = FALSE; alphaindex[i] = 5; alphachar[5] = c;} 
      else if (c=='h' || c=='H')
        {alphaflag[i] = FALSE; alphaindex[i] = 6; alphachar[6] = c;} 
      else if (c=='i' || c=='I')
        {alphaflag[i] = FALSE; alphaindex[i] = 7; alphachar[7] = c;} 
      else if (c=='k' || c=='K')
        {alphaflag[i] = FALSE; alphaindex[i] = 8; alphachar[8] = c;} 
      else if (c=='l' || c=='L')
        {alphaflag[i] = FALSE; alphaindex[i] = 9; alphachar[9] = c;} 
      else if (c=='m' || c=='M')
        {alphaflag[i] = FALSE; alphaindex[i] = 10; alphachar[10] = c;} 
      else if (c=='n' || c=='N')
        {alphaflag[i] = FALSE; alphaindex[i] = 11; alphachar[11] = c;}
      else if (c=='p' || c=='P')
        {alphaflag[i] = FALSE; alphaindex[i] = 12; alphachar[12] = c;} 
      else if (c=='q' || c=='Q')
        {alphaflag[i] = FALSE; alphaindex[i] = 13; alphachar[13] = c;} 
      else if (c=='r' || c=='R')
        {alphaflag[i] = FALSE; alphaindex[i] = 14; alphachar[14] = c;} 
      else if (c=='s' || c=='S')
        {alphaflag[i] = FALSE; alphaindex[i] = 15; alphachar[15] = c;} 
      else if (c=='t' || c=='T')
        {alphaflag[i] = FALSE; alphaindex[i] = 16; alphachar[16] = c;} 
      else if (c=='v' || c=='V')
        {alphaflag[i] = FALSE; alphaindex[i] = 17; alphachar[17] = c;} 
      else if (c=='w' || c=='W')
        {alphaflag[i] = FALSE; alphaindex[i] = 18; alphachar[18] = c;} 
      else if (c=='y' || c=='Y')
        {alphaflag[i] = FALSE; alphaindex[i] = 19; alphachar[19] = c;} 
      else
        {alphaflag[i] = TRUE; alphaindex[i] = 20;}
     }

   palpha->alphaindex = alphaindex;
   palpha->alphaflag = alphaflag;
   palpha->alphachar = alphachar;

   return (palpha);
  }

/*-------------------------------------------------------------(NA4alpha)---*/

AlphaPtr NA4alpha (void)

  {
   AlphaPtr palpha;
   Int4Ptr alphaindex;
   BoolPtr alphaflag;
   CharPtr alphachar;
   Int4 i;
   Char c;

   palpha = (AlphaPtr) MemNew (sizeof(Alpha));

   palpha->alphabet = NA4;
   palpha->alphasize = 4;
   palpha->lnalphasize = LN4;

   alphaindex = (Int4Ptr) MemNew (CHAR_SET * sizeof(Int4));
   alphaflag = (BoolPtr) MemNew (CHAR_SET * sizeof(Boolean));
   alphachar = (CharPtr) MemNew (palpha->alphasize * sizeof(Char));

   for (i=0; i<128; i++)
     {
      c = (Char) i;
      if (c=='a' || c=='A')
        {alphaflag[i] = FALSE; alphaindex[i] = 0; alphachar[0] = c;}
      else if (c=='c' || c=='C')
        {alphaflag[i] = FALSE; alphaindex[i] = 1; alphachar[1] = c;}
      else if (c=='g' || c=='G')
        {alphaflag[i] = FALSE; alphaindex[i] = 2; alphachar[2] = c;} 
      else if (c=='t' || c=='T')
        {alphaflag[i] = FALSE; alphaindex[i] = 3; alphachar[3] = c;} 
      else if (c=='u' || c=='U')
        {alphaflag[i] = FALSE; alphaindex[i] = 3;}
      else 
        {alphaflag[i] = TRUE; alphaindex[i] = 4;}
     }

   palpha->alphaindex = alphaindex;
   palpha->alphaflag = alphaflag;
   palpha->alphachar = alphachar;

   return (palpha);
  }

/*------------------------------------------------------------(AlphaFree)---*/

void AlphaFree (AlphaPtr palpha)

  {
   if (!palpha) return;

   MemFree (palpha->alphaindex);
   MemFree (palpha->alphaflag);
   MemFree (palpha->alphachar);
   MemFree (palpha);

   return;
 }

/*------------------------------------------------------------(AlphaCopy)---*/

AlphaPtr AlphaCopy (AlphaPtr palpha)

  {
   AlphaPtr pbeta;
   Int2 i;

   if (!palpha) return((AlphaPtr) NULL);

   pbeta = (AlphaPtr) MemNew(sizeof(Alpha));
   pbeta->alphabet = palpha->alphabet;
   pbeta->alphasize = palpha->alphasize;
   pbeta->lnalphasize = palpha->lnalphasize;

   pbeta->alphaindex = (Int4Ptr) MemNew (CHAR_SET * sizeof(Int4));
   pbeta->alphaflag = (BoolPtr) MemNew (CHAR_SET * sizeof(Boolean));
   pbeta->alphachar = (CharPtr) MemNew (palpha->alphasize * sizeof(Char));

   for (i=0; i<CHAR_SET; i++)
     {
      pbeta->alphaindex[i] = palpha->alphaindex[i];
      pbeta->alphaflag[i] = palpha->alphaflag[i];
     }

   for (i=0; i<palpha->alphasize; i++)
     {
      pbeta->alphachar[i] = palpha->alphachar[i];
     }

   return(pbeta);
  }

/*--------------------------------------------------------(SegParamsFree)---*/

void SegParamsFree(SegParamsPtr sparamsp)

  {
   if (!sparamsp) return;
   AlphaFree(sparamsp->palpha);
   MemFree(sparamsp);
   return;
  }

/*--------------------------------------------------------------(openwin)---*/

SequencePtr openwin(SequencePtr parent, Int4 start, Int4 length)

  {
   SequencePtr win;

   if (start<0 || length<0 || start+length>parent->length)
     {
      return((SequencePtr) NULL);
     }

   win = (SequencePtr) MemNew(sizeof(Sequence));

/*---                                          ---[set links, up and down]---*/

   win->parent = parent;
   win->palpha = parent->palpha;

/*---                          ---[install the local copy of the sequence]---*/

   win->start = start;
   win->length = length;
#if 0                                                    /* Hi Warren! */
   win->seq = (CharPtr) MemNew(sizeof(Char)*length + 1);
   memcpy(win->seq, (parent->seq)+start, length);
   win->seq[length] = '\0';
#else
	win->seq = parent->seq + start;
#endif

        win->bogus = 0;
	win->punctuation = FALSE;

	win->entropy = -2.;
	win->state = (Int4Ptr) NULL;
	win->composition = (Int4Ptr) NULL;

	stateon(win);

	return win;
}

/*------------------------------------------------------------(shiftwin1)---*/

Int4 shiftwin1(SequencePtr win)

{
	Int4 j, length;
	Int4Ptr comp;
        Int4Ptr alphaindex;
        BoolPtr alphaflag;

	length = win->length;
	comp = win->composition;
        alphaindex = win->palpha->alphaindex;
        alphaflag = win->palpha->alphaflag;

	if ((++win->start + length) > win->parent->length) {
		--win->start;
		return FALSE;
	}

	if (!alphaflag[j = win->seq[0]])
		decrementsv(win->state, comp[alphaindex[j]]--);
        else win->bogus--;

	j = win->seq[length];
	++win->seq;

	if (!alphaflag[j])
		incrementsv(win->state, comp[alphaindex[j]]++);
        else win->bogus++;

	if (win->entropy > -2.)
		win->entropy = entropy(win->state);

	return TRUE;
}

/*-------------------------------------------------------------(closewin)---*/

void closewin(SequencePtr win)

  {
   if (win==NULL) return;

   if (win->state!=NULL)       MemFree(win->state);
   if (win->composition!=NULL) MemFree(win->composition);

   MemFree(win);
   return;
  }

/*---------------------------------------------------------------(compon)---*/

void compon(SequencePtr win)

{
	Int4Ptr comp;
	Int4 letter;
	CharPtr seq, seqmax;
        Int4Ptr alphaindex;
        BoolPtr alphaflag;
        Int4 alphasize;

        alphasize = win->palpha->alphasize;
        alphaindex = win->palpha->alphaindex;
        alphaflag = win->palpha->alphaflag;

	win->composition = comp =
                (Int4Ptr) MemNew(alphasize*sizeof(Int4));
	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		letter = *seq++;
		if (!alphaflag[letter])
			comp[alphaindex[letter]]++;
                else win->bogus++;
	}

	return;
}

/*------------------------------------------------------------(state_cmp)---*/

static int LIBCALLBACK state_cmp(VoidPtr s1, VoidPtr s2)

{
	int *np1, *np2;

	np1 = (int *) s1;
	np2 = (int *) s2;

	return (*np2 - *np1);
}

/*--------------------------------------------------------------(stateon)---*/

void stateon(SequencePtr win)

{
	Int4 letter, nel, c;
        Int4 alphasize;

        alphasize = win->palpha->alphasize;

	if (win->composition == NULL)
		compon(win);

	win->state = (Int4Ptr) MemNew((alphasize+1)*sizeof(win->state[0]));

	for (letter = nel = 0; letter < alphasize; ++letter) {
		if ((c = win->composition[letter]) == 0)
			continue;
		win->state[nel++] = c;
	}
	for (letter = nel; letter < alphasize+1; ++letter)
		win->state[letter] = 0;

	HeapSort(win->state, nel, sizeof(win->state[0]), state_cmp);

	return;
}

/*----------------------------------------------------------------(enton)---*/

void enton(SequencePtr win)

  {
   if (win->state==NULL) {stateon(win);}

   win->entropy = entropy(win->state);

   return;
  }

/*--------------------------------------------------------------(entropy)---*/

FloatHi entropy(Int4Ptr sv)

  {FloatHi ent;
   Int4 i, total;

   total = 0;
   for (i=0; sv[i]!=0; i++)
     {
      total += sv[i];
     }
   if (total==0) return(0.);

   ent = 0.0;
   for (i=0; sv[i]!=0; i++)
     {
      ent += ((FloatHi)sv[i])*log(((FloatHi)sv[i])/(double)total)/LN2;
     }

   ent = fabs(ent/(FloatHi)total);

   return(ent);
  }


/*----------------------------------------------------------(decrementsv)---*/

void decrementsv(Int4Ptr sv, Int4 class)

{
	Int4	svi;

	while ((svi = *sv++) != 0) {
		if (svi == class && *sv < class) {
			sv[-1] = svi - 1;
			break;
		}
	}
}

/*----------------------------------------------------------(incrementsv)---*/

void incrementsv(Int4Ptr sv, Int4 class)

{
	for (;;) {
		if (*sv++ == class) {
			sv[-1]++;
			break;
		}
	}
}

/*----------------------------------------------------------------(upper)---*/

extern void upper(CharPtr string, size_t len)

{
	CharPtr stringmax;
        Char c;

	for (stringmax = string + len; string < stringmax; ++string)
		if (islower(c = *string))
			*string = toupper(c);
}

/*----------------------------------------------------------------(lower)---*/

extern void lower(CharPtr string, size_t len)

{
	CharPtr stringmax;
        Char c;

	for (stringmax = string + len; string < stringmax; ++string)
		if (isupper(c = *string))
			*string = tolower(c);
}

/*--------------------------------------------------------------------------*/
