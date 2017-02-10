/*  edutil.c
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
* File Name:  edutil.c
*
* Author:  James Ostell
*   
* Version Creation Date: 2/4/94
*
* $Revision: 6.56 $
*
* File Description:  Sequence editing utilities
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: edutil.c,v $
* Revision 6.56  2006/04/04 18:00:47  kans
* SeqLocAddEx properly returns value to &last argument, makes SeqLocMix from DeltaSeqsToSeqLocs
*
* Revision 6.55  2006/03/30 19:50:15  kans
* DeltaSeqsToSeqLocs calls SeqLocAddEx for efficient list usage
*
* Revision 6.54  2006/02/07 13:41:29  bollin
* added function AdjustFeatureForGapChange, which changes a feature to accommodate
* a change in the length of a gap
*
* Revision 6.53  2005/12/12 14:12:54  bollin
* BioseqCopyEx was not correctly handling copying the data contents of a
* delta sequence
*
* Revision 6.52  2005/09/22 19:21:34  bollin
* In the sequence editor, if the user inserts Ns into a gap of known length,
* the gap length will be increased instead of creating two gaps on either side
* with N sequence characters in the middle.
*
* Revision 6.51  2005/09/13 15:21:57  bollin
* fixed bug when inserting characters inside a gap that was incorrectly setting
* the lengths of the split gap
*
* Revision 6.50  2005/09/13 14:14:31  bollin
* fixed bug that was preventing the removal of gaps of length 1
*
* Revision 6.49  2005/07/15 19:01:37  kans
* minor fixes for Xcode warnings
*
* Revision 6.48  2005/05/02 14:20:02  bollin
* when inserting gaps, adjust coding region locations to not include gaps.
* when removing gaps, if a feature location has intervals that stop and start
* again at the point where the gap was removed, connect the intervals.
*
* Revision 6.47  2005/04/28 20:10:31  bollin
* added new function AdjustFeaturesForInsertion which is called by BioseqInsert
* and also by a new function in sequin3.c for converting a raw bioseq to a delta
* and inserting gaps
*
* Revision 6.46  2005/04/06 19:33:15  bollin
* made it possible to insert and remove gaps from delta sequences
*
* Revision 6.45  2005/03/18 20:51:10  bollin
* only change frame when CDS location has been changed, change anticodon locations
* and code breaks when locations have just been shifted
*
* Revision 6.44  2005/03/08 21:14:44  bollin
* strand argument in SeqLocCopyRegion is Seq_strand_minus when features
* should be reverse-complemented, does not actually indicate the strand to
* which a feature should be copied
*
* Revision 6.43  2005/02/28 16:53:40  bollin
* corrected Unix compiler warnings
*
* Revision 6.42  2005/02/28 16:08:35  bollin
* added utilities for editing delta sequences
*
* Revision 6.41  2005/01/24 17:00:58  bollin
* only change frames, fix code break locations, and fix anticodon locations
* when feature location is changed in SeqFeatDelete
*
* Revision 6.40  2004/11/17 21:19:18  lavr
* AffectedFeatFree() to return NULL on afp == NULL
*
* Revision 6.39  2004/10/08 16:04:16  bollin
* added ability to check when an action will remove a feature
*
* Revision 6.38  2004/10/08 15:19:07  bollin
* do not set partial flag when deleting from bioseq location in feature
*
* Revision 6.37  2004/09/29 18:49:57  bollin
* fixed bugs in sequence editing, can now undo a nucleotide deletion that
* removes an entire feature location (feature will be restored)
*
* Revision 6.36  2004/09/23 14:59:51  bollin
* moved functions that depend on functions that depend on BLAST functions
* into seqpanel.c, made function scalled by those functions extern
*
* Revision 6.35  2004/09/22 20:12:27  bollin
* fixed error in deleting sequence location for point features
*
* Revision 6.34  2004/09/22 18:20:32  bollin
* added functions for playing and unplaying a sequence editor action to translate
* a CDS
*
* Revision 6.33  2004/09/07 14:52:29  bollin
* when deleting location from a feature, adjust frame if deleting from 5' end of 5' partial feature.
*
* Revision 6.32  2004/08/24 13:16:57  bollin
* do not free list of product features taken from ObjectMgrDataPtr
*
* Revision 6.31  2004/08/06 19:56:20  bollin
* allow deletion from the end of a sequence
*
* Revision 6.30  2004/08/05 18:15:02  bollin
* when maintaining partials during feature drag, use partial in orig_loc
* instead of current feature location
*
* Revision 6.29  2004/08/05 18:07:03  bollin
* maintain partials for features when dragging or sliding intervals
*
* Revision 6.28  2004/07/30 18:46:55  bollin
* added function for reordering intervals after they have been dragged by
* the sequence editor
*
* Revision 6.27  2004/07/30 13:34:50  bollin
* in SeqLocCopyRegion, when copying from the minus strand to a non-minus-strand,
* be sure to set the strand.
*
* Revision 6.26  2004/07/28 20:06:19  bollin
* added journaling for undo/redo of dragged sequence location changes
*
* Revision 6.25  2004/07/28 15:22:15  bollin
* moved functions for moving feature locations around to edutil.c from
* seqpanel.c
*
* Revision 6.24  2004/07/27 19:46:42  bollin
* fixed errors in feature location adjustment when deleting nucleotides
* with new sequence editor
*
* Revision 6.23  2004/07/22 16:08:20  bazhin
* Changes to parse gaps of unknown lengths (like "gap(unk100)")
* within location strings.
*
* Revision 6.22  2004/07/12 12:29:45  bollin
* moved new sequence editor editing functions here
*
* Revision 6.21  2003/11/03 19:37:42  bollin
* SegLocToPartsEx now handles SEQLOC_PNT as well as SEQLOC_INT
*
* Revision 6.20  2003/06/03 20:25:34  kans
* SeqLocReplaceID works on bonds if both ends bonded to the same Seq-id
*
* Revision 6.19  2003/02/10 22:57:45  kans
* added BioseqCopyEx, which takes a BioseqPtr instead of a SeqIdPtr for the source
*
* Revision 6.18  2002/07/26 20:15:55  kans
* BioseqInsert can do feature indexed collection of features to adjust
*
* Revision 6.17  2002/07/17 15:39:40  kans
* BioseqInsert calls Nlm_BSAdd, need to figure out when not to call
*
* Revision 6.16  2002/07/11 17:45:53  kans
* BioseqInsert does not call Nlm_BSAdd due to a bug in that code
*
* Revision 6.15  2002/07/02 13:23:42  kans
* added SeqLocDeleteEx
*
* Revision 6.14  2001/06/01 18:07:20  kans
* changes to SeqLocAdd to allow one plus and one unknown strand to be accepted
*
* Revision 6.13  2001/02/23 21:30:09  shkeda
* Fixed SeqLocAdd: Int-fuzz pointers should be set to NULL after IntFuzzFree
*
* Revision 6.12  2001/02/23 01:26:07  ostell
* Added support to BioseqDelete() for delta seqs
*
* Revision 6.11  2000/10/31 17:11:06  kans
* SeqLocReplaceID was handling SEQLOC_PACKED_PNT incorrectly
*
* Revision 6.10  1999/12/20 20:47:12  kans
* oldscope test was wrong everywhere
*
* Revision 6.9  1999/12/15 20:52:16  kans
* added IndexedSeqFeatsCopy if SeqMgrFeaturesAreIndexed
*
* Revision 6.8  1999/12/07 20:32:13  kans
* for most editing functions, if BioseqFind failed, temporarily clear scope/try again/reset scope
*
* Revision 6.7  1999/11/19 19:54:19  kans
* SeqLocAdd checks for NULL slp before dereferencing
*
* Revision 6.6  1998/09/03 20:43:52  kans
* added delta bioseq support to BioseqCopy
*
* Revision 6.5  1998/06/22 20:00:46  kans
* DelFeat was a bit too agressive when there were multiple feature tables
*
* Revision 6.4  1998/06/17 21:50:11  kans
* fixed unix compiler warnings, including 64-bit SGI
*
* Revision 6.3  1997/11/10 19:40:48  bazhin
* Fixed incorrect comment for ISAGappedSeqLoc() function.
*
* Revision 6.2  1997/10/24 19:16:17  bazhin
* Added three easy functions GapToSeqLoc(...), ISAGappedSeqLoc(...)
* and GappedSeqLocsToDeltaSeqs(...) for processing "gap(...)" tokens
* in CONTIG line.
*
* Revision 6.1  1997/10/10 20:18:02  ostell
* removed tab character from SeqLitTag for DeltaSeqsToSeqLoc
*
* Revision 6.0  1997/08/25 18:05:24  madden
* Revision changed to 6.0
*
* Revision 5.10  1997/07/25 20:34:51  kans
* added SegLocToPartsEx
*
* Revision 5.9  1997/06/19 18:37:30  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.8  1996/12/20 17:59:34  kans
* SeqLocCopyRegion already reversed order for Seq_strand_minus, so no need
* to reverse it again (JO + JK)
*
 * Revision 5.7  1996/10/21  18:56:19  ostell
 * made SegLocToParts accept a complicated Seq-loc argument
 *
 * Revision 5.6  1996/10/09  17:27:34  chappey
 * *** empty log message ***
 *
 * Revision 5.5  1996/10/09  16:34:59  chappey
 * added SeqLocReplaceID() that replaces the Seq-Id of a Seq-Loc
 *
 * Revision 5.4  1996/07/15  14:43:51  epstein
 * change SeqLocAdd() so that it merges identical SEQLOC_PNTs
 *
 * Revision 5.3  1996/06/12  18:29:41  epstein
 * move SeqLocIntNew() and SeqLocPntNew() from edutil to sequtil
 *
 * Revision 5.1  1996/06/10  15:07:17  epstein
 * replace make_seq_loc() with SeqLocIntNew() and make_pnt_loc with SeqLocPntNew()
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 4.10  1996/03/19  19:45:24  kans
 * fix of SegLocToParts (JO)
 *
 * Revision 4.9  1996/03/12  22:14:22  ostell
 * added SeqLocToParts()
 *
 * Revision 4.7  1996/02/19  19:58:05  ostell
 * added support for Code-break and tRNA.anticodon
 *
 * Revision 4.6  1996/01/30  16:24:04  ostell
 * changed name of SeqLocPack() to SeqLocPackage()
 *
 * Revision 4.5  1996/01/29  22:03:52  ostell
 * revised SeqLocAdd
 * added SeqLocPack
 *
 * Revision 4.4  1996/01/10  22:25:25  ostell
 * added SeqLocIntNew()
 *
 * Revision 4.3  1995/12/29  21:31:44  ostell
 * added mapping functions between delta seq and seq loc, for editing utilities
 *
 * Revision 4.2  1995/12/21  02:35:50  ostell
 * changed call for BSAdd
 *
 * Revision 4.1  1995/11/15  20:40:20  ostell
 * fixed SeqLocCopyPart so it correctly handles SEQLOC_NULL in segmented
 * records
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 1.22  1995/05/15  21:46:05  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

#include <edutil.h>
#include <explore.h>
#include <sqnutils.h>
#include <objfdef.h>
#include <gather.h>

/*****************************************************************************
*
*   SeqLocPackage(head)
*     head is a chain of 1 or more SeqLocs connected by slp->next
*     Assumes was built by SeqLocAdd to remove redundancy
*     Frees the last element if it is a NULL.
*     If more than one element left, then packages the chain into a SEQLOC_MIX,
*       or SEQLOC_PACKED_INT as appropriate
*     returns pointer to the head of the resulting single SeqLoc
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SeqLocPackage (SeqLocPtr head)
{
	SeqLocPtr newhead = NULL, tmp, prev;
	Boolean packed_int = TRUE;
	Int4 ctr = 0;

	if (head == NULL) return head;

	prev = NULL;    /* remove trailing NULL */
	for (tmp = head; tmp->next != NULL; tmp = tmp->next)
		prev = tmp;

	if (tmp->choice == SEQLOC_NULL)
	{
		SeqLocFree(tmp);
		if (prev != NULL)
			prev->next = NULL;
		else
			return NULL;   /* nothing left */
	}

	for (tmp = head; tmp != NULL; tmp = tmp->next)
	{
		ctr++;
		if (tmp->choice != SEQLOC_INT)
			packed_int = FALSE;
	}

	if (ctr == 1)
		return head;

	newhead = ValNodeNew(NULL);
	if (packed_int)
		newhead->choice = SEQLOC_PACKED_INT;
	else
		newhead->choice = SEQLOC_MIX;
	newhead->data.ptrvalue = head;

	return newhead;
}

/*****************************************************************************
*
*   SeqLocAdd(headptr, slp, merge, do_copy)
*   	creates a linked list of SeqLocs.
*       returns a pointer to the last SeqLoc in the chain
*       if (merge)
*   	  deletes double NULLs or Nulls at start (application must delete at stop)
*         merges adjacent intervals on the same strand
*       if (do_copy)
*   	  Makes copies of incoming SeqLocs
*         if incoming is merged, deletes the incoming SeqLoc
*
*****************************************************************************/
static SeqLocPtr LIBCALL SeqLocAddEx (SeqLocPtr PNTR head, SeqLocPtr PNTR lastp, SeqLocPtr slp, Boolean merge, Boolean do_copy)
{
	SeqLocPtr tmp, last = NULL, retval = NULL;
	Boolean merged = FALSE;   /* intervals were merged */

	if (slp == NULL) return NULL;

    if (lastp != NULL) {
        last = *lastp;
    } else if (head != NULL && *head != NULL)
	{
		for (tmp = *head; tmp != NULL; tmp = tmp->next)
		{
			last = tmp;
		}
	}

	if ((slp->choice == SEQLOC_NULL) && (merge))  /* no null at start, or two in a row */
	{
		if (last == NULL)  /* first one */
		{
			merged = TRUE;
			goto ret;
		}
		if (last->choice == SEQLOC_NULL)  /* double NULL */
		{
			merged = TRUE;
			goto ret;
		}
	}

	if ((last != NULL) && (merge))     /* check for merging intervals */
	{
		if ((last->choice == SEQLOC_INT) && (slp->choice == SEQLOC_INT))
		{
			SeqIntPtr sip1, sip2;
			Boolean samestrand;
			Uint1 strand = Seq_strand_unknown;

			sip1 = (SeqIntPtr)(last->data.ptrvalue);
			sip2 = (SeqIntPtr)(slp->data.ptrvalue);
			samestrand = FALSE;
			if ((sip1->strand == sip2->strand) ||
				(sip1->strand == Seq_strand_unknown && sip2->strand != Seq_strand_minus) ||
          		(sip1->strand == Seq_strand_unknown && sip2->strand != Seq_strand_minus)) {
				samestrand = TRUE;
				if (sip1->strand == Seq_strand_minus || sip1->strand == Seq_strand_minus) {
					strand = Seq_strand_minus;
				} else if (sip1->strand == Seq_strand_plus || sip1->strand == Seq_strand_plus) {
					strand = Seq_strand_plus;
				} else {
					strand = Seq_strand_unknown;
				}
          	}
			if (samestrand && (SeqIdForSameBioseq(sip1->id, sip2->id)))
			{
				if (strand == Seq_strand_minus)
				{
					if (sip1->from == (sip2->to + 1))  /* they are adjacent */
					{
						sip1->from = sip2->from;
						sip1->if_from = IntFuzzFree(sip1->if_from);
						if (sip2->if_from != NULL)   /* copy the fuzz */
						{
							if (do_copy)
								sip1->if_from = (IntFuzzPtr)AsnIoMemCopy((Pointer)(sip2->if_from),
								    (AsnReadFunc)IntFuzzAsnRead, (AsnWriteFunc)IntFuzzAsnWrite);
							else
							{
								sip1->if_from = sip2->if_from;
								sip2->if_from = NULL;
							}
							sip1->strand = strand;
						}
						merged = TRUE;
					}
				}
				else
				{
					if (sip1->to == (sip2->from - 1))  /* they are adjacent */
					{
						sip1->to = sip2->to;
						sip1->if_to = IntFuzzFree(sip1->if_to);
						if (sip2->if_to != NULL)   /* copy the fuzz */
						{
							if (do_copy)
								sip1->if_to = (IntFuzzPtr)AsnIoMemCopy((Pointer)(sip2->if_to),
								    (AsnReadFunc)IntFuzzAsnRead, (AsnWriteFunc)IntFuzzAsnWrite);
							else
							{
								sip1->if_to = sip2->if_to;
								sip2->if_to = NULL;
							}
							sip1->strand = strand;
						}
						merged = TRUE;
					}
				}
			}
		} else if ((last->choice == SEQLOC_PNT) && (slp->choice == SEQLOC_PNT))
		{
			SeqPntPtr sip1, sip2;

			sip1 = (SeqPntPtr)(last->data.ptrvalue);
			sip2 = (SeqPntPtr)(slp->data.ptrvalue);
			if ((sip1->strand == sip2->strand) && sip1->point == sip2->point && (SeqIdForSameBioseq(sip1->id, sip2->id)))
			{
				sip1->fuzz = IntFuzzFree(sip1->fuzz);
				if (sip2->fuzz != NULL)   /* copy the fuzz */
				{
					if (do_copy)
						sip1->fuzz = (IntFuzzPtr)AsnIoMemCopy((Pointer)(sip2->fuzz),
						    (AsnReadFunc)IntFuzzAsnRead, (AsnWriteFunc)IntFuzzAsnWrite);
					else
					{
						sip1->fuzz = sip2->fuzz;
						sip2->fuzz = NULL;
					}
				}
				merged = TRUE;
			}
		}
	}

ret:
	if (! merged)  /* then have to add a new one */
	{
		if (do_copy)
			tmp = (SeqLocPtr)AsnIoMemCopy((Pointer)slp, (AsnReadFunc)SeqLocAsnRead, (AsnWriteFunc)SeqLocAsnWrite);
		else
			tmp = slp;

		if (tmp != NULL) {
			tmp->next = NULL;
		}

		if (last != NULL) {
			last->next = tmp;
		} else if (head != NULL) {
			*head = tmp;
		}
		last = tmp;
		retval = tmp;
	}
	else
	{
		retval = last;
		if (! do_copy)   /* got to free it here */
			SeqLocFree(slp);
	}
	if (lastp != NULL) {
	    *lastp = last;
	}
		
	return retval;
}

NLM_EXTERN SeqLocPtr LIBCALL SeqLocAdd (SeqLocPtr PNTR head, SeqLocPtr slp, Boolean merge, Boolean do_copy)
{
	SeqLocPtr tmp, last;

	if (slp == NULL) return NULL;

	last = NULL;
	if (* head != NULL)
	{
		for (tmp = *head; tmp != NULL; tmp = tmp->next)
		{
			last = tmp;
		}
	}
	return SeqLocAddEx (head, &last, slp, merge, do_copy);
}

/*****************************************************************************
*
*   SegLocToParts(BioseqPtr seg, SeqLocPtr slp)
*   	seg must be a segmented Bioseq
*       slp must be a SeqLoc on it
*       function maps slp to the components of seg
*       returns a new SeqLocPtr
*       does not delete slp
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SegLocToPartsEx (BioseqPtr seg, SeqLocPtr slp, Boolean nullsBetween)
{
	SeqLocPtr newloc = NULL, tmp, tmp2, tmp3, next, curr;
	ValNode thead;
	SeqIdPtr sip, tsip;
	Int4 left_end, right_end, tlen, tstart;
	SeqIntPtr sintp;
	Boolean split, notFirst = FALSE;

	if ((seg == NULL) || (slp == NULL)) return newloc;
	if (seg->repr != Seq_repr_seg) return newloc;

	sip = SeqLocId(slp);
	if (sip == NULL) return newloc;
	if (! SeqIdIn(sip, seg->id)) return newloc;

	MemSet(&thead, 0, sizeof(ValNode));
	thead.choice = SEQLOC_MIX;
	thead.data.ptrvalue = seg->seq_ext;

	curr = NULL;
	while ((curr = SeqLocFindNext(slp, curr)) != NULL)
	{
		left_end = 0;
		tmp = NULL;
		while ((tmp = SeqLocFindNext(&thead, tmp)) != NULL)
		{
			tlen = SeqLocLen(tmp);
			if (tlen > 0)
			{
				right_end = left_end + tlen - 1;
				tsip = SeqLocId(tmp);
				tstart = SeqLocStart(tmp);
				tmp2 = SeqLocCopyRegion(tsip, curr, seg, left_end, right_end, SeqLocStrand(tmp),
					&split);
				while (tmp2 != NULL)
				{
				  next = tmp2->next;
				  tmp2->next = NULL;
				  if (tmp2->choice == SEQLOC_INT)
				  {
				    if (nullsBetween  && notFirst) {
				      tmp3 = ValNodeNew (NULL);
				      if (tmp3 != NULL) {
				        tmp3->choice = SEQLOC_NULL;
				        SeqLocAdd (&newloc, tmp3, TRUE, FALSE);
				      }
				    }
				    notFirst = TRUE;
				    sintp = (SeqIntPtr)(tmp2->data.ptrvalue);
				    sintp->from += tstart;
				    sintp->to += tstart;
				    SeqLocAdd(&newloc, tmp2, TRUE, FALSE);
				  }
                                  else if (tmp2->choice == SEQLOC_PNT)
                                  {
				    if (nullsBetween  && notFirst) {
				      tmp3 = ValNodeNew (NULL);
				      if (tmp3 != NULL) {
				        tmp3->choice = SEQLOC_NULL;
				        SeqLocAdd (&newloc, tmp3, TRUE, FALSE);
				      }
				    }
				    notFirst = TRUE;
                                    SeqLocAdd (&newloc, tmp2, TRUE, FALSE);
                                  }
				  tmp2 = next;
				}
				left_end = right_end + 1;
			}
		}
	}

	if (newloc != NULL)
		newloc = SeqLocPackage(newloc);
	return newloc;
}

NLM_EXTERN SeqLocPtr LIBCALL SegLocToParts (BioseqPtr seg, SeqLocPtr slp)

{
	return SegLocToPartsEx (seg, slp, FALSE);
}

static CharPtr seqlitdbtag = "SeqLit";
static CharPtr unkseqlitdbtag = "UnkSeqLit";
/*****************************************************************************
*
*   ISADeltaSeqsToSeqLoc(slp)
*   	returns Index (> 0) if this (one) SeqLoc was converted from a Delta Seq by
*         DeltaSeqsToSeqLocs() by looking for the special Dbtag name
*
*****************************************************************************/
NLM_EXTERN Int4 LIBCALL ISADeltaSeqsToSeqLoc (SeqLocPtr slp)
{
	SeqIdPtr sip;
	Int4 retval = 0;

	if (slp == NULL) return retval;
	sip = SeqLocId(slp);
	if (sip == NULL) return retval;

	if (sip->choice != SEQID_GENERAL) return retval;

	if (! StringCmp(seqlitdbtag, ((DbtagPtr)(sip->data.ptrvalue))->db) ||
	    ! StringCmp(unkseqlitdbtag, ((DbtagPtr)(sip->data.ptrvalue))->db))
		retval = (((DbtagPtr)(sip->data.ptrvalue))->tag->id);

	return retval;
}

/*****************************************************************************
*
*   DeltaSeqsToSeqLocs(dsp)
*   	converts a chain of delta seqs to seqlocs
*   	each SeqLit is converted to SeqLoc of type Int with a SeqId of type
*          Dbtag where db="Seq\tLit" and objectId.id which is the index of the
*          element in the delta seq chain where 1 is the first one.
*   	Returned SeqLoc is of type "mix" and must be freed by caller.
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL DeltaSeqsToSeqLocs (DeltaSeqPtr dsp)
{
	SeqLocPtr head = NULL, thead = NULL, last = NULL;
	DeltaSeqPtr curr;
	SeqInt si;
	Dbtag db;
	ObjectId oi;
	ValNode vn, vn2;

	MemSet(&vn, 0, sizeof(ValNode));
	MemSet(&vn2, 0, sizeof(ValNode));
	MemSet(&si, 0, sizeof(SeqInt));
	MemSet(&db, 0, sizeof(Dbtag));
	MemSet(&oi, 0, sizeof(ObjectId));
	vn.choice = SEQLOC_INT;
	vn.data.ptrvalue = &si;
	si.id = &vn2;
	vn2.choice = SEQID_GENERAL;
	vn2.data.ptrvalue = &db;
	db.db = seqlitdbtag;
	db.tag = &oi;
	oi.id = 1;

	
	
	for (curr = dsp; curr != NULL; curr = curr->next)
	{
		if (curr->choice == 1)   /* a SeqLoc */
			SeqLocAddEx (&thead, &last, (SeqLocPtr)(curr->data.ptrvalue), TRUE, TRUE);
		else
		{
			si.to = ((SeqLitPtr) (curr->data.ptrvalue))->length - 1;
			SeqLocAddEx (&thead, &last, &vn, TRUE, TRUE); 
		}
		oi.id++;
	}

	head = SeqLocPackage(thead);
	return head;
}

/*****************************************************************************
* GOHERE
*   SeqLocsToDeltaSeqs(dsp, slp)
*   	converts a chain of seqlocs	generated by DeltaSeqToSeqLocs() back into
*         delta seqs. dsp is the original chain of DeltaSeqs, which is required
*         to convert the delta seqs back.
*
*****************************************************************************/
NLM_EXTERN DeltaSeqPtr LIBCALL SeqLocsToDeltaSeqs (DeltaSeqPtr dsp, SeqLocPtr slp)
{
	DeltaSeqPtr dhead=NULL, dcurr=NULL, dtmp;
	SeqLocPtr scurr;
	Int4 ctr, index, strt, stp;
	SeqIdPtr sip;
	Uint1 strand, newcode;
	SeqLitPtr slitp, slitp_new;
	SeqPortPtr spps;
	ByteStorePtr bsp;
	Int2 residue;
	ValNode vn;

	if ((dsp == NULL) || (slp == NULL))
		return dhead;

	vn.choice = SEQLOC_MIX;
	vn.next = NULL;
	vn.data.ptrvalue = slp;
	scurr = NULL;
	while ((scurr = SeqLocFindNext(&vn, scurr)) != NULL)
	{
		dcurr = ValNodeNew(dhead);
		if (dhead == NULL)
			dhead = dcurr;

		index = ISADeltaSeqsToSeqLoc(scurr);

		if (index == 0)   /* just a SeqLoc */
		{
			dcurr->choice = 1;
			dcurr->data.ptrvalue = NULL;
			dcurr->data.ptrvalue = AsnIoMemCopy((Pointer)scurr, (AsnReadFunc)SeqLocAsnRead, (AsnWriteFunc)SeqLocAsnWrite);

		}
		else                                 /* convert to a delta seq */
		{
			dcurr->choice = 2;
			sip = SeqLocId(scurr);
			dtmp = dsp;
			for (ctr = 1; ctr < index; ctr++)
				dtmp = dtmp->next;

			if (dtmp->choice != 2)   /* wups */
			{
				ErrPostEx(SEV_ERROR,0,0,"Wrong type in SeqLocsToDeltaSeqs");
				dhead = DeltaSeqFree(dhead);
				return dhead;
			}
			slitp = (SeqLitPtr)(dtmp->data.ptrvalue);

			strt = SeqLocStart(scurr);
			stp = SeqLocStop(scurr);
			strand = SeqLocStrand(scurr);

			if ((strt == 0) && (stp == (slitp->length - 1)) && (strand != Seq_strand_minus))  /* no change */
			{
				dcurr->data.ptrvalue = AsnIoMemCopy((Pointer)slitp, (AsnReadFunc)SeqLitAsnRead, (AsnWriteFunc)SeqLitAsnWrite);
			}
			else   /* got to copy part of it */
			{
				switch (slitp->seq_data_type)
				{
					case Seq_code_iupacna:
					case Seq_code_iupacaa:
					case Seq_code_ncbi8na:
					case Seq_code_ncbi8aa:
					case Seq_code_ncbieaa:
					case Seq_code_ncbistdaa:
					case Seq_code_iupacaa3:
						newcode = slitp->seq_data_type;     /* one byte codes.. fine */
						break;
					case Seq_code_ncbipna:
						ErrPostEx(SEV_ERROR,0,0,"Converting from P residue codes");
						newcode = Seq_code_ncbieaa;
						break;
					case Seq_code_ncbipaa:
						ErrPostEx(SEV_ERROR,0,0,"Converting from P residue codes");
					case Seq_code_ncbi2na:
					case Seq_code_ncbi4na:
						newcode = Seq_code_iupacna;
						break;
					default:
						ErrPostEx(SEV_FATAL,0,0,"Unrecognized residue code [%d] in SeqLocsToDeltaSeqs",
							(int)(slitp->seq_data_type));
						return DeltaSeqFree(dhead);
				}
	 			spps = MemNew(sizeof(SeqPort));
				SeqPortSetUpFields (spps, strt, stp, strand, newcode);
				SeqPortSetUpAlphabet(spps, slitp->seq_data_type, newcode);
				spps->bp = slitp->seq_data;
				slitp_new = SeqLitNew();
				dcurr->data.ptrvalue = slitp_new;
				slitp_new->seq_data_type = newcode;
				slitp_new->length = (stp - strt + 1);
				bsp = BSNew(slitp_new->length);
				slitp_new->seq_data = bsp;
				SeqPortSeek(spps, 0, SEEK_SET);
				BSSeek(bsp, 0, SEEK_SET);
			    while (stp >= strt)
				{
					residue = SeqPortGetResidue(spps);
					BSPutByte(bsp, residue);
					strt++;
				}
				SeqPortFree(spps);
			}

		}

	}
	return dhead;
}
/*****************************************************************************
*
*   BioseqDelete (target, from, to, do_feat, do_split)
*      Deletes the region of sequence between from-to, inclusive, on the
*        Bioseq whose SeqId is target.
*      If do_feat, the feature table is updated to reflect the deletion
*        using SeqEntryDelFeat()
*      If do_split, the features across the deleted region are split into
*        two intervals on either side. If not, the feature is just shortened.
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqDelete (SeqIdPtr target, Int4 from, Int4 to, Boolean do_feat, Boolean do_split)
{
	Boolean retval = FALSE;
	BioseqPtr bsp;
	SeqLocPtr tmp, head;
	Int4 len, deleted;
	Int4 totlen, templen, tfrom, tto, diff1, diff2;
	SeqLocPtr slp, tloc, newhead, prev;
	ValNode vn;
	SeqInt si;
	SeqLocPtr PNTR newheadptr;
	SeqFeatPtr sfpcurr, sfpnext, sfpprev;
	Int2 dropped;
	SeqEntryPtr oldscope;
	DeltaSeqPtr tdsp;

	bsp = BioseqFind(target);
	if (bsp == NULL) {
		oldscope = SeqEntrySetScope (NULL);
		if (oldscope != NULL) {
			bsp = BioseqFind(target);
			SeqEntrySetScope (oldscope);
		}
	}
	if (bsp == NULL) return retval;

	if ((from < 0) || (from >= bsp->length) || (to < 0) ||
		(to >= bsp->length) || (from > to)) return retval;

	if (do_feat)
		SeqEntryDelFeat(NULL, target, from, to, do_split);

	len = to - from + 1;
	           /* if actual sequence present */

	if ((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const))
	{
		if (ISA_na(bsp->mol))
		{
			if (bsp->seq_data_type != Seq_code_iupacna)  /* need 1 byte/base */
				BioseqRawConvert(bsp, Seq_code_iupacna);
		}
		else
		{
			if (bsp->seq_data_type != Seq_code_ncbieaa)
				BioseqRawConvert(bsp, Seq_code_ncbieaa);
		}

		BSSeek(bsp->seq_data, from, SEEK_SET);
		deleted = BSDelete(bsp->seq_data, len);
		if (deleted != len)  /* error */
			ErrPost(CTX_NCBIOBJ, 1, "Delete of %ld residues failed", len);
		else
			retval = TRUE;
	}

			   /* update segmented sequence */
	if ((bsp->repr == Seq_repr_seg) || (bsp->repr == Seq_repr_delta))
	{
		head = ValNodeNew(NULL);  /* allocate to facilitate SeqLocFree */
		head->choice = SEQLOC_MIX;   /* make a SeqLoc out of the extension */
		if (bsp->repr == Seq_repr_seg)
			head->data.ptrvalue = bsp->seq_ext;
		else
		{
			tdsp = (DeltaSeqPtr)(bsp->seq_ext);
			head->data.ptrvalue = DeltaSeqsToSeqLocs(tdsp);
		}
		
		newhead = NULL;
		newheadptr = &newhead;

		tloc = &vn;
		MemSet((Pointer)tloc, 0, sizeof(ValNode));
		MemSet((Pointer)&si, 0, sizeof(SeqInt));
		tloc->choice = SEQLOC_INT;
		tloc->data.ptrvalue = (Pointer)(&si);
		
		slp = NULL;
		totlen = 0;
		while ((slp = SeqLocFindNext(head, slp)) != NULL)
		{
			templen = SeqLocLen(slp);
		    tfrom = SeqLocStart(slp);
			tto = SeqLocStop(slp);
			
			if (((totlen + templen - 1) < from) ||   /* before cut */
				(totlen > to))						  /* after cut */
				tmp = SeqLocAdd(newheadptr, slp, TRUE, TRUE); /* add whole SeqLoc */
			else                            
			{
				retval = 1;    /* will modify or drop interval */
		 		diff1 = from - totlen;        /* partial beginning? */
				diff2 = (templen + totlen - 1) - to;  /* partial end? */
				si.id = SeqLocId(slp);
				si.strand = SeqLocStrand(slp);
				
				if (diff1 > 0)	  /* partial start */
				{
					if (si.strand != Seq_strand_minus)
					{
					   si.from = tfrom;
					   si.to = tfrom + diff1 - 1;
					}
					else
					{
						si.from = tto - diff1 + 1;
						si.to = tto;
					}
					tmp = SeqLocAdd(newheadptr, tloc, TRUE, TRUE);
				}

				if (diff2 > 0)    /* partial end */
				{
					if (si.strand != Seq_strand_minus)
					{
					   si.from = tto - diff2 + 1;
					   si.to = tto;
					}
					else
					{
						si.from = tfrom;
						si.to = tfrom + diff2 - 1;
					}
					tmp = SeqLocAdd(newheadptr, tloc, TRUE, TRUE);
				}
				
			}
			totlen += templen;
		}

		prev = NULL;
		for (tmp = newhead; tmp != NULL; tmp = tmp->next)
		{
			if (tmp->next == NULL)   /* last one */
			{
				if (tmp->choice == SEQLOC_NULL)
				{
					if (prev != NULL)
						prev->next = NULL;
					else				  /* only a NULL left */
					{
						newhead = NULL;
					}
					MemFree(tmp);
					break;
				}
			}
			prev = tmp;
		}

		if (bsp->repr == Seq_repr_seg)
			bsp->seq_ext = newhead;
		else
		{
			bsp->seq_ext = SeqLocsToDeltaSeqs(tdsp, newhead);
			DeltaSeqSetFree(tdsp);
			SeqLocSetFree(newhead);
		}
		SeqLocFree(head);
		retval = TRUE;
	}

	if (bsp->repr == Seq_repr_map)      /* map bioseq */
	{
		sfpprev = NULL;
		sfpnext = NULL;
		sfpcurr = (SeqFeatPtr)(bsp->seq_ext);
		bsp->seq_ext = NULL;
		for (; sfpcurr != NULL; sfpcurr = sfpnext)
		{
			sfpnext = sfpcurr->next;
			dropped = SeqFeatDelete(sfpcurr, target, from, to, TRUE);
			if (dropped == 2)   /* completely gone */
			{
				SeqFeatFree(sfpcurr);
			}
			else
			{
				if (sfpprev == NULL)
					bsp->seq_ext = (Pointer)sfpcurr;
				else
					sfpprev->next = sfpcurr;
				sfpcurr->next = NULL;
				sfpprev = sfpcurr;
			}
		}
		retval = TRUE;
	}

	if (bsp->repr == Seq_repr_virtual)
		retval = TRUE;                 /* nothing to do */

	if (retval)
		bsp->length -= len;
	return retval;
}


/*****************************************************************************
*
*   BioseqOverwrite (target, pos, residue, seqcode)
*      Overwrites the residue at pos with residue in the
*        Bioseq whose SeqId is target.
*      residue is iupacna for DNA or ncbieaa for protein
*      target MUST be a raw Bioseq right now
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqOverwrite (SeqIdPtr target, Int4 pos, Uint1 residue)
{
	BioseqPtr bsp;
	Boolean retval = FALSE;
	SeqEntryPtr oldscope;


	bsp = BioseqFind(target);
	if (bsp == NULL) {
		oldscope = SeqEntrySetScope (NULL);
		if (oldscope != NULL) {
			bsp = BioseqFind(target);
			SeqEntrySetScope (oldscope);
		}
	}
	if (bsp == NULL) return retval;

	if ((pos < 0) || (pos >= bsp->length)) return retval;
	if (bsp->repr != Seq_repr_raw) return retval;

	if (ISA_na(bsp->mol))
	{
		if (bsp->seq_data_type != Seq_code_iupacna)  /* need 1 byte/base */
			BioseqRawConvert(bsp, Seq_code_iupacna);
	}
	else
	{
		if (bsp->seq_data_type != Seq_code_ncbieaa)
			BioseqRawConvert(bsp, Seq_code_ncbieaa);
	}

	BSSeek(bsp->seq_data, pos, SEEK_SET);
	BSPutByte(bsp->seq_data, (Int2)(TO_UPPER(residue)));
	retval = TRUE;

	return retval;
}


/*****************************************************************************
*
*   SeqInsertByLoc (target, offset, fragment)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqInsertByLoc (SeqIdPtr target, Int4 offset, SeqLocPtr fragment)
{
	return TRUE;
}


/*****************************************************************************
*
*   SeqDeleteByLoc (slp, do_feat, do_split)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqDeleteByLoc (SeqLocPtr slp, Boolean do_feat, Boolean do_split)
{
	SeqLocPtr tmp;
	Boolean retval = FALSE;
	Int2 numloc, i = 0, ctr, pick, totloc;
	SeqLocPtr PNTR locs, PNTR tlocs, PNTR theorder;
	BioseqPtr bsp;
	Int4 tstart, tstop;

	if (slp == NULL) return retval;

	numloc = 0;
	totloc = 0;
	locs = NULL;
	tmp = NULL;

	while ((tmp = SeqLocFindNext(slp, tmp)) != NULL)
	{
		switch (tmp->choice)
		{
			case SEQLOC_INT:
			case SEQLOC_PNT:
				if (BioseqFind(SeqLocId(tmp)) != NULL)
				{
					if (numloc == totloc)
					{
						tlocs = locs;
						locs = (SeqLocPtr PNTR)(MemNew((totloc+20) * sizeof(SeqLocPtr)));
						MemCopy(locs, tlocs, (size_t)(totloc * sizeof(SeqLocPtr)));
						MemFree(tlocs);
						totloc += 20;
					}
					locs[numloc] = tmp;
					numloc++;
				}
				break;
			default:
				Message(MSG_ERROR, "Unsupported Seqloc [%d] in SeqDeleteByLoc",
					(int)(tmp->choice));
				break;

		}
	}

	if (! numloc) return retval;

	              
				/***********************************************************
				*
				*   first gather all the seqlocs, grouped by Bioseq, and
				*   ordered from end to beginning. They must be ordered
				*   before the underlying Bioseq is changed.
				*
				***********************************************************/

	retval = TRUE;

	bsp = NULL;
	theorder = (SeqLocPtr PNTR)MemNew((sizeof(SeqLocPtr) * numloc));
	for (ctr = 0; ctr < numloc; ctr++)
	{
		pick = -1;   /* flag none found */
		if (bsp != NULL)
		{
			for (i = 0; i < numloc; i++)
			{
				if (locs[i] != NULL)
				{
				  	if (SeqIdIn(SeqLocId(locs[i]), bsp->id))
					{
						pick = i;
						i++;
						break;
					}
				}
			}
			if (pick < 0)
				bsp = NULL;   /* no more locs on this bioseq */
		}

		if (bsp == NULL)  /* have to find a new bioseq */
		{
			for (i = 0; i < numloc; i++)
			{
				if (locs[i] != NULL)
				{
					bsp = BioseqFind(SeqLocId(locs[i]));
					pick = i;
					i++;
					break;
				}
			}
		}

		while (i < numloc)
		{
			if (SeqLocOrder(locs[pick], locs[i], bsp) == (-1)) /* it's after */
				pick = i;
			i++;
		}

		theorder[ctr] = locs[pick];
		locs[pick] = NULL;
	}

	MemFree(locs);   /* finished with original list */

				/*************************************************************
				*
				*   Now do the actual deletions
				*
				*************************************************************/


	for (ctr = 0; ctr < numloc; ctr++)
	{
		tstart = SeqLocStart(theorder[ctr]);
		tstop = SeqLocStop(theorder[ctr]);
		BioseqDelete(SeqLocId(theorder[ctr]), tstart, tstop, do_feat, do_split);
	}

	MemFree(theorder);

	return retval;
}


/*****************************************************************************
*
*   SeqFeatDelete()
*     0 = no changes made to location or product
*     1 = changes made but feature still has some location
*     2 = all of sfp->location in deleted interval
*
*   if (merge)
*      1) correct numbers > to by subtraction
*      2) do not split intervals spanning the deletion
*   else
*      1) do not change numbers > to
*      2) split intervals which span the deletions
*
*****************************************************************************/
NLM_EXTERN Int2 LIBCALL SeqFeatDelete (SeqFeatPtr sfp, SeqIdPtr target, Int4 from, Int4 to, Boolean merge)
{
	ValNode      vn;
	SeqLocPtr    tloc;
	SeqInt       si;
	Boolean      changed = FALSE, tmpbool = FALSE;
	CdRegionPtr  crp;
	CodeBreakPtr cbp, prevcbp, nextcbp;
	RnaRefPtr    rrp;
	tRNAPtr      trp;
	Boolean      partial5, partial3;
	Uint1        strand;
	BioseqPtr    bsp;
	Int4         new_frame;

	tloc = &vn;
	MemSet((Pointer)tloc, 0, sizeof(ValNode));
	MemSet((Pointer)&si, 0, sizeof(SeqInt));
	tloc->choice = SEQLOC_INT;
	tloc->data.ptrvalue = (Pointer)(&si);
	si.id = target;
	si.from = from;
	si.to = to;

        CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
        strand = SeqLocStrand (sfp->location);
        bsp = BioseqFindFromSeqLoc (sfp->location);
	sfp->location = SeqLocDelete(sfp->location, target, from, to, merge, &changed);

	sfp->product = SeqLocDelete(sfp->product, target, from, to, merge, &changed);

	if (sfp->location == NULL)
		return 2;
	
	switch (sfp->data.choice)
	{
		case SEQFEAT_CDREGION:   /* cdregion */
			crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
      if (changed)
      {
		  	/* adjust frame */
			  if ((strand == Seq_strand_minus && bsp != NULL && to == bsp->length - 1 && partial5)
			    || (strand != Seq_strand_minus && from == 0 && partial5))
			  {
  	      if (crp->frame == 0)
			    {
			      crp->frame = 1;
			    }
			    new_frame = crp->frame - ((to - from + 1) % 3);
			    if (new_frame < 1)
			    {
			  	  new_frame += 3;
			    }
          crp->frame = new_frame;
			  }      
      }
			/* fix code_break locations */
			prevcbp = NULL;
			for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
			{
				nextcbp = cbp->next;
				cbp->loc = SeqLocDelete(cbp->loc, target, from, to, merge, &tmpbool);
				if (cbp->loc == NULL)
				{
					if (prevcbp != NULL)
						prevcbp->next = nextcbp;
					else
						crp->code_break = nextcbp;
					cbp->next = NULL;
					CodeBreakFree(cbp);
				}
				else
					prevcbp = cbp;
			}
			break;
		case SEQFEAT_RNA:
			rrp = (RnaRefPtr)(sfp->data.value.ptrvalue);
			if (rrp->ext.choice == 2)   /* tRNA */
			{
				trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
				if (trp->anticodon != NULL)
				{
					trp->anticodon = SeqLocDelete(trp->anticodon, target, from, to, merge, &tmpbool);
				}
			}
			break;
		default:
			break;
	}
			
	if (changed)
	{
		sfp->partial = TRUE;
		return 1;
	}
	else
		return 0;
}

/*****************************************************************************
*
*   SeqLocDelete()
*   	returns altered head or NULL if nothing left.
*   sets changed=TRUE if all or part of loc is deleted
*   does NOT set changed if location coordinates are only moved
*   if (merge) then corrects coordinates upstream of to
*   else
*     splits intervals covering from-to, does not correct upstream of to
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SeqLocDeleteEx (SeqLocPtr head, SeqIdPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed, BoolPtr partial5, BoolPtr partial3)
{
	SeqIntPtr sip, sip2;
	SeqPntPtr spp;
	PackSeqPntPtr pspp, pspp2;
	SeqBondPtr sbp;
	SeqIdPtr sidp;
	SeqLocPtr slp, tmp, prev, next, thead;
	Int4 diff, numpnt, i, tpos;
	BioseqPtr bsp;
	Boolean part5, part3, first;

	if ((head == NULL) || (target == NULL))
		return head;

	head->next = NULL;   /* caller maintains chains */
	diff = to - from + 1;
	
    switch (head->choice)
    {
        case SEQLOC_BOND:   /* bond -- 2 seqs */
			sbp = (SeqBondPtr)(head->data.ptrvalue);
			spp = sbp->a;
			if (SeqIdForSameBioseq(spp->id, target))
			{
				if (spp->point >= from)
				{
					if (spp->point <= to)   /* delete it */
					{
					    *changed = TRUE;
						sbp->a = SeqPntFree(spp);
					}
					else if (merge)
						spp->point -= diff;
				}
			}
			spp = sbp->b;
			if (spp != NULL)
			{
				if (SeqIdForSameBioseq(spp->id, target))
				{
					if (spp->point >= from)
					{
						if (spp->point <= to)   /* delete it */
						{
						    *changed = TRUE;
							sbp->b = SeqPntFree(spp);
						}
						else if (merge)
							spp->point -= diff;
					}
				}
			}
			if (sbp->a == NULL)
			{
				if (sbp->b != NULL)   /* only a required */
				{
					sbp->a = sbp->b;
					sbp->b = NULL;
				}
				else
				{
					head = SeqLocFree(head);
				}
			}
			break;
        case SEQLOC_FEAT:   /* feat -- can't track yet */
        case SEQLOC_NULL:    /* NULL */
        case SEQLOC_EMPTY:    /* empty */
			break;
        case SEQLOC_WHOLE:    /* whole */
			sidp = (SeqIdPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(sidp, target))
			{
				bsp = BioseqFind(target);
				if (bsp != NULL)           /* split it */
				{
					if ((from == 0) && (to >= (bsp->length - 1)))
					{					   /* complete delete */
						head = SeqLocFree(head);
						*changed = TRUE;
						break;
					}

					if (! merge)   /* split it up */
					{
						SeqIdFree(sidp);
						head->choice = SEQLOC_PACKED_INT;
						head->data.ptrvalue = NULL;
						slp = NULL;
						if (from != 0)
						{
							sip = SeqIntNew();
							sip->from = 0;
							sip->to = from - 1;
							sip->id = SeqIdDup(target);
							slp = ValNodeNew(NULL);
							slp->choice = SEQLOC_INT;
							slp->data.ptrvalue = sip;
							head->data.ptrvalue = slp;
							*changed = TRUE;
						}
						if (to < (bsp->length - 1))
						{
							sip = SeqIntNew();
							sip->from = to + 1;
							sip->to = bsp->length - 1;
							sip->id = SeqIdDup(target);
							tmp = ValNodeNew(NULL);
							tmp->choice = SEQLOC_INT;
							tmp->data.ptrvalue = sip;
							if (slp != NULL)
								slp->next = tmp;
							else
								head->data.ptrvalue = tmp;
							*changed = TRUE;
						}

					}
				}
			}
			break;
        case SEQLOC_MIX:    /* mix -- more than one seq */
        case SEQLOC_EQUIV:    /* equiv -- ditto */
        case SEQLOC_PACKED_INT:    /* packed int */
			prev = NULL;
			thead = NULL;
			part5 = FALSE;
			part3 = FALSE;
			first = TRUE;
			for (slp = (SeqLocPtr)(head->data.ptrvalue); slp != NULL; slp = next)
			{
				next = slp->next;
				tmp = SeqLocDeleteEx (slp, target, from, to, merge, changed, &part5, &part3);
				if (first) {
					if (partial5 != NULL) {
						*partial5 = part5;
					}
				}
				first = FALSE;
				if (tmp != NULL)
				{
					if (prev != NULL)
					{
						if ((merge) && (prev->choice == SEQLOC_INT) && (tmp->choice == SEQLOC_INT))
						{
							sip = (SeqIntPtr)(prev->data.ptrvalue);
							sip2 = (SeqIntPtr)(tmp->data.ptrvalue);

							if (SeqIdForSameBioseq(sip->id, sip2->id))
							{
							         /* merge intervals? */
								if ((sip->strand == Seq_strand_minus) &&
									(sip2->strand == Seq_strand_minus))
								{
									if (sip->from == (sip2->to + 1))
									{
										sip->from = sip2->from;
										sip->if_from = sip2->if_from;
										sip2->if_from = NULL;
										tmp = SeqLocFree(tmp);
									}
								}
								else if((sip->strand != Seq_strand_minus) &&
									(sip2->strand != Seq_strand_minus))
								{
									if (sip->to == (sip2->from - 1))
									{
										sip->to = sip2->to;
										sip->if_to = sip2->if_to;
										sip2->if_to = NULL;
										tmp = SeqLocFree(tmp);
									}
								}
							}
						}
						else if ((prev->choice == SEQLOC_NULL) && (tmp->choice == SEQLOC_NULL))
						{
							tmp = SeqLocFree(tmp);
							*changed = TRUE;
						}
					}
					else if (tmp->choice == SEQLOC_NULL)
					{
						tmp = SeqLocFree(tmp);
						*changed = TRUE;
					}

					if (tmp != NULL)   /* still have one? */
					{
						if (prev != NULL)
							prev->next = tmp;
						else
							thead = tmp;
						prev = tmp;
					}
				}
				else
					*changed = TRUE;
			}
			if (partial3 != NULL) {
				*partial3 = part3;
			}
			if (prev != NULL)
			{
				if (prev->choice == SEQLOC_NULL)  /* ends with NULL */
				{
					prev = NULL;
					for (slp = thead; slp->next != NULL; slp = slp->next)
						prev = slp;
					if (prev != NULL)
					{
						prev->next = NULL;
						SeqLocFree(slp);
					}
					else
					{
						thead = SeqLocFree(thead);
					}
					*changed = TRUE;
				}
			}
			head->data.ptrvalue = thead;
			if (thead == NULL)
				head = SeqLocFree(head);
            break;
        case SEQLOC_INT:    /* int */
			sip = (SeqIntPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(sip->id, target))
			{
				if (sip->to < from)  /* completely before cut */
					break;

								     /* completely contained in cut */
				if ((sip->from >= from) && (sip->to <= to))
				{
					head = SeqLocFree(head);
					*changed = TRUE;
					break;
				}

				if (sip->from > to)  /* completely past cut */
				{
					if (merge)
					{
						sip->from -= diff;
						sip->to -= diff;
					}
					break;
				}
									/* overlap here */

				if (sip->to > to)
				{
					if (merge)
						sip->to -= diff;
				}
				else                /* to inside cut, so partial delete */
				{
					sip->to = from - 1;
					*changed = TRUE;
					if (partial3 != NULL) {
						*partial3 = TRUE;
					}
				}
		
				if (sip->from >= from)   /* from inside cut, partial del */
				{
					*changed = TRUE;
					sip->from = to + 1;
					if (merge)
						sip->from -= diff;
					if (partial5 != NULL) {
						*partial5 = TRUE;
					}
				}

				if (merge)
					break;

						   /* interval spans cut.. only in non-merge */
				           /* have to split */

				if ((sip->from < from) && (sip->to > to))
				{
					*changed = TRUE;
					head->choice = SEQLOC_PACKED_INT;
					head->data.ptrvalue = NULL;
					tmp = ValNodeNew(NULL);
					tmp->choice = SEQLOC_INT;
					tmp->data.ptrvalue = sip;

					sip2 = SeqIntNew();
					sip2->from = to + 1;
					sip2->to = sip->to;
					sip2->strand = sip->strand;
					sip2->if_to = sip->if_to;
					sip2->id = SeqIdDup(target);
					slp = ValNodeNew(NULL);
					slp->choice = SEQLOC_INT;
					slp->data.ptrvalue = sip2;

					sip->if_to = NULL;
					sip->to = from - 1;

					if (sip->strand == Seq_strand_minus)
					{
						head->data.ptrvalue = slp;
						slp->next = tmp;
					}
					else
					{
						head->data.ptrvalue = tmp;
						tmp->next = slp;
					}

				}

			}
            break;
        case SEQLOC_PNT:    /* pnt */
			spp = (SeqPntPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(spp->id, target))
			{
				if ((spp->point >= from) && (spp->point <= to))
				{
					head = SeqLocFree(head);
					*changed = TRUE;
				}
				else if (spp->point > to)
				{
					if (merge)
						spp->point -= diff;
				}
			}
            break;
        case SEQLOC_PACKED_PNT:    /* packed pnt */
			pspp = (PackSeqPntPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(pspp->id, target))
			{
				numpnt = PackSeqPntNum(pspp);
				pspp2 = PackSeqPntNew();
				head->data.ptrvalue = pspp2;
				for (i = 0; i < numpnt; i++)
				{
					tpos = PackSeqPntGet(pspp, i);
					if (tpos < from)
						PackSeqPntPut(pspp2, tpos);
					else
					{
						if (tpos > to)
						{
							if (merge)
								tpos -= diff;
							PackSeqPntPut(pspp2, tpos);
						}
						else
							*changed = TRUE;
					}
				}
				pspp2->id = pspp->id;
				pspp->id = NULL;
				pspp2->fuzz = pspp->fuzz;
				pspp->fuzz = NULL;
				pspp2->strand = pspp->strand;
				PackSeqPntFree(pspp);
				numpnt = PackSeqPntNum(pspp2);
				if (! numpnt)
					head = SeqLocFree(head);

			}
            break;
        default:
            break;
    }

	return head;
}

NLM_EXTERN SeqLocPtr LIBCALL SeqLocDelete (SeqLocPtr head, SeqIdPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed)

{
  return SeqLocDeleteEx (head, target, from, to, merge, changed, NULL, NULL);
}

typedef struct delstruct {
	SeqIdPtr sip;
	Int4 from, to;
	Boolean merge;
} DelStruct, PNTR DelStructPtr;

NLM_EXTERN void DelFeat (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);

NLM_EXTERN void DelFeat (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	DelStructPtr dsp;
	BioseqPtr bsp;
	BioseqSetPtr bssp;
	SeqAnnotPtr sap, nextsap;
	SeqFeatPtr sfp, nextsfp;
	Pointer PNTR prevsap, PNTR prevsfp;

	dsp = (DelStructPtr)data;
	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		sap = bsp->annot;
		prevsap = (Pointer PNTR) &(bsp->annot);
	}
	else
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		sap = bssp->annot;
		prevsap = (Pointer PNTR) &(bssp->annot);
	}

	while (sap != NULL)
	{
		nextsap = sap->next;
		if (sap->type == 1)   /* feature table */
		{
			sfp = (SeqFeatPtr) sap->data;
			prevsfp = (Pointer PNTR) &(sap->data);
			while (sfp != NULL)
			{
				nextsfp = sfp->next;
				if (SeqFeatDelete(sfp, dsp->sip, dsp->from, dsp->to, dsp->merge) == 2)
				{
					/* location completely gone */
					*(prevsfp) = sfp->next;
					sfp->next = NULL;
					SeqFeatFree(sfp);
				} else {
					prevsfp = (Pointer PNTR) &(sfp->next);
				}
				sfp = nextsfp;
			}
		}

		if (sap->data == NULL)  /* all features deleted */
		{
			*(prevsap) = sap->next;
			sap->next = NULL;
			SeqAnnotFree (sap);
		} else {
			prevsap = (Pointer PNTR) &(sap->next);
		}

		sap = nextsap;
	}

	return;
}

/*****************************************************************************
*
*   SeqEntryDelFeat(sep, id, from, to, do_split)
*   	Deletes or truncates features on Bioseq (id) in the range
*       from-to, inclusive
*       
*		Moves features > to left to account for decrease in length
*       if do_split, breaks intervals across the deletion
*       else just reduces their size
*
*       If sep == NULL, then calls SeqEntryFind(id) to set scope to look
*       for features.
*   
*****************************************************************************/
NLM_EXTERN Boolean	LIBCALL SeqEntryDelFeat (SeqEntryPtr sep, SeqIdPtr sip, Int4 from, Int4 to, Boolean do_split)
{

	DelStruct ds;

	if (sip == NULL)
		return FALSE;

	if (sep == NULL)
		sep	= SeqEntryFind(sip);

	if (sep == NULL) return FALSE;

	ds.sip = sip;
	ds.from = from;
	ds.to = to;
	if (do_split)
		ds.merge = FALSE;
	else
		ds.merge = TRUE;

	SeqEntryExplore(sep, (Pointer)(&ds), DelFeat);

	return TRUE;
}

/*****************************************************************************
*
*   DescrToFeatures(sep)
*   	Moves all Seqdescr to features in sep where possible
*
*****************************************************************************/

static DeltaSeqPtr CopyDeltaSeqPtrChain (DeltaSeqPtr dsp)
{
  DeltaSeqPtr new_chain = NULL;
  SeqLocPtr   slp_orig, slp_new;
  SeqLitPtr   slip_orig, slip_new;
  
  while (dsp != NULL) {
    if (dsp->choice == 1) {
      slp_orig = (SeqLocPtr) dsp->data.ptrvalue;
      slp_new = AsnIoMemCopy (slp_orig, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
      ValNodeAddPointer (&new_chain, 1, slp_new);
    }
    else if (dsp->choice ==2) 
    {
      slip_orig = (SeqLitPtr) dsp->data.ptrvalue;
      slip_new = AsnIoMemCopy(slip_orig, (AsnReadFunc) SeqLitAsnRead, (AsnWriteFunc) SeqLitAsnWrite);
      ValNodeAddPointer (&new_chain, 2, slip_new);
    }
    dsp = dsp->next;
  }
  
  return new_chain;  
}

/*****************************************************************************
*
*   BioseqCopy(newid, sourceid, from, to, strand, do_feat)
*      Creates a new Bioseq from sourceid in the range from-to inclusive.
*      If strand==Seq_strand_minus, reverse complements the sequence in
*        the copy and (if do_feat) corrects the feature table appropriately.
*      Names new Bioseq as newid, if not NULL
*        else Creates seqid.local = "Clipboard" if newid is NULL
*      If do_feat == TRUE copies appropriate region of feature table from
*        sourceid to new copy using SeqFeatsCopy().
*
*****************************************************************************/
NLM_EXTERN BioseqPtr LIBCALL BioseqCopyEx (SeqIdPtr newid, BioseqPtr oldbsp, Int4 from, Int4 to,
                               Uint1 strand, Boolean do_feat)
{
	BioseqPtr newbsp=NULL, tmpbsp;
	SeqPortPtr spp=NULL;
	ByteStorePtr bsp;
	Uint1 seqtype;
	ValNodePtr tmp;
	ObjectIdPtr oid;
	Int4 len, i;
	Int2 residue;
	ValNode fake;
	SeqLocPtr the_segs, head, curr;
	Boolean handled = FALSE, split;
	SeqFeatPtr sfp, newsfp, lastsfp;
	DeltaSeqPtr dsp;
	SeqEntryPtr oldscope;


	if ((oldbsp == NULL) || (from < 0)) return FALSE;

	len = to - from + 1;
	if (len <= 0) return NULL;

	newbsp = BioseqNew();
	if (newid != NULL)
		newbsp->id = SeqIdDup(newid);
	else
	{
		tmp = ValNodeNew(NULL);
		tmp->choice = SEQID_LOCAL;
		oid = ObjectIdNew();
		tmp->data.ptrvalue = (Pointer)oid;
		oid->str = StringSave("Clipboard");
		tmpbsp = BioseqFind(tmp);   /* old clipboard present? */
		if (tmpbsp == NULL) {
			oldscope = SeqEntrySetScope (NULL);
			if (oldscope != NULL) {
				tmpbsp = BioseqFind(tmp);
				SeqEntrySetScope (oldscope);
			}
		}
		if (tmpbsp != NULL)
			BioseqFree(tmpbsp);
		newbsp->id = tmp;
	}

	newbsp->repr = oldbsp->repr;
	newbsp->mol = oldbsp->mol;
	newbsp->length = len;
	newbsp->seq_ext_type = oldbsp->seq_ext_type;

	if (newbsp->repr == Seq_repr_virtual)
		handled = TRUE;               /* no more to do */

	if ((newbsp->repr == Seq_repr_raw) ||
		(newbsp->repr == Seq_repr_const))
	{
		if (ISA_aa(newbsp->mol))
		{
			seqtype = Seq_code_ncbieaa;
		}
		else
		{
			seqtype = Seq_code_iupacna;
		}
		newbsp->seq_data_type = seqtype;
		bsp = BSNew(len);
		if (bsp == NULL) goto erret;

		newbsp->seq_data = bsp;
		spp = SeqPortNew(oldbsp, from, to, strand, seqtype);
		if (spp == NULL) goto erret;

		for (i = 0; i < len; i++)
		{
			residue = SeqPortGetResidue(spp);
			if (! IS_residue(residue)) goto erret;
			BSPutByte(bsp, residue);
		}

		SeqPortFree(spp);
		handled = TRUE;
	}

	if ((newbsp->repr == Seq_repr_seg) ||
		(newbsp->repr == Seq_repr_ref) ||
		(newbsp->repr == Seq_repr_delta))
	{
        if (newbsp->repr == Seq_repr_seg)  /* segmented */
		{
			fake.choice = SEQLOC_MIX;   /* make SEQUENCE OF Seq-loc, into one */
			fake.data.ptrvalue = oldbsp->seq_ext;
			fake.next = NULL;
			the_segs = (SeqLocPtr)&fake;
			head = SeqLocCopyPart (the_segs, from, to, strand, FALSE, NULL, NULL);
		}
		else if (newbsp->repr == Seq_repr_ref)  /* reference: is a Seq-loc */
		{
			head = SeqLocCopyPart ((SeqLocPtr)(oldbsp->seq_ext), from, to,
			                         strand, TRUE, NULL, NULL);
		}
		else if (newbsp->repr == Seq_repr_delta)
		{
			dsp = (DeltaSeqPtr)(oldbsp->seq_ext);  /* real data is here */
			
			head = CopyDeltaSeqPtrChain (dsp);
		}

        newbsp->seq_ext = (Pointer)head;
		handled = TRUE;
	}

	if (newbsp->repr == Seq_repr_map)
	{
		lastsfp = NULL;
		for (sfp = (SeqFeatPtr)(oldbsp->seq_ext); sfp != NULL; sfp = sfp->next)
		{
			split = FALSE;
			curr = SeqLocCopyRegion(newbsp->id, sfp->location, oldbsp, from, to, strand, &split);
			if (curr != NULL)   /* got one */
			{
				newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
				SeqLocFree(newsfp->location);
				newsfp->location = curr;
				if (split)
					newsfp->partial = TRUE;
				if (lastsfp == NULL)  /* first one */
					newbsp->seq_ext = (Pointer)newsfp;
				else
					lastsfp->next = newsfp;
				lastsfp = newsfp;
			}
		}
		handled = TRUE;
	}
	

	if (! handled) goto erret;

	               /* get descriptors */
	               /* get features */

	if (do_feat)
		SeqFeatsCopy (newbsp, oldbsp, from, to, strand);

	return newbsp;

erret:
	BioseqFree(newbsp);
	SeqPortFree(spp);
	return NULL;
}

NLM_EXTERN BioseqPtr LIBCALL BioseqCopy (SeqIdPtr newid, SeqIdPtr sourceid, Int4 from, Int4 to,
                               Uint1 strand, Boolean do_feat)
{
	BioseqPtr oldbsp;
	SeqEntryPtr oldscope;

	if ((sourceid == NULL) || (from < 0)) return FALSE;

	oldbsp = BioseqFind(sourceid);
	if (oldbsp == NULL) {
		oldscope = SeqEntrySetScope (NULL);
		if (oldscope != NULL) {
			oldbsp = BioseqFind(sourceid);
			SeqEntrySetScope (oldscope);
		}
	}
	if (oldbsp == NULL) return NULL;

	return BioseqCopyEx (newid, oldbsp, from, to, strand, do_feat);
}

/*****************************************************************************
*
*	SeqLocCopyPart (the_segs, from, to, strand, group, first_segp, last_segp)
*      cuts out from the_segs the part from offset from to offset to
*      reverse complements resulting seqloc if strand == Seq_strand_minus
*      if (group) puts resulting intervals into a new Seq-loc (of type
*        PACKED_INT if no SEQLOC_NULL, else SEQLOC_MIX).
*      Currently this always makes intervals or nulls. Is really for segmented and
*        reference sequence extensions
*      If first_segp and last_segp are not NULL, then they are filled in with the
*        ordinal number of the source segments that remain in the copy, based
*        on SeqLocFindNext, where 1 is the first one. Thus if the third and
*        fourth segments were copied, first is 3 and last is 4. If the
*        location was reverse complemented, first is 4 and last is 3.
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SeqLocCopyPart (SeqLocPtr the_segs, Int4 from, Int4 to, Uint1 strand,
					Boolean group, Int2Ptr first_segp, Int2Ptr last_segp)
{
	SeqLocPtr currseg, newhead, head, prev, curr, last;
	Int2 numloc, first_seg = 0, last_seg = 0, seg_ctr = 0;
	Int4 oldpos, tlen, tfrom, tto, tstart, tstop, xfrom, xto;
	Uint1 tstrand;
	SeqIdPtr tid;
        SeqIntPtr sip;
	Boolean done, started, wasa_null, hada_null;
	BioseqPtr bsp;

	if (the_segs == NULL) return NULL;
	if ((from < 0) || (to < 0)) return NULL;

   currseg = NULL;
   oldpos = 0;	   /* position in old sequence */
   done = FALSE;
   started = FALSE;
   head = NULL;
   prev = NULL;
   numloc = 0;
   wasa_null = FALSE;
   hada_null = FALSE;
   while ((oldpos <= to) && ((currseg = SeqLocFindNext(the_segs, currseg)) != NULL))
   {
		seg_ctr++;
   		tlen = SeqLocLen(currseg);
   		tid = SeqLocId(currseg);
   		if (tlen < 0) {
   			bsp = BioseqLockById (tid);  /* only necessary for locations of type WHOLE */
   			tlen = SeqLocLen (currseg);
   			BioseqUnlock (bsp);
   		}
	   	tstrand = SeqLocStrand(currseg);
   		tfrom = SeqLocStart(currseg);
	   	tto = SeqLocStop(currseg);

	   	if (! started)
   		{
			wasa_null = FALSE;
   			if (((oldpos + tlen - 1) >= from) &&
   				(currseg->choice != SEQLOC_NULL))
	   		{
   				tstart = from - oldpos;
   				started = TRUE;
				first_seg = seg_ctr;
	   		}
   			else
   				tstart = -1;
	   	}
   		else
		{
			if (currseg->choice == SEQLOC_NULL)
			{
				wasa_null = TRUE;
				tstart = -1;  /* skip it till later */
			}
			else
				tstart = 0;
		}

	   	if (tstart >= 0)   /* have a start */
   		{
   			if ((oldpos + tlen - 1) >= to)
	   		{
   				done = TRUE;   /* hit the end */
   				tstop = ((oldpos + tlen - 1) - to);
	   		}
   			else
   				tstop = 0;

	   		if (tstrand == Seq_strand_minus)
   			{
   				xfrom = tfrom + tstop;
   				xto = tto - tstart;
	   		}
   			else
   			{
   				xfrom = tfrom + tstart;
	   			xto = tto - tstop;
   			}

   			sip = SeqIntNew();
	   		sip->id = SeqIdDup(tid);
   			sip->strand = tstrand;
   			sip->from = xfrom;
	   		sip->to = xto;
			if (wasa_null)  /* previous SEQLOC_NULL */
			{
				curr = ValNodeAddInt(&head, SEQLOC_NULL, 0);
				numloc++;
				wasa_null = FALSE;
				hada_null = TRUE;
			}
   			curr = ValNodeAddPointer(&head, SEQLOC_INT, (Pointer)sip);
   			numloc++;
			last_seg = seg_ctr;
	   	}

   		oldpos += tlen;
   }

   if (strand == Seq_strand_minus)  /* reverse order and complement */
   {
	   	newhead = NULL;
   		last = NULL;
	   	while (head != NULL)
   		{
   			prev = NULL;
	   		for (curr = head; curr->next != NULL; curr = curr->next)
   				prev = curr;
   			if (prev != NULL)
   				prev->next = NULL;
	   		else
   				head = NULL;

   			if (newhead == NULL)
   				newhead = curr;
	   		else
   				last->next = curr;
   			last = curr;
			if (curr->choice == SEQLOC_INT)
			{
				sip = (SeqIntPtr)(curr->data.ptrvalue);
				sip->strand = StrandCmp(sip->strand);
			}
	   	}

   		head = newhead;
		seg_ctr = last_seg;
		last_seg = first_seg;
		first_seg = seg_ctr;
   }

   if ((numloc) && (group))
   {
   		curr = ValNodeNew(NULL);
		if (hada_null)
			curr->choice = SEQLOC_MIX;
		else
			curr->choice = SEQLOC_PACKED_INT;
	   	curr->data.ptrvalue = (Pointer)head;
   		head = curr;
   }

   if (first_segp != NULL)
	 *first_segp = first_seg;
   if (last_segp != NULL)
	 *last_segp = last_seg;

   return head;
}

/*****************************************************************************
*
*   SeqFeatCopy(new, old, from, to, strand)
*
*****************************************************************************/
static Int2 LIBCALL IndexedSeqFeatsCopy (BioseqPtr newbsp, BioseqPtr oldbsp, Int4 from, Int4 to, Uint1 strand)

{
	Int2 ctr=0;
	SeqFeatPtr sfp, last=NULL, newsfp;
	SeqInt si;
	ValNode vn;
	ValNodePtr region;
	SeqLocPtr newloc;
	Boolean split = FALSE;
	SeqAnnotPtr sap = NULL, saptmp;
	CdRegionPtr crp;
	CodeBreakPtr cbp, prevcbp, nextcbp;
	RnaRefPtr rrp;
	tRNAPtr trp;
	SeqMgrFeatContext fcontext;

	region = &vn;
	vn.choice = SEQLOC_INT;
	vn.data.ptrvalue = (Pointer)(&si);
	si.from = from;
	si.to = to;
	si.id = oldbsp->id;
	si.if_from = NULL;
	si.if_to = NULL;

	sfp = NULL;
	while ((sfp = SeqMgrGetNextFeature (oldbsp, sfp, 0, 0, &fcontext)) != NULL)
	{
		/* can exit once past rightmost limit */
		if (fcontext.left > to) return ctr;

		if (fcontext.right >= from && fcontext.left <= to) {

			split = FALSE;
			newloc = SeqLocCopyRegion(newbsp->id, sfp->location, oldbsp, from, to, strand, &split);
			if (newloc != NULL)   /* got one */
			{
				newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
				SeqLocFree(newsfp->location);
				newsfp->location = newloc;
				if (split)
					newsfp->partial = TRUE;
				if (last == NULL)  /* first one */
				{
					sap = SeqAnnotNew();
					if (newbsp->annot == NULL)
						newbsp->annot = sap;
					else
					{
						for (saptmp = newbsp->annot; saptmp->next != NULL; saptmp = saptmp->next)
							continue;
						saptmp->next = sap;
					}
					sap->type = 1;   /* feature table */
					sap->data = (Pointer)newsfp;
				}
				else
					last->next = newsfp;
				last = newsfp;

				switch (newsfp->data.choice)
				{
					case SEQFEAT_CDREGION:   /* cdregion */
						crp = (CdRegionPtr)(newsfp->data.value.ptrvalue);
						prevcbp = NULL;
						for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
						{
							nextcbp = cbp->next;
							cbp->loc = SeqLocCopyRegion(newbsp->id, cbp->loc, oldbsp, from, to, strand, &split);
							if (cbp->loc == NULL)
							{
								if (prevcbp != NULL)
									prevcbp->next = nextcbp;
								else
									crp->code_break = nextcbp;
								cbp->next = NULL;
								CodeBreakFree(cbp);
							}
							else
								prevcbp = cbp;
						}
						break;
					case SEQFEAT_RNA:
						rrp = (RnaRefPtr)(newsfp->data.value.ptrvalue);
						if (rrp->ext.choice == 2)   /* tRNA */
						{
							trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
							if (trp->anticodon != NULL)
							{
								trp->anticodon = SeqLocCopyRegion(newbsp->id, trp->anticodon, oldbsp, from, to, strand, &split);
							}
						}
						break;
					default:
						break;
				}
			}
		}
		
	}
	return ctr;
}

NLM_EXTERN Int2 LIBCALL SeqFeatsCopy (BioseqPtr newbsp, BioseqPtr oldbsp, Int4 from, Int4 to, Uint1 strand)
{
	Int2 ctr=0;
	BioseqContextPtr bcp = NULL;
	SeqFeatPtr sfp, last=NULL, newsfp;
	SeqInt si;
	ValNode vn;
	ValNodePtr region;
	SeqLocPtr newloc;
	Boolean split = FALSE;
	SeqAnnotPtr sap = NULL, saptmp;
	CdRegionPtr crp;
	CodeBreakPtr cbp, prevcbp, nextcbp;
	RnaRefPtr rrp;
	tRNAPtr trp;
	Uint2 entityID;

	if (oldbsp == NULL) return ctr;

	entityID = ObjMgrGetEntityIDForPointer (oldbsp);
	if (entityID > 0 && SeqMgrFeaturesAreIndexed (entityID)) {
		/* indexed version should be much faster */
		return IndexedSeqFeatsCopy (newbsp, oldbsp, from, to, strand);
	}

	bcp = BioseqContextNew(oldbsp);
	if (bcp == NULL) return ctr;

	region = &vn;
	vn.choice = SEQLOC_INT;
	vn.data.ptrvalue = (Pointer)(&si);
	si.from = from;
	si.to = to;
	si.id = oldbsp->id;
	si.if_from = NULL;
	si.if_to = NULL;

	sfp = NULL;
	while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 0)) != NULL)
	{
		split = FALSE;
		newloc = SeqLocCopyRegion(newbsp->id, sfp->location, oldbsp, from, to, strand, &split);
		if (newloc != NULL)   /* got one */
		{
			newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
			SeqLocFree(newsfp->location);
			newsfp->location = newloc;
			if (split)
				newsfp->partial = TRUE;
			if (last == NULL)  /* first one */
			{
				sap = SeqAnnotNew();
				if (newbsp->annot == NULL)
					newbsp->annot = sap;
				else
				{
					for (saptmp = newbsp->annot; saptmp->next != NULL; saptmp = saptmp->next)
						continue;
					saptmp->next = sap;
				}
				sap->type = 1;   /* feature table */
				sap->data = (Pointer)newsfp;
			}
			else
				last->next = newsfp;
			last = newsfp;

			switch (newsfp->data.choice)
			{
				case SEQFEAT_CDREGION:   /* cdregion */
					crp = (CdRegionPtr)(newsfp->data.value.ptrvalue);
					prevcbp = NULL;
					for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
					{
						nextcbp = cbp->next;
						cbp->loc = SeqLocCopyRegion(newbsp->id, cbp->loc, oldbsp, from, to, strand, &split);
						if (cbp->loc == NULL)
						{
							if (prevcbp != NULL)
								prevcbp->next = nextcbp;
							else
								crp->code_break = nextcbp;
							cbp->next = NULL;
							CodeBreakFree(cbp);
						}
						else
							prevcbp = cbp;
					}
					break;
				case SEQFEAT_RNA:
					rrp = (RnaRefPtr)(newsfp->data.value.ptrvalue);
					if (rrp->ext.choice == 2)   /* tRNA */
					{
						trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
						if (trp->anticodon != NULL)
						{
							trp->anticodon = SeqLocCopyRegion(newbsp->id, trp->anticodon, oldbsp, from, to, strand, &split);
						}
					}
					break;
				default:
					break;
			}
		}
		
	}
	BioseqContextFree (bcp);
	return ctr;
}


NLM_EXTERN SeqLocPtr LIBCALL SeqLocCopyRegion(SeqIdPtr newid, SeqLocPtr head, BioseqPtr oldbsp,
	Int4 from, Int4 to, Uint1 strand, BoolPtr split)
{
	SeqLocPtr newhead = NULL, tmp, slp, prev, next, thead;
	SeqIntPtr sip, sip2;
	SeqPntPtr spp, spp2;
	PackSeqPntPtr pspp, pspp2;
	SeqBondPtr sbp, sbp2;
	SeqIdPtr sidp, oldids;
	Int4 numpnt, i, tpos, len, intcnt, othercnt;
	Boolean dropped_one;
	IntFuzzPtr ifp;
	ValNode vn;
	
	if ((head == NULL) || (oldbsp == NULL)) return NULL;

	oldids = oldbsp->id;
	len = to - from + 1;
	switch (head->choice)
	{
        case SEQLOC_BOND:   /* bond -- 2 seqs */
			sbp2 = NULL;
			sbp = (SeqBondPtr)(head->data.ptrvalue);
			vn.choice = SEQLOC_PNT;
			vn.data.ptrvalue = sbp->a;
			vn.next = NULL;
			tmp = SeqLocCopyRegion(newid, (SeqLocPtr)(&vn), oldbsp, from, to, strand, split);
			if (tmp != NULL)
			{
			 	sbp2 = SeqBondNew();
				sbp2->a = (SeqPntPtr)(tmp->data.ptrvalue);
				MemFree(tmp);
			}
			if (sbp->b != NULL)
			{
				vn.data.ptrvalue = sbp->b;
				tmp = SeqLocCopyRegion(newid, (SeqLocPtr)(&vn), oldbsp, from, to, strand, split);
				if (tmp != NULL)
				{
					if (sbp2 == NULL)
					{
					 	sbp2 = SeqBondNew();
						sbp2->a = (SeqPntPtr)(tmp->data.ptrvalue);
					}
					else
						sbp2->b = (SeqPntPtr)(tmp->data.ptrvalue);
					MemFree(tmp);
				}
			}
			if (sbp2 != NULL)
			{
				newhead = ValNodeNew(NULL);
				newhead->choice = SEQLOC_BOND;
				newhead->data.ptrvalue = sbp2;
				if ((sbp->b != NULL) && (sbp2->b == NULL))
					*split = TRUE;
			}
			break;
        case SEQLOC_FEAT:   /* feat -- can't track yet */
        case SEQLOC_NULL:    /* NULL */
        case SEQLOC_EMPTY:    /* empty */
			break;
        case SEQLOC_WHOLE:    /* whole */
			sidp = (SeqIdPtr)(head->data.ptrvalue);
			if (SeqIdIn(sidp, oldids))
			{
				if ((from != 0) || (to != (oldbsp->length - 1)))
				{
					*split = TRUE;
				}
				newhead = ValNodeNew(NULL);
				sip2 = SeqIntNew();
				sip2->id = SeqIdDup(newid);
				sip2->from = 0;
				sip2->to = to - from;
				newhead->choice = SEQLOC_INT;
				newhead->data.ptrvalue = (Pointer)sip2;
				if (strand == Seq_strand_minus)
				{
				  sip2->strand = Seq_strand_minus;					
				}
				else if (sip2->strand == Seq_strand_minus)
				{
				  sip2->strand = strand;
				}
			}
			break;
        case SEQLOC_EQUIV:    /* does it stay equiv? */
        case SEQLOC_MIX:    /* mix -- more than one seq */
        case SEQLOC_PACKED_INT:    /* packed int */
			prev = NULL;
			thead = NULL;
			dropped_one = FALSE;
			for (slp = (SeqLocPtr)(head->data.ptrvalue); slp != NULL; slp = next)
			{
				next = slp->next;
				tmp = SeqLocCopyRegion(newid, slp, oldbsp, from, to, strand, split);
				if (tmp != NULL)
				{
					if (prev != NULL)
					{
						if ((prev->choice == SEQLOC_INT) && (tmp->choice == SEQLOC_INT))
						{
							sip = (SeqIntPtr)(prev->data.ptrvalue);
							sip2 = (SeqIntPtr)(tmp->data.ptrvalue);

							if ((sip->strand == Seq_strand_minus) &&
								(sip2->strand == Seq_strand_minus))
							{
								if (sip->from == (sip2->to + 1))
								{
									sip->from = sip2->from;
									sip->if_from = sip2->if_from;
									sip2->if_from = NULL;
									tmp = SeqLocFree(tmp);
								}
							}
							else if((sip->strand != Seq_strand_minus) &&
								(sip2->strand != Seq_strand_minus))
							{
								if (sip->to == (sip2->from - 1))
								{
									sip->to = sip2->to;
									sip->if_to = sip2->if_to;
									sip2->if_to = NULL;
									tmp = SeqLocFree(tmp);
								}
							}
						}
						else if ((prev->choice == SEQLOC_NULL) && (tmp->choice == SEQLOC_NULL))
						{
							tmp = SeqLocFree(tmp);
							dropped_one = TRUE;
						}
					}
					else if (tmp->choice == SEQLOC_NULL)
					{
						tmp = SeqLocFree(tmp);
						dropped_one = TRUE;
					}

					if (tmp != NULL)   /* still have one? */
					{
						if (prev != NULL)
							prev->next = tmp;
						else
							thead = tmp;
						prev = tmp;
					}
					else
						dropped_one = TRUE;
				}
				else
					dropped_one = TRUE;
			}
			if (prev != NULL)
			{
				if (prev->choice == SEQLOC_NULL)  /* ends with NULL */
				{
					prev = NULL;
					for (slp = thead; slp->next != NULL; slp = slp->next)
						prev = slp;
					if (prev != NULL)
					{
						prev->next = NULL;
						SeqLocFree(slp);
					}
					else
					{
						thead = SeqLocFree(thead);
					}
					dropped_one = TRUE;
				}
			}
			if (thead != NULL)
			{
				if (dropped_one)
					*split = TRUE;
				intcnt = 0;
				othercnt = 0;
				for (slp = thead; slp != NULL; slp = slp->next)
				{
					if (slp->choice == SEQLOC_INT)
						intcnt++;
					else
						othercnt++;
				}
				if ((intcnt + othercnt) > 1)
				{
					newhead = ValNodeNew(NULL);
					if (head->choice == SEQLOC_EQUIV)
						newhead->choice = SEQLOC_EQUIV;
					else
					{
						if (othercnt == 0)
							newhead->choice = SEQLOC_PACKED_INT;
						else
							newhead->choice = SEQLOC_MIX;
					}

					newhead->data.ptrvalue = (Pointer)thead;
				}
				else                 /* only one SeqLoc left */
					newhead = thead;

			}
            break;
        case SEQLOC_INT:    /* int */
			sip = (SeqIntPtr)(head->data.ptrvalue);
			if (SeqIdIn(sip->id, oldids))
			{
				if (sip->to < from)  /* completely before cut */
					break;
				if (sip->from > to)  /* completely after cut */
					break;

				sip2 = SeqIntNew();
				sip2->id = SeqIdDup(newid);
				sip2->strand = sip->strand;

				if (sip->to > to)
				{
					sip2->to = to;
					*split = TRUE;
					ifp = IntFuzzNew();
					ifp->choice = 4;   /* lim */
					ifp->a = 1;        /* greater than */
					sip2->if_to = ifp;
				}
				else
				{
					sip2->to = sip->to;
					if (sip->if_to != NULL)
					{
						ifp = IntFuzzNew();
						MemCopy((Pointer)ifp, (Pointer)(sip->if_to), sizeof(IntFuzz));
						sip2->if_to = ifp;
					}
				}

				if (sip->from < from)
				{
					sip2->from = from;
					*split = TRUE;
					ifp = IntFuzzNew();
					ifp->choice = 4;   /* lim */
					ifp->a = 2;        /* less than */
					sip2->if_from = ifp;
				}
				else
				{
					sip2->from = sip->from;
					if (sip->if_from != NULL)
					{
						ifp = IntFuzzNew();
						MemCopy((Pointer)ifp, (Pointer)(sip->if_from), sizeof(IntFuzz));
						sip2->if_from = ifp;
					}
				}
									  /* set to region coordinates */
				sip2->from -= from;
				sip2->to -= from;
				IntFuzzClip(sip2->if_from, from, to, strand, split);
				IntFuzzClip(sip2->if_to, from, to, strand, split);

				if (strand == Seq_strand_minus)  /* rev comp */
				{
					sip2->strand = StrandCmp(sip2->strand);
					tpos = len - sip2->from - 1;
					sip2->from = len - sip2->to - 1;
					sip2->to = tpos;
					      /* IntFuzz already complemented by IntFuzzClip */
					      /* just switch order */
					ifp = sip2->if_from;
					sip2->if_from = sip2->if_to;
					sip2->if_to = ifp;
				}

				newhead = ValNodeNew(NULL);
				newhead->choice = SEQLOC_INT;
				newhead->data.ptrvalue = (Pointer)sip2;
			}
            break;
        case SEQLOC_PNT:    /* pnt */
			spp = (SeqPntPtr)(head->data.ptrvalue);
			if (SeqIdIn(spp->id, oldids))
			{
				if ((spp->point >= from) && (spp->point <= to))
				{
					spp2 = SeqPntNew();
					spp2->id = SeqIdDup(newid);
					spp2->point = spp->point - from;
					spp2->strand = spp->strand;
					if (spp->fuzz != NULL)
					{
						ifp = IntFuzzNew();
						spp2->fuzz = ifp;
						MemCopy((Pointer)ifp, (Pointer)spp->fuzz, sizeof(IntFuzz));
						IntFuzzClip(ifp, from, to, strand, split);
					}
					if (strand == Seq_strand_minus)
					{
						spp2->point = len - spp2->point - 1;
						spp2->strand = StrandCmp(spp->strand);
					}
					else if (spp2->strand == Seq_strand_minus)
					{
					  spp2->strand = strand;
					}
					newhead = ValNodeNew(NULL);
					newhead->choice = SEQLOC_PNT;
					newhead->data.ptrvalue = (Pointer)spp2;
				}
			}
            break;
        case SEQLOC_PACKED_PNT:    /* packed pnt */
			pspp = (PackSeqPntPtr)(head->data.ptrvalue);
			if (SeqIdIn(pspp->id, oldids))
			{
				numpnt = PackSeqPntNum(pspp);
				pspp2 = PackSeqPntNew();
				pspp2->strand = pspp->strand;
				intcnt = 0;	     /* use for included points */
				othercnt = 0;	 /* use for exclued points */
				for (i = 0; i < numpnt; i++)
				{
					tpos = PackSeqPntGet(pspp, i);
					if ((tpos < from) || (tpos > to))
					{
						othercnt++;
					}
					else
					{
						intcnt++;
						PackSeqPntPut(pspp2, tpos - from);
					}
				}
				if (! intcnt)  /* no points in region */
				{
					PackSeqPntFree(pspp2);
					break;
				}
				if (othercnt)
					*split = TRUE;
				if (pspp->fuzz != NULL)
				{
					ifp = IntFuzzNew();
					MemCopy((Pointer)ifp, (Pointer)(pspp->fuzz), sizeof(IntFuzz));
				}
				else
					ifp = NULL;

				if (strand == Seq_strand_minus)  /* rev comp */
				{
					IntFuzzClip(ifp, from, to, strand, split);
					pspp = pspp2;
					pspp2 = PackSeqPntNew();
					pspp2->strand = StrandCmp(pspp->strand);
					numpnt = PackSeqPntNum(pspp);
					numpnt--;
					for (i = numpnt; i >= 0; i--)	 /* reverse order */
					{
						tpos = PackSeqPntGet(pspp, i);
						PackSeqPntPut(pspp2, (len - tpos - 1));
					}
					PackSeqPntFree(pspp);
				}
				else if (pspp2->strand == Seq_strand_minus)
				{
					pspp2->strand = strand;
				}
				pspp2->id = SeqIdDup(newid);
				pspp2->fuzz = ifp;

				newhead = ValNodeNew(NULL);
				newhead->choice = SEQLOC_PACKED_PNT;
				newhead->data.ptrvalue = (Pointer)pspp2;

			}
            break;
        default:
            break;

	}
	return newhead;
}

/*****************************************************************************
*
*   IntFuzzClip()
*   	returns TRUE if clipped range values
*       in all cases, adjusts and/or complements IntFuzz
*       Designed for IntFuzz on SeqLocs
*
*****************************************************************************/
NLM_EXTERN void LIBCALL IntFuzzClip(IntFuzzPtr ifp, Int4 from, Int4 to, Uint1 strand, BoolPtr split)
{
	Int4 len, tmp;

	if (ifp == NULL) return;
	len = to - from + 1;
	switch (ifp->choice)
	{
		case 1:      /* plus/minus - no changes */
		case 3:      /* percent - no changes */
			break;
		case 2:      /* range */
			if (ifp->a > to)     /* max */
			{
				*split = TRUE;
				ifp->a = to;
			}
			if (ifp->a < from) 
			{
				*split = TRUE;
				ifp->a = from;
			}
			if (ifp->b > to)     /* min */
			{
				*split = TRUE;
				ifp->b = to;
			}
			if (ifp->b < from) 
			{
				*split = TRUE;
				ifp->b = from;
			}
			ifp->a -= from;     /* adjust to window */
			ifp->b -= to;
			if (strand == Seq_strand_minus)
			{
				tmp = len - ifp->a;   /* reverse/complement */
				ifp->a = len - ifp->b;
				ifp->b = tmp;
			}
			break;
		case 4:     /* lim */
			if (strand == Seq_strand_minus)  /* reverse/complement */
			{
				switch (ifp->a)
				{
					case 1:    /* greater than */
						ifp->a = 2;
						break;
					case 2:    /* less than */
						ifp->a = 1;
						break;
					case 3:    /* to right of residue */
						ifp->a = 4;
						break;
					case 4:    /* to left of residue */
						ifp->a = 3;
						break;
					default:
						break;
				}
			}
			break;
	}
	return;
}

extern void 
AdjustFeaturesForInsertion 
(BioseqPtr tobsp, 
 SeqIdPtr  to_id,
 Int4 pos, 
 Int4 len, 
 Boolean do_split)
{
  Uint2             entityID;
  SeqFeatPtr        sfp;
  CdRegionPtr       crp;
  CodeBreakPtr      cbp, prevcbp, nextcbp;
  RnaRefPtr         rrp;
  tRNAPtr           trp;
	SeqMgrFeatContext fcontext;
	ValNodePtr        prods, vnp;
	BioseqContextPtr  bcp;
  
  if (tobsp == NULL || to_id == NULL)
  {
    return;
  }
  
	entityID = ObjMgrGetEntityIDForPointer (tobsp);
	if (entityID > 0 && SeqMgrFeaturesAreIndexed (entityID)) {
    sfp = NULL;
		while ((sfp = SeqMgrGetNextFeature (tobsp, sfp, 0, 0, &fcontext)) != NULL)
		{
			sfp->location = SeqLocInsert (sfp->location, to_id,pos, len, do_split, NULL);
			switch (sfp->data.choice)
			{
				case SEQFEAT_CDREGION:   /* cdregion */
					crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
				  prevcbp = NULL;
				  for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
					{
						nextcbp = cbp->next;
						cbp->loc = SeqLocInsert (cbp->loc, to_id,pos, len, do_split, NULL);
						if (cbp->loc == NULL)
						{
							if (prevcbp != NULL)
								prevcbp->next = nextcbp;
							else
								crp->code_break = nextcbp;
							cbp->next = NULL;
							CodeBreakFree (cbp);
						}
						else
							prevcbp = cbp;
					}
					break;
				case SEQFEAT_RNA:
					rrp = (RnaRefPtr)(sfp->data.value.ptrvalue);
					if (rrp->ext.choice == 2)   /* tRNA */
					{
						trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
						if (trp->anticodon != NULL)
						{
							trp->anticodon = SeqLocInsert (trp->anticodon, to_id,pos, len, do_split, NULL);
						}
					}
					break;
				default:
					break;
			}
		}

		/* adjust features pointing by product */
		prods = SeqMgrGetSfpProductList (tobsp);
		for (vnp = prods; vnp != NULL; vnp = vnp->next) {
			sfp = (SeqFeatPtr) vnp->data.ptrvalue;
			if (sfp == NULL) continue;
			sfp->product = SeqLocInsert (sfp->product, to_id,pos, len, do_split, NULL);
		}

	} else {
		bcp = BioseqContextNew(tobsp);
		sfp = NULL;
	  /* adjust features pointing by location */
		while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 0)) != NULL)
		{
			sfp->location = SeqLocInsert(sfp->location, to_id,pos, len, do_split, NULL);
			switch (sfp->data.choice)
			{
				case SEQFEAT_CDREGION:   /* cdregion */
					crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
					prevcbp = NULL;
					for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
					{
						nextcbp = cbp->next;
						cbp->loc = SeqLocInsert(cbp->loc, to_id,pos, len, do_split, NULL);
						if (cbp->loc == NULL)
						{
							if (prevcbp != NULL)
								prevcbp->next = nextcbp;
							else
								crp->code_break = nextcbp;
							cbp->next = NULL;
							CodeBreakFree(cbp);
						}
						else
							prevcbp = cbp;
					}
					break;
				case SEQFEAT_RNA:
					rrp = (RnaRefPtr)(sfp->data.value.ptrvalue);
					if (rrp->ext.choice == 2)   /* tRNA */
					{
						trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
						if (trp->anticodon != NULL)
						{
							trp->anticodon = SeqLocInsert(trp->anticodon, to_id,pos, len, do_split, NULL);
						}
					}
					break;
				default:
					break;
			}
		}

		sfp = NULL;
	  /* adjust features pointing by product */
		while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 1)) != NULL)
			sfp->product = SeqLocInsert(sfp->product, to_id,pos, len, do_split, NULL);
		BioseqContextFree(bcp);
	}
}

/*****************************************************************************
*
* BioseqInsert (from_id, from, to, strand, to_id, pos, from_feat, to_feat,
*                                                                  do_split)
*   	Inserts a copy the region "from"-"to" on "strand" of the Bioseq
*          identified by "from_id" into the Bioseq identified by "to_id" 
*          before "pos".
*       if from_feat = TRUE, copies the feature table from "from" and updates
*          to locations to point to the proper residues in "to_id"
*       If to_feat = TRUE, updates feature table on "to_id" as well.
*          if do_split == TRUE, then splits features in "to_id" (to_feat must
*             be TRUE as well). Otherwise expands features at insertion.
*
*       All operations are copies. "frombsp" is unchanged.
*       Insert will only occur between certain Bioseq.repr classes as below
*
*   From Bioseq.repr                      To Bioseq.repr
*   
*					      virtual       raw      segmented        map
*					   +---------------------------------------------------
*	         virtual   |   length	    inst	  SeqLoc		 length
*					   +---------------------------------------------------
*				 raw   |   error        copy      SeqLoc         error
*					   +---------------------------------------------------
*		   segmented   |   error        inst      SeqLoc*        error
*					   +---------------------------------------------------
*				 map   |   error        inst*     SeqLoc         copy
*					   +---------------------------------------------------
*
*   length = changes length of "to" by length of "from"
*   error  = insertion not allowed
*   inst   = "from" instantiated as residues ("N" or "X" for virtual "from")
*   inst*  = as above, but a restriction map can instantiate other bases
*            than "N" for known restriction recognition sites.
*   copy   = copy of "from" inserted into "to"
*   SeqLoc = a SeqLoc added to "to" which points to "from". No copy of residues.
*   SeqLoc* = as above, but note that "to" points to "from" directly, not
*             what "from" itself may point to.
*   
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioseqInsert (SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos,
			Boolean from_feat, Boolean to_feat, Boolean do_split)
{
	BioseqPtr tobsp, frombsp;
	Int4 len, i, ctr, tlen;
	Boolean from_type, to_type;
	Uint1 seqtype;
	SeqAnnotPtr sap, newsap;
	SeqFeatPtr sfp, newsfp, prevsfp, sfphead = NULL;
	BioseqContextPtr bcp;
	Boolean handled = FALSE;
	SeqPortPtr spp;
	Int2 residue;
	Boolean split, added = FALSE, do_bsadd = TRUE;
	SeqLocPtr newloc, curr, head, tloc, xloc, yloc, fake;
	SeqIntPtr sip;
	CdRegionPtr crp;
	CodeBreakPtr cbp, prevcbp, nextcbp;
	RnaRefPtr rrp;
	tRNAPtr trp;
	SeqEntryPtr oldscope;

	if ((from_id == NULL) || (to_id == NULL)) return FALSE;

	tobsp = BioseqFind(to_id);
	if (tobsp == NULL) {
		oldscope = SeqEntrySetScope (NULL);
		if (oldscope != NULL) {
			tobsp = BioseqFind(to_id);
			SeqEntrySetScope (oldscope);
		}
	}
	if (tobsp == NULL) return FALSE;

	len = BioseqGetLen(tobsp);

	if (pos == LAST_RESIDUE)
		pos = len - 1;
	else if (pos == APPEND_RESIDUE) {
		pos = len;
	}

	if ((pos < 0) || (pos > len)) return FALSE;

	frombsp = BioseqFind(from_id);
	if (frombsp == NULL) {
		oldscope = SeqEntrySetScope (NULL);
		if (oldscope != NULL) {
			frombsp = BioseqFind(from_id);
			SeqEntrySetScope (oldscope);
		}
	}
	if (frombsp == NULL) return FALSE;
	
	from_type = ISA_na(frombsp->mol);
	to_type = ISA_na(tobsp->mol);

	if (from_type != to_type) return FALSE;

	len = BioseqGetLen(frombsp);
	if (to == LAST_RESIDUE)
		to = len - 1;
	
	if ((from < 0) || (to >= len)) return FALSE;

	len = to - from + 1;

	if (tobsp->repr == Seq_repr_virtual)
	{
		if (frombsp->repr != Seq_repr_virtual)
			return FALSE;

		handled = TRUE;                    /* just length and features */
	}

 	if ((tobsp->repr == Seq_repr_raw) || (tobsp->repr == Seq_repr_const))
	{
		if (ISA_na(tobsp->mol))
		{
			seqtype = Seq_code_iupacna;
		}
		else
		{
			seqtype = Seq_code_ncbieaa;
		}

		if (tobsp->seq_data_type != seqtype)
			BioseqRawConvert(tobsp, seqtype);
		BSSeek(tobsp->seq_data, pos, SEEK_SET);
		if (do_bsadd) {
			Nlm_BSAdd(tobsp->seq_data, len, FALSE);
		}

		i = 0;

		spp = SeqPortNew(frombsp, from, to, strand, seqtype);
		while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{
			if (! IS_residue(residue))
			{
				ErrPost(CTX_NCBIOBJ, 1, "Non-residue in BioseqInsert [%d]",
					(int)residue);
			}
			else
			{
				BSPutByte(tobsp->seq_data, residue);
				i++;
			}
		}
		SeqPortFree(spp);

		if (i != len)
		{
			ErrPost(CTX_NCBIOBJ, 1, "Tried to insert %ld residues but %ld went in",
				len, i);
			return FALSE;
		}

		handled = TRUE;
	}

	if ((tobsp->repr == Seq_repr_seg) || (tobsp->repr == Seq_repr_ref))
	{
		sip = SeqIntNew();
		sip->id = SeqIdDup(from_id);
		sip->from = from;
		sip->to = to;
		sip->strand = strand;
		tloc = ValNodeNew(NULL);
		tloc->choice = SEQLOC_INT;
		tloc->data.ptrvalue = (Pointer)sip;
		head = NULL;
		if (tobsp->repr == Seq_repr_seg)
		{
			fake = ValNodeNew(NULL);
			fake->choice = SEQLOC_MIX;
			fake->data.ptrvalue = (Pointer)(tobsp->seq_ext);
		}
		else
			fake = (SeqLocPtr)(tobsp->seq_ext);
		curr = NULL;
		ctr = 0;
		while ((curr = SeqLocFindNext(fake, curr)) != NULL)
		{
			if ((! added) && (ctr == pos))
			{
				newloc = SeqLocAdd(&head, tloc, TRUE, TRUE);
				added = TRUE;
			}
			tlen = SeqLocLen(curr);
			if ((! added) && ((ctr + tlen) > pos))  /* split interval */
			{
				yloc = NULL;
				xloc = SeqLocAdd(&yloc, curr, TRUE, TRUE);
				i = (pos - ctr) + SeqLocStart(curr);
			    newloc = SeqLocInsert(xloc, SeqLocId(xloc), i, 0, TRUE, NULL);
				xloc = newloc;
				yloc = newloc->next;
				SeqLocAdd(&head, xloc, TRUE, TRUE);
				SeqLocAdd(&head, tloc, TRUE, TRUE);
				SeqLocAdd(&head, yloc, TRUE, TRUE);
				SeqLocFree(xloc);
				SeqLocFree(yloc);
				added = TRUE;
			}
			else
				newloc = SeqLocAdd(&head, curr, TRUE, TRUE);
			ctr += tlen;
		}
		if ((! added) && (ctr == pos))
		{
			newloc = SeqLocAdd(&head, tloc, TRUE, TRUE);
			added = TRUE;
		}
		SeqLocFree(tloc);
		SeqLocFree(fake);
		if (tobsp->repr == Seq_repr_seg)
		{
			tobsp->seq_ext = (Pointer)head;
		}
		else
		{
			tobsp->seq_ext = SeqLocPackage(head);
		}
		handled = TRUE;
	}

	if (tobsp->repr == Seq_repr_map)
	{
		if (! ((frombsp->repr == Seq_repr_map) || (frombsp->repr == Seq_repr_virtual)))
			return FALSE;

		prevsfp = NULL;
		for (sfp = (SeqFeatPtr)(tobsp->seq_ext); sfp != NULL; sfp = sfp->next)
		{
			sfp->location = SeqLocInsert(sfp->location, to_id, pos, len, TRUE, NULL);
			prevsfp = sfp;
		}

		if (frombsp->repr == Seq_repr_map)
		{
			for (sfp = (SeqFeatPtr)(frombsp->seq_ext); sfp != NULL; sfp = sfp->next)
			{
				split = FALSE;
				newloc = SeqLocCopyRegion(to_id, sfp->location, frombsp, from, to, strand, &split);
				if (newloc != NULL)   /* got one */
				{
					newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
					SeqLocFree(newsfp->location);
					newsfp->location = newloc;
					if (split)
						newsfp->partial = TRUE;
					
					if (prevsfp == NULL)
						tobsp->seq_ext = (Pointer)newsfp;
					else
						prevsfp->next = newsfp;
					prevsfp = newsfp;
		
					newsfp->location = SeqLocInsert(newsfp->location, to_id, 0,
				                                          pos, TRUE, to_id);
				}
			}
		}
		handled = TRUE;
	}

	if (! handled) return FALSE;

	tobsp->length += len;

	if (to_feat)		     /* fix up sourceid Bioseq feature table(s) */
	{
    AdjustFeaturesForInsertion (tobsp, to_id, pos, len, do_split);
	}

	if (from_feat)				/* add source Bioseq features to sourceid */
	{
		bcp = BioseqContextNew(frombsp);
		sfp = NULL;					/* NOTE: should make NEW feature table */
		prevsfp = NULL;
	                            /* is there an old feature table to use? */
		for (newsap = tobsp->annot; newsap != NULL; newsap = newsap->next)
		{
			if (newsap->type == 1)  /* feature table */
				break;
		}
		if (newsap != NULL)
		{							/* create a new one if necessary */
			for (prevsfp = (SeqFeatPtr)(newsap->data); prevsfp != NULL;
			                                            prevsfp = prevsfp->next)
			{
				if (prevsfp->next == NULL)
					break;
			}
		}
		                                     /* get features by location */
		while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 0)) != NULL)
		{									/* copy all old features */
			split = FALSE;
			newloc = SeqLocCopyRegion(to_id, sfp->location, frombsp, from, to, strand, &split);
			if (newloc != NULL)   /* got one */
			{
				newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
				SeqLocFree(newsfp->location);
				newsfp->location = newloc;

				if (split)
					newsfp->partial = TRUE;

				if (prevsfp == NULL)
					sfphead = newsfp;
				else
					prevsfp->next = newsfp;
				prevsfp = newsfp;

				newsfp->location = SeqLocInsert(newsfp->location, to_id, 0,
				                                          pos, TRUE, to_id);
				switch (newsfp->data.choice)
				{
					case SEQFEAT_CDREGION:   /* cdregion */
						crp = (CdRegionPtr)(newsfp->data.value.ptrvalue);
						prevcbp = NULL;
						for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
						{
							nextcbp = cbp->next;
							cbp->loc = SeqLocCopyRegion(to_id, cbp->loc, frombsp, from, to, strand, &split);
							if (cbp->loc == NULL)
							{
								if (prevcbp != NULL)
									prevcbp->next = nextcbp;
								else
									crp->code_break = nextcbp;
								cbp->next = NULL;
								CodeBreakFree(cbp);
							}
							else
							{
								cbp->loc = SeqLocInsert(cbp->loc, to_id, 0,
				                                          pos, TRUE, to_id);
								prevcbp = cbp;
							}
						}
						break;
					case SEQFEAT_RNA:
						rrp = (RnaRefPtr)(newsfp->data.value.ptrvalue);
						if (rrp->ext.choice == 2)   /* tRNA */
						{
							trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
							if (trp->anticodon != NULL)
							{
								trp->anticodon = SeqLocCopyRegion(to_id, trp->anticodon, frombsp, from, to, strand, &split);
								trp->anticodon = SeqLocInsert(trp->anticodon, to_id, 0,
				                                          pos, TRUE, to_id);
							}
						}
						break;
					default:
						break;
				}
			}
		}

		sfp = NULL;
								/* get features by product */
		while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 1)) != NULL)
		{									/* copy all old features */
			split = FALSE;
			newloc = SeqLocCopyRegion(to_id, sfp->product, frombsp, from, to, strand, &split);
			if (newloc != NULL)   /* got one */
			{
				newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
				SeqLocFree(newsfp->product);
				newsfp->product = newloc;
				if (split)
					newsfp->partial = TRUE;

				if (prevsfp == NULL)
					sfphead = newsfp;
				else
					prevsfp->next = newsfp;
				prevsfp = newsfp;

				newsfp->product = SeqLocInsert(newsfp->product, to_id, 0, pos,
				                                                TRUE, to_id);
			}
		}
		BioseqContextFree(bcp);
	

		if (sfphead != NULL)    /* orphan chain of seqfeats to attach */
		{
			if (newsap == NULL)
			{
				for (sap = tobsp->annot; sap != NULL; sap = sap->next)
				{
					if (sap->next == NULL)
						break;
				}
				newsap = SeqAnnotNew();
				newsap->type = 1;
				if (sap == NULL)
					tobsp->annot = newsap;
				else
					sap->next = newsap;
			}

			newsap->data = (Pointer)sfphead;
		}
	}

	return TRUE;
}

/*****************************************************************************
*
*   SeqLocInsert()
*       alters "head" by insert "len" residues before "pos" in any SeqLoc
*         on the Bioseq "target"
*       all SeqLocs not on "target" are unaltered
*       for SeqLocs on "target"
*          all SeqLocs before "pos" are unaltered
*          all SeqLocs >= "pos" are incremented by "len"
*          all SeqLocs spanning "pos"
*             if "split" == TRUE, are split into two SeqLocs, one to the
*               left of the insertion, the other to right
*             if "split" != TRUE, the SeqLoc is increased in length to cover
*               the insertion
*   	returns altered head or NULL if nothing left.
*       if ("newid" != NULL) replaces "target" with "newid" whether the
*          SeqLoc is altered on not.
*
*       Usage hints:
*          1) To update a feature location on "target" when 10 residues of
*               sequence have been inserted before position 5
*          SeqFeatPtr->location = SeqLocInsert ( SeqFeatPtr->location ,
*                "target", 5, 10, TRUE, NULL);  [for some feature types
*                      you may want "split" equal FALSE]
*          2) To insert the complete feature table from "source" into a
*                different Bioseq "dest" before position 20 in "dest"
*          SFP->location = SeqLocInsert(SFP->location, "source", 0, 20,
*                FALSE, "dest");
*   
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SeqLocInsert (SeqLocPtr head, SeqIdPtr target, Int4 pos, Int4 len,
                                               Boolean split, SeqIdPtr newid)
{
	SeqIntPtr sip, sip2;
	SeqPntPtr spp;
	PackSeqPntPtr pspp, pspp2;
	SeqBondPtr sbp;
	SeqLocPtr slp, tmp, prev, next, thead, tmp2;
	Int4 diff, numpnt, i, tpos;
	Uint1 oldchoice;
	ValNode vn;
	SeqIdPtr sidp;

	if ((head == NULL) || (target == NULL))
		return head;

	head->next = NULL;   /* caller maintains chains */

	diff = len;

    switch (head->choice)
    {
        case SEQLOC_BOND:   /* bond -- 2 seqs */
			vn.next = NULL;
			vn.choice = SEQLOC_PNT;

			sbp = (SeqBondPtr)(head->data.ptrvalue);
			vn.data.ptrvalue = (Pointer)(sbp->a);
			SeqLocInsert(&vn, target, pos, len, split, newid);
			sbp->a = (SeqPntPtr)(vn.data.ptrvalue);
			if (sbp->b != NULL)
			{
				vn.data.ptrvalue = (Pointer)(sbp->b);
				SeqLocInsert(&vn, target, pos, len, split, newid);
				sbp->b = (SeqPntPtr)(vn.data.ptrvalue);
			}
			break;
        case SEQLOC_FEAT:   /* feat -- can't track yet */
        case SEQLOC_NULL:    /* NULL */
			break;
        case SEQLOC_EMPTY:    /* empty */
        case SEQLOC_WHOLE:    /* whole */
			if (newid != NULL)
			{
				sidp = (SeqIdPtr)(head->data.ptrvalue);
				if (SeqIdForSameBioseq(sidp, target))
				{
					SeqIdFree(sidp);
					sidp = SeqIdDup(newid);
					head->data.ptrvalue = (Pointer)sidp;
				}
			}
			break;
        case SEQLOC_MIX:    /* mix -- more than one seq */
        case SEQLOC_EQUIV:    /* equiv -- ditto */
        case SEQLOC_PACKED_INT:    /* packed int */
			prev = NULL;
			thead = NULL;
			for (slp = (SeqLocPtr)(head->data.ptrvalue); slp != NULL; slp = next)
			{
				next = slp->next;
				oldchoice = slp->choice;
				tmp = SeqLocInsert(slp, target, pos, len, split, newid);
				if (tmp != NULL)
				{
					if ((head->choice != SEQLOC_EQUIV) &&
						(oldchoice != tmp->choice))  /* split interval? */
					{
						if ((oldchoice == SEQLOC_INT) &&
							(tmp->choice == SEQLOC_PACKED_INT))
						{
							tmp2 = tmp;
							tmp = (SeqLocPtr)(tmp2->data.ptrvalue);
							MemFree(tmp2);
							while (tmp->next != NULL)
							{
								if (prev != NULL)
									prev->next = tmp;
								else
									thead = tmp;
								prev = tmp;
								tmp = tmp->next;
							}
						}
					}
					if (prev != NULL)
						prev->next = tmp;
					else
						thead = tmp;
					prev = tmp;
				}
			}
			head->data.ptrvalue = thead;
			if (thead == NULL)
				head = SeqLocFree(head);
            break;
        case SEQLOC_INT:    /* int */
			sip = (SeqIntPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(sip->id, target))
			{
				if (newid != NULL)   /* change id? */
				{
					SeqIdFree(sip->id);
					sip->id = SeqIdDup(newid);
				}

				if (sip->to < pos)  /* completely before insertion */
				{
					break;
				}

				if ((! split) || (sip->from >= pos))  /* interval unbroken */
				{
					if (sip->from >= pos)
						sip->from += len;
					sip->to += len;
					break;
				}

						                          /* split interval */
				sip2 = SeqIntNew();
				slp = ValNodeNew(NULL);
				slp->choice = SEQLOC_INT;
				slp->data.ptrvalue = (Pointer)sip2;
				sip2->strand = sip->strand;
				sip2->id = SeqIdDup(sip->id);

				sip2->to = sip->to + len;
				sip2->from = pos + len;
				sip2->if_to = sip->if_to;
				sip->if_to = NULL;
				sip->to = pos - 1;
				head->next = slp;

				if (sip->strand == Seq_strand_minus)  /* reverse order */
				{
					head->data.ptrvalue = (Pointer)sip2;
					slp->data.ptrvalue = (Pointer)sip;
				}

				thead = head;   /* make split interval into PACKED_INT */
				head = ValNodeNew(NULL);
				head->choice = SEQLOC_PACKED_INT;
				head->data.ptrvalue = thead;

			}
            break;
        case SEQLOC_PNT:    /* pnt */
			spp = (SeqPntPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(spp->id, target))
			{
				if (newid != NULL)   /* change id? */
				{
					SeqIdFree(spp->id);
					spp->id = SeqIdDup(newid);
				}

				if (spp->point >= pos)
					spp->point += len;
			}
            break;
        case SEQLOC_PACKED_PNT:    /* packed pnt */
			pspp = (PackSeqPntPtr)(head->data.ptrvalue);
			if (SeqIdForSameBioseq(pspp->id, target))
			{
				if (newid != NULL)   /* change id? */
				{
					SeqIdFree(pspp->id);
					pspp->id = SeqIdDup(newid);
				}

				numpnt = PackSeqPntNum(pspp);
				pspp2 = PackSeqPntNew();
				head->data.ptrvalue = pspp2;
				for (i = 0; i < numpnt; i++)
				{
					tpos = PackSeqPntGet(pspp, i);
					if (tpos >= pos)
						tpos += len;
					PackSeqPntPut(pspp2, tpos);
				}
				pspp2->id = pspp->id;
				pspp->id = NULL;
				pspp2->fuzz = pspp->fuzz;
				pspp->fuzz = NULL;
				pspp2->strand = pspp->strand;
				PackSeqPntFree(pspp);
			}
            break;
        default:
            break;
    }

	if (head == NULL)
		ErrPost(CTX_NCBIOBJ, 1, "SeqLocInsert: lost a SeqLoc");

	return head;
}

/*****************************************************************************
*
*   SeqLocSubtract (SeqLocPtr head, SeqLocPtr piece)
*   	Deletes piece from head.
*       head may be changed.
*       returns the changed head.
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SeqLocSubtract (SeqLocPtr head, SeqLocPtr piece)
{
	SeqLocPtr slp = NULL;
	SeqIdPtr sip;
	Int4 from, to;
	Boolean changed = FALSE;

	if ((head == NULL) || (piece == NULL))
		return NULL;

	while ((slp = SeqLocFindNext(piece, slp)) != NULL)
	{
		sip = SeqLocId(slp);
		from = SeqLocStart(slp);
		to = SeqLocStop(slp);
		head = SeqLocDelete(head, sip, from, to, FALSE, &changed);
	}

	return head;
}

/********************************************************************
*
* SeqLocReplaceID
*   replaces the Seq-Id in a Seq-Loc (slp) with a new Seq-Id (new_sip)
*
**********************************************************************/
NLM_EXTERN SeqLocPtr SeqLocReplaceID (SeqLocPtr slp, SeqIdPtr new_sip)
{
  SeqLocPtr        curr;
  PackSeqPntPtr    pspp;
  SeqIntPtr        target_sit;
  SeqBondPtr       sbp;
  SeqPntPtr        spp;

  switch (slp->choice) {
     case SEQLOC_PACKED_INT :
     case SEQLOC_MIX :
     case SEQLOC_EQUIV :
        curr = NULL;
        while ((curr = SeqLocFindNext (slp, curr)) != NULL) {
           curr = SeqLocReplaceID (curr, new_sip);
        }
        break;
     case SEQLOC_PACKED_PNT :
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          SeqIdFree (pspp->id);
          pspp->id = SeqIdDup (new_sip);
        }
        break;
     case SEQLOC_EMPTY :
     case SEQLOC_WHOLE :
        SeqIdFree ((SeqIdPtr) slp->data.ptrvalue);
        slp->data.ptrvalue = (Pointer) SeqIdDup (new_sip);
        break;
     case SEQLOC_INT :
        target_sit = (SeqIntPtr) slp->data.ptrvalue;
        SeqIdFree (target_sit->id);
        target_sit->id = SeqIdDup (new_sip);
        break;
     case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        SeqIdFree(spp->id);
        spp->id = SeqIdDup(new_sip);
        break;
     case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp == NULL || sbp->a == NULL || sbp->b == NULL) break;
        /* only do this if both ends bonded to same Seq-id */
        if (SeqIdMatch (sbp->a->id, sbp->b->id)) {
          spp = sbp->a;
          SeqIdFree(spp->id);
          spp->id = SeqIdDup(new_sip);
          spp = sbp->b;
          SeqIdFree(spp->id);
          spp->id = SeqIdDup(new_sip);
        }
        break;
     default :
        break;
  }
  return slp;
}

/**********************************************************
 *
 *   NLM_EXTERN SeqLocPtr LIBCALL GapToSeqLoc(range):
 *
 *      Gets the size of gap and constructs SeqLoc block with
 *   $(seqlitdbtag) value as Dbtag.db and Dbtag.tag.id = 0.
 *
 **********************************************************/
NLM_EXTERN SeqLocPtr LIBCALL GapToSeqLoc(Int4 range)
{
    SeqLocPtr slp;
    SeqIntPtr sip;
    SeqIdPtr  sidp;
    DbtagPtr  dp;

    if(range < 0)
        return(NULL);

    slp = ValNodeNew(NULL);
    if(range == 0)
    {
        slp->choice = SEQLOC_NULL;
        slp->data.ptrvalue = NULL;
        slp->next = NULL;
        return(slp);
    }

    dp = DbtagNew();
    dp->db = StringSave(seqlitdbtag);
    dp->tag = ObjectIdNew();
    dp->tag->id = 0;
    dp->tag->str = NULL;

    sidp = ValNodeNew(NULL);
    sidp->choice = SEQID_GENERAL;
    sidp->data.ptrvalue = dp;

    sip = SeqIntNew();
    sip->from = 0;
    sip->to = range - 1;
    sip->id = sidp;

    slp->choice = SEQLOC_INT;
    slp->data.ptrvalue = sip;

    return(slp);
}

/**********************************************************
 *
 *   NLM_EXTERN SeqLocPtr LIBCALL GapToSeqLocEx(range, unknown):
 *
 *      Gets the size of gap and constructs SeqLoc block with
 *   $(seqlitdbtag) value as Dbtag.db and Dbtag.tag.id = 0.
 *
 **********************************************************/
NLM_EXTERN SeqLocPtr LIBCALL GapToSeqLocEx(Int4 range, Boolean unknown)
{
    SeqLocPtr slp;
    SeqIntPtr sip;
    SeqIdPtr  sidp;
    DbtagPtr  dp;

    if(range < 0)
        return(NULL);

    slp = ValNodeNew(NULL);
    if(range == 0)
    {
        slp->choice = SEQLOC_NULL;
        slp->data.ptrvalue = NULL;
        slp->next = NULL;
        return(slp);
    }

    dp = DbtagNew();
    if(unknown == FALSE)
        dp->db = StringSave(seqlitdbtag);
    else
        dp->db = StringSave(unkseqlitdbtag);
    dp->tag = ObjectIdNew();
    dp->tag->id = 0;
    dp->tag->str = NULL;

    sidp = ValNodeNew(NULL);
    sidp->choice = SEQID_GENERAL;
    sidp->data.ptrvalue = dp;

    sip = SeqIntNew();
    sip->from = 0;
    sip->to = range - 1;
    sip->id = sidp;

    slp->choice = SEQLOC_INT;
    slp->data.ptrvalue = sip;

    return(slp);
}

/**********************************************************
 *
 *   NLM_EXTERN Boolean LIBCALL ISAGappedSeqLoc(slp):
 *
 *      Looks at a single SeqLoc item. If it has the SeqId
 *   of type GENERAL with Dbtag.db == $(seqlitdbtag) and
 *   Dbtag.tag.id == 0, then returns TRUE, otherwise
 *   returns FALSE.
 *
 **********************************************************/
NLM_EXTERN Boolean LIBCALL ISAGappedSeqLoc(SeqLocPtr slp)
{
    SeqIdPtr sip;
    DbtagPtr dp;

    if(slp == NULL)
        return(FALSE);

    sip = SeqLocId(slp);
    if(sip == NULL || sip->choice != SEQID_GENERAL)
        return(FALSE);

    dp = (DbtagPtr) sip->data.ptrvalue;
    if(dp == NULL || dp->db == NULL || dp->tag == NULL)
        return(FALSE);

    if((StringCmp(seqlitdbtag, dp->db) == 0 ||
        StringCmp(unkseqlitdbtag, dp->db) == 0) && dp->tag->id == 0)
        return(TRUE);

    return(FALSE);
}

/**********************************************************
 *
 *   NLM_EXTERN DeltaSeqPtr LIBCALL GappedSeqLocsToDeltaSeqs(slp):
 *
 *      This functions is used only in the case, if ISAGappedSeqLoc()
 *   has returned TRUE.
 *      Converts SeqLoc set to the sequence of DeltaSeqs.
 *   Gbtag'ed SeqLocs it turns into SeqLits with the only "length"
 *   element. The regular SeqLocs saves as they are. Returns
 *   obtained DeltaSeq.
 *
 **********************************************************/
NLM_EXTERN DeltaSeqPtr LIBCALL GappedSeqLocsToDeltaSeqs(SeqLocPtr slp)
{
    DeltaSeqPtr res;
    DeltaSeqPtr dsp;
    SeqIntPtr   sip;
    SeqLitPtr   slip;
    SeqIdPtr    id;
    DbtagPtr    dp;

    dsp = ValNodeNew(NULL);
    dsp->next = NULL;
    dsp->choice = 0;
    res = dsp;
    for(; slp != NULL; slp = slp->next)
    {
        if(ISAGappedSeqLoc(slp) != FALSE)
        {
            dsp->next = ValNodeNew(NULL);
            dsp = dsp->next;
            sip = slp->data.ptrvalue;
            slip = SeqLitNew();
            slip->length = sip->to - sip->from + 1;
            dsp->choice = 2;
            dsp->data.ptrvalue = slip;
            id = SeqLocId(slp);
            if(id != NULL)
            {
                dp = (DbtagPtr) id->data.ptrvalue;
                if(dp != NULL && dp->db != NULL &&
                   StringCmp(unkseqlitdbtag, dp->db) == 0)
                {
                    slip->fuzz = IntFuzzNew();
                    slip->fuzz->choice = 4;
                }
            }
        }
        else
        {
            dsp->next = ValNodeNew(NULL);
            dsp = dsp->next;
            dsp->choice = 1;
            dsp->data.ptrvalue = AsnIoMemCopy((Pointer) slp,
                                              (AsnReadFunc) SeqLocAsnRead,
                                              (AsnWriteFunc) SeqLocAsnWrite);
        }
    }
    dsp = res->next;
    MemFree(res);
    return(dsp);
}

/* This structure and the functions following it are used to track the prior locations
 * of features that were affected by the removal of nucleotides, so that they may be
 * returned to their original status in an undo.
 */
typedef struct affectedfeat
{
  SeqFeatPtr feat_before;
  SeqFeatPtr feat_after;
} AffectedFeatData, PNTR AffectedFeatPtr;

static AffectedFeatPtr AffectedFeatNew (void)
{
  AffectedFeatPtr afp;
  
  afp = (AffectedFeatPtr) MemNew (sizeof (AffectedFeatData));
  if (afp != NULL)
  {
    afp->feat_before = NULL;
    afp->feat_after = NULL;
  }
  return afp;
}

static AffectedFeatPtr AffectedFeatFree (AffectedFeatPtr afp)
{
  if (afp == NULL) return NULL;
  afp->feat_before = SeqFeatFree (afp->feat_before);
  afp->feat_after = SeqFeatFree (afp->feat_after);
  afp = MemFree (afp);
  return NULL;
}

static ValNodePtr SeqEdJournalAffectedFeatsFree (ValNodePtr vnp)
{
  if (vnp == NULL) return NULL;
  vnp->next = SeqEdJournalAffectedFeatsFree (vnp->next);
  vnp->data.ptrvalue = AffectedFeatFree ((AffectedFeatPtr) (vnp->data.ptrvalue));
  ValNodeFree (vnp);
  return NULL;
}

static Boolean SeqEdRecreateDeletedFeats (SeqEdJournalPtr sejp)
{
  ValNodePtr      vnp;
  AffectedFeatPtr afp = NULL;
  Boolean         recreated_feats = FALSE;
  SeqEntryPtr     sep = NULL;
  SeqFeatPtr      sfp;

  for (vnp = sejp->affected_feats; vnp != NULL && afp == NULL; vnp = vnp->next)
  {
    if (vnp->choice == 1 || vnp->data.ptrvalue == NULL) continue;
    afp = (AffectedFeatPtr) vnp->data.ptrvalue;
    if (afp->feat_after == NULL && afp->feat_before != NULL)
    {
      vnp->choice = 1;
      if (sep == NULL)
      {
        sep = SeqMgrGetSeqEntryForData (sejp->bsp);        
        if (sep == NULL) return FALSE;
      }
      sfp = CreateNewFeature (sep, NULL, afp->feat_before->data.choice, afp->feat_before);
      afp->feat_before = NULL;
      recreated_feats = TRUE;
    }
  } 
  return recreated_feats;
}


/* This section of code deals with inserting new characters into a Bioseq and adjusting the
 * locations of the affected features.   It is adapted from code from SeqLocInsert.
 */
 
NLM_EXTERN void SeqEdInsertAdjustCdRgn 
(SeqFeatPtr sfp,
 BioseqPtr  bsp,
 Int4       insert_pos,
 Int4       len,
 Boolean    do_split)
{
  CdRegionPtr  crp;
  CodeBreakPtr prevcbp, cbp, nextcbp;
  
  if (sfp == NULL || bsp == NULL) return;
  crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
  if (crp == NULL) return;
  
  prevcbp = NULL;
  for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
  {
    nextcbp = cbp->next;
    cbp->loc = SeqEdSeqLocInsert (cbp->loc, bsp, insert_pos, len, do_split, NULL);
    if (cbp->loc == NULL)
    {
      if (prevcbp != NULL)
        prevcbp->next = nextcbp;
      else
        crp->code_break = nextcbp;
      cbp->next = NULL;
      CodeBreakFree (cbp);
    }
    else
    {
      prevcbp = cbp;
    }
  }
}

NLM_EXTERN void SeqEdInsertAdjustRNA 
(SeqFeatPtr sfp,
 BioseqPtr  bsp,
 Int4       insert_pos,
 Int4       len,
 Boolean    do_split)
{
  RnaRefPtr rrp;
  tRNAPtr   trp;
    
  if (sfp == NULL || bsp == NULL) return;
  rrp = (RnaRefPtr)(sfp->data.value.ptrvalue);
  if (rrp == NULL) return;
  if (rrp->ext.choice == 2)   /* tRNA */
  {
    trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
    if (trp->anticodon != NULL)
    {
      trp->anticodon = SeqEdSeqLocInsert (trp->anticodon, bsp, insert_pos, len, do_split, NULL);
    }
  }
}

/*****************************************************************************
*
*   SeqEdSeqLocInsert()
*       alters "head" by insert "len" residues before "pos" in any SeqLoc
*         on the Bioseq "target"
*       all SeqLocs not on "target" are unaltered
*       for SeqLocs on "target"
*          all SeqLocs before "pos" are unaltered
*          all SeqLocs >= "pos" are incremented by "len"
*          all SeqLocs spanning "pos"
*             if "split" == TRUE, are split into two SeqLocs, one to the
*               left of the insertion, the other to right
*             if "split" != TRUE, the SeqLoc is increased in length to cover
*               the insertion
*   	returns altered head or NULL if nothing left.
*       if ("newid" != NULL) replaces "target" with "newid" whether the
*          SeqLoc is altered on not.
*
*       Usage hints:
*          1) To update a feature location on "target" when 10 residues of
*               sequence have been inserted before position 5
*          SeqFeatPtr->location = SeqLocInsert ( SeqFeatPtr->location ,
*                "target", 5, 10, TRUE, NULL);  [for some feature types
*                      you may want "split" equal FALSE]
*          2) To insert the complete feature table from "source" into a
*                different Bioseq "dest" before position 20 in "dest"
*          SFP->location = SeqLocInsert(SFP->location, "source", 0, 20,
*                FALSE, "dest");
*   
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr LIBCALL SeqEdSeqLocInsert (SeqLocPtr head, BioseqPtr target, Int4 pos, Int4 len,
                                               Boolean split, SeqIdPtr newid)
{
	SeqIntPtr sip, sip2;
	SeqPntPtr spp;
	PackSeqPntPtr pspp, pspp2;
	SeqBondPtr sbp;
	SeqLocPtr slp, tmp, prev, next, thead, tmp2;
	Int4 diff, numpnt, i, tpos;
	Uint1 oldchoice;
	ValNode vn;
	SeqIdPtr sidp;

	if ((head == NULL) || (target == NULL))
		return head;

	head->next = NULL;   /* caller maintains chains */

	diff = len;

    switch (head->choice)
    {
        case SEQLOC_BOND:   /* bond -- 2 seqs */
			vn.next = NULL;
			vn.choice = SEQLOC_PNT;

			sbp = (SeqBondPtr)(head->data.ptrvalue);
			vn.data.ptrvalue = (Pointer)(sbp->a);
			SeqEdSeqLocInsert(&vn, target, pos, len, split, newid);
			sbp->a = (SeqPntPtr)(vn.data.ptrvalue);
			if (sbp->b != NULL)
			{
				vn.data.ptrvalue = (Pointer)(sbp->b);
				SeqEdSeqLocInsert(&vn, target, pos, len, split, newid);
				sbp->b = (SeqPntPtr)(vn.data.ptrvalue);
			}
			break;
        case SEQLOC_FEAT:   /* feat -- can't track yet */
        case SEQLOC_NULL:    /* NULL */
			break;
        case SEQLOC_EMPTY:    /* empty */
        case SEQLOC_WHOLE:    /* whole */
			if (newid != NULL)
			{
				sidp = (SeqIdPtr)(head->data.ptrvalue);
				if ( SeqIdIn(sidp, target->id))
				{
					SeqIdFree(sidp);
					sidp = SeqIdDup(newid);
					head->data.ptrvalue = (Pointer)sidp;
				}
			}
			break;
        case SEQLOC_MIX:    /* mix -- more than one seq */
        case SEQLOC_EQUIV:    /* equiv -- ditto */
        case SEQLOC_PACKED_INT:    /* packed int */
			prev = NULL;
			thead = NULL;
			for (slp = (SeqLocPtr)(head->data.ptrvalue); slp != NULL; slp = next)
			{
				next = slp->next;
				oldchoice = slp->choice;
				tmp = SeqEdSeqLocInsert(slp, target, pos, len, split, newid);
				if (tmp != NULL)
				{
					if ((head->choice != SEQLOC_EQUIV) &&
						(oldchoice != tmp->choice))  /* split interval? */
					{
						if ((oldchoice == SEQLOC_INT) &&
							(tmp->choice == SEQLOC_PACKED_INT))
						{
							tmp2 = tmp;
							tmp = (SeqLocPtr)(tmp2->data.ptrvalue);
							MemFree(tmp2);
							while (tmp->next != NULL)
							{
								if (prev != NULL)
									prev->next = tmp;
								else
									thead = tmp;
								prev = tmp;
								tmp = tmp->next;
							}
						}
					}
					if (prev != NULL)
						prev->next = tmp;
					else
						thead = tmp;
					prev = tmp;
				}
			}
			head->data.ptrvalue = thead;
			if (thead == NULL)
				head = SeqLocFree(head);
            break;
        case SEQLOC_INT:    /* int */
			sip = (SeqIntPtr)(head->data.ptrvalue);
			if (SeqIdIn(sip->id, target->id))
			{
				if (newid != NULL)   /* change id? */
				{
					SeqIdFree(sip->id);
					sip->id = SeqIdDup(newid);
				}

				if (sip->to < pos)  /* completely before insertion */
				{
					break;
				}

				if ((! split) || (sip->from >= pos))  /* interval unbroken */
				{
					if (sip->from >= pos)
						sip->from += len;
					sip->to += len;
					break;
				}

						                          /* split interval */
				sip2 = SeqIntNew();
				slp = ValNodeNew(NULL);
				slp->choice = SEQLOC_INT;
				slp->data.ptrvalue = (Pointer)sip2;
				sip2->strand = sip->strand;
				sip2->id = SeqIdDup(sip->id);

				sip2->to = sip->to + len;
				sip2->from = pos + len;
				sip2->if_to = sip->if_to;
				sip->if_to = NULL;
				sip->to = pos - 1;
				head->next = slp;

				if (sip->strand == Seq_strand_minus)  /* reverse order */
				{
					head->data.ptrvalue = (Pointer)sip2;
					slp->data.ptrvalue = (Pointer)sip;
				}

				thead = head;   /* make split interval into PACKED_INT */
				head = ValNodeNew(NULL);
				head->choice = SEQLOC_PACKED_INT;
				head->data.ptrvalue = thead;

			}
            break;
        case SEQLOC_PNT:    /* pnt */
			spp = (SeqPntPtr)(head->data.ptrvalue);
			if (SeqIdIn(spp->id, target->id))
			{
				if (newid != NULL)   /* change id? */
				{
					SeqIdFree(spp->id);
					spp->id = SeqIdDup(newid);
				}

				if (spp->point >= pos)
					spp->point += len;
			}
            break;
        case SEQLOC_PACKED_PNT:    /* packed pnt */
			pspp = (PackSeqPntPtr)(head->data.ptrvalue);
			if (SeqIdIn(pspp->id, target->id))
			{
				if (newid != NULL)   /* change id? */
				{
					SeqIdFree(pspp->id);
					pspp->id = SeqIdDup(newid);
				}

				numpnt = PackSeqPntNum(pspp);
				pspp2 = PackSeqPntNew();
				head->data.ptrvalue = pspp2;
				for (i = 0; i < numpnt; i++)
				{
					tpos = PackSeqPntGet(pspp, i);
					if (tpos >= pos)
						tpos += len;
					PackSeqPntPut(pspp2, tpos);
				}
				pspp2->id = pspp->id;
				pspp->id = NULL;
				pspp2->fuzz = pspp->fuzz;
				pspp->fuzz = NULL;
				pspp2->strand = pspp->strand;
				PackSeqPntFree(pspp);
			}
            break;
        default:
            break;
    }

	if (head == NULL)
		ErrPost(CTX_NCBIOBJ, 1, "SeqEdSeqLocInsert: lost a SeqLoc");

	return head;
}

static SeqLocPtr 
SeqEdDeleteFromSeqLocBond (SeqLocPtr head, SeqIdPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed)
{
  SeqBondPtr sbp;
  SeqPntPtr  spp;
  Int4       diff;

  if (head == NULL || target == NULL || head->choice != SEQLOC_BOND) return NULL;
  sbp = (SeqBondPtr)(head->data.ptrvalue);
  spp = sbp->a;
  diff = to - from + 1;

  if (SeqIdForSameBioseq(spp->id, target))
  {
	if (spp->point >= from)
	{
	  if (spp->point <= to)   /* delete it */
	  {
	    *changed = TRUE;
		sbp->a = SeqPntFree(spp);
	  }
      else if (merge)
	  {
		spp->point -= diff;
	  }
	}
  }
  spp = sbp->b;
  if (spp != NULL)
  {
	if (SeqIdForSameBioseq(spp->id, target))
	{
	  if (spp->point >= from)
	  {
		if (spp->point <= to)   /* delete it */
		{
		  *changed = TRUE;
		  sbp->b = SeqPntFree(spp);
		}
		else if (merge)
		{
		  spp->point -= diff;
		}
      }
	}
  }
  if (sbp->a == NULL)
  {
	if (sbp->b != NULL)   /* only a required */
	{
      sbp->a = sbp->b;
      sbp->b = NULL;
	}
	else
	{
	  head = SeqLocFree(head);
	}
  }
  return head;
}

static SeqLocPtr DeleteFromSeqLocWhole 
(SeqLocPtr head, BioseqPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed)
{
  SeqIdPtr  sidp;
  SeqIntPtr sip;
  SeqLocPtr slp, tmp;

  if (head == NULL || target == NULL || head->choice != SEQLOC_WHOLE) return NULL;

  sidp = (SeqIdPtr)(head->data.ptrvalue);

  if ( SeqIdIn(sidp, target->id))
  {
	if ((from == 0) && (to >= (target->length - 1)))
	{					   /* complete delete */
	  head = SeqLocFree(head);
      *changed = TRUE;
	  return head;
	}

	if (! merge)   /* split it up */
	{
      SeqIdFree(sidp);
      head->choice = SEQLOC_PACKED_INT;
      head->data.ptrvalue = NULL;
      slp = NULL;
      if (from != 0)
      {
		sip = SeqIntNew();
		sip->from = 0;
		sip->to = from - 1;
		sip->id = SeqIdDup(target->id);
		slp = ValNodeNew(NULL);
		slp->choice = SEQLOC_INT;
		slp->data.ptrvalue = sip;
		head->data.ptrvalue = slp;
		*changed = TRUE;
	  }
      if (to < (target->length - 1))
	  {
		sip = SeqIntNew();
		sip->from = to + 1;
		sip->to = target->length - 1;
		sip->id = SeqIdDup(target->id);
		tmp = ValNodeNew(NULL);
		tmp->choice = SEQLOC_INT;
		tmp->data.ptrvalue = sip;
		if (slp != NULL)
		  slp->next = tmp;
		else
		  head->data.ptrvalue = tmp;
		*changed = TRUE;
	  }
	}
  }
  return head;
}

static SeqLocPtr SeqEdDeleteFromSeqLocPackedInt 
(SeqLocPtr head, BioseqPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed, BoolPtr partial5, BoolPtr partial3)
{
  Boolean   part5, part3, first;
  SeqLocPtr slp, tmp, prev, next, thead;
  SeqIntPtr sip, sip2;

  if (head == NULL || target == NULL) return NULL;
  if (head->choice != SEQLOC_MIX && head->choice != SEQLOC_EQUIV && head->choice != SEQLOC_PACKED_INT) 
	return NULL;
  prev = NULL;
  thead = NULL;
  part5 = FALSE;
  part3 = FALSE;
  first = TRUE;
  for (slp = (SeqLocPtr)(head->data.ptrvalue); slp != NULL; slp = next)
  {
	next = slp->next;
	tmp = SeqEdSeqLocDelete (slp, target, from, to, merge, changed, &part5, &part3);
	if (first) 
	{
	  if (partial5 != NULL) 
	  {
		*partial5 = part5;
      }
    }
    first = FALSE;
    if (tmp != NULL)
    {
      if (prev != NULL)
        {
          if ((merge) && (prev->choice == SEQLOC_INT) && (tmp->choice == SEQLOC_INT))
          {
            sip = (SeqIntPtr)(prev->data.ptrvalue);
            sip2 = (SeqIntPtr)(tmp->data.ptrvalue);

            if (SeqIdForSameBioseq(sip->id, sip2->id))
            {
              /* merge intervals? */
              if ((sip->strand == Seq_strand_minus) &&
                  (sip2->strand == Seq_strand_minus))
              {
                if (sip->from == (sip2->to + 1))
                {
                  sip->from = sip2->from;
                  sip->if_from = sip2->if_from;
                  sip2->if_from = NULL;
                  tmp = SeqLocFree(tmp);
                }
              }
              else if((sip->strand != Seq_strand_minus) &&
                      (sip2->strand != Seq_strand_minus))
              {
                if (sip->to == (sip2->from - 1))
                {
                  sip->to = sip2->to;
                  sip->if_to = sip2->if_to;
                  sip2->if_to = NULL;
                  tmp = SeqLocFree(tmp);
                }
              }
            }
          }
          else if ((prev->choice == SEQLOC_NULL) && (tmp->choice == SEQLOC_NULL))
          {
            tmp = SeqLocFree(tmp);
            *changed = TRUE;
          }
        }
        else if (tmp->choice == SEQLOC_NULL)
        {
          tmp = SeqLocFree(tmp);
          *changed = TRUE;
        }

        if (tmp != NULL)   /* still have one? */
        {
          if (prev != NULL)
            prev->next = tmp;
          else
            thead = tmp;
          prev = tmp;
        }
      }
      else
	  {
        *changed = TRUE;
	  }
    }
  if (partial3 != NULL) 
  {
    *partial3 = part3;
  }
  if (prev != NULL)
  {
    if (prev->choice == SEQLOC_NULL)  /* ends with NULL */
    {
      prev = NULL;
      for (slp = thead; slp->next != NULL; slp = slp->next)
      {
        prev = slp;
	  }
      if (prev != NULL)
      {
        prev->next = NULL;
        SeqLocFree(slp);
      }
      else
      {
        thead = SeqLocFree(thead);
      }
      *changed = TRUE;
    }
  }
  head->data.ptrvalue = thead;
  if (thead == NULL)
    head = SeqLocFree(head);
  return head;
}

static SeqLocPtr SeqEdDeleteFromSeqLocInt (SeqLocPtr head, BioseqPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed, BoolPtr partial5, BoolPtr partial3)
{
  Int4      diff;
  SeqIntPtr sip, sip2;
  SeqLocPtr slp, tmp;

  if (head == NULL || target == NULL || head->choice != SEQLOC_INT) return NULL;

  sip = (SeqIntPtr)(head->data.ptrvalue);
  if ( !SeqIdIn(sip->id, target->id)) return head;

  diff = to - from + 1;

  if (sip->to < from)  /* completely before cut */
	return head;

  /* completely contained in cut */
  if ((sip->from >= from) && (sip->to <= to))
  {
    head = SeqLocFree(head);
  	*changed = TRUE;
    return head;
  }

  if (sip->from > to)  /* completely past cut */
  {
    sip->from -= diff;
    sip->to -= diff;
    return head;
  }
  /* overlap here */
  if (sip->to > to)
  {
    sip->to -= diff;
  }
  else                /* to inside cut, so partial delete */
  {
	sip->to = from - 1;
	*changed = TRUE;
	if (partial3 != NULL) 
	{
	  *partial3 = TRUE;
	}
  }
		
  if (sip->from >= from)   /* from inside cut, partial del */
  {
    *changed = TRUE;
    sip->from = to + 1;
    sip->from -= diff;
    if (partial5 != NULL)
    {
      *partial5 = TRUE;
    }

    if (merge)
      return head;
	
	/* interval spans cut.. only in non-merge */
    /* have to split */

	if ((sip->from < from) && (sip->to > to))
	{
      *changed = TRUE;
      head->choice = SEQLOC_PACKED_INT;
      head->data.ptrvalue = NULL;
      tmp = ValNodeNew(NULL);
      tmp->choice = SEQLOC_INT;
      tmp->data.ptrvalue = sip;

      sip2 = SeqIntNew();
      sip2->from = to + 1;
      sip2->to = sip->to;
      sip2->strand = sip->strand;
      sip2->if_to = sip->if_to;
      sip2->id = SeqIdDup(target->id);
      slp = ValNodeNew(NULL);
      slp->choice = SEQLOC_INT;
      slp->data.ptrvalue = sip2;

      sip->if_to = NULL;
      sip->to = from - 1;

      if (sip->strand == Seq_strand_minus)
      {
        head->data.ptrvalue = slp;
        slp->next = tmp;
      }
      else
      {
        head->data.ptrvalue = tmp;
        tmp->next = slp;
      }

    }
  }
  return head;
}

static SeqLocPtr SeqEdDeleteFromSeqLocPackedPnt 
(SeqLocPtr head, BioseqPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed)
{
  PackSeqPntPtr pspp, pspp2;
  Int4          i, diff, numpnt, tpos;

  if (head == NULL || target == NULL || head->choice != SEQLOC_PACKED_PNT) return NULL;

  pspp = (PackSeqPntPtr)(head->data.ptrvalue);
  if (!SeqIdIn (pspp->id, target->id)) return head;

  diff = to - from + 1;

  numpnt = PackSeqPntNum(pspp);
  pspp2 = PackSeqPntNew();
  head->data.ptrvalue = pspp2;
  for (i = 0; i < numpnt; i++)
  {
    tpos = PackSeqPntGet(pspp, i);
    if (tpos < from)
	{
      PackSeqPntPut(pspp2, tpos);
	}
    else
    {
      if (tpos > to)
      {
        if (merge)
		{
          tpos -= diff;
		}
        PackSeqPntPut(pspp2, tpos);
      }
      else
	  {
        *changed = TRUE;
	  }
    }
  }
  pspp2->id = pspp->id;
  pspp->id = NULL;
  pspp2->fuzz = pspp->fuzz;
  pspp->fuzz = NULL;
  pspp2->strand = pspp->strand;
  PackSeqPntFree(pspp);
  numpnt = PackSeqPntNum(pspp2);
  if (! numpnt)
  {
    head = SeqLocFree(head);
  }
  return head;
}


/*****************************************************************************
*
*   SeqEdSeqLocDelete()
*   	returns altered head or NULL if nothing left.
*   sets changed=TRUE if all or part of loc is deleted
*   does NOT set changed if location coordinates are only moved
*   if (merge) then corrects coordinates upstream of to
*   else
*     splits intervals covering from-to, does not correct upstream of to
*
*****************************************************************************/
NLM_EXTERN SeqLocPtr SeqEdSeqLocDelete (SeqLocPtr head, BioseqPtr target, Int4 from, Int4 to, Boolean merge, BoolPtr changed, BoolPtr partial5, BoolPtr partial3)
{
	SeqPntPtr spp;
	Int4 diff;

	if ((head == NULL) || (target == NULL))
		return head;

	head->next = NULL;   /* caller maintains chains */
	diff = to - from + 1;
	
  switch (head->choice)
  {
    case SEQLOC_BOND:   /* bond -- 2 seqs */
      head = SeqEdDeleteFromSeqLocBond (head, target->id, from, to, merge, changed);
			break;
    case SEQLOC_FEAT:   /* feat -- can't track yet */
    case SEQLOC_NULL:    /* NULL */
    case SEQLOC_EMPTY:    /* empty */
			break;
    case SEQLOC_WHOLE:    /* whole */
      head = DeleteFromSeqLocWhole (head, target, from, to, merge, changed);
			break;
    case SEQLOC_MIX:    /* mix -- more than one seq */
    case SEQLOC_EQUIV:    /* equiv -- ditto */
    case SEQLOC_PACKED_INT:    /* packed int */
      head = SeqEdDeleteFromSeqLocPackedInt (head, target, from, to, merge, changed, partial5, partial3);
      break;
    case SEQLOC_INT:    /* int */
      head = SeqEdDeleteFromSeqLocInt (head, target, from, to, merge, changed, partial5, partial3);
      break;
    case SEQLOC_PNT:    /* pnt */
			spp = (SeqPntPtr)(head->data.ptrvalue);
			if (SeqIdIn (spp->id, target->id)) 
			{
				if ((spp->point >= from) && (spp->point <= to))
				{
					head = SeqLocFree(head);
					*changed = TRUE;
				}
				else if (spp->point > to)
				{
						spp->point -= diff;
				}
			}
      break;
    case SEQLOC_PACKED_PNT:    /* packed pnt */
      head = SeqEdDeleteFromSeqLocPackedPnt (head, target, from, to, merge, changed);
      break;
    default:
      break;
  }

	return head;
}

NLM_EXTERN SeqFeatPtr 
SeqEdGetNextFeature 
(BioseqPtr bsp, 
 SeqFeatPtr curr,
 Uint1 seqFeatChoice,
 Uint1 featDefChoice,
 SeqMgrFeatContext PNTR context,
 Boolean byLabel,
 Boolean byLocusTag,
 Uint2         entityID)

{
  SMFeatItemPtr PNTR  array = NULL;
  BioseqExtraPtr      bspextra;
  Uint4               i;
  SMFeatItemPtr       item;
  Int4                num = 0;
  ObjMgrDataPtr       omdp;
  ObjMgrPtr           omp;
  Uint1               seqfeattype;

  if (context == NULL) return NULL;

  /* if curr is NULL, initialize context fields (in user's stack) */


  if (curr == NULL) {
    if (bsp == NULL) return NULL;
    omdp = (ObjMgrDataPtr) bsp->omdp;
    if (omdp == NULL) 
    {
      omp = ObjMgrWriteLock ();
      omdp = ObjMgrFindByData (omp, bsp);
      ObjMgrUnlock ();
      bsp->omdp = (Pointer) omdp;	
    }
    if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return NULL;

    context->omdp = (Pointer) omdp;
    context->index = 0;
  }

  omdp = (ObjMgrDataPtr) context->omdp;
  if (omdp == NULL) return NULL;
  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return NULL;
  if (byLocusTag) {
    array = bspextra->genesByLocusTag;
    num = bspextra->numgenes;
  } else if (byLabel) {
    array = bspextra->featsByLabel;
    num = bspextra->numfeats;
  } else {
    array = bspextra->featsByPos;
    num = bspextra->numfeats;
  }
  if (array == NULL || num < 1) return NULL;

  i = context->index;

  /* now look for next appropriate feature */

  while (i < num) {
    item = array [i];
    if (item != NULL) {
      curr = item->sfp;
      i++;
      if (curr != NULL) {
        seqfeattype = curr->data.choice;
        if ((seqFeatChoice == 0 || seqfeattype == seqFeatChoice) &&
            (featDefChoice == 0 || item->subtype == featDefChoice) &&
            (! item->ignore)) {
          context->entityID = entityID;
          context->itemID = item->itemID;
          context->sfp = curr;
          context->sap = item->sap;
          context->bsp = item->bsp;
          context->label = item->label;
          context->left = item->left;
          context->right = item->right;
          context->dnaStop = item->dnaStop;
          context->partialL = item->partialL;
          context->partialR = item->partialR;
          context->farloc = item->farloc;
          context->strand = item->strand;
          context->seqfeattype = seqfeattype;
          context->featdeftype = item->subtype;
          context->numivals = item->numivals;
          context->ivals = item->ivals;
          context->userdata = NULL;
          context->omdp = (Pointer) omdp;
          if (byLocusTag) {
            context->index = i;
          } else if (byLabel) {
            context->index = i;
          } else {
            context->index = item->index + 1;
          }
          return curr;
        }
      }
    }
  }

  return NULL;
}

static void ReindexExtendedFeatures (SeqEdJournalPtr sejp)
{
  ValNodePtr        vnp;
  AffectedFeatPtr   afp;
  SeqFeatPtr        affected_sfp, real_sfp;
  SeqMgrFeatContext fcontext;
  
  for (vnp = sejp->affected_feats; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == 1 && vnp->data.ptrvalue != NULL)
    {
      afp = (AffectedFeatPtr) vnp->data.ptrvalue;
      affected_sfp = afp->feat_after;
      if (affected_sfp != NULL)
      {
        real_sfp = SeqMgrGetDesiredFeature (sejp->entityID, sejp->bsp, affected_sfp->idx.itemID, 0, NULL, &fcontext);
        SeqEdReindexFeature (real_sfp, sejp->bsp);
      }
    }
  }
}


static Boolean DoesSeqFeatMatch (SeqFeatPtr a, SeqFeatPtr b)
{
  if (a == b) return TRUE;
  if (a == NULL || b == NULL) return FALSE;
  
  if (a->data.choice != b->data.choice) return FALSE;
  if (SeqLocCompare (a->location, b->location) != SLC_A_EQ_B)
  {
    return FALSE;
  }
  return TRUE;
}


static void SeqEdInsertAdjustFeat (SeqFeatPtr sfp, SeqEdJournalPtr sejp, Int4 insert_point)
{
  ValNodePtr      vnp;
  AffectedFeatPtr afp = NULL;
  SeqLocPtr       tmp_loc;
  Boolean         split_mode;
  
  if (sfp == NULL || sejp == NULL)
  {
    return;
  }
  
  for (vnp = sejp->affected_feats; vnp != NULL && afp == NULL; vnp = vnp->next)
  {
    afp = (AffectedFeatPtr) vnp->data.ptrvalue;
    if (afp != NULL && DoesSeqFeatMatch (afp->feat_after, sfp))
    {
      vnp->choice = 1;
    }
    else
    {
      afp = NULL;
    }
  }

  /* if we're inserting a gap and the feature is a coding region, need to split location
   * regardless of mode */
  split_mode = sejp->spliteditmode;
  if (sejp->action == eSeqEdInsertGap || sejp->action == eSeqEdDeleteGap
      && sfp->data.choice == SEQFEAT_CDREGION)
  {
    split_mode = TRUE;
  }
  
  if (afp != NULL)
  {
    tmp_loc = sfp->location;
    sfp->location = afp->feat_before->location;
    afp->feat_before->location = tmp_loc;
  }
  else
  {
    sfp->location = SeqEdSeqLocInsert (sfp->location, sejp->bsp, insert_point,
                                       sejp->num_chars, split_mode, NULL);
  }
	switch (sfp->data.choice)
	{
    case SEQFEAT_CDREGION:   /* cdregion */
      SeqEdInsertAdjustCdRgn (sfp, sejp->bsp, insert_point, sejp->num_chars,
                              split_mode);
      break;
    case SEQFEAT_RNA:
      SeqEdInsertAdjustRNA (sfp, sejp->bsp, insert_point, sejp->num_chars,
                            split_mode);
      break;
    default:
      break;
  }
}

static Boolean IsDeltaSeqGap (DeltaSeqPtr dsp)
{
  SeqLitPtr   slip;
  if (dsp == NULL || dsp->choice != 2 || dsp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  slip = (SeqLitPtr) (dsp->data.ptrvalue);
  if (slip->seq_data == NULL)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static Boolean IsDeltaSeqUnknownGap (DeltaSeqPtr dsp)
{
  SeqLitPtr   slip;
  if (dsp == NULL || dsp->choice != 2 || dsp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  slip = (SeqLitPtr) (dsp->data.ptrvalue);
  if (slip->seq_data == NULL && slip->fuzz != NULL && slip->fuzz->choice == 4)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static DeltaSeqPtr GetDeltaSeqForOffset (BioseqPtr bsp, Int4 offset, Int4Ptr seqstart)
{
  Int4        curr_pos = 0;
  Boolean     found = FALSE;
  SeqLocPtr   slp;
  SeqLitPtr   slip = NULL;
  DeltaSeqPtr dsp;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL
      || offset < 0)
  {
    return NULL;
  }
  
  if (seqstart != NULL)
  {
    *seqstart = 0;
  }
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL && !found)
  {
    if (dsp->data.ptrvalue == NULL) continue;
		if (dsp->choice == 1) 
		{  /* SeqLoc */
		  slp = (SeqLocPtr)(dsp->data.ptrvalue);
      curr_pos += SeqLocLen (slp);
		}
		else if (dsp->choice == 2)
		{
		  slip = (SeqLitPtr) (dsp->data.ptrvalue);
		  curr_pos += slip->length;
		}
		if (curr_pos > offset
		    || (curr_pos == offset 
		      && (dsp->next == NULL || ! IsDeltaSeqGap (dsp))))
		{
		  found = TRUE;
		}
		else
		{
		  if (seqstart != NULL)
		  {
		    *seqstart = curr_pos;
		  }
		  dsp=dsp->next;
		}
  }

  return dsp;
}

static Boolean 
SeqEdInsertByteStore 
(ByteStorePtr seq_data, 
 Int4         insert_point, 
 CharPtr      char_data,
 Int4         num_chars,
 Uint1        moltype)
{
  Char         ch;
  Int4         i;
  
  if (seq_data == NULL || insert_point < 0 || char_data == NULL || num_chars < 1)
  {
    return FALSE;
  }
  BSSeek(seq_data, insert_point, SEEK_SET);
  Nlm_BSAdd(seq_data, num_chars, FALSE);
  BSSeek(seq_data, insert_point, SEEK_SET);
  for (i = 0; i < num_chars; i++)
  {
    ch = TO_UPPER (char_data [i]); 
    if ( ISA_na (moltype) ) {
      if (ch == 'U') ch = 'T';
      if (ch == 'X') ch = 'N';
      if ( StringChr ("EFIJLOPQXZ-.*", ch) == NULL ) { 
        BSPutByte ( seq_data, (Int2) ch );
      }
    } 
    else 
    {
      if ( StringChr("JO-.", ch) == NULL ) { 
        BSPutByte ( seq_data, (Int2) ch );
      }
    } 
  } 
  return TRUE;
}

static Boolean SeqEdInsertRaw (SeqEdJournalPtr sejp, Int4 insert_point)
{
  Boolean      rval;

  if (sejp == NULL || sejp->bsp == NULL || sejp->bsp->repr != Seq_repr_raw 
      || sejp->char_data == NULL || sejp->num_chars == 0 || insert_point < 0) 
  {
    return FALSE;    
  }

  rval = SeqEdInsertByteStore (sejp->bsp->seq_data, insert_point, 
                               sejp->char_data, sejp->num_chars, sejp->moltype);

  if (rval)
  {
    sejp->bsp->length += sejp->num_chars;
  }
  return rval;
}

static Boolean
SeqEdInsertIntoDeltaGap 
(DeltaSeqPtr dsp,
 SeqEdJournalPtr sejp, 
 Int4 insert_point)
{
  SeqLitPtr    slip, slip_data, slip_second_gap;
  Boolean      rval = FALSE;
  DeltaSeqPtr  dsp_data, dsp_second_gap;
  IntFuzzPtr   ifp = NULL;
  
  if (dsp == NULL || dsp->choice != 2 || dsp->data.ptrvalue == NULL)
  {
    return rval;
  }
  slip = (SeqLitPtr) dsp->data.ptrvalue;
  if (slip->seq_data != NULL)
  {
    return rval;
  }
  
  if (slip->fuzz != NULL && slip->fuzz->choice == 4)
  {
    ifp = IntFuzzNew ();
    ifp->choice = 4;
  }
  
  /* split the gap in two and create a new DeltaSeqPtr in the middle */
  slip_data = SeqLitNew ();
  slip_data->seq_data_type = Seq_code_iupacna;
  slip_data->seq_data = BSNew (sejp->num_chars);
  rval = SeqEdInsertByteStore (slip_data->seq_data, 0, 
                               sejp->char_data, sejp->num_chars, sejp->moltype);
  if (rval)
  {
    slip_data->length = sejp->num_chars;
    /* create second gap */
    slip_second_gap = SeqLitNew ();
    slip_second_gap->length = slip->length - insert_point;
    slip_second_gap->fuzz = ifp;
    /* truncate first gap */
    slip->length = insert_point;
    dsp_data = ValNodeNew (NULL);
    dsp_data->choice = 2;
    dsp_data->data.ptrvalue = slip_data;
    dsp_second_gap = ValNodeNew (NULL);
    dsp_second_gap->choice = 2;
    dsp_second_gap->data.ptrvalue = slip_second_gap;
    dsp_second_gap->next = dsp->next;
    dsp_data->next = dsp_second_gap;
    dsp->next = dsp_data;
  }
  return rval;
}

static Boolean IsInsertAllNs (SeqEdJournalPtr sejp)
{
  Int4 k;
  
  if (sejp == NULL || sejp->char_data == NULL || sejp->num_chars < 1)
  {
    return FALSE;
  }
  
  for (k = 0; k < sejp->num_chars; k++)
  {
    if (TO_LOWER (sejp->char_data [k]) != 'n')
    {
      return FALSE;
    } 
  }
  return TRUE;
}

static Boolean SeqEdInsertDelta (SeqEdJournalPtr sejp, Int4 insert_point)
{
  DeltaSeqPtr  dsp;
  SeqLitPtr    slip;
  Int4         seqstart = 0;
  ByteStorePtr bs_new;
  Boolean      rval;
  
  if (sejp == NULL || sejp->bsp == NULL || sejp->bsp->repr != Seq_repr_delta 
      || sejp->bsp->seq_ext_type != 4
      || sejp->char_data == NULL || sejp->num_chars == 0
      || insert_point < 0) 
  {
    return FALSE;    
  }

  dsp =  GetDeltaSeqForOffset (sejp->bsp, insert_point, &seqstart);
  
  if (dsp == NULL || dsp->choice != 2 || dsp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  
  slip = (SeqLitPtr) dsp->data.ptrvalue;
  insert_point -= seqstart;

  if (IsDeltaSeqGap (dsp))
  {
    if (slip->fuzz == NULL && IsInsertAllNs (sejp))
    {
      slip->length += sejp->num_chars;
      rval = TRUE;
    }
    else
    {
      rval = SeqEdInsertIntoDeltaGap (dsp, sejp, insert_point);    
    }
  }
  else
  {
    if (slip->seq_data_type != Seq_code_iupacna)
    {
      bs_new = BSConvertSeq(slip->seq_data, Seq_code_iupacna, 
                            slip->seq_data_type, 
                            slip->length);
      slip->seq_data_type = Seq_code_iupacna;
      slip->seq_data = bs_new;
    }

    rval = SeqEdInsertByteStore (slip->seq_data, insert_point, 
                                 sejp->char_data, sejp->num_chars,
                                 sejp->moltype);
    if (rval)
    {
      slip->length += sejp->num_chars;
    }
  }

  if (rval)
  {
    sejp->bsp->length += sejp->num_chars;
  }
  return rval;
}

static Boolean
SeqEdInsertGap (SeqEdJournalPtr sejp, Int4 insert_point)
{
  DeltaSeqPtr  dsp, dsp_gap, dsp_after;
  Int4         seqstart = 0;
  SeqLitPtr    slip, slip_before, slip_gap, slip_after;
  ByteStorePtr bs_new;
  
  if (sejp == NULL || sejp->bsp == NULL || sejp->bsp->repr != Seq_repr_delta 
      || sejp->bsp->seq_ext_type != 4
      || sejp->char_data == NULL || sejp->num_chars == 0
      || insert_point < 0) 
  {
    return FALSE;    
  }

  dsp =  GetDeltaSeqForOffset (sejp->bsp, insert_point, &seqstart);  
  
  if (dsp == NULL || dsp->choice != 2 || dsp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  
  slip_gap = SeqLitNew ();
  slip_gap->seq_data_type = 0;
  slip_gap->seq_data = NULL;
  slip_gap->length = sejp->num_chars;
  if (sejp->unknown_gap)
  {
    slip_gap->fuzz = IntFuzzNew ();
    slip_gap->fuzz->choice = 4;
  }
  
  slip = (SeqLitPtr) (dsp->data.ptrvalue);
  
  /* make insert_point relative to start of this SeqLit */
  insert_point -= seqstart;
  
  if (insert_point == 0)
  {
    /* insert gap before */
    dsp_after = ValNodeNew (NULL);
    dsp_after->choice = 2;
    dsp_after->data.ptrvalue = slip;
    dsp_after->next = dsp->next;
    dsp->next = dsp_after;
    dsp->data.ptrvalue = slip_gap;
  }
  else if (insert_point == slip->length)
  {
    /* insert gap after */
    dsp_after = ValNodeNew (NULL);
    dsp_after->choice = 2;
    dsp_after->data.ptrvalue = slip_gap;
    dsp_after->next = dsp->next;
    dsp->next = dsp_after;
  }
  else if (IsDeltaSeqUnknownGap (dsp))
  {
    /* can't insert gap inside gap of unknown length */
    slip_gap = SeqLitFree (slip_gap);
    return FALSE;
  }
  else if (IsDeltaSeqGap (dsp) && !sejp->unknown_gap)
  {
    slip_gap = SeqLitFree (slip_gap);
    slip->length += sejp->num_chars;
  }
  else
  {
    slip_before = SeqLitNew ();
    slip_before->seq_data_type = Seq_code_iupacna;
    slip_before->length = insert_point;
    
    slip_after = SeqLitNew ();
    slip_after->seq_data_type = Seq_code_iupacna;
    slip_after->length = slip->length - insert_point;
    
    if (slip->seq_data != NULL)
    {
      if (slip->seq_data_type != Seq_code_iupacna)
      {
        bs_new = BSConvertSeq(slip->seq_data, Seq_code_iupacna, 
                              slip->seq_data_type, 
                              slip->length);
        slip->seq_data_type = Seq_code_iupacna;
        slip->seq_data = bs_new;
      }
      slip_before->seq_data = BSNew (slip_before->length);
      slip_after->seq_data = BSNew (slip_after->length);
      
      BSSeek(slip->seq_data, 0, SEEK_SET);
      BSInsertFromBS (slip_before->seq_data, slip->seq_data, slip_before->length);
      BSInsertFromBS (slip_after->seq_data, slip->seq_data, slip_after->length);
    }
   
    dsp_after = ValNodeNew (NULL);
    dsp_after->choice = 2;
    dsp_after->data.ptrvalue = slip_after;
    dsp_after->next = dsp->next;
    
    dsp_gap = ValNodeNew (NULL);
    dsp_gap->choice = 2;
    dsp_gap->data.ptrvalue = slip_gap;
    dsp_gap->next = dsp_after;
    
    dsp->data.ptrvalue = slip_before;
    dsp->next = dsp_gap;
    slip = SeqLitFree (slip);
  }
  
  sejp->bsp->length += sejp->num_chars;

  return TRUE;
}

NLM_EXTERN Boolean 
SeqEdInsert (SeqEdJournalPtr sejp)
{
  Int4              len;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  ValNodePtr        prods, vnp;
  BioseqContextPtr  bcp;
  Int4              insert_point;
  Boolean           recreated_feats = FALSE;
  Boolean           rval = FALSE;
  
  if (sejp == NULL || sejp->bsp == NULL 
      || sejp->char_data == NULL || sejp->num_chars == 0) 
  {
    return FALSE;    
  }
  
  len = BioseqGetLen(sejp->bsp);
  insert_point = sejp->offset;

  if (insert_point == LAST_RESIDUE)
  {
    insert_point = len - 1;
  }
  else if (insert_point == APPEND_RESIDUE) 
  {
    insert_point = len;
  }

  if ((insert_point < 0) || (insert_point > len)) return FALSE;

  if (sejp->action == eSeqEdInsertGap || sejp->action == eSeqEdDeleteGap)
  {
    rval = SeqEdInsertGap (sejp, insert_point);
  }
  else if (sejp->bsp->repr == Seq_repr_raw)
  {
    rval = SeqEdInsertRaw (sejp, insert_point);
  }
  else if (sejp->bsp->repr == Seq_repr_delta)
  {
    rval = SeqEdInsertDelta (sejp, insert_point);
  }
  
  if (!rval)
  {
    return FALSE;
  }
  
  /* fix features */
  if (sejp->entityID > 0 && SeqMgrFeaturesAreIndexed (sejp->entityID)) 
  {
    sfp = NULL;
    while ((sfp = SeqEdGetNextFeature (sejp->bsp, sfp, 0, 0, &fcontext, FALSE, FALSE, sejp->entityID)) != NULL)
    {
      SeqEdInsertAdjustFeat (sfp, sejp, insert_point);
	  }

	  /* adjust features pointing by product */
	  prods = SeqMgrGetSfpProductList (sejp->bsp);
	  for (vnp = prods; vnp != NULL; vnp = vnp->next) 
	  {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (sfp == NULL) continue;
      sfp->product = SeqEdSeqLocInsert (sfp->product, sejp->bsp, insert_point, sejp->num_chars, sejp->spliteditmode, NULL);
	  }
  } else {
    bcp = BioseqContextNew(sejp->bsp);
    sfp = NULL;
    /* adjust features pointing by location */
    while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 0)) != NULL)
    {
      SeqEdInsertAdjustFeat (sfp, sejp, insert_point);
    }
    sfp = NULL;
	  /* adjust features pointing by product */
	  while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 1)) != NULL)
	  {
      sfp->product = SeqEdSeqLocInsert (sfp->product, sejp->bsp, insert_point, sejp->num_chars, sejp->spliteditmode, NULL);
	  }
	  BioseqContextFree(bcp);
  }

	recreated_feats = SeqEdRecreateDeletedFeats (sejp);
  
  if (recreated_feats)
  {
    SeqMgrIndexFeatures (sejp->entityID, NULL);
  }
  else
  {
    SeqEdReindexAffectedFeatures (sejp->offset, sejp->num_chars, 
                                    sejp->spliteditmode, sejp->bsp);
    ReindexExtendedFeatures (sejp);    
  }
  sejp->affected_feats = SeqEdJournalAffectedFeatsFree (sejp->affected_feats);
  return TRUE;
}


/* This section contains code for deleting from sequences and feature locations, adapted from
 * that found in edutil.c */

/*****************************************************************************
*
*   SeqEdSeqFeatDelete()
*     0 = no changes made to location or product
*     1 = changes made but feature still has some location
*     2 = all of sfp->location in deleted interval
*
*   if (merge)
*      1) correct numbers > to by subtraction
*      2) do not split intervals spanning the deletion
*   else
*      1) do not change numbers > to
*      2) split intervals which span the deletions
*
*****************************************************************************/
static Int2 LIBCALL SeqEdSeqFeatDelete (SeqFeatPtr sfp, BioseqPtr target, Int4 from, Int4 to, Boolean merge)
{
	ValNode vn;
	SeqLocPtr tloc;
	SeqInt si;
	Boolean changed = FALSE, tmpbool = FALSE;
	CdRegionPtr crp;
	CodeBreakPtr cbp, prevcbp, nextcbp;
	RnaRefPtr rrp;
	tRNAPtr trp;

	tloc = &vn;
	MemSet((Pointer)tloc, 0, sizeof(ValNode));
	MemSet((Pointer)&si, 0, sizeof(SeqInt));
	tloc->choice = SEQLOC_INT;
	tloc->data.ptrvalue = (Pointer)(&si);
	si.id = target->id;
	si.from = from;
	si.to = to;

	sfp->location = SeqEdSeqLocDelete (sfp->location, target, from, to, merge, &changed, NULL, NULL);
	sfp->product = SeqEdSeqLocDelete(sfp->product, target, from, to, merge, &changed, NULL, NULL);

	if (sfp->location == NULL)
		return 2;

	switch (sfp->data.choice)
	{
		case SEQFEAT_CDREGION:   /* cdregion */
			crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
			prevcbp = NULL;
			for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp)
			{
				nextcbp = cbp->next;
				cbp->loc = SeqEdSeqLocDelete(cbp->loc, target, from, to, merge, &tmpbool, NULL, NULL);
				if (cbp->loc == NULL)
				{
					if (prevcbp != NULL)
						prevcbp->next = nextcbp;
					else
						crp->code_break = nextcbp;
					cbp->next = NULL;
					CodeBreakFree(cbp);
				}
				else
					prevcbp = cbp;
			}
			break;
		case SEQFEAT_RNA:
			rrp = (RnaRefPtr)(sfp->data.value.ptrvalue);
			if (rrp->ext.choice == 2)   /* tRNA */
			{
				trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
				if (trp->anticodon != NULL)
				{
					trp->anticodon = SeqEdSeqLocDelete(trp->anticodon, target, from, to, merge, &tmpbool, NULL, NULL);
				}
			}
			break;
		default:
			break;
	}
			
	if (changed)
	{
		return 1;
	}
	else
		return 0;
}

/*
static Boolean SeqEdDeleteFromDeltaSeq (DeltaSeqPtr dsp, Int4 from, Int4 to)
{
  ByteStorePtr    bs_new;
  SeqLitPtr       slip;

  if (dsp == NULL || dsp->choice != 2 || dsp->data.ptrvalue == NULL)
  {
    return FALSE;
  }

  slip = (SeqLitPtr) dsp->data.ptrvalue;
  
  if (from < 0 || to > slip->length)
  {
    return FALSE;
  }
  if (to < 0)
  {
    to = slip->length - 1;
  }

  if (! IsDeltaSeqGap (dsp))
  {
    if (slip->seq_data_type != Seq_code_iupacna)
    {
      bs_new = BSConvertSeq(slip->seq_data, Seq_code_iupacna, 
                            slip->seq_data_type, 
                            slip->length);
      slip->seq_data_type = Seq_code_iupacna;
      slip->seq_data = bs_new;
    }
    BSSeek(slip->seq_data, from, SEEK_SET);
    Nlm_BSDelete (slip->seq_data, to - from + 1);    
  }
  slip->length -= (to - from + 1);
  
  return TRUE;
}
*/

static void DeleteFromSeqLit (SeqLitPtr slip, Int4 from, Int4 to)
{
  ByteStorePtr    bs_new;

  if (slip == NULL)
  {
    return;
  }
  if (from < 0)
  {
    from = 0;
  }
  
  if (to > slip->length - 1 || to < 0)
  {
    to = slip->length - 1;
  }
  
  if (slip->seq_data != NULL)
  {
    if (slip->seq_data_type != Seq_code_iupacna)
    {
      bs_new = BSConvertSeq(slip->seq_data, Seq_code_iupacna, 
                            slip->seq_data_type, 
                            slip->length);
      slip->seq_data_type = Seq_code_iupacna;
      slip->seq_data = bs_new;
    }
    BSSeek(slip->seq_data, from, SEEK_SET);
    Nlm_BSDelete (slip->seq_data, to - from + 1);    
  }
  slip->length -= (to - from + 1);
}

static Boolean SeqEdDeleteFromDeltaBsp (BioseqPtr bsp, Int4 from, Int4 to)
{
  Boolean         retval = FALSE;
  DeltaSeqPtr     dsp, dsp_next, prev_dsp;
  SeqLitPtr       slip;
  Int4            curr_pos = 0;
  Int4            del_to, del_from;
  Int4            piece_len;
  SeqLocPtr       slp;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL)
  {
    return retval;
  }
    
  prev_dsp = NULL;  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL && curr_pos <= to)
  {
    dsp_next = dsp->next;
    piece_len = 0;
    /* remove empty dsps */
    if (dsp->data.ptrvalue == NULL)
    {
      /* skip */
      prev_dsp = dsp;
    }
		else if (dsp->choice == 1) 
		{  /* SeqLoc */
		  slp = (SeqLocPtr)(dsp->data.ptrvalue);
      piece_len = SeqLocLen (slp);
      prev_dsp = dsp;
		}
		else if (dsp->choice == 2)
		{
		  slip = (SeqLitPtr) (dsp->data.ptrvalue);
		  piece_len = slip->length;
		  if (curr_pos + piece_len > from)
		  {
		    if (from > curr_pos)
		    {
		      del_from = from - curr_pos;
		    }
		    else
		    {
		      del_from = 0;
		    }
		    
		    if (to - curr_pos < slip->length - 1)
		    {
		      del_to = to - curr_pos;
		    }
		    else
		    {
		      del_to = slip->length - 1;
		    }
		    DeleteFromSeqLit (slip, del_from, del_to);
		    
		    /* remove empty delta seq parts */
		    if (slip->length == 0)
		    {
		      if (prev_dsp == NULL)
		      {
		        bsp->seq_ext = dsp->next;
		      }
		      else
		      {
		        prev_dsp->next = dsp->next;
		      }
		      dsp->next = NULL;
		      slip = SeqLitFree (slip);
		      dsp = ValNodeFree (dsp);
		    }
		    else
		    {
		      prev_dsp = dsp;
		    }
		  }
		  else
		  {
		    prev_dsp = dsp;
		  }
		}
		curr_pos += piece_len;
		dsp = dsp_next;
  }
  return TRUE;
}

static Boolean SeqEdDeleteFromSegOrDeltaBsp (BioseqPtr bsp, Int4 from, Int4 to)
{
  SeqLocPtr       tmp, head;
  DeltaSeqPtr     tdsp;
  SeqLocPtr PNTR  newheadptr;
  Int4            totlen, templen, tfrom, tto, diff1, diff2;
  SeqLocPtr       slp, tloc, newhead, prev;
  Boolean         retval = FALSE;
  SeqInt          si;
  ValNode         vn;

  if (bsp == NULL) return retval;
  if (bsp->repr != Seq_repr_seg && bsp->repr != Seq_repr_delta) return retval;

  head = ValNodeNew(NULL);  /* allocate to facilitate SeqLocFree */
  head->choice = SEQLOC_MIX;   /* make a SeqLoc out of the extension */
  if (bsp->repr == Seq_repr_seg)
    head->data.ptrvalue = bsp->seq_ext;
  else
  {
    tdsp = (DeltaSeqPtr)(bsp->seq_ext);
    head->data.ptrvalue = DeltaSeqsToSeqLocs(tdsp);
  }
		
  newhead = NULL;
  newheadptr = &newhead;

  tloc = &vn;
  MemSet((Pointer)tloc, 0, sizeof(ValNode));
  MemSet((Pointer)&si, 0, sizeof(SeqInt));
  tloc->choice = SEQLOC_INT;
  tloc->data.ptrvalue = (Pointer)(&si);
		
  slp = NULL;
  totlen = 0;
  while ((slp = SeqLocFindNext(head, slp)) != NULL)
  {
    templen = SeqLocLen(slp);
    tfrom = SeqLocStart(slp);
    tto = SeqLocStop(slp);
			
    if (((totlen + templen - 1) < from) ||   /* before cut */
	  (totlen > to))						  /* after cut */
	{
      tmp = SeqLocAdd(newheadptr, slp, TRUE, TRUE); /* add whole SeqLoc */
	}
	else                            
    {
      retval = TRUE;    /* will modify or drop interval */
	  diff1 = from - totlen;        /* partial beginning? */
	  diff2 = (templen + totlen - 1) - to;  /* partial end? */
	  si.id = SeqLocId(slp);
	  si.strand = SeqLocStrand(slp);
				
	  if (diff1 > 0)	  /* partial start */
	  {
        if (si.strand != Seq_strand_minus)
	    {
	      si.from = tfrom;
		  si.to = tfrom + diff1 - 1;
		}
	    else
        {
		  si.from = tto - diff1 + 1;
	      si.to = tto;
        }
        tmp = SeqLocAdd(newheadptr, tloc, TRUE, TRUE);
	  }

	  if (diff2 > 0)    /* partial end */
	  {
		if (si.strand != Seq_strand_minus)
		{
	      si.from = tto - diff2 + 1;
          si.to = tto;
		}
		else
		{
          si.from = tfrom;
          si.to = tfrom + diff2 - 1;
        }
		tmp = SeqLocAdd(newheadptr, tloc, TRUE, TRUE);
	  }			
	}
	totlen += templen;
  }

  prev = NULL;
  for (tmp = newhead; tmp != NULL; tmp = tmp->next)
  {
	if (tmp->next == NULL)   /* last one */
	{
	  if (tmp->choice == SEQLOC_NULL)
	  {
		if (prev != NULL)
		  prev->next = NULL;
		else				  /* only a NULL left */
		{
		  newhead = NULL;
		}
		MemFree(tmp);
		break;
	  }
	}
	prev = tmp;
  }

  if (bsp->repr == Seq_repr_seg)
	bsp->seq_ext = newhead;
  else
  {
    bsp->seq_ext = SeqLocsToDeltaSeqs(tdsp, newhead);
	DeltaSeqSetFree(tdsp);
	SeqLocSetFree(newhead);
  }
  SeqLocFree(head);
  return TRUE;
}

static Boolean SeqEdDeleteFromMapBioseq (BioseqPtr bsp, Int4 from, Int4 to)
{
  SeqFeatPtr sfpcurr, sfpnext, sfpprev;
  Int2 dropped;

  if (bsp == NULL || bsp->repr != Seq_repr_map) return FALSE;

  sfpprev = NULL;
  sfpnext = NULL;
  sfpcurr = (SeqFeatPtr)(bsp->seq_ext);
  bsp->seq_ext = NULL;
  for (; sfpcurr != NULL; sfpcurr = sfpnext)
  {
	sfpnext = sfpcurr->next;
	dropped = SeqEdSeqFeatDelete(sfpcurr, bsp, from, to, TRUE);
	if (dropped == 2)   /* completely gone */
	{
	  SeqFeatFree(sfpcurr);
	}
	else
	{
	  if (sfpprev == NULL)
		bsp->seq_ext = (Pointer)sfpcurr;
	  else
		sfpprev->next = sfpcurr;
	  sfpcurr->next = NULL;
	  sfpprev = sfpcurr;
	}
  }
  return TRUE;
}

static SeqLocPtr FreeSeqLocList (SeqLocPtr slp)
{
  if (slp == NULL)
  {
    return NULL;
  }
  slp->next = SeqLocFree (slp->next);
  slp = SeqLocFree (slp);
  return slp;
}

static Boolean ReStitchLocation (Int4 delete_point, SeqFeatPtr sfp)
{
  Int4      this_start, this_stop, next_start, next_stop;
  SeqLocPtr this_slp, next_slp, loc_list = NULL, tmp_slp, last_slp = NULL, tmp_next;
  SeqIdPtr  this_id, next_id;
  Boolean   merged = FALSE;
  Uint2     this_strand, next_strand;
  
  if (sfp->location == NULL)
  {
    return FALSE;
  }
  
  this_start = SeqLocStart (sfp->location);
  this_stop = SeqLocStop (sfp->location);
  if (delete_point <= this_start || delete_point >= this_stop)
  {
    return FALSE;
  }
  
  this_slp = SeqLocFindNext (sfp->location, NULL);
  if (this_slp == NULL)
  {
    return FALSE;
  }
  next_slp = SeqLocFindNext (sfp->location, this_slp);
  
  while (next_slp != NULL)
  {
    this_start = SeqLocStart (this_slp);
    this_stop = SeqLocStop (this_slp);
    this_id = SeqLocId (this_slp);
    this_strand = SeqLocStrand (this_slp);
    next_start = SeqLocStart (next_slp);
    next_stop = SeqLocStop (next_slp);
    next_id = SeqLocId (next_slp);
    next_strand = SeqLocStrand (next_slp);
    if (this_stop + 1 == next_start 
        && next_start == delete_point
        && SeqIdComp (this_id, next_id) == SIC_YES
        && this_strand == next_strand)
    {
      tmp_slp = SeqLocIntNew (this_start, next_stop, this_strand, this_id);
      next_slp = SeqLocFindNext (sfp->location, next_slp);
      merged = TRUE;
    }
    else
    {
      tmp_next = this_slp->next;
      this_slp->next = NULL;
      tmp_slp = SeqLocCopy (this_slp);
      this_slp->next = tmp_next;
    }
    if (tmp_slp != NULL)
    {
      if (last_slp == NULL)
      {
        loc_list = tmp_slp;
      }
      else
      {
        last_slp->next = tmp_slp;
      }
      last_slp = tmp_slp;
    }
    
    this_slp = next_slp;
    if (this_slp != NULL)
    {
      next_slp = SeqLocFindNext (sfp->location, this_slp);
    }
  }
  if (merged && loc_list != NULL)
  {
    if (this_slp != NULL)
    {
      this_start = SeqLocStart (this_slp);
      this_stop = SeqLocStop (this_slp);
      tmp_next = this_slp->next;
      this_slp->next = NULL;
      tmp_slp = SeqLocCopy (this_slp);
      this_slp->next = tmp_next;
      if (last_slp == NULL)
      {
        loc_list = tmp_slp;
      }
      else
      {
        last_slp->next = tmp_slp;
      }
    }
    if (loc_list->next == NULL)
    {
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = loc_list;
    }
    else
    {
      /* already mix, just need to replace list */
      sfp->location->data.ptrvalue = FreeSeqLocList (sfp->location->data.ptrvalue);
      sfp->location->data.ptrvalue = loc_list;
    }
    return TRUE;
  }
  else
  {
    loc_list = FreeSeqLocList (loc_list);
    return FALSE;
  }
}

/* ideally, this should take a SeqJournalEntry and perform the deletion.
 * We will always be deleting a contiguous section of characters.
 * This function will only delete from the specified Bioseq, so there should 
 * be no need to call BioseqFind (which is expensive).
 */
NLM_EXTERN Boolean SeqEdDeleteFromBsp (SeqEdJournalPtr sejp, BoolPtr pfeats_deleted)
{
  Boolean           retval = FALSE;
  Boolean           feats_altered = FALSE;
  Int4              deleted;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  BioseqContextPtr  bcp;
  Int2              feat_change;
  Int2              feats_deleted = FALSE;
  SeqFeatPtr        tmp_sfp;
  AffectedFeatPtr   afp;
  Boolean           merge_mode;
  Boolean           location_restitched = FALSE;

  if (sejp == NULL || sejp->bsp == NULL || sejp->offset < 0 || sejp->offset >= sejp->bsp->length
	  || sejp->offset + sejp->num_chars + 1 < 0 || sejp->offset + sejp->num_chars > sejp->bsp->length
	  || sejp->num_chars < 1)
  {
    return retval;
  }
  
  if (sejp->affected_feats != NULL)
  {
    sejp->affected_feats = SeqEdJournalAffectedFeatsFree (sejp->affected_feats);
  }
  
  /* fix features */
  if (sejp->entityID > 0 && SeqMgrFeaturesAreIndexed (sejp->entityID)) {
    sfp = NULL;
    while ((sfp = SeqEdGetNextFeature (sejp->bsp, sfp, 0, 0, &fcontext, FALSE, FALSE, sejp->entityID)) != NULL)
    {
      if ((sejp->offset <= fcontext.left && sejp->offset + sejp->num_chars >= fcontext.left)
          || (sejp->offset >= fcontext.left && sejp->offset + sejp->num_chars <= fcontext.right)
          || (sejp->offset <= fcontext.right && sejp->offset + sejp->num_chars >= fcontext.right))
      {
        tmp_sfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);        
      }
      else
      {
        tmp_sfp = NULL;
      }
      /* if we're deleting a gap and the feature is a coding region, merge location
       * by default */
      merge_mode = sejp->spliteditmode;
      if (sejp->action == eSeqEdInsertGap || sejp->action == eSeqEdDeleteGap
          && sfp->data.choice == SEQFEAT_CDREGION)
      {
        merge_mode = TRUE;
      }

      feat_change = SeqEdSeqFeatDelete (sfp, sejp->bsp, sejp->offset, 
                                        sejp->offset + sejp->num_chars - 1,
                                        sejp->spliteditmode);
                                        
      if (feat_change == 0 || feat_change == 1)
      {
        if (ReStitchLocation (sejp->offset, sfp))
        {
          feat_change = 1;
          location_restitched = TRUE;
        }
      }
      
      if (feat_change > 0)
	    {
	      if (feat_change == 2)
	      {
	        /* remove from index and SeqAnnot */
	        sfp->idx.deleteme = TRUE;
	        feats_deleted = TRUE;
	      }

	      afp = AffectedFeatNew ();
	      if (afp != NULL)
	      {
	        afp->feat_before = tmp_sfp;
	        if (feat_change != 2)
	        {
	          afp->feat_after = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
	          if (afp->feat_after != NULL)
	          {
	            afp->feat_after->idx.itemID = sfp->idx.itemID;	          
	          }	          
	        }
	      }
	      ValNodeAddPointer (&sejp->affected_feats, 0, afp);
	      feats_altered = TRUE;
	    }
	    else
	    {
	      SeqFeatFree (tmp_sfp);
	    }
	  }

  } else {
    bcp = BioseqContextNew(sejp->bsp);
    sfp = NULL;
    /* adjust features pointing by location */
    while ((sfp = BioseqContextGetSeqFeat(bcp, 0, sfp, NULL, 0)) != NULL)
    {
      tmp_sfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
      /* if we're deleting a gap and the feature is a coding region, merge location
       * by default */
      merge_mode = sejp->spliteditmode;
      if (sejp->action == eSeqEdInsertGap || sejp->action == eSeqEdDeleteGap
          && sfp->data.choice == SEQFEAT_CDREGION)
      {
        merge_mode = TRUE;
      }      
      feat_change = SeqEdSeqFeatDelete (sfp, sejp->bsp, sejp->offset, 
                                        sejp->offset + sejp->num_chars - 1, 
                                        sejp->spliteditmode);

      if (feat_change == 0 || feat_change == 1)
      {
        if (ReStitchLocation (sejp->offset, sfp))
        {
          feat_change = 1;
          location_restitched = TRUE;
        }
      }

      if (feat_change > 0)
	    {
	      if (feat_change == 2)
	      {
	        /* remove from index and SeqAnnot */
	        sfp->idx.deleteme = TRUE;
	        feats_deleted = TRUE;
	      }
	      afp = AffectedFeatNew ();
	      if (afp != NULL)
	      {
	        afp->feat_before = tmp_sfp;
	        afp->feat_after = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
	      }
	      ValNodeAddPointer (&sejp->affected_feats, 0, afp);
	      feats_altered = TRUE;
	    }
	    else
	    {
	      SeqFeatFree (tmp_sfp);
	    }
    }
	  BioseqContextFree(bcp);
  }

  /* now delete nucleotides from bioseq */
  switch (sejp->bsp->repr)
  {
    case Seq_repr_raw:
    case Seq_repr_const:
      /* if actual sequence present */
      if (ISA_na(sejp->bsp->mol))
      {
        if (sejp->bsp->seq_data_type != Seq_code_iupacna)  /* need 1 byte/base */
          BioseqRawConvert(sejp->bsp, Seq_code_iupacna);
	    }
      else
      {
        if (sejp->bsp->seq_data_type != Seq_code_ncbieaa)
          BioseqRawConvert(sejp->bsp, Seq_code_ncbieaa);
      }

      BSSeek(sejp->bsp->seq_data, sejp->offset, SEEK_SET);
      deleted = BSDelete(sejp->bsp->seq_data, sejp->num_chars);
      if (deleted != sejp->num_chars)  /* error */
        ErrPost(CTX_NCBIOBJ, 1, "Delete of %ld residues failed", sejp->num_chars);
      else
        retval = TRUE;
      break;
    case Seq_repr_seg:
      /* update segmented sequence */
      retval = SeqEdDeleteFromSegOrDeltaBsp (sejp->bsp, sejp->offset, sejp->offset + sejp->num_chars - 1);
      break;
    case Seq_repr_delta:
	    /* update delta sequence */
      retval = SeqEdDeleteFromDeltaBsp (sejp->bsp, sejp->offset, sejp->offset + sejp->num_chars - 1);
      break;
    case Seq_repr_map:
      /* map bioseq */
      retval = SeqEdDeleteFromMapBioseq (sejp->bsp, sejp->offset, sejp->offset + sejp->num_chars - 1);
      break;
    case Seq_repr_virtual:
      retval = TRUE;                 /* nothing to do */
      break;
  }

  if (retval)
    sejp->bsp->length -= sejp->num_chars;
  
  if (feats_deleted)
  {
    DeleteMarkedObjects (sejp->entityID, 0, NULL);
    SeqMgrIndexFeatures (sejp->entityID, NULL);
  }
  else if (location_restitched)
  {
    SeqMgrIndexFeatures (sejp->entityID, NULL);
  }
  else
  {
    SeqEdReindexAffectedFeatures (sejp->offset, 0 - sejp->num_chars, 
                                  sejp->spliteditmode, sejp->bsp);
    
  }
  
  if (pfeats_deleted != NULL)
  {
    *pfeats_deleted = feats_deleted;
  }

  return retval;
}

/* this function will indicate whether the interval on the Bioseq specified contains
 * any gaps of unknown length.
 */

/*
static Boolean DoesIntervalContainUnknownGap (BioseqPtr bsp, Int4 from, Int4 to)
{
  DeltaSeqPtr from_dsp, to_dsp, this_dsp;
  Int4        from_start = 0, to_start = 0;
  Boolean     unknown_gap = FALSE;
  
  if (bsp == NULL || from < 0 || from >= bsp->length || to < 0 || to >= bsp->length)
  {
    return FALSE;
  }
  
  from_dsp = GetDeltaSeqForOffset (bsp, from, &from_start);
  to_dsp = GetDeltaSeqForOffset (bsp, to, &to_start);
  
  this_dsp = from_dsp;
  while (!unknown_gap && this_dsp != NULL && (to_dsp == NULL || this_dsp != to_dsp->next))
  {
    unknown_gap = IsDeltaSeqUnknownGap (this_dsp);
    this_dsp = this_dsp->next;
  }
  
  return unknown_gap;
}
*

/* This section of code deals with editing the sequence by inserting and removing characters.
 * Functions are needed to change the indices for the affected features so that they will
 * display properly.
 */
static void SeqEdFixExtraIndex
(SMFeatItemPtr PNTR array,
 Int4               num,
 Int4               shift_start,
 Int4               shift_amt,
 Boolean            split,
 BioseqPtr          bsp)
{
  SMFeatItemPtr       item;
  Int4                i = 0, j, k, n;
  Int4Ptr             newivals;
	
  if (array == NULL || num < 1 || bsp == NULL) return;  
  while (i < num) {
    item = array [i];
    i++;
    if (item != NULL) {
      if (item->right >= shift_start)
      {
      	if (item->left > shift_start 
      	    || (shift_amt > 0 && item->left == shift_start))
      	{
      	  /* move left and right indexed endpoints */
      	  item->left += shift_amt;
      	  if (item->left < 0)
      	  {
      	    item->left = 0;
      	  }
      	  item->right += shift_amt;
      	  /* move all ivals */
      	  for (j = 0; j < item->numivals; j++)
      	  {
      	  	item->ivals [2 * j] += shift_amt;
      	  	if (item->ivals [2 * j] < 0)
      	  	{
      	  	  item->ivals [2 * j] = 0;
      	  	}
      	  	item->ivals [2 * j + 1] += shift_amt;
      	  	if (item->ivals [2 * j + 1] < 0)
      	  	{
      	  	  item->ivals [2 * j + 1] = 0;
      	  	}
      	  }
      	}
      	else
      	{
      	  item->right += shift_amt;
          for (j = 0; j < item->numivals; j++)
          {
          	if (item->ivals [2 * j] < shift_start && item->ivals[2 * j + 1] < shift_start)
          	{
          	  /* upstream - we may safely ignore */
          	}
            else if ((item->ivals [2 * j] > shift_start && item->ivals [2 * j + 1] > shift_start)
                    || (shift_amt > 0 && item->ivals [2 * j] >= shift_start 
                                      && item->ivals [2 * j + 1] >= shift_start))
            {
              /* downstream - shift both endpoints */
          	  item->ivals [2 * j] += shift_amt;
          	  item->ivals [2 * j + 1] += shift_amt;
            }
      	    else if (split)
      	    {
      	      /* create a new list of ivals */
      	  	  newivals = (Int4Ptr) MemNew (sizeof (Int4) * (item->numivals + 1) * 2);
      	  	  /* copy all ivals up to j into new list */
      	  	  for (k = 0; k < j; k++)
      	  	  {
      	  	    newivals [2 * k] = item->ivals [2 * k];
      	  	    newivals [2 * k + 1] = item->ivals [2 * k + 1];
      	  	  }
      	  	  /* create two intervals using split */
      	  	  if (item->ivals [2 * j] < item->ivals [2 * j + 1])
      	  	  {
      	  	    /* plus strand */
      	  	    newivals [2 * k] = item->ivals [2 * j];
      	  	    newivals [2 * k + 1] = shift_start - 1;
      	  	    k++;
      	  	    newivals [2 * k] = shift_start + shift_amt;
      	  	    newivals [2 * k + 1] = item->ivals [2 * j + 1] + shift_amt;
      	  	    k++;
      	  	  }
      	  	  else
      	  	  {
      	  	    /* minus strand */
      	  	    newivals [2 * k] = item->ivals [2 * j] + shift_amt;
      	  	    newivals [2 * k + 1] = shift_start + shift_amt;
      	  	    k++;
      	  	    newivals [2 * k] = shift_start - 1;
      	  	    newivals [2 * k + 1] = item->ivals [2 * j + 1];
      	  	    k++;      	  	  	
      	  	  }
      	  	  /* copy remaining intervals (they will be shifted later in the loop */
      	  	  n = j + 1;
      	  	  while (n < item->numivals)
      	  	  {
      	  	    newivals[2 * k] = item->ivals [2 * n];
      	  	    newivals[2 * k + 1] = item->ivals [2 * n + 1];
      	  	    k++;
      	  	    n++;
      	   	  }
      	  	  MemFree (item->ivals);
      	  	  item->ivals = newivals;
      	  	  item->numivals ++;  
      	  	  /* increment j so that we will not re-increment the second interval */
      	  	  j++;	  	    	    	     
      	    }
      	    else
      	    {
      	      /* move only downstream endpoint */
      	      if (item->ivals [2 * j] > shift_start
      	          || (shift_amt > 0 && item->ivals [2 * j] == shift_start))
      	      {
      	      	item->ivals [2 * j] += shift_amt;
      	      	if (item->ivals [2 * j] < 0)
      	      	{
      	      	  item->ivals [2 * j] = 0;
      	      	}
      	      }
      	      else
      	      {
      	      	item->ivals [2 * j + 1] += shift_amt;
      	      	if (item->ivals [2 * j + 1] < 0)
      	      	{
      	      	  item->ivals [2 * j + 1] = 0;
      	      	}
      	      }
      	    }
          }
        }
      }
    }
  }
}

NLM_EXTERN void SeqEdReindexAffectedFeatures (Int4 shift_start, Int4 shift_amt, 
                                          Boolean split, BioseqPtr bsp)
{
  ObjMgrDataPtr       omdp;
  BioseqExtraPtr      bspextra;
  ObjMgrPtr           omp;

  if (bsp == NULL) return;
  
  omdp = (ObjMgrDataPtr) bsp->omdp;
  if (omdp == NULL) 
  {
    omp = ObjMgrWriteLock ();
    omdp = ObjMgrFindByData (omp, bsp);
    ObjMgrUnlock ();
    bsp->omdp = (Pointer) omdp;	
  }
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;
  
  SeqEdFixExtraIndex (bspextra->featsByPos, bspextra->numfeats,
                      shift_start, shift_amt, split, bsp);
}

NLM_EXTERN void SeqEdReindexFeature (SeqFeatPtr sfp, BioseqPtr bsp)
{
  ObjMgrDataPtr       omdp;
  BioseqExtraPtr      bspextra;
  ObjMgrPtr           omp;
  Int4                i;
  SeqLocPtr           this_slp;
  SMFeatItemPtr       item = NULL;
  Int4                numivals;
  Int4                start, stop;
  Int4                left, right;

  if (sfp == NULL || bsp == NULL) return;
  omdp = (ObjMgrDataPtr) bsp->omdp;
  if (omdp == NULL) 
  {
    omp = ObjMgrWriteLock ();
    omdp = ObjMgrFindByData (omp, bsp);
    ObjMgrUnlock ();
    bsp->omdp = (Pointer) omdp;	
  }
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;

  for (i = 0; i < bspextra->numfeats; i++)
  {
    item = bspextra->featsByPos [i];
    if (item != NULL && item->itemID == sfp->idx.itemID) 
    {
      /* first, find out how many intervals we have, so we can make sure our ivals
       * array is the right size */
      for (this_slp = SeqLocFindNext (sfp->location, NULL), numivals = 0;
           this_slp != NULL;
           this_slp = this_slp->next, numivals ++)
      {
      	
      }
      if (numivals != item->numivals)
      {
      	item->ivals = MemFree (item->ivals);
      	item->ivals = (Int4Ptr) MemNew (2 * numivals * sizeof (Int4));
      	if (item->ivals == NULL) return;
      	item->numivals = numivals;
      }
      
      /* now populate the ivals */
      
      left = -1;
      right = -1;
      for (this_slp = SeqLocFindNext (sfp->location, NULL), numivals = 0;
           this_slp != NULL;
           this_slp = this_slp->next, numivals ++)
      {
      	start = GetOffsetInBioseq (this_slp, bsp, SEQLOC_START);
      	stop = GetOffsetInBioseq (this_slp, bsp, SEQLOC_STOP);
      	item->ivals [2 * numivals] = start;
      	item->ivals [2 * numivals + 1] = stop;
      	if (left == -1 || start < left)
      	{
      	  left = start;
      	}
      	if (stop < left)
      	{
      	  left = stop;
      	}
      	if (right == -1 || right < start)
      	{
      	  right = start;
      	}
      	if (right < stop)
      	{
      	  right = stop;
      	}
      }
      item->left = left;
      item->right = right;      
    }
    else
    {
      item = NULL;
    }
  }
}


/* This function will repair any problems with the interval order that
 * moving the feature interval around may have caused.
 */
NLM_EXTERN void SeqEdRepairIntervalOrder (SeqFeatPtr sfp, BioseqPtr bsp)
{
  Boolean   hasNulls;
  SeqLocPtr gslp;
  Boolean   noLeft, noRight;
  
  hasNulls = LocationHasNullsBetween (sfp->location);
  gslp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, hasNulls);
  if (gslp != NULL) 
  {
    CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = gslp;
    if (bsp->repr == Seq_repr_seg) 
    {
      gslp = SegLocToParts (bsp, sfp->location);
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = gslp;
    }
    FreeAllFuzz (sfp->location);
    SetSeqLocPartial (sfp->location, noLeft, noRight);
  }
}

/* This function recursively frees a list of SeqEdJournalPtr, working in the next direction,
 * and fixes the prev pointer for the previous entry in the SeqEdJournalPtr list (if there is one).
 */
NLM_EXTERN void SeqEdJournalFree (SeqEdJournalPtr sejp)
{
  SeqEdJournalPtr prev;
  
  if (sejp == NULL) return;
  SeqEdJournalFree (sejp->next);
  sejp->slp = SeqLocFree (sejp->slp);
  MemFree (sejp->char_data);
  sejp->affected_feats = SeqEdJournalAffectedFeatsFree (sejp->affected_feats);
  prev = sejp->prev;
  if (prev != NULL)
    prev->next = NULL;
  MemFree (sejp);
}

NLM_EXTERN SeqEdJournalPtr SeqEdJournalNewSeqEdit 
(ESeqEdJournalAction action,
 Int4                offset,
 Int4                num_chars,
 CharPtr             char_data,
 Boolean             spliteditmode,
 BioseqPtr           bsp,
 Uint2               moltype,
 Uint2               entityID)
{
  SeqEdJournalPtr sejp;
  
  if (num_chars == 0) return NULL;
  sejp = (SeqEdJournalPtr) MemNew (sizeof (SeqEdJournalData));
  if (sejp == NULL) return NULL;
  sejp->action = action;
  sejp->offset = offset;
  sejp->num_chars = num_chars;
  sejp->spliteditmode = spliteditmode;
  sejp->affected_feats = NULL;
  sejp->sfp = NULL;
  sejp->slp = NULL;
  sejp->bsp = bsp;
  sejp->moltype = moltype;
  sejp->entityID = entityID;
  sejp->char_data = MemNew (sejp->num_chars + 1);
  if (char_data != NULL)
  {
    StringCpy (sejp->char_data, char_data);   	
  }
  sejp->prev = NULL;
  sejp->next = NULL;
  return sejp;
}

NLM_EXTERN SeqEdJournalPtr SeqEdJournalNewFeatEdit
(ESeqEdJournalAction action,
 SeqFeatPtr          sfp,
 SeqLocPtr           slp,
 BioseqPtr           bsp,
 Uint2               moltype,
 Uint2               entityID)
{
  SeqEdJournalPtr sejp;
  
  if (sfp == NULL || slp == NULL) return NULL;
  sejp = (SeqEdJournalPtr) MemNew (sizeof (SeqEdJournalData));
  if (sejp == NULL) return NULL;
  sejp->action = action;
  sejp->offset = 0;
  sejp->num_chars = 0;
  sejp->spliteditmode = FALSE;
  sejp->sfp = sfp;
  sejp->slp = slp;
  sejp->bsp = bsp;
  sejp->affected_feats = NULL;
  sejp->moltype = moltype;
  sejp->entityID = entityID;
  sejp->char_data = NULL;
  sejp->prev = NULL;
  sejp->next = NULL;
  return sejp;	
}

/* This section of code contains functions used by the new sequence editor for moving feature
 * intervals.
 */
static Boolean SeqEdAdjustFeatureInterval
(SeqLocPtr slp, Int4 change, EMoveType move_type, Int4 interval_offset, BioseqPtr bsp)
{
  SeqIntPtr sint;
  SeqPntPtr spp;
  SeqLocPtr this_slp;
  Boolean   rval = FALSE;
  
  if (slp == NULL || bsp == NULL) return rval;

  if (slp->choice == SEQLOC_INT)
  {
  	if (interval_offset != 0)
  	{
  	  return rval;
  	}
    sint = (SeqIntPtr)slp->data.ptrvalue;
    switch (move_type)
    {
      case eLeftEnd:
        if (sint->from + change < sint->to 
            && sint->from + change > -1
            && sint->from + change < bsp->length)
        {
          sint->from += change;
          rval = TRUE;
        }
        break;
      case eRightEnd:
        if (sint->to + change > sint->from
            && sint->to + change > -1
            && sint->to + change < bsp->length)
        {
      	  sint->to += change;
      	  rval = TRUE;
        }
        break;
      case eSlide:
        if (sint->from + change > -1 && sint->from + change < bsp->length
          && sint->to + change > -1 && sint->to + change < bsp->length)
        {
          sint->from += change;
          sint->to += change;
          rval = TRUE;
        }
    }
  }
  else if (slp->choice == SEQLOC_PNT)
  {
    if (interval_offset != 0)
    {
      return rval;
    }
    spp = (SeqPntPtr)(slp->data.ptrvalue);
    if (spp->point + change > -1 && spp->point + change < bsp->length)
    {
      spp->point += change;
      rval = TRUE;
    }
  }
  else
  {
    for (this_slp = SeqLocFindNext (slp, NULL);
         this_slp != NULL && interval_offset > 0;
         this_slp = SeqLocFindNext (slp, this_slp), interval_offset --)
    {}
    if (this_slp != NULL && interval_offset == 0)
    {
      rval = SeqEdAdjustFeatureInterval (this_slp, change, move_type, interval_offset, bsp);
    }
  }
  return rval;
}


NLM_EXTERN Boolean SeqEdGetNthIntervalEndPoints 
(SeqLocPtr slp, Int4 n, Int4Ptr left, Int4Ptr right)
{
  Boolean rval = FALSE;
  SeqIntPtr sintp;
  SeqPntPtr spp;
  SeqLocPtr this_slp;
  
  if (slp == NULL || left == NULL || right == NULL || n < 0) return FALSE;
  switch (slp->choice)
  {
  	case SEQLOC_INT:
  	  if (n == 0)
  	  {
  	    sintp = (SeqIntPtr) slp->data.ptrvalue;
  	  	*left = sintp->from;
  	  	*right = sintp->to;
  	  	rval = TRUE;
  	  }
  	  break;
  	case SEQLOC_PNT:
  	  if (n == 0)
  	  {
  	  	spp = (SeqPntPtr) slp->data.ptrvalue;
  	  	*left = spp->point;
  	  	*right = spp->point;
  	  	rval = TRUE;
  	  }
  	  break;
  	default:
      for (this_slp = SeqLocFindNext (slp, NULL);
           this_slp != NULL && n > 0;
           this_slp = SeqLocFindNext (slp, this_slp), n --)
      {}
      if (this_slp != NULL && n == 0)
      {
        rval = SeqEdGetNthIntervalEndPoints (this_slp, n, left, right);
      }
      break;
  }
  return rval;
}

static void 
SeqEdFixFeatureIndexForFeatureLocAdjust 
(BioseqPtr  bsp,
 SeqFeatPtr sfp, 
 Int4       change, 
 Int4       move_type,
 Int4       interval_offset)
{
  ObjMgrDataPtr       omdp;
  BioseqExtraPtr      bspextra;
  ObjMgrPtr           omp;
  SMFeatItemPtr       item;
  Int4                i, j;
  Int4                left, right;
  
  if (bsp == NULL || sfp == NULL) return;
  
  omdp = (ObjMgrDataPtr) bsp->omdp;
  if (omdp == NULL) 
  {
    omp = ObjMgrWriteLock ();
    omdp = ObjMgrFindByData (omp, bsp);
    ObjMgrUnlock ();
    bsp->omdp = (Pointer) omdp;	
  }
  if (omdp == NULL || omdp->datatype != OBJ_BIOSEQ) return;

  bspextra = (BioseqExtraPtr) omdp->extradata;
  if (bspextra == NULL) return;

  if (! SeqEdGetNthIntervalEndPoints (sfp->location, interval_offset, &left, &right))
  {
  	return;
  }
  
  i = 0;
  while (i < bspextra->numfeats) {
    item = bspextra->featsByPos [i];
    i++;
    if (item != NULL && item->itemID == sfp->idx.itemID) 
    {
      if (interval_offset >= item->numivals || interval_offset < 0) return;
      if (item->ivals [ 2 * interval_offset] < item->ivals [2 * interval_offset + 1])
      {
      	item->ivals [2 * interval_offset] = left;
      	item->ivals [2 * interval_offset + 1] = right;
      }
      else
      {
      	item->ivals [2 * interval_offset + 1] = left;
      	item->ivals [2 * interval_offset] = right;
      }
      /* correct item left and right values */ 
      if (item->ivals [0] > item->ivals [1])
      {
      	item->right = item->ivals [0];
      	item->left = item->ivals [1];
      }
      else
      {
      	item->left = item->ivals [0];
      	item->right = item->ivals [1];
      }
      for (j = 1; j < item->numivals; j++)
      {
      	if (item->left > item->ivals[2 * j])
      	{
      	  item->left = item->ivals [2 * j];
      	}
      	if (item->left > item->ivals [2 * j + 1])
      	{
      	  item->left = item->ivals [2 * j + 1];
      	}
      	if (item->right < item->ivals [2 * j])
      	{
      	  item->right = item->ivals [2 * j];
      	}
      	if (item->right < item->ivals [2 * j + 1])
      	{
      	  item->right = item->ivals [ 2 * j + 1];
      	}
      }
    }
  }
}


NLM_EXTERN void SeqEdFeatureAdjust
(SeqFeatPtr sfp,
 SeqLocPtr  orig_loc,
 Int4       change,
 EMoveType  move_type,
 Int4       interval_offset,
 BioseqPtr  bsp)
{
  SeqLocPtr new_loc;
  Boolean   partial3, partial5;

  if (sfp == NULL || bsp == NULL)
  {
  	return;
  }
  
  CheckSeqLocForPartial (orig_loc, &partial5, &partial3);
  new_loc = SeqLocMerge (bsp, orig_loc, NULL, FALSE, FALSE, FALSE);
  if (new_loc == NULL)
  {
  	return;
  }
  SetSeqLocPartial (new_loc, partial5, partial3);
  
  if (SeqEdAdjustFeatureInterval (new_loc, change, move_type, interval_offset, bsp))
  {
    SeqLocFree (sfp->location);
    sfp->location = new_loc;

    /* need to reindex feature */
    SeqEdFixFeatureIndexForFeatureLocAdjust (bsp, sfp, change, move_type, interval_offset);
  }
}


NLM_EXTERN void 
AdjustFeatureForGapChange 
(SeqFeatPtr sfp,
 BioseqPtr  bsp, 
 Int4       offset, 
 Int4       len_diff)
{
  if (sfp == NULL || bsp == NULL || offset < 0 || len_diff == 0)
  {
    return;
  }
  
  if (len_diff > 0)
  {
    SeqEdSeqFeatDelete (sfp, bsp, offset, offset + len_diff - 1, TRUE);
  }
  else
  {
    sfp->location = SeqEdSeqLocInsert (sfp->location, bsp, offset, -len_diff, FALSE, NULL);
    if (sfp->data.choice == SEQFEAT_CDREGION)
	  {
      SeqEdInsertAdjustCdRgn (sfp, bsp, offset, -len_diff, FALSE);
	  }
	  else if (sfp->data.choice == SEQFEAT_RNA)
	  {
      SeqEdInsertAdjustRNA (sfp, bsp, offset, -len_diff, FALSE);
	  }
  }    
}




