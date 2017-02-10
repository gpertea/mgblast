static char const rcsid[] = "$Id: blfmtutl.c,v 1.17 2006/04/26 12:42:36 madden Exp $";

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

File name: blfmtutl.c

Author: Tom Madden

Contents: Utilities for BLAST formatting

******************************************************************************/
/*
* $Revision: 
* $Log: blfmtutl.c,v $
* Revision 1.17  2006/04/26 12:42:36  madden
* BlastSetUserErrorString and BlastDeleteUserErrorString moved from blastool.c to blfmtutl.c
*
* Revision 1.16  2006/04/07 19:46:59  coulouri
* correction to previous commit
*
* Revision 1.15  2006/04/07 18:38:19  coulouri
* bump version
*
* Revision 1.14  2006/01/24 18:38:47  papadopo
* from Mike Gertz: Fixed a typo in a name in a format string: Aravaind -> Aravind
*
* Revision 1.13  2005/12/29 19:55:04  madden
* Added functions to print tabular output
*
* Revision 1.12  2005/11/22 13:44:24  coulouri
* bump version
*
* Revision 1.11  2005/10/17 12:47:30  camacho
* From Alejandro Schaffer: Updated reference for compositional adjustment
*
* Revision 1.10  2005/08/05 12:10:48  coulouri
* bump version
*
* Revision 1.9  2005/07/25 12:48:39  camacho
* Updated reference for compositional adjustment
*
* Revision 1.8  2005/06/05 02:54:41  coulouri
* bump date
*
* Revision 1.7  2005/05/20 15:28:03  coulouri
* bump date
*
* Revision 1.6  2005/05/16 17:42:19  papadopo
* From Alejandro Schaffer: Print references for composition-based statistics
* and for compositional score matrix adjustment, if either method was used.
*
* Revision 1.5  2005/05/08 13:32:52  coulouri
* bump version to 2.2.11
*
* Revision 1.4  2004/10/19 15:28:59  coulouri
* bump version and date
*
* Revision 1.3  2004/10/04 17:54:14  madden
* BlastPrintVersionInfo[Ex] now takes const char* as arg for program
*
* Revision 1.2  2004/07/22 15:18:45  jianye
* correct blast paper url
*
* Revision 1.1  2004/06/30 12:31:15  madden
* Structures and prototypes for blast formatting utilities
*
*/

#include <ncbi.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <readdb.h>
#include <ncbithr.h>
#include <txalign.h>
#include <blfmtutl.h>


/* the version of BLAST. */
#define BLAST_ENGINE_VERSION "2.2.14"
#define BLAST_RELEASE_DATE "Apr-09-2006"

#define BUFFER_LENGTH 255

/*
	adds the new string to the buffer, separating by a tilde.
	Checks the size of the buffer for FormatBlastParameters and
	allocates longer replacement if needed.
*/

Boolean LIBCALL
add_string_to_bufferEx(CharPtr buffer, CharPtr *old, Int2Ptr old_length, Boolean add_tilde)

{
	CharPtr new, ptr;
	Int2 length, new_length;

	length = (StringLen(*old));

	if((StringLen(buffer)+length+3) > *old_length)
	{
		new_length = *old_length + 255;
		new = MemNew(new_length*sizeof(Char));
		if (*old_length > 0 && *old != NULL)
		{
			MemCpy(new, *old, *old_length);
			*old = MemFree(*old);
		}
		*old = new;
		*old_length = new_length;
	}

	ptr = *old;
	ptr += length;
	if (add_tilde)
	{
		*ptr = '~';
		ptr++;
	}

	while (*buffer != NULLB)
	{
		*ptr = *buffer;
		buffer++; ptr++;
	}

	return TRUE;
}

Boolean LIBCALL
add_string_to_buffer(CharPtr buffer, CharPtr *old, Int2Ptr old_length)

{
	return add_string_to_bufferEx(buffer, old, old_length, TRUE);
}

/*
	Print the buffer, adding newlines where tildes are found.
*/

Boolean LIBCALL
PrintTildeSepLines(CharPtr buffer, Int4 line_length, FILE *outfp)

{
	if (outfp == NULL || buffer == NULL)
		return FALSE;

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	while (*buffer != NULLB)
	{
		if (*buffer != '~')
			ff_AddChar(*buffer);
		else
			NewContLine();
		buffer++;
	}
	ff_EndPrint();

	return TRUE;
}

/*
	Print the Karlin-Altschul parameters.

	if gapped is TRUE, then slightly different formatting is used.
*/

Boolean LIBCALL
PrintKAParameters(Nlm_FloatHi Lambda, Nlm_FloatHi K, Nlm_FloatHi H, Int4 line_length, FILE *outfp, Boolean gapped)

{
	return PrintKAParametersExtra(Lambda, K, H, 0.0, line_length, outfp, gapped);
}

Boolean LIBCALL
PrintKAParametersExtra(Nlm_FloatHi Lambda, Nlm_FloatHi K, Nlm_FloatHi H, Nlm_FloatHi C, Int4 line_length, FILE *outfp, Boolean gapped)

{
	Char buffer[BUFFER_LENGTH];

	if (outfp == NULL)
		return FALSE;

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	if (gapped)
	{
		ff_AddString("Gapped");
		NewContLine();
	}
	
	if (C == 0.0)
		ff_AddString("Lambda     K      H");
	else
		ff_AddString("Lambda     K      H      C");
	NewContLine();
	sprintf(buffer, "%#8.3g ", Lambda);
	ff_AddString(buffer);
	sprintf(buffer, "%#8.3g ", K);
	ff_AddString(buffer);
	sprintf(buffer, "%#8.3g ", H);
	ff_AddString(buffer);
	if (C != 0.0)
	{
		sprintf(buffer, "%#8.3g ", C);
		ff_AddString(buffer);
	}
	NewContLine();
	ff_EndPrint();

	return TRUE;

}


TxDfDbInfoPtr LIBCALL 
TxDfDbInfoNew (TxDfDbInfoPtr old)

{
	TxDfDbInfoPtr dbinfo;
	dbinfo = MemNew(sizeof(TxDfDbInfo));
	if (old)
		old->next = dbinfo;
	return dbinfo;
}

TxDfDbInfoPtr LIBCALL 
TxDfDbInfoDestruct (TxDfDbInfoPtr dbinfo)

{
	TxDfDbInfoPtr next;

	if (dbinfo == NULL)
		return NULL;

	while (dbinfo)
	{
		dbinfo->name = MemFree(dbinfo->name);
		dbinfo->definition = MemFree(dbinfo->definition);
		dbinfo->date = MemFree(dbinfo->date);
		next = dbinfo->next;
		dbinfo = MemFree(dbinfo);
		dbinfo = next;
	}

	return dbinfo;
}

Boolean LIBCALL
PrintDbReport(TxDfDbInfoPtr dbinfo, Int4 line_length, FILE *outfp)

{

	if (dbinfo == NULL || outfp == NULL)
		return FALSE;

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(2, 2, line_length, NULL);

	if (dbinfo->subset == FALSE)
	{
		ff_AddString("Database: ");
		ff_AddString(dbinfo->definition);
		NewContLine();
		ff_AddString("  Posted date:  ");
		ff_AddString(dbinfo->date);
		NewContLine();
		ff_AddString("Number of letters in database: "); 
		ff_AddString(Nlm_Int8tostr((Int8) dbinfo->total_length, 1));
		NewContLine();
		ff_AddString("Number of sequences in database:  ");
		ff_AddString(Ltostr((long) dbinfo->number_seqs, 1));
		NewContLine();
	}
	else
	{
		ff_AddString("Subset of the database(s) listed below");
		NewContLine();
		ff_AddString("   Number of letters searched: "); 
		ff_AddString(Nlm_Int8tostr((Int8) dbinfo->total_length, 1));
		NewContLine();
		ff_AddString("   Number of sequences searched:  ");
		ff_AddString(Ltostr((long) dbinfo->number_seqs, 1));
		NewContLine();
	}
	ff_EndPrint();

	return TRUE;
}

/*
	Prints an acknowledgement of the Blast Query, in the standard
	BLAST format.
*/


Boolean LIBCALL
AcknowledgeBlastQuery(BioseqPtr bsp, Int4 line_length, FILE *outfp, Boolean believe_query, Boolean html)

{
	Char buffer[BUFFER_LENGTH];

	if (bsp == NULL || outfp == NULL)
		return FALSE;
	
	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	if (html)
		ff_AddString("<b>Query=</b> ");
	else
		ff_AddString("Query= ");
	if (bsp->id && (bsp->id->choice != SEQID_LOCAL || believe_query))
	{
		SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
		if (StringNCmp(buffer, "lcl|", 4) == 0)
			ff_AddString(buffer+4);
		else
			ff_AddString(buffer);
		ff_AddChar(' ');
	}
	ff_AddString(BioseqGetTitle(bsp));
	NewContLine();
	TabToColumn(10);
	ff_AddChar('(');
	ff_AddString(Ltostr((long) BioseqGetLen(bsp), 1));
	ff_AddString(" letters)");
	NewContLine();
        ff_EndPrint();

        return TRUE;
}

/*
	return the version of BLAST as a char. string.
*/
CharPtr LIBCALL
BlastGetReleaseDate (void)

{
	return BLAST_RELEASE_DATE;
}


/*
	return the version of BLAST as a char. string.
*/
CharPtr LIBCALL
BlastGetVersionNumber (void)

{
	return BLAST_ENGINE_VERSION;
}

Boolean BlastPrintVersionInfo (const char* program, Boolean html, FILE *outfp)

{
	return BlastPrintVersionInfoEx(program, html, BlastGetVersionNumber(), BlastGetReleaseDate(), outfp);
}

Boolean BlastPrintVersionInfoEx (const char* program, Boolean html, CharPtr version, CharPtr date, FILE *outfp)

{
	CharPtr ret_buffer;


	if (outfp == NULL)
		return FALSE;

	ret_buffer = StringSave(program);
	Nlm_StrUpper(ret_buffer);
	if (html)
		fprintf(outfp, "<b>%s %s [%s]</b>\n", ret_buffer, version, date);
	else
		fprintf(outfp, "%s %s [%s]\n", ret_buffer, version, date);
	ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/* 
	Returns a reference for the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/

CharPtr LIBCALL
BlastGetReference(Boolean html)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
		add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9254694&dopt=Citation\">Reference</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
		add_string_to_bufferEx("Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch&auml;ffer, ", &ret_buffer, &ret_buffer_length, TRUE);
	} else
		add_string_to_bufferEx("Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"Gapped BLAST and PSI-BLAST: a new generation of protein database search", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("programs\",  Nucleic Acids Res. 25:3389-3402.", &ret_buffer, &ret_buffer_length, TRUE);
	
	return ret_buffer;
}

Boolean LIBCALL
MegaBlastPrintReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	if (outfp == NULL)
		return FALSE;
	
	if (html) {
           add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=10890397&dopt=Citation\">Reference</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
           add_string_to_bufferEx("Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000),", &ret_buffer, &ret_buffer_length, TRUE);
	} else
           add_string_to_bufferEx("Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"A greedy algorithm for aligning DNA sequences\", ", 
                               &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("J Comput Biol 2000; 7(1-2):203-14.", 
                               &ret_buffer, &ret_buffer_length, TRUE);
	
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);
	return TRUE;
}

Boolean LIBCALL
BlastPrintReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;
	
        ret_buffer = BlastGetReference(html);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/* 
	Returns a reference for the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/


/* 
	Returns a reference for composition-based statistics to use
        in the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/

CharPtr LIBCALL
CBStatisticsGetReference(Boolean html, Boolean firstRound, Boolean moreRounds)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
        if (firstRound) {
	  if (html) {
	    add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=11452024&dopt=Citation\">Reference for composition-based statistics</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
	    add_string_to_bufferEx("Sch&auml;ffer, Alejandro A., L. Aravind, Thomas L. Madden, ", &ret_buffer, &ret_buffer_length, TRUE);
	} else
	  add_string_to_bufferEx("Reference for composition-based statistics:", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Schaffer, Alejandro A., L. Aravind, Thomas L. Madden,", &ret_buffer, &ret_buffer_length, TRUE);
	}
	else {
	  if (html) {
	    add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=11452024&dopt=Citation\">Reference for composition-based statistics </a></b>", &ret_buffer, &ret_buffer_length, TRUE);
	    add_string_to_bufferEx("starting in round 2:", &ret_buffer, &ret_buffer_length, TRUE);

	    add_string_to_bufferEx("Sch&auml;ffer, Alejandro A., L. Aravind, Thomas L. Madden, ", &ret_buffer, &ret_buffer_length, TRUE);
	  } else {
	    add_string_to_bufferEx("Reference for composition-based statistics starting in round 2:", &ret_buffer, &ret_buffer_length, TRUE);
	    add_string_to_bufferEx("Schaffer, Alejandro A., L. Aravind, Thomas L. Madden,", &ret_buffer, &ret_buffer_length, TRUE);
	  }
	}
	add_string_to_bufferEx("Sergei Shavirin, John L. Spouge, Yuri I. Wolf,  ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Eugene V. Koonin, and Stephen F. Altschul (2001), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"Improving the accuracy of PSI-BLAST protein database searches with ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("composition-based statistics and other refinements\",  Nucleic Acids Res. 29:2994-3005.", &ret_buffer, &ret_buffer_length, TRUE);
	return ret_buffer;
}

/*print the reference for composition-based statistics when they are used*/
Boolean LIBCALL
CBStatisticsPrintReference(Boolean html, Int4 line_length, 
			   Boolean firstRound, Boolean moreRounds, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;

	if (!(firstRound || moreRounds))
	  return FALSE;
	
        ret_buffer = CBStatisticsGetReference(html,firstRound, moreRounds);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/* 
	Returns a reference for the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/

CharPtr LIBCALL
CAdjustmentGetReference(Boolean html)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	if (html) {
	  add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=16218944&dopt=Citation\">Reference for compositional score matrix adjustment</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Altschul, Stephen F., John C. Wootton, E. Michael Gertz, Richa Agarwala,", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Aleksandr Morgulis, Alejandro A. Sch&auml;ffer, and Yi-Kuo Yu (2005) \"Protein database", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("searches using compositionally adjusted substitution matrices\", FEBS J. 272:5101-5109.", &ret_buffer, &ret_buffer_length, TRUE);	
	}
	else {
	  add_string_to_bufferEx("Reference for compositional score matrix adjustment: Altschul, Stephen F., ", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("John C. Wootton, E. Michael Gertz, Richa Agarwala, Aleksandr Morgulis,", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("Alejandro A. Schaffer, and Yi-Kuo Yu (2005) \"Protein database searches", &ret_buffer, &ret_buffer_length, TRUE);
	  add_string_to_bufferEx("using compositionally adjusted substitution matrices\", FEBS J. 272:5101-5109.", &ret_buffer, &ret_buffer_length, TRUE);	
	}
	return ret_buffer;
}

/*print the reference for composition-based statistics when they are used*/
Boolean LIBCALL
CAdjustmentPrintReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;

        ret_buffer = CAdjustmentGetReference(html);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/* 
	Returns a reference for the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/



CharPtr LIBCALL
BlastGetPhiReference(Boolean html)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
		add_string_to_bufferEx("<b><a http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9705509&dopt=Citation\">Reference</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
		add_string_to_bufferEx("Zhang, Zheng, Alejandro A. Sch&auml;ffer, Webb Miller, Thomas L. Madden, ", &ret_buffer, &ret_buffer_length, TRUE);
	} else
		add_string_to_bufferEx("Reference: Zhang, Zheng, Alejandro A. Schaffer, Webb Miller, Thomas L. Madden, ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("David J. Lipman, Eugene V. Koonin, and Stephen F. Altschul (1998), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"Protein sequence similarity searches using patterns as seeds\", ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Nucleic Acids Res. 26:3986-3990.", &ret_buffer, &ret_buffer_length, TRUE);
	
	return ret_buffer;
}

Boolean LIBCALL
BlastPrintPhiReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;
	
        ret_buffer = BlastGetPhiReference(html);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/*
        Counts the number of SeqAligns present.
*/

static Int4
GetSeqAlignCount(SeqAlignPtr sap)

{
        Int4 count = 0;
        SeqIdPtr last_id=NULL, id;

        while (sap)
        {
                id = TxGetSubjectIdFromSeqAlign(sap);
                if (last_id)
                {
                        if(SeqIdComp(id, last_id) != SIC_YES)
                                count++;
                }
                else
                {
                        count = 1;
                }
                last_id = id;
                sap = sap->next;
        }

        return count;

}

/*
        Duplicates a SeqAlignPtr, up to the number of unique
        records specified.
*/

static SeqAlignPtr
GetPrivateSeqAlign(SeqAlignPtr sap, Int4 number, Int4Ptr number_returned)

{
        Int4 count=0;
        SeqIdPtr last_id=NULL, id;
        SeqAlignPtr new_head=NULL, var;

        last_id = TxGetSubjectIdFromSeqAlign(sap);

        while (count<number && sap)
        {
                count++;
                while (sap)
                {
                        id = TxGetSubjectIdFromSeqAlign(sap);
                        if(SeqIdComp(id, last_id) != SIC_YES)
                        {
                                last_id = id;
                                break;
                        }
                        if (new_head == NULL)
                        {
                                new_head = AsnIoMemCopy(sap, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);
                                var = new_head;
                        }
                        else
                        {
                                var->next = AsnIoMemCopy(sap, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);
                                var = var->next;
                        }
                        last_id = id;
                        sap = sap->next;
                }
        }

        *number_returned = count;

        return new_head;
}

/*
        Duplicate a SeqAlignPtr, keeping on the number of unique db
        hits specified.
*/

BlastPruneSapStructPtr LIBCALL
BlastPruneHitsFromSeqAlign(SeqAlignPtr sap, Int4 number, BlastPruneSapStructPtr prune)

{
        if (prune == NULL)
        {
                prune = MemNew(sizeof(BlastPruneSapStruct));
        }
        else
        {
                if (prune->number == number)
                        return prune;
                if (prune->allocated)
                        prune->sap = SeqAlignSetFree(prune->sap);
                prune->sap = NULL;
                prune->allocated = FALSE;
                prune->original_number = 0;
                prune->number = 0;
        }

        prune->original_number = GetSeqAlignCount(sap);

        if (prune->original_number < number)
        {
                prune->number = prune->original_number;
                prune->sap = sap;
                prune->allocated = FALSE;
        }
        else
        {
                prune->sap = GetPrivateSeqAlign(sap, number, &(prune->number));
                prune->allocated = TRUE;
        }

        return prune;
}

BlastPruneSapStructPtr LIBCALL
BlastPruneSapStructDestruct(BlastPruneSapStructPtr prune)

{
        if (prune == NULL)
                return NULL;

        if (prune->allocated)
        {
                prune->sap = SeqAlignSetFree(prune->sap);
        }
        prune = MemFree(prune);

        return prune;
}


void PrintTabularOutputHeader(CharPtr blast_database, BioseqPtr query_bsp,
                              SeqLocPtr query_slp, CharPtr blast_program,
                              Int4 iteration, Boolean believe_query,
                              FILE *outfp)
{
   Char buffer[BUFFER_LENGTH+1];
   Boolean unlock_bioseq = FALSE;

   asn2ff_set_output(outfp, NULL);
   
   ff_StartPrint(0, 0, BUFFER_LENGTH, NULL);

   if (blast_program) {
      CharPtr program = StringSave(blast_program);
      Nlm_StrUpper(program);
      sprintf(buffer, "# %s %s [%s]", program, BlastGetVersionNumber(),
              BlastGetReleaseDate());
      MemFree(program);
      ff_AddString(buffer);
      NewContLine();
   }

   if (iteration > 0) {
      ff_AddString("# Iteration: ");
      ff_AddString(Ltostr((long) iteration, 1));
      NewContLine();
   }

   if (query_bsp || query_slp) {
      CharPtr title;
      const CharPtr str = "# Query: ";
      Int4 string_length = StrLen(str);

      ff_AddString(str);

      if (!query_bsp) {
         Int4 num_queries = ValNodeLen(query_slp);
         if (num_queries > 1) {
            /* Multiple queries: just print the number, without deflines. */
            sprintf(buffer, "%ld sequences", (long)num_queries);
            ff_AddString(buffer);
         } else {
            query_bsp = BioseqLockById(SeqLocId(query_slp));
            unlock_bioseq = TRUE;
         }
      }
      if (query_bsp) {
         if (query_bsp->id && believe_query) {
            SeqIdWrite(query_bsp->id, buffer, PRINTID_FASTA_LONG, 
                       BUFFER_LENGTH);
            if (StringNCmp(buffer, "lcl|", 4) == 0) {
               ff_AddString(buffer+4);
            } else {
               ff_AddString(buffer);
            }
            string_length += StrLen(buffer);
            ff_AddChar(' ');
            string_length++; /* to account for the space above. */
         }

         if ((title = BioseqGetTitle(query_bsp)) != NULL) { 
            /* We do this to keep the entire title on one line 
               (of length BUFFER_LENGTH). */
            StrNCpy(buffer, title, BUFFER_LENGTH - string_length);
            buffer[BUFFER_LENGTH - string_length] = NULLB;
            ff_AddString(buffer);
         }

         if (unlock_bioseq)
            BioseqUnlock(query_bsp);
      }
      NewContLine();
   }
   if (blast_database) {
      ff_AddString("# Database: ");
      ff_AddString(blast_database);
      NewContLine();
   }
   if (getenv("PRINT_SEQUENCES")) {
         ff_AddString("# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score, query seq., subject seq.");
   } else {
         ff_AddString("# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score");
   }

   ff_EndPrint();
}

static Int4
BlastBioseqGetNumIdentical(BioseqPtr q_bsp, BioseqPtr s_bsp, Int4 q_start,
                     Int4 s_start, Int4 length,
                     Uint1 q_strand, Uint1 s_strand)
{
   SeqLocPtr q_slp, s_slp;
   SeqPortPtr q_spp, s_spp;
   Int4 i, ident = 0;
   Uint1 q_res, s_res;

   if (!q_bsp || !s_bsp)
      return 0;

   q_slp = SeqLocIntNew(q_start, q_start+length-1, q_strand, q_bsp->id);
   s_slp = SeqLocIntNew(s_start, s_start+length-1, s_strand, s_bsp->id);
   if (ISA_na(q_bsp->mol))
      q_spp = SeqPortNewByLoc(q_slp, Seq_code_ncbi4na);
   else
      q_spp = SeqPortNewByLoc(q_slp, Seq_code_ncbistdaa);
   if (ISA_na(s_bsp->mol))
      s_spp = SeqPortNewByLoc(s_slp, Seq_code_ncbi4na);
   else
      s_spp = SeqPortNewByLoc(s_slp, Seq_code_ncbistdaa);

   for (i=0; i<length; i++) {
      while ((q_res = SeqPortGetResidue(q_spp)) != SEQPORT_EOF &&
             !IS_residue(q_res));
      while ((s_res = SeqPortGetResidue(s_spp)) != SEQPORT_EOF &&
             !IS_residue(s_res));
      if (q_res == SEQPORT_EOF || s_res == SEQPORT_EOF)
         break;
      else if (q_res == s_res)
         ident++;
   }

   SeqLocFree(q_slp);
   SeqLocFree(s_slp);
   SeqPortFree(q_spp);
   SeqPortFree(s_spp);

   return ident;
}
/* 
   Function to print results in tab-delimited format, given a SeqAlign list.
   q_shift and s_shift are the offsets in query and subject in case of a
   subsequence search 
*/
void BlastPrintTabulatedResults(SeqAlignPtr seqalign, BioseqPtr query_bsp,
                                SeqLocPtr query_slp, Int4 num_alignments, 
                                CharPtr blast_program, Boolean is_ungapped, 
                                Boolean believe_query, Int4 q_shift, 
                                Int4 s_shift, FILE *fp,
                                Boolean print_query_info)
{
   BlastPrintTabulatedResultsEx(seqalign, query_bsp, query_slp, num_alignments,
                                blast_program, is_ungapped, believe_query,
                                q_shift, s_shift, fp, NULL, print_query_info);
}

void BlastPrintTabulatedResultsEx(SeqAlignPtr seqalign, BioseqPtr query_bsp,
                                SeqLocPtr query_slp, Int4 num_alignments, 
                                CharPtr blast_program, Boolean is_ungapped, 
                                Boolean believe_query, Int4 q_shift, 
                                Int4 s_shift, FILE *fp, 
                                int *num_formatted, Boolean print_query_info)
{
   BlastPrintTabularResults(seqalign, query_bsp, query_slp, num_alignments,
      blast_program, is_ungapped, FALSE, believe_query,
      q_shift, s_shift, fp, num_formatted, print_query_info);
}

void BlastPrintTabularResults(SeqAlignPtr seqalign, BioseqPtr query_bsp,
        SeqLocPtr query_slp, Int4 num_alignments, CharPtr blast_program, 
        Boolean is_ungapped, Boolean is_ooframe, Boolean believe_query, 
        Int4 q_shift, Int4 s_shift, FILE *fp, int *num_formatted, 
        Boolean print_query_info)
{
   SeqAlignPtr sap, sap_tmp = NULL;
   FloatHi perc_ident, bit_score, evalue;
   Int4 numseg, num_gap_opens, num_mismatches, num_ident, score;
   Int4 number, align_length, index, i;
   Int4 q_start, q_end, s_start, s_end;
   Char bit_score_buff[10];
   CharPtr eval_buff;
   Boolean is_translated;
   SeqIdPtr query_id, old_query_id = NULL, subject_id, old_subject_id = NULL;
   BioseqPtr subject_bsp=NULL;
   Char query_buffer[BUFFER_LENGTH+1], subject_buffer[BUFFER_LENGTH+1];
   DenseSegPtr dsp;
   StdSegPtr ssp = NULL;
   DenseDiagPtr ddp = NULL;
   AlignSumPtr asp = NULL;
   CharPtr defline, title;
   SeqLocPtr slp;
   Int4 alignments_count;

   is_translated = (StringCmp(blast_program, "blastn") &&
                    StringCmp(blast_program, "blastp"));
   
   if (is_translated) {
      asp = MemNew(sizeof(AlignSum));
      asp->matrix = load_default_matrix();
      asp->is_aa = TRUE;
      asp->ooframe = is_ooframe;
   }

   if (is_ungapped)
      sap_tmp = SeqAlignNew();

   slp = query_slp;
   if (query_bsp)
      query_id = query_bsp->id;

   /* Evalue buffer is dynamically allocated to avoid compiler warnings 
      in calls to ScoreAndEvalueToBuffers. */
   eval_buff = Malloc(10);

   for (sap = seqalign; sap; sap = sap->next) {
      if (query_slp)
         query_id = TxGetQueryIdFromSeqAlign(sap);
      if (SeqIdComp(query_id, old_query_id) != SIC_YES) {
         if (old_query_id && num_formatted)
            (*num_formatted)++;
         alignments_count = num_alignments;
         /* New query: find the corresponding SeqLoc */
         while (slp && SeqIdComp(query_id, SeqLocId(slp)) != SIC_YES)
            slp = slp->next;
         if (slp != NULL) {
            query_id = old_query_id = SeqLocId(slp);
            /* Print new query information */
            if (print_query_info)
               PrintTabularOutputHeader(NULL, NULL, slp, NULL, 0, 
                                        believe_query, fp);
         } else if (query_bsp)
            old_query_id = query_bsp->id;
         defline = (CharPtr) Malloc(BUFFER_LENGTH+1);
         SeqIdWrite(query_id, defline, PRINTID_FASTA_LONG, BUFFER_LENGTH);
         if (StringNCmp(defline, "lcl|", 4))
            StringCpy(query_buffer, defline);
         else if (!believe_query) {
            if (slp) {
               BioseqUnlock(query_bsp);
               query_bsp = BioseqLockById(query_id);
            }
            if ((title = StringSave(BioseqGetTitle(query_bsp))) != NULL) {
               defline = MemFree(defline);
               defline = StringTokMT(title, " ", &title);
               StringNCpy_0(query_buffer, defline, BUFFER_LENGTH);
               defline = MemFree(defline);
            } else
               StringCpy(query_buffer, defline+4);
            defline = MemFree(defline);
         } else
            StringCpy(query_buffer, defline+4);
      } else
         query_id = old_query_id;      

      subject_id = TxGetSubjectIdFromSeqAlign(sap);

      if (SeqIdComp(subject_id, old_subject_id) != SIC_YES) {
         /* New subject sequence has been found in the seqalign list */
         if (--alignments_count < 0)
            continue;
         BioseqUnlock(subject_bsp);
         subject_bsp = BioseqLockById(subject_id);
      
         if (!subject_bsp || !subject_bsp->id)
            continue;
         if (subject_bsp->id->choice != SEQID_GENERAL ||
             StringCmp(((DbtagPtr)subject_id->data.ptrvalue)->db, "BL_ORD_ID")) {
            defline = (CharPtr) Malloc(BUFFER_LENGTH+1);
            SeqIdWrite(subject_bsp->id, defline, PRINTID_FASTA_LONG, BUFFER_LENGTH);
            if (StringNCmp(defline, "lcl|", 4))
               StringCpy(subject_buffer, defline);
            else
               StringCpy(subject_buffer, defline+4);
         } else {
            defline = StringSave(BioseqGetTitle(subject_bsp));
            defline = StringTokMT(defline, " \t", &title);
            StringCpy(subject_buffer, defline);
         }
         defline = MemFree(defline);
      }
      
      perc_ident = 0;
      align_length = 0;
      num_gap_opens = 0;
      num_mismatches = 0;

      GetScoreAndEvalue(sap, &score, &bit_score, &evalue, &number);

      /* Do not allow knocking off digit in evalue buffer, so parsers are 
         not confused. */
      ScoreAndEvalueToBuffers(bit_score, evalue, 
                              bit_score_buff, &eval_buff, 0);

      /* Loop on segments within this seqalign (in ungapped case) */
      while (TRUE) {
         if (sap->segtype == SAS_DENSEG) {
            Boolean get_num_ident = TRUE;
            dsp = (DenseSegPtr) sap->segs;
            numseg = dsp->numseg;
            /* Query Bioseq is needed for calculating number of identities.
               NB: even if number of identities is already filled in the 
               seqalign score list, that is not enough here, because we need to
               know number of identities in each segment in order to calculate
               number of mismatches correctly. */
            if (!query_bsp) {
               query_bsp = BioseqLockById(query_id);
            }

            for (i=0; i<numseg; i++) {
               align_length += dsp->lens[i];
               if (dsp->starts[2*i] != -1 && dsp->starts[2*i+1] != -1) {
                  if (get_num_ident) {
                     num_ident = BlastBioseqGetNumIdentical(query_bsp, subject_bsp, 
                                    dsp->starts[2*i], dsp->starts[2*i+1], 
                                    dsp->lens[i], dsp->strands[2*i], 
                                    dsp->strands[2*i+1]);
                     perc_ident += num_ident;
                     num_mismatches += dsp->lens[i] - num_ident;
                  }
               } else {
                  num_gap_opens++;
               }
            }
            perc_ident = perc_ident / align_length * 100;

            if (dsp->strands[0] != dsp->strands[1]) {
               q_start = dsp->starts[2*numseg-2] + 1;
               q_end = dsp->starts[0] + dsp->lens[0];
               s_end = dsp->starts[1] + 1;
               s_start = dsp->starts[2*numseg-1] + dsp->lens[numseg-1];
            } else {
               q_start = dsp->starts[0] + 1;
               q_end = dsp->starts[2*numseg-2] + dsp->lens[numseg-1];
               s_start = dsp->starts[1] + 1;
               s_end = dsp->starts[2*numseg-1] + dsp->lens[numseg-1];
            }
         } else if (sap->segtype == SAS_STD) {
            if (!ssp)
               ssp = (StdSegPtr) sap->segs;
            
            if (is_ungapped) {
               sap_tmp->segtype = SAS_STD;
               sap_tmp->segs = ssp;
               GetScoreAndEvalue(sap_tmp, &score, &bit_score, &evalue, &number);
               ScoreAndEvalueToBuffers(bit_score, evalue, 
                                       bit_score_buff, &eval_buff, 0);
               find_score_in_align(sap_tmp, 1, asp);
            } else
               find_score_in_align(sap, 1, asp);
            
            if (asp->m_frame < 0)
               q_start = SeqLocStop(ssp->loc) + 1;
            else
               q_start = SeqLocStart(ssp->loc) + 1;
            
            if (asp->t_frame < 0)
               s_start = SeqLocStop(ssp->loc->next) + 1;
            else
               s_start = SeqLocStart(ssp->loc->next) + 1;
            
            if (!is_ungapped) {
               for (index=1; ssp->next; index++)
                  ssp = ssp->next;
               num_gap_opens = index / 2;
            } else 
               num_gap_opens = 0;

            if (asp->m_frame < 0)
               q_end = SeqLocStart(ssp->loc) + 1;
            else
               q_end = SeqLocStop(ssp->loc) + 1;
            
            if (asp->t_frame < 0)
               s_end = SeqLocStart(ssp->loc->next) + 1;
            else
               s_end = SeqLocStop(ssp->loc->next) + 1;
            
            align_length = asp->totlen;
            num_mismatches = asp->totlen - asp->gaps - asp->identical;
            perc_ident = ((FloatHi) 100*asp->identical)/ (asp->totlen);
         } else if (sap->segtype == SAS_DENDIAG) {
            if (!ddp)
               ddp = (DenseDiagPtr) sap->segs;
            sap_tmp->segtype = SAS_DENDIAG;
            sap_tmp->segs = ddp;
            GetScoreAndEvalue(sap_tmp, &score, &bit_score, &evalue, &number);
            ScoreAndEvalueToBuffers(bit_score, evalue, 
                                    bit_score_buff, &eval_buff, 0);

            align_length = ddp->len;
            if (ddp->strands[0] == Seq_strand_minus) {
               q_start = ddp->starts[0] + align_length;
               q_end = ddp->starts[0] + 1;
            } else {
               q_start = ddp->starts[0] + 1;
               q_end = ddp->starts[0] + align_length;
            }

            if (ddp->strands[1] == Seq_strand_minus) {
               s_start = ddp->starts[1] + align_length;
               s_end = ddp->starts[1] + 1;
            } else {
               s_start = ddp->starts[1] + 1;
               s_end = ddp->starts[1] + align_length;
            }
            num_gap_opens = 0;
            /* Query Bioseq is needed for calculating number of identities.
               NB: even if number of identities is already filled in the 
               seqalign score list, that is not enough here, because we need to
               know number of identities in each segment in order to calculate
               number of mismatches correctly. */
            if (!query_bsp) {
               query_bsp = BioseqLockById(query_id);
            }

            num_ident = BlastBioseqGetNumIdentical(query_bsp, subject_bsp, 
                           ddp->starts[0], ddp->starts[1], align_length, 
                           ddp->strands[0], ddp->strands[1]);
            num_mismatches = align_length - num_ident;
            perc_ident = ((FloatHi)num_ident) / align_length * 100;
         }
         if (!is_translated) {
            /* Adjust coordinates if query and/or subject is a subsequence */
            q_start += q_shift;
            q_end += q_shift;
            s_start += s_shift;
            s_end += s_shift;
         }
         
         if (perc_ident >= 99.995 && perc_ident < 100.00)
            perc_ident = 99.99;
         
         fprintf(fp, 
                 "%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
                 query_buffer, subject_buffer, perc_ident, align_length, 
                 num_mismatches, num_gap_opens, q_start, 
                 q_end, s_start, s_end, eval_buff, bit_score_buff);
         old_subject_id = subject_id;
         if (sap->segtype == SAS_DENSEG)
            break;
         else if (sap->segtype == SAS_DENDIAG) {
            if ((ddp = ddp->next) == NULL)
               break;
         } else if (sap->segtype == SAS_STD) {
            if ((ssp = ssp->next) == NULL)
               break;
         }
      }
   }

   eval_buff = MemFree(eval_buff);

   if (is_ungapped)
      sap_tmp = MemFree(sap_tmp);

   if (is_translated) {
      free_default_matrix(asp->matrix);
      MemFree(asp);
   }

   BioseqUnlock(subject_bsp);
   if (query_slp)
      BioseqUnlock(query_bsp);
}



/* Mutex for assignment of db seqs to search. */
TNlmMutex err_message_mutex=NULL;

#define BLAST_ERROR_BULEN 50
/*
	The following functions fill a the Error user string with
	text to identify BLAST and the entry being worked on.
	The SeqIdPtr is used to make a FASTA id, which is appended
	to string.

	A Uint1 is returned, which allows Nlm_ErrUserDelete to delete
	this error string when it's done.
*/

Uint1
BlastSetUserErrorString(CharPtr string, SeqIdPtr sip, Boolean use_id)

{
	BioseqPtr bsp;
	Char buffer[2*BLAST_ERROR_BULEN+1], textid[BLAST_ERROR_BULEN+1];
	CharPtr buf_start, ptr, title;
	Int2 length=0, index;
	Uint1 retval=0;

	buffer[0] = NULLB;
	ptr = buf_start = &buffer[0];

	if (string)
		StringNCpy_0(ptr, string, BLAST_ERROR_BULEN);

	if (sip != NULL)
	{
	    bsp = BioseqLockById(sip);
	    if(bsp)
	    {
		if (use_id)
			sip = bsp->id;
		else
			title = BioseqGetTitle(bsp);
	    }

	    if (string)
	    {
	    	length = StringLen(string);
	    	if (length > BLAST_ERROR_BULEN)
			length = BLAST_ERROR_BULEN;
	    }

	    ptr += length;

	    if (use_id)
	    {
    	    	SeqIdWrite(sip, textid, PRINTID_FASTA_LONG, BLAST_ERROR_BULEN-1);
	    	StringNCpy_0(ptr, textid, BLAST_ERROR_BULEN-1);
	    }
	    else if (title)
	    {
		for (index=0; index<BLAST_ERROR_BULEN-1; index++)
		{
			if (title[index] == NULLB || title[index] == ' ')
			{
				break;
			}
			*ptr = title[index];
			ptr++;
		}
		*ptr = NULLB;
	    }
	    BioseqUnlock(bsp);
	    StringCpy(ptr+StringLen(ptr), ":");
	}
	NlmMutexLockEx(&err_message_mutex);
	retval = Nlm_ErrUserInstall (buf_start, 0);
	NlmMutexUnlock(err_message_mutex);

	return retval;
}

void
BlastDeleteUserErrorString(Uint1 err_id)

{
	NlmMutexLockEx(&err_message_mutex);
	Nlm_ErrUserDelete(err_id);
	NlmMutexUnlock(err_message_mutex);
	return;
}

