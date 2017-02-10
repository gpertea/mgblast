/*------------------------------------------------------------------------------
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
* ------------------------------------------------------------------------------
*
* File Name:  nrpdbsrv.c
*
* File Description: Non-redundant-PDB-Set WWW-server 
*
* Created by: Yo Matsuo
* Version Creation Date: October 20, 1998
------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Header Files Included                                                      */
/*----------------------------------------------------------------------------*/
#include <ncbi.h>
#include <ncbiwww.h>
#include <ncbistr.h>
#include <miscutils.h>
#include <ncbimisc.h>
#include <sys/resource.h>

/*----------------------------------------------------------------------------*/
/* Definitions                                                                */
/*----------------------------------------------------------------------------*/
/* #define NRTABLE	"/net/ncbi/ftp/mmdb/nrtable/nrpdb.latest" */
#define NRTABLE	"nrpdb.latest"
#define NRPDBHOMEPAGE "nrpdb.html"
#define NRPDBSRVPATH "nrpdbsrv.cgi"
#define ENTREZSTRUC     "http://www.ncbi.nlm.nih.gov/Structure/mmdb/mmdbsrv?uid="
#define THIS_PROGRAM_IS "NRPDBSRV"
#define MAXNCHAINS 75000
#define NNRLEVEL 4
#define CPUTIME_MAX 120
#undef	NR_DEBUG

/*----------------------------------------------------------------------------*/
/* global variables                                                           */
/*----------------------------------------------------------------------------*/
/* data items from the NR PDB table (see below for detais) */
static Int4		nchains; /* No. of chains */
static Char		pdbcode[MAXNCHAINS][8];
static Char		chainid[MAXNCHAINS];
static Int4		mmdbid[MAXNCHAINS];
static Int4		groupid[NNRLEVEL][MAXNCHAINS];
static Int4		rankingroup[NNRLEVEL][MAXNCHAINS];
static Int4		repornot[NNRLEVEL][MAXNCHAINS];
static FloatHi 	punk[MAXNCHAINS];
static FloatHi 	picr[MAXNCHAINS];
static FloatHi 	pmr[MAXNCHAINS];
static FloatHi 	pics[MAXNCHAINS];
static FloatHi 	resol[MAXNCHAINS];
static Int4 	nch[MAXNCHAINS];
static Int4 	nht[MAXNCHAINS];
static Int4 	nhtt[MAXNCHAINS];
static Int4 	nres[MAXNCHAINS];
static Char 	expmet[MAXNCHAINS];
static Char 	accornot[MAXNCHAINS];

/* no. of members in each group */
static Int4		nmem[NNRLEVEL][MAXNCHAINS];

/* work space in listing representatives */
static Int4		chainbuffer[MAXNCHAINS];

/* describe each non-redundancy level */
static Char		chNRlevel[NNRLEVEL][40] =
				{"10e-7", "10e-40", "10e-80", "non-identical"};

/* describe each attribute */
/* attribute - description in the table */
/**/
static Char	chPdbcode1[] = "PDB";
static Char chPdbcode2[] = "PDB code";
/**/
static Char	chChainid1[] = "Ch";
static Char chChainid2[] = "Chain ID";
/**/
static Char chRankingroup1[] = "Rank";
static Char chRankingroup2[] = "Rank in this sequence-similar group";
/**/
static Char	chUnk1[] = "Unk";
static Char chUnk2[] = "Percentage of Unknown residues";
/**/
static Char	chIcr1[] = "Icr";
static Char chIcr2[] = "Percentage of Incomplete residues";
/**/
static Char	chMr1[] = "Mr";
static Char chMr2[] = "Percentage of Missing residues";
/**/
static Char	chIcs1[] = "Ics";
static Char chIcs2[] = "Percentage of Residues with Incomplete side-chain";
/**/
static Char chResol1[] = "Resol";
static Char chResol2[] = "Crystallographic Resolution (0.0 for NMR entry)";
/**/
static Char chNch1[] = "Nch";
static Char chNch2[] = "No. of Chains (subunits) in the PDB entry";
/**/
static Char chNht1[] = "Nht";
static Char chNht2[] = "No. of Heterogens (except for water)";
/**/
static Char chNhtt1[] = "Nhtt";
static Char chNhtt2[] = "No. of Types of Heterogens (except for water)";
/**/
static Char chNres1[] = "Nres";
static Char chNres2[] = "No. of Residues";
/**/
static Char chExpmet1[] = "Exp";
static Char chExpmet2[] = "Experimental methods: X = X-ray, N = NMR, M = theoretical";
/**/
static Char chNmem1[] = "Nmem";
static Char chNmem2[] = "No. of sequence-similar members in its group";

/*---------------------------------------*/
/* print HTML Header                     */
/*---------------------------------------*/
static void beginHTML()
{
	fprintf(stdout, "Content-type: text/html\n\n");
	fprintf(stdout, "<HTML>\n");
	fprintf(stdout, "<HEAD>\n");
	fprintf(stdout, "<TITLE>Non-redundant PDB chain set</TITLE>\n");
	fprintf(stdout, "</HEAD>\n");
	fprintf(stdout, "<BODY bgcolor=#FFFFFF>\n");
}

/*---------------------------------------*/
/* print HTML Footer                     */
/*---------------------------------------*/
static void endHTML()
{
	fprintf(stdout, "<BODY/>\n");
	fprintf(stdout, "<HTML/>\n");
}

/*---------------------------------------*/
/* trivial routine when an erro occurred */
/*---------------------------------------*/
static void EarlyOut()
{
	fflush(stdout);
	exit(1);
}

/*--------------------------------------------------------*/
/* Read in the NR set table                               */
/*--------------------------------------------------------*/
static void ReadNrTable()
{
	int 		i, j, k; /* counters */
	FILE 		*fp; /* file pointer for NR PDB set */
	Char 		line[10000]; /* work space to read in string */
	Char 		cWork[10][8]; /* work space to read in string */
	Int4 		iWork[20]; /* work space for reading in integer */
	FloatHi 	rWork[20]; /* work space for reading in floating point number */

/* This is the format of the NR set table.
   A line starting with "#" is a comment line, and should be skipped.
#---------------------------------------------------------------------
# Non-redundant PDB chain set
#---------------------------------------------------------------------
#
# 1: PDB code
# 2: Chain ID
# 3: MMDB ID
#
# 4: Group ID (BLAST pvalue 10e-7)
# 5: Rank (BLAST pvalue 10e-7)
# 6: Representative (1) or not (0) (BLAST pvalue 10e-7)
#
# 7: Group ID (BLAST pvalue 10e-40)
# 8: Rank (BLAST pvalue 10e-40)
# 9: Representative (1) or not (0) (BLAST pvalue 10e-40)
#
# A: Group ID (BLAST pvalue 10e-80)
# B: Rank (BLAST pvalue 10e-80)
# C: Representative (1) or not (0) (BLAST pvalue 10e-40)
#
# D: Group ID (Non-identical sequences)
# E: Rank (Non-identical sequence)
# F: Representative (1) or not (0) (Non-identical sequence)
#
# G: Percentage of Unknown residues
# H: Percentage of Incomplete residues
# I: Percentage of Missing residues
# J: Percentage of Incomplete side-chain residues
# K: Resolution (0.0 if NMR)
# L: No. of chains (subunits) in the PDB entry
# M: No. of heterogens in the PDB entry
# N: No. of different heterogen types in the PDB entry
# O: No. of residues in the chain
# P: Method of coordinate determination (X: X-ray, N: NMR, M: theoretical)
# Q: Acceptable (a) in structural quality or not (n)
#
#---------------------------------------------------------------------------------------------------------------------------
# 1  2       3     4    5 6     7    8 9     A    B C     D    E F    G      H      I      J     K     L   M   N     O P Q
#---------------------------------------------------------------------------------------------------------------------------
1L58      1460     1    1 0     1    1 0     1    1 0  2206    1 1   0.00   0.00   0.00   0.00  1.65   1   1   1   164 X a
237L      7320     1    2 0     1    2 0     1    2 0  3773    1 1   0.00   0.00   0.00   0.00  1.70   1   3   2   164 X a
*/
	if ((fp = FileOpen(NRTABLE, "r")) == NULL)
	{
		fprintf(stdout, "Content-type: text/html\n\n");
		fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
		fprintf(stdout, "<H3>Failed to read data.</H3>\n");
		EarlyOut();
	}

	nchains = 0;
	while (NULL != fgets(line, 500, fp))
	{
		/* skip the comment line */
		if (line[0] == '#') continue;

		/*---------------------------------------------------------
			PDB code (columns 1 - 4)
			chain ID (column 6)
			domain ID or number (column 8) (if to be printed)
		---------------------------------------------------------*/
		/* PDB code */
		for (i = 0; i < 4; ++i) pdbcode[nchains][i] = line[i];
		pdbcode[nchains][4] = '\0';

		/* chain ID */
		chainid[nchains] = line[5];

		for (i = 0; ; ++i)
		{
			line[i] = line[i + 9];
			if (line[i] == '\0') break;
		}

		sscanf(line, "%d%d%d%d%d%d%d%d%d%d%d%d%d%lf%lf%lf%lf%lf%d%d%d%d%s%s",
			&(iWork[0]), &(iWork[1]), &(iWork[2]), &(iWork[3]), &(iWork[4]), &(iWork[5]),
			&(iWork[6]), &(iWork[7]), &(iWork[8]), &(iWork[9]), &(iWork[10]), &(iWork[11]),
			&(iWork[12]),
			&(rWork[0]), &(rWork[1]), &(rWork[2]), &(rWork[3]), &(rWork[4]),
			&(iWork[13]), &(iWork[14]), &(iWork[15]), &(iWork[16]),
			cWork[0], cWork[1]
			);

		/* MMDB ID */
		mmdbid[nchains] = iWork[0];

		/* group ID, rank in the group, and representative or not */
		/* p-value cutoff 10e-7 */
		groupid[0][nchains] = iWork[1];
		rankingroup[0][nchains] = iWork[2];
		repornot[0][nchains] = iWork[3];

		/* group ID, rank in the group, and representative or not */
		/* p-value cutoff 10e-40 */
		groupid[1][nchains] = iWork[4];
		rankingroup[1][nchains] = iWork[5];
		repornot[1][nchains] = iWork[6];

		/* group ID, rank in the group, and representative or not */
		/* p-value cutoff 10e-80 */
		groupid[2][nchains] = iWork[7];
		rankingroup[2][nchains] = iWork[8];
		repornot[2][nchains] = iWork[9];

		/* group ID, rank in the group, and representative or not */
		/* 100% identical */
		groupid[3][nchains] = iWork[10];
		rankingroup[3][nchains] = iWork[11];
		repornot[3][nchains] = iWork[12];

		/* Percentage of Unknown residues */
		punk[nchains] = rWork[0];

		/* Percentage of Incomplete residues */
		picr[nchains] = rWork[1];

		/* Percentage of Missing residues */
		pmr[nchains] = rWork[2];

		/* Percentage of Incomplete side-chain residues */
		pics[nchains] = rWork[3];

		/* Resolution (0.0 if NMR) */
		resol[nchains] = rWork[4];

		/* No. of chains (subunits) in the PDB entry */
		nch[nchains] = iWork[13];

		/* No. of heterogens in the PDB entry */
		nht[nchains] = iWork[14];

		/* No. of different heterogen types in the PDB entry */
		nhtt[nchains] = iWork[15];

		/* No. of residues in the chain */
		nres[nchains] = iWork[16];

		/* Method of coordinate determination (X: X-ray, N: NMR, M: theoretical) */
		expmet[nchains] = cWork[0][0];

		/* Acceptable (a) in structural quality or not (n) */
		accornot[nchains] = cWork[1][0];

		++nchains;
	}
}

/*----------------------*/
/* List representatives */
/*----------------------*/
static void ListRepresentatives(Int2 NRlevel, CharPtr tableoption, Int4 iPage, Int4 iPagesize)
{
	Int4 ichain;
	Int4 ic;

	Int4 nrep;
	Int4 irep;
	Int4 nPage;
	Int4 i1, i2;
	Int4 nextPage;

	/* collect representatives */
	nrep = 0;
	for (ichain = 0; ichain < nchains; ++ichain)
		if (repornot[NRlevel][ichain] == 1) chainbuffer[nrep ++] = ichain;

	if (0 == strcmp(tableoption, "long"))
	{
	/* long table */

	/* No. of pages */
	nPage = nrep / iPagesize;
	if (0 < nrep % iPagesize) ++nPage;
	if (iPage > nPage) iPage = nPage;

	/* range of representatives to be shown in this page */
	i1 = iPagesize * (iPage - 1) + 1;
	i2 = iPagesize * iPage;
	if (i2 > nrep) i2 = nrep;

	nextPage = iPage + 1;
	if (nextPage > nPage) nextPage = 1;

	fprintf(stdout, "<A name=\"tabletop\"></A>\n");
	fprintf(stdout, "<B>No. of Representatives:</B> %d<BR><BR>\n", nrep);
	fprintf(stdout, "<B>Page</B> %d / %d<BR>\n", iPage, nPage);
	fprintf(stdout, "<B>Representatives</B> %d - %d\n", i1, i2);
	fflush(stdout);
	fprintf(stdout, "<TABLE border=1 cellpadding=3>\n");
	fprintf(stdout, "<TR>\n");
	fprintf(stdout, "<TH align=left><A href=\"#iorder\"></A></TH>\n");
	fprintf(stdout, "<TH align=left><A href=\"#pdbcode\">%s</A></TH>\n", chPdbcode1);
	fprintf(stdout, "<TH align=left><A href=\"#chainid\">%s</A></TH>\n", chChainid1);
	fprintf(stdout, "<TH align=right><A href=\"#unk\">%s</A></TH>\n", chUnk1);
	fprintf(stdout, "<TH align=right><A href=\"#icr\">%s</A></TH>\n", chIcr1);
	fprintf(stdout, "<TH align=right><A href=\"#mr\">%s</A></TH>\n", chMr1);
	fprintf(stdout, "<TH align=right><A href=\"#ics\">%s</A></TH>\n", chIcs1);
	fprintf(stdout, "<TH align=right><A href=\"#resol\">%s</A></TH>\n", chResol1);
	fprintf(stdout, "<TH align=right><A href=\"#nch\">%s</A></TH>\n", chNch1);
	fprintf(stdout, "<TH align=right><A href=\"#nht\">%s</A></TH>\n", chNht1);
	fprintf(stdout, "<TH align=right><A href=\"#nhtt\">%s</A></TH>\n", chNhtt1);
	fprintf(stdout, "<TH align=right><A href=\"#nres\">%s</A></TH>\n", chNres1);
	fprintf(stdout, "<TH align=right><A href=\"#expmet\">%s</A></TH>\n", chExpmet1);
	fprintf(stdout, "<TH align=right><A href=\"#nmem\">%s</A></TH>\n", chNmem1);
	fprintf(stdout, "</TR>\n");
	for (irep = (i1 - 1); irep <= (i2 - 1); ++irep)
	{
		ichain = chainbuffer[irep];
		fprintf(stdout, "<TR>\n");
		fprintf(stdout, "<TD align=left>%d</TD>\n", irep + 1);
		fprintf(stdout, "<TD align=left><A href=\"%s%d\">%s</A></TD>\n",
			ENTREZSTRUC, mmdbid[ichain], pdbcode[ichain]);
		if (chainid[ichain] == ' ')
			fprintf(stdout, "<TD align=left>&nbsp;</TD>\n");
		else
			fprintf(stdout, "<TD align=left>%c</TD>\n", chainid[ichain]);
		fprintf(stdout, "<TD align=right>%6.1f</TD>\n", punk[ichain]);
		fprintf(stdout, "<TD align=right>%6.1f</TD>\n", picr[ichain]);
		fprintf(stdout, "<TD align=right>%6.1f</TD>\n", pmr[ichain]);
		fprintf(stdout, "<TD align=right>%6.1f</TD>\n", pics[ichain]);
		fprintf(stdout, "<TD align=right>%5.2f</TD>\n", resol[ichain]);
		fprintf(stdout, "<TD align=right>%3d</TD>\n", nch[ichain]);
		fprintf(stdout, "<TD align=right>%3d</TD>\n", nht[ichain]);
		fprintf(stdout, "<TD align=right>%3d</TD>\n", nhtt[ichain]);
		fprintf(stdout, "<TD align=right>%5d</TD>\n", nres[ichain]);
		fprintf(stdout, "<TD align=right>%c</TD>\n", expmet[ichain]);
		fprintf(stdout, "<TD align=right><A href=\"%s?nrlevel=%d&group=Submit&pdbcode=%s&chainid=%c\">%d</A></TD>\n",
			NRPDBSRVPATH, NRlevel, pdbcode[ichain], chainid[ichain], nmem[NRlevel][groupid[NRlevel][ichain]]);
		fprintf(stdout, "</TR>\n");
	}
	fprintf(stdout, "</TABLE>\n");
	fprintf(stdout, "<BR>\n");

	fprintf(stdout, "<FORM METHOD=\"post\" ACTION=\"%s\">\n", NRPDBSRVPATH);
	fprintf(stdout, "<INPUT TYPE=\"submit\" NAME=\"replist\" VALUE=\"Go to\">&nbsp;<B>page</B>\n");
	fprintf(stdout, "<INPUT TYPE=\"text\" NAME=\"page\" value=%d size=5> of %d\n", nextPage, nPage);
	fprintf(stdout, "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n");
	fprintf(stdout, "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n");
	fprintf(stdout, "<B>No. representatives per page</B>");
	fprintf(stdout, "<INPUT TYPE=\"text\" NAME=\"pagesize\" value=%d size=5>\n", iPagesize);
	fprintf(stdout, "<INPUT TYPE=\"hidden\" NAME=\"ltable\" VALUE=\"long\">\n");
	fprintf(stdout, "<INPUT TYPE=\"hidden\" NAME=\"nrlevel\" VALUE=%d>\n", NRlevel);
	fprintf(stdout, "</FORM>\n");

	fprintf(stdout, "<TABLE border=0>\n");
	fprintf(stdout, "<TR><TD><A name=\"pdbcode\"></A><B>%s:</B>\n", chPdbcode1);
	fprintf(stdout, "<TD>%s</TR>\n", chPdbcode2);
	fprintf(stdout, "<TR><TD><A name=\"chainid\"></A><B>%s:</B>\n", chChainid1);
	fprintf(stdout, "<TD>%s</TR>\n", chChainid2);
	fprintf(stdout, "<TR><TD><A name=\"unk\"></A><B>%s:</B>\n", chUnk1);
	fprintf(stdout, "<TD>%s</TR>\n", chUnk2);
	fprintf(stdout, "<TR><TD><A name=\"icr\"></A><B>%s:</B>\n", chIcr1);
	fprintf(stdout, "<TD>%s</TR>\n", chIcr2);
	fprintf(stdout, "<TR><TD><A name=\"mr\"></A><B>%s:</B>\n", chMr1);
	fprintf(stdout, "<TD>%s</TR>\n", chMr2);
	fprintf(stdout, "<TR><TD><A name=\"ics\"></A><B>%s:</B>\n", chIcs1);
	fprintf(stdout, "<TD>%s</TR>\n", chIcs2);
	fprintf(stdout, "<TR><TD><A name=\"resol\"></A><B>%s:</B>\n", chResol1);
	fprintf(stdout, "<TD>%s</TR>\n", chResol2);
	fprintf(stdout, "<TR><TD><A name=\"nch\"></A><B>%s:</B>\n", chNch1);
	fprintf(stdout, "<TD>%s</TR>\n", chNch2);
	fprintf(stdout, "<TR><TD><A name=\"nht\"></A><B>%s:</B>\n", chNht1);
	fprintf(stdout, "<TD>%s</TR>\n", chNht2);
	fprintf(stdout, "<TR><TD><A name=\"nhtt\"></A><B>%s:</B>\n", chNhtt1);
	fprintf(stdout, "<TD>%s</TR>\n", chNhtt2);
	fprintf(stdout, "<TR><TD><A name=\"nres\"></A><B>%s:</B>\n", chNres1);
	fprintf(stdout, "<TD>%s</TR>\n", chNres2);
	fprintf(stdout, "<TR><TD><A name=\"expmet\"></A><B>%s:</B>\n", chExpmet1);
	fprintf(stdout, "<TD>%s</TR>\n", chExpmet2);
	fprintf(stdout, "<TR><TD valign=top><A name=\"nmem\"></A><B>%s:</B>\n", chNmem1);
	fprintf(stdout, "<TD>%s<BR>", chNmem2);
	fprintf(stdout, "Click on this number, then you can see a list of those members.");
	fprintf(stdout, "</TR>\n");
	fprintf(stdout, "</TABLE>\n");
	fprintf(stdout, "<BR>\n");

	fprintf(stdout, "<A href=\"#tabletop\">Top of the Table</A><BR><BR>\n");
	fprintf(stdout, "<A href=\"%s\">NR set homepage</A><BR>\n", NRPDBHOMEPAGE);

	}
	else
	{
	/* short table */
	fprintf(stdout, "<B>No. of Representatives:</B> %d<BR><BR>\n", nrep);
	fflush(stdout);
	fprintf(stdout, "<TABLE border=1>\n");
	ic = 0;
	for (irep = 0; irep < nrep; ++irep)
	{
		ichain = chainbuffer[irep];
		if (ic == 0) fprintf(stdout, "<TR>\n");
/*
		fprintf(stdout, "<TD align=left><A href=\"%s%d\">%s%c</A>&nbsp\;</TD>\n",
			ENTREZSTRUC, mmdbid[ichain], pdbcode[ichain], chainid[ichain]);
*/
		fprintf(stdout, "<TD align=left>%s%c&nbsp;</TD>\n",
			pdbcode[ichain], chainid[ichain]);
		++ic;
		if (ic == 10)
		{
			fprintf(stdout, "</TR>\n");
			ic = 0;
		}
	}
	fprintf(stdout, "</TABLE>\n");
	}
}

/*----------------------------*/
/* Show members of a NR group */
/*----------------------------*/
static void ShowNRGroup(Int2 NRlevel, Int4 querygroupid, Int4 querychain)
{
	Int4 ichain;

	fprintf(stdout, "<A name=\"tabletop\"></A>\n");
	fprintf(stdout, "<TABLE border=1 cellpadding=3>\n");
	fprintf(stdout, "<TR>\n");
	fprintf(stdout, "<TH align=left><A href=\"#rank\">%s</A></TH>\n", chRankingroup1);
	fprintf(stdout, "<TH align=left><A href=\"#pdbcode\">%s</A></TH>\n", chPdbcode1);
	fprintf(stdout, "<TH align=left><A href=\"#chainid\">%s</A></TH>\n", chChainid1);
	fprintf(stdout, "<TH align=right><A href=\"#unk\">%s</A></TH>\n", chUnk1);
	fprintf(stdout, "<TH align=right><A href=\"#icr\">%s</A></TH>\n", chIcr1);
	fprintf(stdout, "<TH align=right><A href=\"#mr\">%s</A></TH>\n", chMr1);
	fprintf(stdout, "<TH align=right><A href=\"#ics\">%s</A></TH>\n", chIcs1);
	fprintf(stdout, "<TH align=right><A href=\"#resol\">%s</A></TH>\n", chResol1);
	fprintf(stdout, "<TH align=right><A href=\"#nch\">%s</A></TH>\n", chNch1);
	fprintf(stdout, "<TH align=right><A href=\"#nht\">%s</A></TH>\n", chNht1);
	fprintf(stdout, "<TH align=right><A href=\"#nhtt\">%s</A></TH>\n", chNhtt1);
	fprintf(stdout, "<TH align=right><A href=\"#nres\">%s</A></TH>\n", chNres1);
	fprintf(stdout, "<TH align=right><A href=\"#expmet\">%s</A></TH>\n", chExpmet1);
	fprintf(stdout, "</TR>\n");
	for (ichain = 0; ichain < nchains; ++ichain)
	{
		if (groupid[NRlevel][ichain] == querygroupid)
		{
			fprintf(stdout, "<TR>\n");
			fprintf(stdout, "<TD align=left>%d", rankingroup[NRlevel][ichain]);
			if (repornot[NRlevel][ichain] == 1) fprintf(stdout, "*");
			if (ichain == querychain) fprintf(stdout, "#");
			fprintf(stdout, "</TD>\n");
			fprintf(stdout, "<TD align=left><A href=\"%s%d\">%s</A></TD>\n",
				ENTREZSTRUC, mmdbid[ichain], pdbcode[ichain]);
			if (chainid[ichain] == ' ')
				fprintf(stdout, "<TD align=left>&nbsp;</TD>\n");
			else
				fprintf(stdout, "<TD align=left>%c</TD>\n", chainid[ichain]);
			fprintf(stdout, "<TD align=right>%6.1f</TD>\n", punk[ichain]);
			fprintf(stdout, "<TD align=right>%6.1f</TD>\n", picr[ichain]);
			fprintf(stdout, "<TD align=right>%6.1f</TD>\n", pmr[ichain]);
			fprintf(stdout, "<TD align=right>%6.1f</TD>\n", pics[ichain]);
			fprintf(stdout, "<TD align=right>%5.2f</TD>\n", resol[ichain]);
			fprintf(stdout, "<TD align=right>%3d</TD>\n", nch[ichain]);
			fprintf(stdout, "<TD align=right>%3d</TD>\n", nht[ichain]);
			fprintf(stdout, "<TD align=right>%3d</TD>\n", nhtt[ichain]);
			fprintf(stdout, "<TD align=right>%5d</TD>\n", nres[ichain]);
			fprintf(stdout, "<TD align=right>%c</TD>\n", expmet[ichain]);
			fprintf(stdout, "</TR>\n");
		}
	}
	fprintf(stdout, "</TABLE>\n");
	fprintf(stdout, "<BR>\n");
	fprintf(stdout, "<TABLE border=0>\n");
	fprintf(stdout, "<TR><TD valign=top><A name=\"rank\"></A><B>%s:</B>\n", chRankingroup1);
	fprintf(stdout, "<TD>%s<BR>", chRankingroup2);
	fprintf(stdout, "* if it is the representative of this group<BR>");
	fprintf(stdout, "# if it is the query</TR>\n");
	fprintf(stdout, "<TR><TD><A name=\"pdbcode\"></A><B>%s:</B>\n", chPdbcode1);
	fprintf(stdout, "<TD>%s</TR>\n", chPdbcode2);
	fprintf(stdout, "<TR><TD><A name=\"chainid\"></A><B>%s:</B>\n", chChainid1);
	fprintf(stdout, "<TD>%s</TR>\n", chChainid2);
	fprintf(stdout, "<TR><TD><A name=\"unk\"></A><B>%s:</B>\n", chUnk1);
	fprintf(stdout, "<TD>%s</TR>\n", chUnk2);
	fprintf(stdout, "<TR><TD><A name=\"icr\"></A><B>%s:</B>\n", chIcr1);
	fprintf(stdout, "<TD>%s</TR>\n", chIcr2);
	fprintf(stdout, "<TR><TD><A name=\"mr\"></A><B>%s:</B>\n", chMr1);
	fprintf(stdout, "<TD>%s</TR>\n", chMr2);
	fprintf(stdout, "<TR><TD><A name=\"ics\"></A><B>%s:</B>\n", chIcs1);
	fprintf(stdout, "<TD>%s</TR>\n", chIcs2);
	fprintf(stdout, "<TR><TD><A name=\"resol\"></A><B>%s:</B>\n", chResol1);
	fprintf(stdout, "<TD>%s</TR>\n", chResol2);
	fprintf(stdout, "<TR><TD><A name=\"nch\"></A><B>%s:</B>\n", chNch1);
	fprintf(stdout, "<TD>%s</TR>\n", chNch2);
	fprintf(stdout, "<TR><TD><A name=\"nht\"></A><B>%s:</B>\n", chNht1);
	fprintf(stdout, "<TD>%s</TR>\n", chNht2);
	fprintf(stdout, "<TR><TD><A name=\"nhtt\"></A><B>%s:</B>\n", chNhtt1);
	fprintf(stdout, "<TD>%s</TR>\n", chNhtt2);
	fprintf(stdout, "<TR><TD><A name=\"nres\"></A><B>%s:</B>\n", chNres1);
	fprintf(stdout, "<TD>%s</TR>\n", chNres2);
	fprintf(stdout, "<TR><TD><A name=\"expmet\"></A><B>%s:</B>\n", chExpmet1);
	fprintf(stdout, "<TD>%s</TR>\n", chExpmet2);
	fprintf(stdout, "</TABLE>\n");
	fprintf(stdout, "<BR>\n");
	fprintf(stdout, "<A href=\"#tabletop\">Top of the Table</A><BR><BR>\n");
	fprintf(stdout, "<A href=\"%s\">NR set homepage</A><BR>\n", NRPDBHOMEPAGE);
}

/******************************************************************************
*******************************************************************************
* NRPDBSRV Main
*******************************************************************************
*******************************************************************************/
Int2 Main ()
{
	int 		i; /* counters */
	FILE 		*fp; /* file pointer for NR PDB set */
	Char 		line[10000]; /* work space to read in string */

	CharPtr		Value; /* to extract from the info struct */
	struct		rlimit     rl; /* data structure for run limit */
	WWWInfoPtr	www_info; /* holds all info accessible */
	Int4		NumLabels  = 0;    /* number of WWW-parameters */
	Int4		indx; /* used in parsing WWW-input */

	Int4		ichain;
	Int2		NRlevel; /* NR level: 0=10e-7, 1=10e-40, 2=10e-80, 3=nonidentical */
	Char		querypdbcode[8];
	Char		querychainid;
	Int4		querygroupid;
	Char		tableoption[8];
	Int4		iPage;
	Int4		iPagesize;
	Int4		iNRlevel;

	/*----------------------------------------------------------------------------*/
	/* this sets up the unix time limit                                           */  
	/*----------------------------------------------------------------------------*/
	getrlimit(RLIMIT_CPU, &rl);
	rl.rlim_max = rl.rlim_cur = CPUTIME_MAX;
	setrlimit(RLIMIT_CPU, &rl);

	/*----------------------------------------------------------------------------*/
	/* check whether parameters have been submitted along with the posting        */
	/*----------------------------------------------------------------------------*/
	if (WWWGetArgs(&www_info) != WWWErrOk)
	{
		beginHTML();
		fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
		fprintf(stdout, "<H3>Failed to process posting - check your get/post syntax.</H3>\n");
		endHTML();
		EarlyOut();
	}
	if ((NumLabels = WWWGetNumEntries(www_info)) == 0)
	{
		beginHTML();
		fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
		fprintf(stdout, "<H3>No input - nothing to report.</H3>\n");
		endHTML();
		EarlyOut();
	}

	/*--------------------------------------------------------*/
	/* Read in the NR set table                               */
	/*--------------------------------------------------------*/
	ReadNrTable();

	/* count the number of members of each group */
	for (iNRlevel = 0; iNRlevel < NNRLEVEL; ++iNRlevel)
	{
		for (i = 0; i < MAXNCHAINS; ++i) nmem[iNRlevel][i] = 0;
		for (ichain = 0; ichain < nchains; ++ichain)
			nmem[iNRlevel][groupid[iNRlevel][ichain]] ++;
	}
	
	/*----------------------------------------------------------------------------*/
	/* extract from the query the level of non-redundancy that should be used     */
	/*----------------------------------------------------------------------------*/
	/* This is a mandatory option */
	if ((indx = WWWFindName(www_info, "nrlevel")) >= 0)
	{
		Value = WWWGetValueByIndex(www_info, indx);
		NRlevel = atoi(Value);
		if (NRlevel >= NNRLEVEL)
		{
			beginHTML();
			fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
			fprintf(stdout, "<H3>Unexpected non-redundancy level submitted.</H3>\n");
			endHTML();
			EarlyOut();
		}
	}
	else
	{
		beginHTML();
		fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
		fprintf(stdout, "<H3>NON-REDUNDANCY LEVEL UNKNOWN</H3>\n");
		endHTML();
		EarlyOut();
	}

	/*----------------------------------------------*/
	/*----------------------------------------------*/
	/* Action requested                             */
	/*----------------------------------------------*/
	/*----------------------------------------------*/
	/*----------------------*/
	/* List representatives */
	/*----------------------*/
	if ((indx = WWWFindName(www_info, "replist")) >= 0)
	{
		Value = WWWGetValueByIndex(www_info, indx);

		/* short table or long table? */
		if ((indx = WWWFindName(www_info, "ltable")) >= 0)
		{
			Value = WWWGetValueByIndex(www_info, indx);
			sprintf(tableoption, "%s", Value);
		}
		else
		{
			beginHTML();
			fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
			fprintf(stdout, "<H3>TABLE TYPE UNKNOWN</H3>\n");
			endHTML();
			EarlyOut();
		}

		/* Page to be shown */
		if ((indx = WWWFindName(www_info, "page")) >= 0)
		{
			Value = WWWGetValueByIndex(www_info, indx);
			iPage = atoi(Value);
		}
		else
		{
			beginHTML();
			fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
			fprintf(stdout, "<H3>PAGE NOT SPECIFIED</H3>\n");
			endHTML();
			EarlyOut();
		}

		/* Page size */
		if ((indx = WWWFindName(www_info, "pagesize")) >= 0)
		{
			Value = WWWGetValueByIndex(www_info, indx);
			iPagesize = atoi(Value);
		}
		else
		{
			beginHTML();
			fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
			fprintf(stdout, "<H3>PAGE SIZE NOT SPECIFIED</H3>\n");
			endHTML();
			EarlyOut();
		}

		/* send the page */
		beginHTML();
		fprintf(stdout, "<H2>List of Representatives</H2>\n");
		fprintf(stdout, "<B>Non-redundancy level: </B>\n");
		fprintf(stdout, "%s<BR>", chNRlevel[NRlevel]);
		fflush(stdout);
		ListRepresentatives(NRlevel, tableoption, iPage, iPagesize);
		endHTML();
	}
	/*----------------*/
	/* Show the group */
	/*----------------*/
	else if ((indx = WWWFindName(www_info, "group")) >= 0)
	{
		Value = WWWGetValueByIndex(www_info, indx);

		/*----------*/
		/* PDB code */
		/*----------*/
		querypdbcode[0] = '\0';
		if ((indx = WWWFindName(www_info, "pdbcode")) >= 0)
		{
			Value = WWWGetValueByIndex(www_info, indx);
			sprintf(querypdbcode, "%s", Value);
		}
		if (strlen(querypdbcode) != 4)
		{
			beginHTML();
			fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
			fprintf(stdout, "<H3>Invalid PDB code</H3>\n");
			endHTML();
			EarlyOut();
		}
		for (i = 0; i < 4; ++i)
			querypdbcode[i] = toupper(querypdbcode[i]);

		/*----------*/
		/* chain ID */
		/*----------*/
		querychainid = '\0';
		if ((indx = WWWFindName(www_info, "chainid")) >= 0)
		{
			Value = WWWGetValueByIndex(www_info, indx);
			querychainid = Value[0];
		}
		if (querychainid == '\0') querychainid = ' ';
		querychainid = toupper(querychainid);

/*
#ifdef NR_DEBUG
beginHTML();
fprintf(stdout, "<B>Query:</B> %s %c<BR>\n", querypdbcode, querychainid);
endHTML();
#endif
*/

		/*-------------*/
		/* NR group ID */
		/*-------------*/
		querygroupid = -1;
		for (ichain = 0; ichain < nchains; ++ichain)
		{
			if (0 == strcmp(querypdbcode, pdbcode[ichain]) &&
				querychainid == chainid[ichain])
			{
				querygroupid = groupid[NRlevel][ichain];
				break;
			}
		}
		if (querygroupid < 0)
		{
			beginHTML();
			fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
			fprintf(stdout, "<H3>The chain does not exit in this NR set.</H3>\n");
			endHTML();
			EarlyOut();
		}

		beginHTML();
		fprintf(stdout, "<H2>Cluster of Sequence-Similar Chains</H2>\n");
		fprintf(stdout, "<B>Query:</B> <A href=\"%s%d\">%s %c</A><BR>\n",
			ENTREZSTRUC, mmdbid[ichain], querypdbcode, querychainid);
		fprintf(stdout, "<B>No. of chains in the cluster:</B> %d<BR>\n",
			nmem[NRlevel][groupid[NRlevel][ichain]]);
		fprintf(stdout, "<B>Non-redundancy level: </B>\n");
		fprintf(stdout, "%s<BR>", chNRlevel[NRlevel]);
		fprintf(stdout, "<BR>");
		fflush(stdout);
		ShowNRGroup(NRlevel, querygroupid, ichain);
		endHTML();
	}
	else
	{
		beginHTML();
		fprintf(stdout, "<H2>%s</H2>\n", THIS_PROGRAM_IS);
		fprintf(stdout, "<H3>Undefined action requested.</H3>\n");
		endHTML();
		EarlyOut();
	}

	return(0);
}
