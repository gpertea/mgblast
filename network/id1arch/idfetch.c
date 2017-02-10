/*  idfetch.c
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
 * Author Karl Sirotkin
 *
 $Log: idfetch.c,v $
 Revision 1.38  2005/05/16 23:18:34  vysokolo
 Added features 'HPRD' and 'STS' to the key '-F'.

 Revision 1.37  2005/04/13 14:38:12  kans
 prototype for TryGetGi, send NORMAL_STYLE to SeqEntryToGnbk again

 Revision 1.36  2004/10/19 21:51:29  vysokolo
 Bug fix of -s key

 Revision 1.35  2004/10/12 21:39:28  vysokolo
 Added intervals of accessions like: "ABC_000123-ABC_000456"

 Revision 1.34  2004/10/04 19:30:25  vysokolo
 The "strcasecmp" replaced by "StringICmp"

 Revision 1.33  2004/09/30 17:59:26  vysokolo
 Added key -F to enable features by name.

 Revision 1.32  2004/05/25 18:41:35  kans
 removed obsolete STREAM_SEQ_PORT_FIRST flag

 Revision 1.31  2004/02/18 22:18:45  yaschenk
 adding recognition of gnl|sat_name|ent seqids

 Revision 1.30  2004/02/03 21:25:16  yaschenk
 relaxing ranges for -g and -e

 Revision 1.29  2003/12/17 20:35:38  kans
 initialize status, send NORMAL_STYLE to SeqEntrytoGnbk instead of 0 (also fixed in asn2gnbk), pass lookup flags

 Revision 1.28  2003/11/19 16:35:19  yaschenk
 relaxing ranges for -g and -c

 Revision 1.27  2003/03/28 18:48:39  yaschenk
 tuning ObjMgr, adding STREAM_SEQ_PORT_FIRST to SeqEntryToGnbk

 Revision 1.26  2003/01/29 23:08:19  yaschenk
 fixing FASTA on far pointers

 Revision 1.25  2003/01/21 22:27:23  kans
 new CstType parameter for flatfile generator

 Revision 1.24  2002/12/30 22:36:53  yaschenk
 optimizing..

 Revision 1.23  2002/11/07 17:21:55  yaschenk
 switching ID1 to new displatcher

 Revision 1.22  2002/07/23 19:31:43  butanaev
 Filtered out gi -1

 Revision 1.21  2001/11/02 14:24:44  kans
 made Fasta style SeqId args multi-line for Mac window

 Revision 1.20  2001/11/02 12:36:20  kans
 now using public Entrez2 server

 Revision 1.19  2001/09/28 15:56:04  kans
 look for extra and title fields in Entrez2 docsum

 Revision 1.18  2001/09/10 21:09:36  kans
 changed to use new Entrez2DocsumDataPtr - still need to get example of field_name keys

 Revision 1.17  2001/02/12 21:57:11  butanaev
 Made 3 retries to EntrezSynchronousQuery() when the NULL is returned.

 Revision 1.16  2001/02/08 16:13:46  yaschenk
 fixing wrong check for missing version in _PIR and SP

 Revision 1.15  2000/10/06 22:59:44  yaschenk
 strncpy not setting \0 bug

 Revision 1.14  2000/08/10 15:17:38  butanaev
 Updated -t 7 mode: strings like 'gi|3|emb|A00003.1|A00003' retreived from
 Entrez2DocsumPtr->caption.

 Revision 1.13  2000/08/03 17:01:23  kans
 included ni_lib.h for Mac, removed Mac compiler warnings

 Revision 1.12  2000/08/02 16:55:28  yaschenk
 increasing buffer size to 1000

 Revision 1.11  2000/08/02 16:17:00  butanaev
 Added:
 -t 7 - to retrieve Entrez DocSums
 -Q filename - to read Entrez query from the file

 Revision 1.9  2000/07/13 16:46:54  yaschenk
 adding ObjMgrFreeCache(0) to avoid hitting the limit in ObrMgr

 Revision 1.2  2000/06/01 18:05:35  butanaev
 Fixed numerous bugs with control flow...

 Revision 1.1  2000/06/01 16:48:22  butanaev
 New functionality:
 -G parameter, which previously accepted the list of gi's,
 now accepts gi,accession,accession.version,fasta seqid,
 which can be mixed.

 -q parameter generates the list out of Entrez query
   when -q is used -d has special meaning:
   -d n - run query against Nucleotide database
   -d p - run query against Protein database

 -n parameter limits the output to the list of gi's

 Revision 1.6  2000/05/24 17:30:42  yaschenk
 make parameter list look better

 Revision 1.5  2000/05/23 15:42:08  yaschenk
 adding quality score display

 Revision 1.4  2000/03/31 18:35:58  yaschenk
 Adding Jonathan's logic for FF and FASTA

 Revision 1.3  2000/03/30 20:43:51  yaschenk
 adding AsnIoReset between Entries

 Revision 1.2  1999/11/02 18:27:43  yaschenk
 adding -G parameter to idfetch

 Revision 1.1  1998/12/28 17:56:29  yaschenk
 preparing idfetch to go to production

 Revision 1.1  1997/05/29 14:34:07  sirotkin
 syncing sampson from mutant for procs. taking source from sampson. this is now current

 * Revision 4.0  1995/07/26  13:55:55  ostell
 * force revision to 4.0
 *
 * Revision 1.3  1995/06/21  14:14:29  kans
 * replaced asn2ff_entrez with SeqEntryToFlat
 *
 * Revision 1.2  1995/05/17  17:59:15  epstein
 * add RCS log revision history
 *
 * Revision 1.1  94/08/11  13:26:31  ostell
 * Initial revision
 *
 * Revision 1.3  1993/12/02  10:12:41  kans
 * Includes <ncbi.h> instead of <sys/types.h>
 *
 * Revision 1.2  93/11/24  13:25:56  sirotkin
 * First working version
 *
 * Revision 1.1  93/11/23  16:01:51  sirotkin
 * Initial revision
 *
 revised by OStell for public use.
 *
 * Modified by Eugene Yaschenko for ID1 Server
 *
 */
#include <ncbi.h>
#include <objsset.h>
#include <accid1.h>

#if 0
#include <asn2ff.h>
#else
#include <asn2gnbk.h>
#endif

#include <tofasta.h>
#include <ni_types.h>
#include <ni_lib.h>
#include <sqnutils.h>
#include <sequtil.h>
#include <ffprint.h>
#include <ent2api.h>

static Boolean ProcessOneDocSum (Int4 num, Int4Ptr uids);
static void EntrezQuery(char *query);

static Int4 BEGetUidsFromQuery(CharPtr query, Uint4Ptr PNTR uids,
                               Boolean is_na, Boolean count_only);
static Boolean IdFetch_func1(CharPtr data, Int2 maxplex);
static Boolean IdFetch_func(Int4 gi,CharPtr db, Int4 ent,Int2 maxplex);

static int TryGetGi(int choice, char *accession, char *name, int version);

Args myargs[] = {
	{"Filename for output ","stdout", NULL,NULL,FALSE,'o',ARG_FILE_OUT, 0.0,0,NULL},
	{"Output type:\t\
1=text asn.1\n\t\t\t\
2=Binary asn.1\n\t\t\t\
3=Genbank (Seq-entry only)\n\t\t\t\
4=genpept (Seq-entry only)\n\t\t\t\
5=fasta (table for history)\n\t\t\t\
6=quality scores (Seq-entry only)\n\t\t\t\
7=Entrez DocSums\n","1", "1", "7", FALSE, 't', ARG_INT, 0.0, 0, NULL } ,
{"Database to use (special meaning for -q flag: n - nucleotide, p - protein)",NULL,NULL,NULL,TRUE,'d',ARG_STRING,0.0,0,NULL},
	{"Entity number (retrieval number) to dump" ,"0",NULL,NULL,TRUE,'e',ARG_INT,0.0,0,NULL},
        {"Type of lookup:\t\
0 - get Seq-entry\n\t\t\t\
1 - get gi state (output to stderr)\n\t\t\t\
2 - get SeqIds\n\t\t\t\
3 - get gi historyn (sequence change only)\n\t\t\t\
4 - get gi revision history (any change to asn.1)\n", "0","0","4",TRUE,'i',ARG_INT,0.0,0,NULL},
	{"GI id for single Entity to dump" ,"0",NULL,NULL,TRUE,'g',ARG_INT,0.0,0,NULL},
	{"File with list of gi's, accessions, accession.version's, fasta seqid's to dump",NULL,NULL,NULL,TRUE,'G',ARG_FILE_IN,0.0,0,NULL},
	{"Max complexity:\t\
0 - get the whole blob\n\t\t\t\
1 - get the bioseq of interest\n\t\t\t\
2 - get the minimal bioseq-set containing the bioseq of interest\n\t\t\t\
3 - get the minimal nuc-prot containing the bioseq of interest\n\t\t\t\
4 - get the minimal pub-set containing the bioseq of interest\n" ,"0",NULL,NULL,TRUE,'c',ARG_INT,0.0,0,NULL},
 	{"flaTtened SeqId, format: \n		\'type(name,accession,release,version)\'\n			as \'5(HUMHBB)\' or \n		type=accession, or \n		type:number ",
		NULL,NULL,NULL,TRUE,'f',ARG_STRING,0.0,0,NULL},
 	{"Fasta style SeqId ENCLOSED IN QUOTES:\n\t\t\t\
lcl|int or str bbs|int bbm|int gb|acc|loc\n\t\t\t\
emb|acc|loc pir|acc|name sp|acc|name\n\t\t\t\
pat|country|patent|seq gi|int dbj|acc|loc\n\t\t\t\
prf|acc|name pdb|entry|chain",
	NULL,NULL,NULL,TRUE,'s',ARG_STRING,0.0,0,NULL},
        {"Log file", NULL,NULL,NULL,TRUE,'l',ARG_FILE_OUT,0.0,0,NULL},
        {"Generate gi list by entrez query", NULL,NULL,NULL,TRUE,'q',ARG_STRING,0.0,0,NULL},
        {"Generate gi list by entrez query", NULL,NULL,NULL,TRUE,'Q',ARG_FILE_IN,0.0,0,NULL},
        {"Output only the list of gis, used with -q", NULL,NULL,NULL,TRUE,'n',ARG_BOOLEAN,0.0,0,NULL},
        {"Add features delimited by ','. Allowed values are: 'CDD', 'SNP', 'SNP_graph', 'MGC', 'HPRD', 'STS'.", NULL,NULL,NULL,TRUE,'F',ARG_STRING,0.0,0,NULL}
};
int Numarg = sizeof(myargs)/sizeof(myargs[0]);

#define MACRO_SETARG(TAG,P) \
   {\
      P = Nlm_WhichArg (TAG, Numarg, myargs);\
      if( P < 0){\
         ErrPost(CTX_NCBIIDRETRIEVE,10,\
         "Program error looking for arg %c\n", TAG);\
         has_trouble = TRUE;\
      }\
   }

static Nlm_Int2 Nlm_WhichArg PROTO(( Nlm_Char which, Nlm_Int2 numargs, Nlm_ArgPtr ap));
static void MyBioseqToFasta(BioseqPtr bsp, Pointer userdata);

static Boolean CreateMaxPlexParam();
static Int4 GetIntervalAccession( const Char* pAccession, Char* pResult);

Int4 giBuffer[1000];
int giBufferPos = 0;

DataVal Val;
Int2 fileoutarg, logarg, outtypearg,maxplexarg,seqidarg,
  giarg, fastaarg, infotypearg, entarg, dbarg, gifilearg,
  entrezqueryarg, entrezqueryfilearg, onlylistarg, showfeatures, maxplex_param;

FILE * fp = NULL;
AsnIoPtr asnout=NULL;

Int2 Main()
{
  Boolean has_trouble = FALSE;
  Int4 entity_spec_count = 0;
  CharPtr outmode;
  Int4 ent = 0;
  Int4 gi = 0;
  FILE * fp_in = NULL;
  SeqIdPtr sip;
  Char tbuf[1024];
  static CharPtr accession = NULL;
  int type_int;

  /* check command line arguments */

  if( ! GetArgs("idfetch.c",Numarg, myargs))
    return 1;


  /********************************************************************
   ****                                                            ****
   ****  Map Args So Can be Accessed in order independent fashion  ****
   ****                                                            ****
   *******************************************************************/

  MACRO_SETARG('o', fileoutarg)
  MACRO_SETARG('t', outtypearg)
  MACRO_SETARG('i', infotypearg)
  MACRO_SETARG('c', maxplexarg)
  MACRO_SETARG('e', entarg)
  MACRO_SETARG('d', dbarg)
  MACRO_SETARG('g', giarg)
  MACRO_SETARG('G', gifilearg)
  MACRO_SETARG('f', seqidarg)
  MACRO_SETARG('l', logarg)
  MACRO_SETARG('s', fastaarg)
  MACRO_SETARG('q', entrezqueryarg)
  MACRO_SETARG('Q', entrezqueryfilearg)
  MACRO_SETARG('n', onlylistarg)
  MACRO_SETARG('F', showfeatures)

  if( !CreateMaxPlexParam())
  {
    has_trouble=TRUE;
    goto FATAL;
  }

  if(! SeqEntryLoad())
    ErrShow();

  EntrezSetProgramName("IDFETCH");
  /* EntrezSetServer("www.ncbi.nlm.nih.gov", 80,
                  "/entrez/utils/entrez2server.fcgi");
  */

  if(myargs[entrezqueryarg].strvalue || myargs[entrezqueryfilearg].strvalue)
  {
    if(! myargs[dbarg].strvalue ||
       (0 != strcmp(myargs[dbarg].strvalue, "n") &&
        0 != strcmp(myargs[dbarg].strvalue, "p")))
    {
      ErrPost(CTX_NCBIIDRETRIEVE, 10, "-q should be with protein (-dp) or nucleotide (-dn)\n");
      exit(1);
    }
  }

  if(myargs[logarg].strvalue != NULL)
  {
    if(! ErrSetLog (myargs[logarg].strvalue))
    {
      ErrShow();
      has_trouble = TRUE;
    }
    else
    {
      ErrSetOpts (ERR_TEE, ERR_LOG_ON);
    }
  }
  if(myargs[infotypearg].intvalue>1
     && (myargs[outtypearg].intvalue == 3
         || myargs[outtypearg].intvalue == 4
         || myargs[outtypearg].intvalue == 6
        ))
  {
    ErrPostEx(SEV_ERROR,0,0,"-t 3,4,6 can be used only with -i 0");
    has_trouble=TRUE;
    goto FATAL;
  }

  if(myargs[outtypearg].intvalue == 7 && ! myargs[dbarg].strvalue)
  {
    ErrPostEx(SEV_ERROR,0,0,"-t 7 can be used only with protein (-dp) or nucleotide (-dn)");
    has_trouble=TRUE;
    goto FATAL;
  }

  if(myargs[outtypearg].intvalue == 7)
    myargs[infotypearg].intvalue = 2;

  if(myargs[entrezqueryarg].strvalue || myargs[entrezqueryfilearg].strvalue)
    entity_spec_count++;

  if(myargs[giarg].intvalue)
    entity_spec_count++;

  if(myargs[seqidarg].strvalue)
    entity_spec_count++;

  if(myargs[fastaarg].strvalue)
    entity_spec_count++;

  if(myargs[gifilearg].strvalue)
  {
    fp_in=FileOpen(myargs[gifilearg].strvalue,"r");
    if(fp_in==NULL)
    {
      ErrPostEx(SEV_ERROR, 0, 0, "couldn't open file %s", myargs[gifilearg].strvalue);
      has_trouble=TRUE;
      goto FATAL;
    }
    entity_spec_count ++;
  }

  if(entity_spec_count != 1)
  {
    ErrPostEx(SEV_ERROR,0,0, "One and only one parameters may be used: -g,-G,-f,-s -q");
    has_trouble=TRUE;
    goto FATAL;
  }

  {  /** tuning ObjMgr **/
	ObjMgrPtr	omp;
        omp=ObjMgrGet();
        omp->maxtemp=500;
        /*omp->maxobj=64000;*/
        omp->autoclean = TRUE;

  }

  if(myargs[infotypearg].intvalue != 1)
  {
    outmode = "w";
    switch(myargs[outtypearg].intvalue)
    {
    case 2:
      outmode = "wb";
    case 1:
      asnout = AsnIoOpen((CharPtr)myargs[fileoutarg].strvalue, outmode);
      if(asnout == NULL)
      {
        ErrPost(CTX_NCBIIDRETRIEVE,10,
                "Could not open %s for asn output\n",
                myargs[fileoutarg].strvalue);
        has_trouble = TRUE;
      }
      break;
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      fp = FileOpen((CharPtr)myargs[fileoutarg].strvalue, outmode);
      if(fp == NULL)
      {
        ErrPost(CTX_NCBIIDRETRIEVE,10,
                "Could not open %s for output\n",
                myargs[fileoutarg].strvalue);
        has_trouble = TRUE;
      }
      break;
    }
  }

  if( has_trouble )
    exit (1);


	if(   fp_in 
	   || myargs[entrezqueryarg].strvalue
           || myargs[entrezqueryfilearg].strvalue
           || myargs[outtypearg].intvalue == 3
           || myargs[outtypearg].intvalue == 4
           || myargs[outtypearg].intvalue == 5
              ) { /*** Statefull mode ***/
		putenv("CONN_STATELESS=FALSE"); 
	}
	else { /*** Stateless mode ***/
		putenv("CONN_STATELESS=TRUE"); 
	}

  if(!ID1BioseqFetchEnable("idfetch",TRUE))
  {
    ErrPost(CTX_NCBIIDRETRIEVE, 20, "Could not open ID1 service");
    exit(1);
  }

  if(myargs[giarg].intvalue)
  {
    gi = myargs[giarg].intvalue;
  }
  else if(myargs[entrezqueryarg].strvalue || myargs[entrezqueryfilearg].strvalue)
  {}
  else if(fp_in)
  {}
  else if(myargs[fastaarg].strvalue != NULL)
  {
    sip = SeqIdParse((CharPtr)myargs[fastaarg].strvalue);
    if(sip == NULL)
    {
      ErrPostEx(SEV_ERROR,0,0,"Couldn't parse [%s]", myargs[fastaarg].strvalue);
      exit(1);
    }
  }
  else
  {
    /*
     "flaTtened SeqId, format:
     type(name,accession,release,version) or type=accession",
     */
    static CharPtr name = NULL, release = NULL, version = NULL, number = NULL;
    CharPtr p;
    static CharPtr PNTR fields [] = {&name, &accession, &release, &number};
    Boolean found_equals = FALSE, found_left = FALSE,
      found_colon = FALSE, flat_seqid_err = FALSE,
    dna_type = FALSE, any_type = FALSE;
    int dex;
    TextSeqIdPtr tsip;

    type_int = atoi(myargs[seqidarg].strvalue);
    for(p = myargs[seqidarg].strvalue; *p; p++ )
    {
      if( *p == '(' || *p == '='  || *p == ':')
      {  /* ) match */

        if( *p == '('  )
        {  /* ) match */
          found_left = TRUE;
          if(p == myargs[seqidarg].strvalue)
          {
            any_type = TRUE;
            ErrPost(CTX_NCBIIDRETRIEVE,10,
                    "Sorry, any type is not allowed for ID1service");
            exit(1);
          }
          else if( p - myargs[seqidarg].strvalue == 1)
          {
            if(*myargs[seqidarg].strvalue == '0')
            {
              dna_type = TRUE;
              ErrPost(CTX_NCBIIDRETRIEVE,10,
                      "Sorry, 0== nucroe3 type is not allowed for ID1service");
              exit(1);
            }
          }
        }
        else if( *p == '=')
        {
          found_equals = TRUE;
          if(p == myargs[seqidarg].strvalue)
          {
            any_type = TRUE;
          }
          else if( p - myargs[seqidarg].strvalue == 1)
          {
            if(*myargs[seqidarg].strvalue == '0')
            {
              dna_type = TRUE;
            }
          }
        }
        else if( *p == ':')
        {
          found_colon = TRUE;
          if(p == myargs[seqidarg].strvalue)
          {
            any_type = TRUE;
          }
          else if( p - myargs[seqidarg].strvalue == 1)
          {
            if(*myargs[seqidarg].strvalue == '0')
            {
              dna_type = TRUE;
            }
          }
        }
        *p = '\0';
        p++;
        break;
      }
    }
    if( found_left)
    {
      for( * (fields[0]) = p, dex = 0; *p && dex < 4; p++)
      {
        if( *p == ',' || *p == ')' )
        {
          *p = '\0';
          dex ++;
          *(fields[dex]) = p + 1;
        }
      }
    }
    else if(found_equals)
    {
      accession = p;
    }
    else if(found_colon)
    {
      number = p;
    }
    else
    {
      ErrPost(CTX_NCBIIDRETRIEVE, 10,
              "id1test: could not find \'(\' or \'=\' or \':\' in flattened seqid=%s",myargs[seqidarg].strvalue);  /* ) match */
      exit(1);
    }
    sip = ValNodeNew(NULL);
    sip->choice = type_int;
    switch(type_int)
    {
    case SEQID_GIBBSQ :
    case SEQID_GIBBMT :
    case SEQID_GI :
      sip->data.intvalue = atoi(number);
      break;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_PIR :
    case SEQID_SWISSPROT :
    case SEQID_OTHER :
    case SEQID_PRF :
      tsip = TextSeqIdNew();
      sip->data.ptrvalue = tsip;
      if(accession)
        if(!*accession)
          accession = NULL;
      if(release)
        if(!*release)
          release = NULL;
      if(name)
        if(!*name)
          name = NULL;
      tsip->name = StringSave(name);
      tsip->accession = StringSave(accession);
      tsip->release = StringSave(release);
      break;
    case SEQID_PATENT :
    case SEQID_GENERAL :
    case SEQID_PDB :
    case SEQID_LOCAL :
    default:
      ErrPost(CTX_NCBIIDRETRIEVE,30,
              "Sorry, this test program does not support %d patent, general, pdb, or local, try id2asn ",
              type_int);
      exit(1);
      break;
    }
  }

  if(! fp_in && ! gi && ! (myargs[entrezqueryarg].strvalue || myargs[entrezqueryfilearg].strvalue))
  {
    Char tacc[64];
    if( sip->data.ptrvalue && GetIntervalAccession( accession, tacc) > 0 && myargs[fastaarg].strvalue == NULL)
    {
      do{ 
      	Int4 tgi = TryGetGi( type_int, tacc, NULL, 0);
	if(tgi <= 0)
	{
	  SeqIdPrint(sip, tbuf, PRINTID_FASTA_SHORT);
	  ErrPostEx(SEV_ERROR, 0, 0, "Couldn't find SeqId [%s]", tbuf);
	  has_trouble=TRUE;
	  break;
	}
	if( !IdFetch_func(tgi,myargs[dbarg].strvalue, myargs[entarg].intvalue,maxplex_param))
	{
	  has_trouble=TRUE;
	  break;
	}

	if(fp) fflush(fp);
      }
      while( GetIntervalAccession( NULL, tacc) > 0);
      if(fp) fflush(fp);
      if( has_trouble ) goto FATAL;
    }
    else
    {
      gi = ID1ArchGIGet (sip);
      if(gi <= 0)
      {
	  if(sip->choice==SEQID_GENERAL){
		  Dbtag *db = (Dbtag*)sip->data.ptrvalue;
		  gi=myargs[entarg].intvalue=db->tag->str?atoi(db->tag->str):db->tag->id;
		  myargs[dbarg].strvalue=StringSave(db->db);
	  } else {
		SeqIdPrint(sip, tbuf, PRINTID_FASTA_SHORT);
		ErrPostEx(SEV_ERROR, 0, 0, "Couldn't find SeqId [%s]", tbuf);
		goto FATAL;	
	  }
      }
    }
  }
  else if(fp_in)
  {
    Char tacc[1024];
    while(fgets(tbuf,sizeof(tbuf)-1,fp_in)){
      if( GetIntervalAccession( tbuf, tacc) > 0 )
      {
	do{ 
	  IdFetch_func1(tacc, maxplex_param); 
	  if(fp) fflush(fp);
	}
	while( GetIntervalAccession( NULL, tacc) > 0);
      }
      else
      {
	IdFetch_func1(tbuf, maxplex_param);
      }	      
      if(fp) fflush(fp);
    }
  }
  else if(myargs[entrezqueryarg].strvalue || myargs[entrezqueryfilearg].strvalue)
  {
    if(myargs[entrezqueryarg].strvalue)
      EntrezQuery(myargs[entrezqueryarg].strvalue);
    else if(myargs[entrezqueryfilearg].strvalue)
    {
      FILE *fp_query=FileOpen(myargs[entrezqueryfilearg].strvalue,"r");
      if(fp_query==NULL)
      {
        ErrPostEx(SEV_ERROR, 0, 0, "couldn't open file %s", myargs[entrezqueryfilearg].strvalue);
        has_trouble=TRUE;
        goto FATAL;
      }
      {
        char buffer[2000];
        while(fgets(buffer, sizeof(buffer) - 1,fp_query))
        {
          int len = strlen(buffer);
          if(buffer[len - 1] == '\n')
            buffer[len - 1] = 0;
          EntrezQuery(buffer);
        }
      }
    }
  }

  if(gi>0 && !IdFetch_func(gi,myargs[dbarg].strvalue, myargs[entarg].intvalue,maxplex_param))
  {
    has_trouble=TRUE;
    goto FATAL;
  }

  if(giBufferPos != 0)
  {
    if(! ProcessOneDocSum(giBufferPos, giBuffer))
    {
      has_trouble=TRUE;
      goto FATAL;
    }
  }

FATAL:
  if(asnout)
    AsnIoClose(asnout);
  if(fp)
    FileClose(fp);
  if(fp_in)
    FileClose(fp_in);

  /*ID1ArchFini();*/

  return(has_trouble?1:0);
}

void EntrezQuery(char *query)
{
  Uint4 *ids;
  int count;
  int i;
  count = BEGetUidsFromQuery(query,
                             &ids,
                             0 == strcmp(myargs[dbarg].strvalue, "n"), FALSE);
  for(i = 0; i < count; ++i)
  {
    if(ids[i] == -1)
      continue;

    if(myargs[onlylistarg].intvalue)
      printf("%d\n", ids[i]);
    else
      IdFetch_func(ids[i],
                   myargs[dbarg].strvalue,
                   myargs[entarg].intvalue,
		   maxplex_param);
  }
}

static void PrintQualScores (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr  bsp;
  FILE       *fp;

  if(IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    fp = (FILE*) data;
    PrintQualityScores (bsp, fp);
  }
}

static int TryGetGi(int choice, char *accession, char *name, int version)
{
  TextSeqIdPtr tsip;
  SeqIdPtr sip;

  sip = ValNodeNew(NULL);
  tsip = TextSeqIdNew();
  sip->choice = choice;
  sip->data.ptrvalue = tsip;
  tsip->name = NULL;
  tsip->accession = NULL;
  tsip->release = NULL;

  if(name != NULL)
    tsip->name = StringSave(name);
  if(accession != NULL)
    tsip->accession = StringSave(accession);
  tsip->version = version;
  return ID1ArchGIGet(sip);
}

static Boolean IdFetch_func1(CharPtr data, Int2 maxplex)
{
  int seqIdParse = 0, hasVersion = 0, accession = 0;
  int len = strlen(data), i;

  for(i = 0; i < len; ++i)
  {
    if(data[i] == '\n')
      data[i] = 0;
    if(data[i] == '|')
      seqIdParse = 1;
    if(data[i] == '.')
      hasVersion = i;
    if(isalpha(data[i]))
      accession = 1;
  }
  if(seqIdParse)
  {
    SeqIdPtr sip;
    Int4 gi,ent=0;
    CharPtr	sat=NULL;

   
    sip = SeqIdParse(data);
    if(sip == NULL)
    {
      ErrPostEx(SEV_ERROR,0,0,"Couldn't parse [%s]", data);
      exit(1);
    }
    gi = ID1ArchGIGet(sip);
    if(gi <= 0)
    {
	if(sip->choice==SEQID_GENERAL){
                Dbtag *db = (Dbtag*)sip->data.ptrvalue;
                gi=ent=db->tag->str?atoi(db->tag->str):db->tag->id;
                sat=StringSave(db->db);
	} else {
	      SeqIdPrint(sip, data, PRINTID_FASTA_SHORT);
	      ErrPostEx(SEV_ERROR, 0, 0, "Couldn't find SeqId [%s]", data);
	      exit(1);
	}
    }
    return IdFetch_func(gi, sat, ent, maxplex);
  }
  else if(hasVersion || accession)
  {
    char acc[100];
    int ver = 0;
    int gi;
    if(hasVersion)
    {
      strncpy(acc, data, hasVersion);
      acc[hasVersion]='\0';
      ver = atoi(data + hasVersion + 1);
    }
    else
      strcpy(acc, data);

    if(ver == 0)
      ver = INT2_MIN;

    if((gi = TryGetGi(SEQID_GENBANK, acc, NULL, ver)) ||
       (gi = TryGetGi(SEQID_OTHER, acc, NULL, ver)) ||
       (gi = TryGetGi(SEQID_GENBANK, NULL, acc, ver)) ||
       (gi = TryGetGi(SEQID_OTHER, NULL, acc, ver)))
      return IdFetch_func(gi,
                          myargs[dbarg].strvalue,
                          myargs[entarg].intvalue,
                          maxplex_param);

    if(ver == INT2_MIN &&
       (
        (gi = TryGetGi(SEQID_PIR, acc, NULL, ver)) ||
        (gi = TryGetGi(SEQID_SWISSPROT, acc, NULL, ver)) ||
        (gi = TryGetGi(SEQID_PIR, NULL, acc, ver)) ||
        (gi = TryGetGi(SEQID_SWISSPROT, NULL, acc, ver))
       )
      )
      return IdFetch_func(gi,
                          myargs[dbarg].strvalue,
                          myargs[entarg].intvalue,
                          maxplex_param);
    return 0;
  }
  else
    return IdFetch_func(atoi(data), NULL, 0, maxplex);
}

static Entrez2ReplyPtr MyEntrezSynchronousQuery(Entrez2RequestPtr e2rq)
{
  int i;
  for(i = 0; i < 3; ++i)
  {
    Entrez2ReplyPtr reply = EntrezSynchronousQuery(e2rq);
    if(reply != NULL)
      return reply;
  }
  return NULL;
}

static Boolean ProcessOneDocSum (Int4 num, Int4Ptr uids)

{
  Entrez2DocsumPtr      dsp;
  Entrez2DocsumListPtr  e2dl;
  Entrez2RequestPtr     e2rq;
  Entrez2ReplyPtr       e2ry;
  Entrez2DocsumDataPtr  e2ddp;
  CharPtr db;
  Boolean result;

  if(num == 0)
    return TRUE;

  db = strcmp(myargs[dbarg].strvalue, "n") == 0 ? "Nucleotide" : "Protein";
  e2rq = EntrezCreateDocSumRequest (db, 0, num, uids, NULL);
  if (e2rq == NULL) return FALSE;

  e2ry = MyEntrezSynchronousQuery (e2rq);
  e2rq = Entrez2RequestFree(e2rq);
  e2dl = EntrezExtractDocsumReply (e2ry);

  if (e2dl == NULL)
    return FALSE;

  result = TRUE;
  for (dsp = e2dl->list; dsp != NULL; dsp = dsp->next)
  {
    char *title = "";
    char *extra = "";
    for (e2ddp = dsp->docsum_data; e2ddp != NULL; e2ddp = e2ddp->next) {
      if (StringHasNoText (e2ddp->field_value)) continue;
      if (StringICmp (e2ddp->field_name, "Title") == 0) {
        title = e2ddp->field_value;
      } else if (StringICmp (e2ddp->field_name, "Extra") == 0) {
        extra = e2ddp->field_value;
      }
    }

    fprintf(fp,">%s %s\n", extra, title);

  }

  Entrez2DocsumListFree (e2dl);
  return result;
}

static Boolean IdFetch_func(Int4 gi,CharPtr db, Int4 ent,Int2 maxplex)
{
  SeqEntryPtr	sep=NULL;
  Int4		status = 0,gi_state;
  SeqIdPtr	sip_ret=NULL;
  SeqId		si={SEQID_GI,0,0};
  ID1SeqHistPtr	ishp=NULL;
  Char		buf[200],user_string[100];
  ErrStrId        utag;
  Boolean		retval=TRUE;
  BioseqPtr	bsp=NULL;
  Uint2		entityID;
  Uint1		group_segs;

  if(myargs[outtypearg].intvalue == 7)
  {
    Boolean result = TRUE;
    if(giBufferPos == sizeof(giBuffer) / sizeof(giBuffer[0]))
    {
      result =  ProcessOneDocSum(giBufferPos, giBuffer);
      giBufferPos = 0;
    }
    giBuffer[giBufferPos++] = gi;
    return result;
  }

  sprintf(user_string,"GI=%d|db=%s|ent=%d|",gi,db?db:"NULL",ent);
  utag=ErrUserInstall(user_string,0);

  switch(myargs[infotypearg].intvalue)
  {
  case 0:
    if(maxplex == 1 && myargs[outtypearg].intvalue != 1 && myargs[outtypearg].intvalue != 2)
    {
      si.data.intvalue=gi;
      if((bsp=BioseqLockById(&si)) != NULL)
      {
        sep = ObjMgrGetChoiceForData (bsp);
        switch(myargs[outtypearg].intvalue)
        {
        case 3:
        case 4:
          if(bsp->repr == Seq_repr_seg)
          {
            entityID = ObjMgrGetEntityIDForChoice (sep);
            sep = GetBestTopParentForData (entityID, bsp);
          }
          break;
        case 5:
          if(ISA_na (bsp->mol))
          {
            group_segs = 0;
            if(bsp->repr == Seq_repr_seg)
            {
              group_segs = 1;
            }
            else if(bsp->repr == Seq_repr_delta)
            {
              group_segs = 3;
            }
          }
          else
          {
            group_segs=0;
          }
          break;
        }
      }
    }
    else
    {
      sep = ID1ArchSeqEntryGet (gi,db,ent,&status,maxplex);
    }
    if(!sep)
    {
      switch(status)
      {
      case 1:
        ErrPostEx(SEV_WARNING,0,0,"Sequence has been withdrawn");
        break;
      case 2:
        ErrPostEx(SEV_WARNING,0,0,"Sequence is not yet available");
        break;
      case 10:
        ErrPostEx(SEV_WARNING,0,0,"GI <%d> is not found",gi );
        break;
      default:
        ErrPostEx(SEV_WARNING,0,0,"Unable to read ASN.1");
        ID1ArchFini();
        ID1ArchInit();
      }
      retval=FALSE;
      goto DONE;
    }
    if(status==3)
      ErrPostEx(SEV_INFO,0,0,"IS_DEAD");
    break;
  case 1:
    gi_state = ID1ArcgGIStateGet(gi);
    break;
  case 2:
    sip_ret = ID1ArchSeqIdsGet(gi,asnout);
    break;
  case 3:
    ishp = ID1ArchGIHistGet(gi,FALSE,asnout);
    break;
  case 4:
    ishp = ID1ArchGIHistGet(gi,TRUE,asnout);
    break;
  }

  if(myargs[infotypearg].intvalue == 1)
  {
    Char	buf[200];
    id_print_gi_state(gi_state,buf,sizeof(buf));
    printf("gi= %d, states: %s\n",gi,buf);
  }
  else
  {
    putenv("CONN_STATELESS=FALSE");  /*** some formats may ask a lot of information ***/
    switch(myargs[outtypearg].intvalue)
    {
    case 1:
    case 2:
      switch(myargs[infotypearg].intvalue)
      {
      case 0:
        SeqEntryAsnWrite(sep, asnout, NULL);
        AsnIoReset(asnout);
        break;
      }
      break;
    case 3:
      AssignIDsInEntity(0,OBJ_SEQENTRY,sep);
      if(!SeqEntryToGnbk(sep,NULL,GENBANK_FMT,ENTREZ_MODE,NORMAL_STYLE,SHOW_CONTIG_FEATURES|ONLY_NEAR_FEATURES,
         LOOKUP_FAR_COMPONENTS|LOOKUP_FAR_LOCATIONS|LOOKUP_FAR_PRODUCTS|LOOKUP_FAR_HISTORY|LOOKUP_FAR_OTHERS,
         0,NULL,fp)){
        ErrPostEx(SEV_WARNING,0,0,
                  "GenBank Format does not exist for this sequence ");
        retval=FALSE;
        goto DONE;
      }
      break;
    case 4:
      AssignIDsInEntity(0,OBJ_SEQENTRY,sep);
      if(!SeqEntryToGnbk(sep,NULL,GENPEPT_FMT,ENTREZ_MODE,NORMAL_STYLE,SHOW_CONTIG_FEATURES|ONLY_NEAR_FEATURES,
        LOOKUP_FAR_COMPONENTS|LOOKUP_FAR_LOCATIONS|LOOKUP_FAR_PRODUCTS|LOOKUP_FAR_HISTORY|LOOKUP_FAR_OTHERS,
        0,NULL,fp))
      {
        ErrPostEx(SEV_WARNING,0,0,
                  "GenPept Format does not exist for this sequence");
        retval=FALSE;
        goto DONE;
      }
      break;
    case 6:
      SeqEntryExplore (sep, (Pointer) fp, PrintQualScores);
      break;
    case 5:
      switch(myargs[infotypearg].intvalue){
      case 0:
        if(bsp){
		MyBioseqToFasta(bsp,(Pointer)fp);
        }
        else
        {
		VisitBioseqsInSep(sep,(Pointer)fp, MyBioseqToFasta);
        }
        break;
      case 2:
        SeqIdWrite(sip_ret, buf, PRINTID_FASTA_LONG, sizeof(buf) - 1);
        fprintf(fp, "%s\n", buf);
        break;

      case 3:
      case 4:
        SeqHistPrintTable(ishp,fp);
        break;
      }
      break;
    }
  }
DONE:
  if(bsp)
  {
    static Uint2  reap_cnt;
    BioseqUnlock(bsp);
#if 0
    reap_cnt++;
    if(reap_cnt > 128){
      ObjMgrFreeCache(0);
      reap_cnt=0;
    }
#endif
  }
  else if(sep)
  {
    SeqEntryFree(sep);
  }
  if(sip_ret)	SeqIdFree(sip_ret);
  if(ishp)	ID1SeqHistFree(ishp);
  ErrUserDelete(utag);
  return retval;
}
/*****************************************************************************
 *
 *   Nlm_WhichArg(ap)
 *     returns array position for a tag
 *
 *****************************************************************************/
static Nlm_Int2 Nlm_WhichArg( Nlm_Char which, Nlm_Int2 numargs, Nlm_ArgPtr ap)
{
  Nlm_Boolean okay = FALSE;
  Nlm_Int2 i;
  Nlm_ArgPtr curarg;
  Nlm_Int2 retval = -1;

  if((ap == NULL) || (numargs == 0) )
    return okay;

  curarg = ap;                        /* set defaults */

  for(i = 0; i < numargs; i++, curarg++)
  {
    if(curarg->tag == which)
    {
      retval = i;
      break;
    }
  }

  return retval;
}

/* This function is interface to the Entrez2 engine. It may be used
   to get list of gis corresponding to the Entrez Boolean string or
   just number of such hits in the Entrez database */

static Int4 BEGetUidsFromQuery(CharPtr query, Uint4Ptr PNTR uids, 
                               Boolean is_na, Boolean count_only)
{
  Entrez2ReplyPtr e2ry;
  Entrez2RequestPtr  e2rq;
  E2ReplyPtr e2rp;
  Int4 count = 0;
  Entrez2BooleanReplyPtr e2br;
  Entrez2IdListPtr e2idlist;

  *uids = NULL;

  e2rq = EntrezCreateBooleanRequest(!count_only, FALSE,
                                    is_na? "Nucleotide" : "Protein",
                                    query, 0, 0, NULL, 0, 0);

  e2ry = MyEntrezSynchronousQuery(e2rq);

  if(e2ry == NULL)
  {
    ErrPostEx(SEV_ERROR, 0, 0, "NULL returned from EntrezSynchronousQuery()");
    return -1;
  }

  if((e2rp = e2ry->reply) == NULL)
  {
    ErrPostEx(SEV_ERROR, 0, 0, "Invalid ASN.1: E2ReplyPtr==NULL");
    return -1;
  }

  switch(e2rp->choice)
  {

  case E2Reply_error:
    ErrPostEx(SEV_ERROR, 0, 0, (CharPtr) e2rp->data.ptrvalue);
    count = -1;
    break;
  case E2Reply_eval_boolean:
    e2br = (Entrez2BooleanReplyPtr) e2rp->data.ptrvalue;
    count = e2br->count;
    if((e2idlist = e2br->uids) != NULL) {
      count = e2idlist->num;
      *uids = MemNew(sizeof(Int4)*count);
      BSSeek((ByteStorePtr) e2idlist->uids, 0, SEEK_SET);
      BSRead((ByteStorePtr) e2idlist->uids, *uids, sizeof(Int4)*count);

    }
    break;
  default:
    ErrPostEx(SEV_ERROR, 0, 0, "Invalid reply type from the server: %d", e2rp->choice);
    count = -1;
    break;
  }
  Entrez2ReplyFree(e2ry);
  Entrez2RequestFree(e2rq);
  return count;
}

#define FASTA_LINE_SIZE         70
#define FASTA_LINES_IN_CHUNK    5000

static void
MyBioseqToFasta(BioseqPtr bsp, Pointer userdata)
{
        SeqPortPtr      spp=NULL;
        Char            buf[2048];
        Char         	str[200];
        ValNodePtr      vnp;
        MolInfoPtr      mip=NULL;
	FILE PNTR fp=(FILE PNTR)userdata;
	Int4    start=0,step=FASTA_LINE_SIZE*FASTA_LINES_IN_CHUNK,stop;

	for(vnp=bsp->descr;vnp;vnp=vnp->next){
                if(vnp->choice==Seq_descr_molinfo){
                        mip = (MolInfoPtr)(vnp->data.ptrvalue);
                }
        }
	SeqIdWrite(bsp->id,str,PRINTID_FASTA_LONG,199);
	if(!CreateDefLine(NULL, bsp, buf, sizeof(buf), mip?mip->tech:0, NULL, NULL)) return;
	fprintf(fp, ">%s %s\n", str,buf);
	while(start < bsp->length){
		stop=start+step-1;
		if(stop >= bsp->length) stop=bsp->length-1;
		spp = SeqPortNew(bsp,start,stop,0, (ISA_na(bsp->mol))?Seq_code_iupacna:Seq_code_ncbieaa);
		if(spp==NULL) return;
		SeqPortSet_do_virtual(spp, TRUE);
		while (FastaSeqLineEx(spp, buf, FASTA_LINE_SIZE, ISA_na(bsp->mol),TRUE)) {
			fprintf(fp,"%s\n", buf);
			buf[0] = '\0';
		}
		if(spp) {
			SeqPortFree(spp);
			spp=NULL;
		}
		start=stop+1;
	}
}

/*
select * from annot_types;
 id          name                           is_private  dependencies
 ----------- ------------------------------ ----------- ------------
           1 SNP                                      0            0
           2 WGS descriptor                           1            0
           3 SNP graph                               -2            1
           4 CDD                                      0            0
           5 MGC                                      0            0
           6 HPRD                                     0            0
           7 STS                                      0            0

1 "SNP"       ffef 5
3 "SNP_graph" ffbf 7
4 "CDD"       ff7f 8
5 "MGC"       feff 9
6 "HPRD"      fdff 10
7 "STS"       fbff 11
*/

Boolean CreateMaxPlexParam()
{
  Char buf[1024];
  Char *ptoken = NULL;
  maxplex_param = 0xfffffff0 | myargs[maxplexarg].intvalue;

  if(myargs[showfeatures].strvalue)
  {
    strncpy( buf, myargs[showfeatures].strvalue, 1024);
    ptoken = strtok( buf, ",");
    while( ptoken)
    {
      if( !StringICmp( ptoken, "CDD"))
      {
	maxplex_param &= 0xffffff7f;
      }
      else if( !StringICmp( ptoken, "SNP"))
      {
	maxplex_param &= 0xffffffef;
      }
      else if( !StringICmp( ptoken, "SNP_graph"))
      {
	maxplex_param &= 0xffffffbf;
      }
      else if( !StringICmp( ptoken, "MGC"))
      {
	maxplex_param &= 0xfffffeff;
      }
      else if( !StringICmp( ptoken, "HPRD"))
      {
	maxplex_param &= 0xfffffdff;
      }
      else if( !StringICmp( ptoken, "STS"))
      {
	maxplex_param &= 0xfffffbff;
      }
      else
      {
	/* Error: unknown feature */
	ErrPostEx(SEV_ERROR,0,0,"Unknown feature type [%s]", ptoken);
	return FALSE;
      }
      ptoken = strtok( NULL, ",");
    }
  }
  return TRUE;
}

static Int4 GetPrefixIndexFromAcession( const Char* pSrc, Char* pDest, Int4* pDiditsLen)
{
  Int4 ret;
  while( *pSrc )
  {
    if( isdigit( *pSrc)) break;
    *pDest++ = *pSrc++;
  }
  *pDest = '\0';
  if( !*pSrc ) return -1;
  ret = atoi( pSrc);
  *pDiditsLen = 0;
  while( isdigit( *pSrc++)) ++(*pDiditsLen);
  return ret;
}

/*************************************************
* GetIntervalAccession parses string as 'ABC_0000123-ABC_0000456'
* and *returns* in 'pResult' strings like 'ABC_0000123', 'ABC_0000124', 'ABC_0000125', ..., 'ABC_0000456'
* iteratively using 'iterator'
* First call 'pAccession' must be a string like 'ABC_0000123-ABC_0000456', and all other times NULL
* returns: 0 if no more iterations
*         -1 if wrong format of input string
*         >0 if next iteration are present
**************************************************/
static Int4 GetIntervalAccession( const Char* pAccession, Char* pResult)
{
  static Int4 iterator=0, Start=0, End=0, DigitsStart=0;
  static Char pStartPrefix[64]={'\0'};
  static Char pFormat[64]={'\0'};
  
  if( pAccession )
  {
    Int4 DigitsEnd=0;
    Char pBuff[64]={'\0'};
    Char pEndPrefix[64];
    Char *pEnd, *pStart;

    iterator = 0;

    strncpy( pBuff, pAccession, 64);
    /*parse interval*/

    /*find dash*/
    pEnd = pStart = pBuff;
    while( *pEnd++ )
    {
      if( *pEnd == '-' ) break;
    }
    if( !*pEnd ) return -1; /*no dash found*/

    *pEnd++ = '\0';
    
    /*find index of start and end in the interval*/
    Start = GetPrefixIndexFromAcession( pStart, pStartPrefix, &DigitsStart);
    End = GetPrefixIndexFromAcession( pEnd, pEndPrefix, &DigitsEnd);

    if( DigitsStart != DigitsEnd || strcmp( pStartPrefix, pEndPrefix) || Start >= End || DigitsStart < 1)
    {
      return -1; /* wrong interval format */
    }
    sprintf( pFormat, "%%s%%0%hdhd", DigitsStart);
  }
  
  if( Start+iterator > End ) return 0; /*stop*/
  
  /* return next accession */
  sprintf( pResult, pFormat, pStartPrefix, Start+iterator);

  return ++iterator;
}

