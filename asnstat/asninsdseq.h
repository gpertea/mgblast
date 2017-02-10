/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asninsdseq.h16";
static AsnType atx[83] = {
  {401, "INSDSeq" ,1,0,0,0,0,0,0,0,NULL,&atx[48],&atx[1],0,&atx[19]} ,
  {0, "locus" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "length" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "strandedness" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[6]} ,
  {0, "moltype" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[7]} ,
  {0, "topology" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[8]} ,
  {0, "division" ,128,5,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[9]} ,
  {0, "update-date" ,128,6,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[10]} ,
  {0, "create-date" ,128,7,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[11]} ,
  {0, "update-release" ,128,8,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[12]} ,
  {0, "create-release" ,128,9,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[13]} ,
  {0, "definition" ,128,10,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[14]} ,
  {0, "primary-accession" ,128,11,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[15]} ,
  {0, "entry-version" ,128,12,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[16]} ,
  {0, "accession-version" ,128,13,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[17]} ,
  {0, "other-seqids" ,128,14,0,1,0,0,0,0,NULL,&atx[20],&atx[18],0,&atx[21]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[19],NULL,0,NULL} ,
  {402, "INSDSeqid" ,1,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[23]} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "secondary-accessions" ,128,15,0,1,0,0,0,0,NULL,&atx[20],&atx[22],0,&atx[24]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[23],NULL,0,NULL} ,
  {403, "INSDSecondary-accn" ,1,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[27]} ,
  {0, "project" ,128,16,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[25]} ,
  {0, "keywords" ,128,17,0,1,0,0,0,0,NULL,&atx[20],&atx[26],0,&atx[28]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[27],NULL,0,NULL} ,
  {404, "INSDKeyword" ,1,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[34]} ,
  {0, "segment" ,128,18,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[29]} ,
  {0, "source" ,128,19,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[30]} ,
  {0, "organism" ,128,20,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[31]} ,
  {0, "taxonomy" ,128,21,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[32]} ,
  {0, "references" ,128,22,0,1,0,0,0,0,NULL,&atx[20],&atx[33],0,&atx[52]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[34],NULL,0,NULL} ,
  {405, "INSDReference" ,1,0,0,0,0,0,0,0,NULL,&atx[48],&atx[35],0,&atx[58]} ,
  {0, "reference" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[36]} ,
  {0, "position" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[37]} ,
  {0, "authors" ,128,2,0,1,0,0,0,0,NULL,&atx[20],&atx[38],0,&atx[40]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[39],NULL,0,NULL} ,
  {407, "INSDAuthor" ,1,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[45]} ,
  {0, "consortium" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[41]} ,
  {0, "title" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[42]} ,
  {0, "journal" ,128,5,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[43]} ,
  {0, "xref" ,128,6,0,1,0,0,0,0,NULL,&atx[49],&atx[44],0,&atx[50]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[45],NULL,0,NULL} ,
  {408, "INSDXref" ,1,0,0,0,0,0,0,0,NULL,&atx[48],&atx[46],0,&atx[63]} ,
  {0, "dbname" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[47]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "pubmed" ,128,7,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[51]} ,
  {0, "remark" ,128,8,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "comment" ,128,23,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[53]} ,
  {0, "primary" ,128,24,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[54]} ,
  {0, "source-db" ,128,25,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[55]} ,
  {0, "database-reference" ,128,26,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[56]} ,
  {0, "feature-table" ,128,27,0,1,0,0,0,0,NULL,&atx[20],&atx[57],0,&atx[79]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[58],NULL,0,NULL} ,
  {406, "INSDFeature" ,1,0,0,0,0,0,0,0,NULL,&atx[48],&atx[59],0,&atx[39]} ,
  {0, "key" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[60]} ,
  {0, "location" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[61]} ,
  {0, "intervals" ,128,2,0,1,0,0,0,0,NULL,&atx[20],&atx[62],0,&atx[71]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[63],NULL,0,NULL} ,
  {409, "INSDInterval" ,1,0,0,0,0,0,0,0,NULL,&atx[48],&atx[64],0,&atx[76]} ,
  {0, "from" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[65]} ,
  {0, "to" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[66]} ,
  {0, "point" ,128,2,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[67]} ,
  {0, "iscomp" ,128,3,0,1,0,0,0,0,NULL,&atx[68],NULL,0,&atx[69]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "interbp" ,128,4,0,1,0,0,0,0,NULL,&atx[68],NULL,0,&atx[70]} ,
  {0, "accession" ,128,5,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "operator" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[72]} ,
  {0, "partial5" ,128,4,0,1,0,0,0,0,NULL,&atx[68],NULL,0,&atx[73]} ,
  {0, "partial3" ,128,5,0,1,0,0,0,0,NULL,&atx[68],NULL,0,&atx[74]} ,
  {0, "quals" ,128,6,0,1,0,0,0,0,NULL,&atx[20],&atx[75],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[76],NULL,0,NULL} ,
  {410, "INSDQualifier" ,1,0,0,0,0,0,0,0,NULL,&atx[48],&atx[77],0,&atx[81]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[78]} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "sequence" ,128,28,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[80]} ,
  {0, "contig" ,128,29,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {411, "INSDSet" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[82],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "INSD-INSDSeq" , "asninsdseq.h16",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = NULL;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module INSD-INSDSeq
*
**************************************************/

#define INSDSEQ &at[0]
#define INSDSEQ_locus &at[1]
#define INSDSEQ_length &at[3]
#define INSDSEQ_strandedness &at[5]
#define INSDSEQ_moltype &at[6]
#define INSDSEQ_topology &at[7]
#define INSDSEQ_division &at[8]
#define INSDSEQ_update_date &at[9]
#define INSDSEQ_create_date &at[10]
#define INSDSEQ_update_release &at[11]
#define INSDSEQ_create_release &at[12]
#define INSDSEQ_definition &at[13]
#define INSDSEQ_primary_accession &at[14]
#define INSDSEQ_entry_version &at[15]
#define INSDSEQ_accession_version &at[16]
#define INSDSEQ_other_seqids &at[17]
#define INSDSEQ_other_seqids_E &at[18]
#define INSDSEQ_secondary_accessions &at[21]
#define INSDSEQ_secondary_accessions_E &at[22]
#define INSDSEQ_project &at[24]
#define INSDSEQ_keywords &at[25]
#define INSDSEQ_keywords_E &at[26]
#define INSDSEQ_segment &at[28]
#define INSDSEQ_source &at[29]
#define INSDSEQ_organism &at[30]
#define INSDSEQ_taxonomy &at[31]
#define INSDSEQ_references &at[32]
#define INSDSEQ_references_E &at[33]
#define INSDSEQ_comment &at[52]
#define INSDSEQ_primary &at[53]
#define INSDSEQ_source_db &at[54]
#define INSDSEQ_database_reference &at[55]
#define INSDSEQ_feature_table &at[56]
#define INSDSEQ_feature_table_E &at[57]
#define INSDSEQ_sequence &at[79]
#define INSDSEQ_contig &at[80]

#define INSDSEQID &at[19]

#define INSDSECONDARY_ACCN &at[23]

#define INSDKEYWORD &at[27]

#define INSDREFERENCE &at[34]
#define INSDREFERENCE_reference &at[35]
#define INSDREFERENCE_position &at[36]
#define INSDREFERENCE_authors &at[37]
#define INSDREFERENCE_authors_E &at[38]
#define INSDREFERENCE_consortium &at[40]
#define INSDREFERENCE_title &at[41]
#define INSDREFERENCE_journal &at[42]
#define INSDREFERENCE_xref &at[43]
#define INSDREFERENCE_xref_E &at[44]
#define INSDREFERENCE_pubmed &at[50]
#define INSDREFERENCE_remark &at[51]

#define INSDFEATURE &at[58]
#define INSDFEATURE_key &at[59]
#define INSDFEATURE_location &at[60]
#define INSDFEATURE_intervals &at[61]
#define INSDFEATURE_intervals_E &at[62]
#define INSDFEATURE_operator &at[71]
#define INSDFEATURE_partial5 &at[72]
#define INSDFEATURE_partial3 &at[73]
#define INSDFEATURE_quals &at[74]
#define INSDFEATURE_quals_E &at[75]

#define INSDAUTHOR &at[39]

#define INSDXREF &at[45]
#define INSDXREF_dbname &at[46]
#define INSDXREF_id &at[47]

#define INSDINTERVAL &at[63]
#define INSDINTERVAL_from &at[64]
#define INSDINTERVAL_to &at[65]
#define INSDINTERVAL_point &at[66]
#define INSDINTERVAL_iscomp &at[67]
#define INSDINTERVAL_interbp &at[69]
#define INSDINTERVAL_accession &at[70]

#define INSDQUALIFIER &at[76]
#define INSDQUALIFIER_name &at[77]
#define INSDQUALIFIER_value &at[78]

#define INSDSET &at[81]
#define INSDSET_E &at[82]
