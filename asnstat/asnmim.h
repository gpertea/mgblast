/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnmim.h";
static AsnValxNode avnx[11] = {
    {20,"none" ,0,0.0,&avnx[1] } ,
    {20,"star" ,1,0.0,&avnx[2] } ,
    {20,"caret" ,2,0.0,&avnx[3] } ,
    {20,"pound" ,3,0.0,&avnx[4] } ,
    {20,"plus" ,4,0.0,&avnx[5] } ,
    {20,"perc" ,5,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[7] } ,
    {20,"citation" ,1,0.0,&avnx[8] } ,
    {20,"book" ,2,0.0,&avnx[9] } ,
    {20,"personal-communication" ,3,0.0,&avnx[10] } ,
    {20,"book-citation" ,4,0.0,NULL } };

static AsnType atx[129] = {
  {401, "Mim-entries" ,1,0,0,0,0,0,0,0,NULL,&atx[13],&atx[1],0,&atx[2]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {402, "Mim-entry" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[3],0,&atx[125]} ,
  {0, "mimNumber" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "mimType" ,128,1,0,0,0,0,0,0,NULL,&atx[6],&avnx[0],0,&atx[7]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "title" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[8]} ,
  {0, "copyright" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[9]} ,
  {0, "symbol" ,128,4,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[10]} ,
  {0, "locus" ,128,5,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[11]} ,
  {0, "synonyms" ,128,6,0,1,0,0,0,0,NULL,&atx[13],&atx[12],0,&atx[14]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "aliases" ,128,7,0,1,0,0,0,0,NULL,&atx[13],&atx[15],0,&atx[16]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "included" ,128,8,0,1,0,0,0,0,NULL,&atx[13],&atx[17],0,&atx[18]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "seeAlso" ,128,9,0,1,0,0,0,0,NULL,&atx[13],&atx[19],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,NULL} ,
  {405, "Mim-cit" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[21],0,&atx[28]} ,
  {0, "number" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[22]} ,
  {0, "author" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[23]} ,
  {0, "others" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[24]} ,
  {0, "year" ,128,3,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "text" ,128,10,0,1,0,0,0,0,NULL,&atx[13],&atx[27],0,&atx[36]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,NULL} ,
  {406, "Mim-text" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[29],0,&atx[44]} ,
  {0, "label" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[30]} ,
  {0, "text" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[31]} ,
  {0, "neighbors" ,128,2,0,1,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {411, "Mim-link" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[33],0,&atx[90]} ,
  {0, "num" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[34]} ,
  {0, "uids" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[35]} ,
  {0, "numRelevant" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "textfields" ,128,11,0,1,0,0,0,0,NULL,&atx[13],&atx[37],0,&atx[38]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,NULL} ,
  {0, "hasSummary" ,128,12,0,1,0,0,0,0,NULL,&atx[39],NULL,0,&atx[40]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "summary" ,128,13,0,1,0,0,0,0,NULL,&atx[13],&atx[41],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,NULL} ,
  {0, "summaryAttribution" ,128,14,0,1,0,0,0,0,NULL,&atx[13],&atx[43],0,&atx[51]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {407, "Mim-edit-item" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[45],0,&atx[56]} ,
  {0, "author" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[46]} ,
  {0, "modDate" ,128,1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {404, "Mim-date" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[48],0,&atx[20]} ,
  {0, "year" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[49]} ,
  {0, "month" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[50]} ,
  {0, "day" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "summaryEditHistory" ,128,15,0,1,0,0,0,0,NULL,&atx[13],&atx[52],0,&atx[53]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {0, "summaryCreationDate" ,128,16,0,1,0,0,0,0,NULL,&atx[44],NULL,0,&atx[54]} ,
  {0, "allelicVariants" ,128,17,0,1,0,0,0,0,NULL,&atx[13],&atx[55],0,&atx[66]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {408, "Mim-allelic-variant" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[57],0,&atx[69]} ,
  {0, "number" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[58]} ,
  {0, "name" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[59]} ,
  {0, "aliases" ,128,2,0,1,0,0,0,0,NULL,&atx[13],&atx[60],0,&atx[61]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "mutation" ,128,3,0,1,0,0,0,0,NULL,&atx[13],&atx[62],0,&atx[63]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,NULL} ,
  {0, "description" ,128,4,0,1,0,0,0,0,NULL,&atx[13],&atx[64],0,&atx[65]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,NULL} ,
  {0, "snpLinks" ,128,5,0,1,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {0, "hasSynopsis" ,128,18,0,1,0,0,0,0,NULL,&atx[39],NULL,0,&atx[67]} ,
  {0, "clinicalSynopsis" ,128,19,0,1,0,0,0,0,NULL,&atx[13],&atx[68],0,&atx[73]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[69],NULL,0,NULL} ,
  {409, "Mim-index-term" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[70],0,&atx[83]} ,
  {0, "key" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[71]} ,
  {0, "terms" ,128,1,0,0,0,0,0,0,NULL,&atx[13],&atx[72],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "synopsisAttribution" ,128,20,0,1,0,0,0,0,NULL,&atx[13],&atx[74],0,&atx[75]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {0, "synopsisEditHistory" ,128,21,0,1,0,0,0,0,NULL,&atx[13],&atx[76],0,&atx[77]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {0, "synopsisCreationDate" ,128,22,0,1,0,0,0,0,NULL,&atx[44],NULL,0,&atx[78]} ,
  {0, "editHistory" ,128,23,0,1,0,0,0,0,NULL,&atx[13],&atx[79],0,&atx[80]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {0, "creationDate" ,128,24,0,1,0,0,0,0,NULL,&atx[44],NULL,0,&atx[81]} ,
  {0, "references" ,128,25,0,1,0,0,0,0,NULL,&atx[13],&atx[82],0,&atx[117]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[83],NULL,0,NULL} ,
  {410, "Mim-reference" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[84],0,&atx[32]} ,
  {0, "number" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[85]} ,
  {0, "origNumber" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[86]} ,
  {0, "type" ,128,2,0,1,0,0,0,0,NULL,&atx[87],&avnx[6],0,&atx[88]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "authors" ,128,3,0,0,0,0,0,0,NULL,&atx[13],&atx[89],0,&atx[93]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[90],NULL,0,NULL} ,
  {412, "Mim-author" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[91],0,&atx[110]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[92]} ,
  {0, "index" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "primaryAuthor" ,128,4,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[94]} ,
  {0, "otherAuthors" ,128,5,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[95]} ,
  {0, "citationTitle" ,128,6,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[96]} ,
  {0, "citationType" ,128,7,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[97]} ,
  {0, "bookTitle" ,128,8,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[98]} ,
  {0, "editors" ,128,9,0,1,0,0,0,0,NULL,&atx[13],&atx[99],0,&atx[100]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[90],NULL,0,NULL} ,
  {0, "volume" ,128,10,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[101]} ,
  {0, "edition" ,128,11,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[102]} ,
  {0, "journal" ,128,12,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[103]} ,
  {0, "series" ,128,13,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[104]} ,
  {0, "publisher" ,128,14,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[105]} ,
  {0, "place" ,128,15,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[106]} ,
  {0, "commNote" ,128,16,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[107]} ,
  {0, "pubDate" ,128,17,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[108]} ,
  {0, "pages" ,128,18,0,1,0,0,0,0,NULL,&atx[13],&atx[109],0,&atx[113]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[110],NULL,0,NULL} ,
  {413, "Mim-page" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[111],0,NULL} ,
  {0, "from" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[112]} ,
  {0, "to" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "miscInfo" ,128,19,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[114]} ,
  {0, "pubmedUID" ,128,20,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[115]} ,
  {0, "ambiguous" ,128,21,0,0,0,0,0,0,NULL,&atx[39],NULL,0,&atx[116]} ,
  {0, "noLink" ,128,22,0,1,0,0,0,0,NULL,&atx[39],NULL,0,NULL} ,
  {0, "attribution" ,128,26,0,1,0,0,0,0,NULL,&atx[13],&atx[118],0,&atx[119]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {0, "numGeneMaps" ,128,27,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[120]} ,
  {0, "medlineLinks" ,128,28,0,1,0,0,0,0,NULL,&atx[32],NULL,0,&atx[121]} ,
  {0, "proteinLinks" ,128,29,0,1,0,0,0,0,NULL,&atx[32],NULL,0,&atx[122]} ,
  {0, "nucleotideLinks" ,128,30,0,1,0,0,0,0,NULL,&atx[32],NULL,0,&atx[123]} ,
  {0, "structureLinks" ,128,31,0,1,0,0,0,0,NULL,&atx[32],NULL,0,&atx[124]} ,
  {0, "genomeLinks" ,128,32,0,1,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {403, "Mim-set" ,1,0,0,0,0,0,0,0,NULL,&atx[25],&atx[126],0,&atx[47]} ,
  {0, "releaseDate" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[127]} ,
  {0, "mimEntries" ,128,1,0,0,0,0,0,0,NULL,&atx[13],&atx[128],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Mim" , "asnmim.h",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Mim
*
**************************************************/

#define MIM_ENTRIES &at[0]
#define MIM_ENTRIES_E &at[1]

#define MIM_ENTRY &at[2]
#define MIM_ENTRY_mimNumber &at[3]
#define MIM_ENTRY_mimType &at[5]
#define MIM_ENTRY_title &at[7]
#define MIM_ENTRY_copyright &at[8]
#define MIM_ENTRY_symbol &at[9]
#define MIM_ENTRY_locus &at[10]
#define MIM_ENTRY_synonyms &at[11]
#define MIM_ENTRY_synonyms_E &at[12]
#define MIM_ENTRY_aliases &at[14]
#define MIM_ENTRY_aliases_E &at[15]
#define MIM_ENTRY_included &at[16]
#define MIM_ENTRY_included_E &at[17]
#define MIM_ENTRY_seeAlso &at[18]
#define MIM_ENTRY_seeAlso_E &at[19]
#define MIM_ENTRY_text &at[26]
#define MIM_ENTRY_text_E &at[27]
#define MIM_ENTRY_textfields &at[36]
#define MIM_ENTRY_textfields_E &at[37]
#define MIM_ENTRY_hasSummary &at[38]
#define MIM_ENTRY_summary &at[40]
#define MIM_ENTRY_summary_E &at[41]
#define MIM_ENTRY_summaryAttribution &at[42]
#define MIM_ENTRY_summaryAttribution_E &at[43]
#define MIM_ENTRY_summaryEditHistory &at[51]
#define MIM_ENTRY_summaryEditHistory_E &at[52]
#define MIM_ENTRY_summaryCreationDate &at[53]
#define MIM_ENTRY_allelicVariants &at[54]
#define MIM_ENTRY_allelicVariants_E &at[55]
#define MIM_ENTRY_hasSynopsis &at[66]
#define MIM_ENTRY_clinicalSynopsis &at[67]
#define MIM_ENTRY_clinicalSynopsis_E &at[68]
#define MIM_ENTRY_synopsisAttribution &at[73]
#define MIM_ENTRY_synopsisAttribution_E &at[74]
#define MIM_ENTRY_synopsisEditHistory &at[75]
#define MIM_ENTRY_synopsisEditHistory_E &at[76]
#define MIM_ENTRY_synopsisCreationDate &at[77]
#define MIM_ENTRY_editHistory &at[78]
#define MIM_ENTRY_editHistory_E &at[79]
#define MIM_ENTRY_creationDate &at[80]
#define MIM_ENTRY_references &at[81]
#define MIM_ENTRY_references_E &at[82]
#define MIM_ENTRY_attribution &at[117]
#define MIM_ENTRY_attribution_E &at[118]
#define MIM_ENTRY_numGeneMaps &at[119]
#define MIM_ENTRY_medlineLinks &at[120]
#define MIM_ENTRY_proteinLinks &at[121]
#define MIM_ENTRY_nucleotideLinks &at[122]
#define MIM_ENTRY_structureLinks &at[123]
#define MIM_ENTRY_genomeLinks &at[124]

#define MIM_SET &at[125]
#define MIM_SET_releaseDate &at[126]
#define MIM_SET_mimEntries &at[127]
#define MIM_SET_mimEntries_E &at[128]

#define MIM_DATE &at[47]
#define MIM_DATE_year &at[48]
#define MIM_DATE_month &at[49]
#define MIM_DATE_day &at[50]

#define MIM_CIT &at[20]
#define MIM_CIT_number &at[21]
#define MIM_CIT_author &at[22]
#define MIM_CIT_others &at[23]
#define MIM_CIT_year &at[24]

#define MIM_TEXT &at[28]
#define MIM_TEXT_label &at[29]
#define MIM_TEXT_text &at[30]
#define MIM_TEXT_neighbors &at[31]

#define MIM_EDIT_ITEM &at[44]
#define MIM_EDIT_ITEM_author &at[45]
#define MIM_EDIT_ITEM_modDate &at[46]

#define MIM_ALLELIC_VARIANT &at[56]
#define MIM_ALLELIC_VARIANT_number &at[57]
#define MIM_ALLELIC_VARIANT_name &at[58]
#define MIM_ALLELIC_VARIANT_aliases &at[59]
#define MIM_ALLELIC_VARIANT_aliases_E &at[60]
#define MIM_ALLELIC_VARIANT_mutation &at[61]
#define MIM_ALLELIC_VARIANT_mutation_E &at[62]
#define MIM_ALLELIC_VARIANT_description &at[63]
#define ALLELIC_VARIANT_description_E &at[64]
#define MIM_ALLELIC_VARIANT_snpLinks &at[65]

#define MIM_INDEX_TERM &at[69]
#define MIM_INDEX_TERM_key &at[70]
#define MIM_INDEX_TERM_terms &at[71]
#define MIM_INDEX_TERM_terms_E &at[72]

#define MIM_REFERENCE &at[83]
#define MIM_REFERENCE_number &at[84]
#define MIM_REFERENCE_origNumber &at[85]
#define MIM_REFERENCE_type &at[86]
#define MIM_REFERENCE_authors &at[88]
#define MIM_REFERENCE_authors_E &at[89]
#define MIM_REFERENCE_primaryAuthor &at[93]
#define MIM_REFERENCE_otherAuthors &at[94]
#define MIM_REFERENCE_citationTitle &at[95]
#define MIM_REFERENCE_citationType &at[96]
#define MIM_REFERENCE_bookTitle &at[97]
#define MIM_REFERENCE_editors &at[98]
#define MIM_REFERENCE_editors_E &at[99]
#define MIM_REFERENCE_volume &at[100]
#define MIM_REFERENCE_edition &at[101]
#define MIM_REFERENCE_journal &at[102]
#define MIM_REFERENCE_series &at[103]
#define MIM_REFERENCE_publisher &at[104]
#define MIM_REFERENCE_place &at[105]
#define MIM_REFERENCE_commNote &at[106]
#define MIM_REFERENCE_pubDate &at[107]
#define MIM_REFERENCE_pages &at[108]
#define MIM_REFERENCE_pages_E &at[109]
#define MIM_REFERENCE_miscInfo &at[113]
#define MIM_REFERENCE_pubmedUID &at[114]
#define MIM_REFERENCE_ambiguous &at[115]
#define MIM_REFERENCE_noLink &at[116]

#define MIM_LINK &at[32]
#define MIM_LINK_num &at[33]
#define MIM_LINK_uids &at[34]
#define MIM_LINK_numRelevant &at[35]

#define MIM_AUTHOR &at[90]
#define MIM_AUTHOR_name &at[91]
#define MIM_AUTHOR_index &at[92]

#define MIM_PAGE &at[110]
#define MIM_PAGE_from &at[111]
#define MIM_PAGE_to &at[112]
