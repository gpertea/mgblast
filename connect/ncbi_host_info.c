/*  $Id: ncbi_host_info.c,v 6.10 2006/03/17 16:41:20 lavr Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   NCBI host info constructor and getters
 *
 */

#include "ncbi_lbsmd.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
/* Not defined on MacOS.9 :-( */
#  define M_PI 3.14159265358979323846
#endif


typedef struct SHostInfoTag {
    const char* env;
    const char* arg;
    const char* val;
    double      pad;    /* for proper 'hinfo' alignment; also as a magic */
} SHOST_Info;


HOST_INFO HINFO_Create(const void* hinfo, size_t hinfo_size, const char* env,
                       const char* arg, const char* val)
{
    SHOST_Info* host_info;
    size_t      size;
    size_t      e_s;
    size_t      a_s;
    size_t      v_s;
    char*       s;

    if (!hinfo)
        return 0;
    e_s = env && *env ? strlen(env) + 1 : 0;
    a_s = arg && *arg ? strlen(arg) + 1 : 0;
    v_s = a_s &&  val ? strlen(val) + 1 : 0;
    size = sizeof(*host_info) + hinfo_size;
    if (!(host_info = (SHOST_Info*) calloc(1, size + e_s + a_s + v_s)))
        return 0;
    memcpy((char*) host_info + sizeof(*host_info), hinfo, hinfo_size);
    s = (char*) host_info + size;
    if (e_s) {
        host_info->env = (const char*) memcpy(s, env, e_s);
        s += e_s;
    }
    if (a_s) {
        host_info->arg = (const char*) memcpy(s, arg, a_s);
        s += a_s;
    }
    if (v_s) {
        host_info->val = (const char*) memcpy(s, val, v_s);
        s += v_s;
    }
    host_info->pad = M_PI;
    return host_info;
}


int HINFO_CpuCount(const HOST_INFO host_info)
{
    if (!host_info || host_info->pad != M_PI)
        return -1;
    return LBSM_HINFO_CpuCount((const char*) host_info + sizeof(*host_info));
}


int HINFO_TaskCount(const HOST_INFO host_info)
{
    if (!host_info || host_info->pad != M_PI)
        return -1;
    return LBSM_HINFO_TaskCount((const char*) host_info + sizeof(*host_info));
}


int/*bool*/ HINFO_LoadAverage(const HOST_INFO host_info, double lavg[2])
{
    if (!host_info || host_info->pad != M_PI)
        return 0;
    return LBSM_HINFO_LoadAverage((const char*) host_info + sizeof(*host_info),
                                  lavg);
}


int/*bool*/ HINFO_Status(const HOST_INFO host_info, double status[2])
{
    if (!host_info || host_info->pad != M_PI)
        return 0;
    return LBSM_HINFO_Status((const char*) host_info + sizeof(*host_info),
                             status);
}


/*ARGSUSED*/
int/*bool*/ HINFO_BLASTParams(const HOST_INFO host_info, unsigned int blast[8])
{
    return 0;
}


const char* HINFO_Environment(const HOST_INFO host_info)
{
    if (!host_info || host_info->pad != M_PI)
        return 0;
    return host_info->env;
}


const char* HINFO_AffinityArgument(const HOST_INFO host_info)
{
    if (!host_info || host_info->pad != M_PI)
        return 0;
    return host_info->arg;
}


const char* HINFO_AffinityArgvalue(const HOST_INFO host_info)
{
    if (!host_info || host_info->pad != M_PI)
        return 0;
    return host_info->val;
}


/*
 * --------------------------------------------------------------------------
 * $Log: ncbi_host_info.c,v $
 * Revision 6.10  2006/03/17 16:41:20  lavr
 * Fix a typo in the last rev's change log
 *
 * Revision 6.9  2006/03/17 16:39:56  lavr
 * Explicit casts to keep compiler happy
 *
 * Revision 6.8  2006/03/06 20:24:20  lavr
 * Added "const" qualifier to all host-infos when passed to getters
 *
 * Revision 6.7  2006/03/05 17:36:52  lavr
 * +HINFO_AffinityArgument, +HINFO_AffinityArgvalue; HINFO_Create modified
 *
 * Revision 6.6  2003/01/17 19:44:46  lavr
 * Reduce dependencies
 *
 * Revision 6.5  2002/10/29 22:19:07  lavr
 * Fix typo in the file description
 *
 * Revision 6.4  2002/10/29 00:31:08  lavr
 * Fixed hinfo overflow from the use of precalculated size
 *
 * Revision 6.3  2002/10/28 21:55:38  lavr
 * LBSM_HINFO introduced for readability to replace plain "const void*"
 *
 * Revision 6.2  2002/10/28 20:49:04  lavr
 * Conditionally define M_PI if it is not already defined by <math.h>
 *
 * Revision 6.1  2002/10/28 20:13:45  lavr
 * Initial revision
 *
 * ==========================================================================
 */
