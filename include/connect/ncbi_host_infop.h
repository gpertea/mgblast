#ifndef CONNECT___NCBI_HOST_INFOP__H
#define CONNECT___NCBI_HOST_INFOP__H

/*  $Id: ncbi_host_infop.h,v 6.3 2006/03/05 17:37:14 lavr Exp $
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
 * Author:  Anton Lavrentiev, Denis Vakatov
 *
 * File Description:
 *   Private API for host information
 *
 */

#include <connect/ncbi_host_info.h>
#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif


HOST_INFO HINFO_Create(const void* hinfo, size_t hinfo_size,
                       const char* env, const char* arg, const char* val);


#ifdef __cplusplus
}
#endif


/*
 * --------------------------------------------------------------------------
 * $Log: ncbi_host_infop.h,v $
 * Revision 6.3  2006/03/05 17:37:14  lavr
 * New proto for HINFO_Create
 *
 * Revision 6.2  2002/10/28 21:55:38  lavr
 * LBSM_HINFO introduced for readability to replace plain "const void*"
 *
 * Revision 6.1  2002/10/28 20:13:45  lavr
 * Initial revision
 *
 * ==========================================================================
 */

#endif /* CONNECT___NCBI_HOST_INFOP__H */
