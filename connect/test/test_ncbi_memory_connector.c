/*  $Id: test_ncbi_memory_connector.c,v 6.5 2006/03/30 17:46:40 lavr Exp $
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
 *   Standard test for the MEMORY CONNECTOR
 *
 */

#include "ncbi_conntest.h"
#include <connect/ncbi_memory_connector.h>
#include <connect/ncbi_util.h>
/* This header must go last */
#include "test_assert.h"


int main(void)
{
    static STimeout timeout/* = 0 */;
    CONNECTOR       connector;
    FILE*           data_file;

    /* Log and data-log streams */
    CORE_SetLOGFILE(stderr, 0/*false*/);
    data_file = fopen("test_ncbi_memory_connector.log", "wb");
    assert(data_file);

     /* Run the tests */
    connector = MEMORY_CreateConnector();
    CONN_TestConnector(connector, &timeout, data_file, fTC_SingleBounceCheck);

    connector = MEMORY_CreateConnectorEx(0);
    CONN_TestConnector(connector, &timeout, data_file, fTC_Everything);

    /* Cleanup and Exit */
    fclose(data_file);
    CORE_SetLOG(0);
    return 0;
}


/*
 * --------------------------------------------------------------------------
 * $Log: test_ncbi_memory_connector.c,v $
 * Revision 6.5  2006/03/30 17:46:40  lavr
 * Adjust for lock-less MEMORY_Connector API
 *
 * Revision 6.4  2002/12/04 16:58:49  lavr
 * Move change log to end
 *
 * Revision 6.3  2002/03/22 19:47:31  lavr
 * Test_assert.h made last among the include files
 *
 * Revision 6.2  2002/02/20 20:53:48  lavr
 * Use xconntest to perform standard tests
 *
 * Revision 6.1  2002/02/20 19:14:40  lavr
 * Initial revision
 *
 * ==========================================================================
 */
