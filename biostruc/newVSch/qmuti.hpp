#ifndef _VSMmdbQmUti_HPP
#define _VSMmdbQmUti_HPP

#include <cgi/cgictx.hpp>
#include "qman.hpp"
#include "hVSMmdbcmd.hpp"
#include "hUtilib.hpp"


void ChkReqId(CCgiContext& ctx, QmReqID reqId);
void ChkGrpId(CCgiContext& ctx, QmReqID ori_gid, QmReqID new_gid);


#endif
